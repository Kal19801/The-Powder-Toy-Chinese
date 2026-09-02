#ifndef EMFIELD_DEBUG
#define EMFIELD_DEBUG 0
#endif
#if EMFIELD_DEBUG

#include "EMSelfTest.h"

#include "gui/game/GameController.h"
#include "gui/game/GameModel.h"
#include "gui/interface/Engine.h"
#include "simulation/EMField.h"
#include "simulation/ElementClasses.h"
#include "simulation/Simulation.h"
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <iostream>

// EM field self test, a debug block (see the EMFIELD_DEBUG guard in
// GameController.cpp): drives the rewritten current system from the main loop
// and asserts its invariants. Compiled out of normal builds.

namespace
{
        int frame = 0;
        int failCount = 0;
        double closedE0 = 0;
        double openE0 = 0;
        // --- superconductor regression state ---
        float scTestY0 = 0, scControlY0 = 0;
        double lambdaPxSmall = 0;
        int emtxCellGi = 0;

        // --- EMADJ continuous add/subtract test state (frames 165..236) ---
        // Reproduces the tool-level 0.2s repeat throttle of
        // EMAdjustTool::ApplyBrush (first application immediate, then every 12
        // frames while held) around the REAL per-frame game loop, so the
        // accumulation runs through ApplyEMProperty -> SyncMaterials exactly
        // like in-game. One job per property/material cell rect.
        struct EmAdjJob
        {
                int property;
                int mode;
                float value;
                int x0, y0, x1, y1; // particle rect of the test material
                const char *what;
        };
        const EmAdjJob emAdjJobs[] = {
                { EMADJP_CONDUCT, EMADJA_ADD, 0.07f, 120, 360, 136, 372, "conduct ADD on a fair conductor" },
                { EMADJP_CONDUCT, EMADJA_ADD, 0.07f, 200, 360, 216, 372, "conduct ADD on a ferromagnetic conductor" },
                { EMADJP_PERM,    EMADJA_ADD, 0.50f, 200, 360, 216, 372, "perm ADD on a ferromagnet" },
                { EMADJP_CONDUCT, EMADJA_SUB, 0.07f, 280, 360, 296, 372, "conduct SUB on a perfect conductor" },
                { EMADJP_MEDIUM,  EMADJA_SUB, 10.0f, 360, 360, 376, 372, "medium SUB on a dielectric" },
                { EMADJP_J,       EMADJA_ADD, 0.10f, 440, 360, 456, 372, "current ADD on a current source" },
                { EMADJP_MAG_DIR, EMADJA_ADD, 90.0f, 520, 360, 536, 372, "mag direction ADD on a magnet" },
                { EMADJP_MAG_STR, EMADJA_ADD, 0.10f, 520, 360, 536, 372, "mag strength ADD on a magnet" },
        };
        int emAdjRepeatCounter = 0;
        bool emAdjStrokeActive = false;

        void EmAdjStrokeTick(EMField *emf)
        {
                // verbatim replica of the EMAdjustTool::ApplyBrush add/subtract
                // branch: apply the whole brush, then wait 12 frames
                if (!emAdjStrokeActive || emAdjRepeatCounter <= 0)
                {
                        for (const auto &job : emAdjJobs)
                        {
                                for (int y = job.y0; y < job.y1; ++y)
                                {
                                        for (int x = job.x0; x < job.x1; ++x)
                                        {
                                                emf->ApplyEMProperty(job.property, job.mode, emf->CellIndex(x, y), job.value);
                                        }
                                }
                        }
                        emAdjStrokeActive = true;
                        emAdjRepeatCounter = 12;
                }
                else
                {
                        emAdjRepeatCounter--;
                }
        }

        void Check(bool cond, const char *what)
        {
                if (cond)
                {
                        std::cout << "[EMSELFTEST] PASS: " << what << std::endl;
                }
                else
                {
                        failCount++;
                        std::cerr << "[EMSELFTEST] FAIL: " << what << std::endl;
                }
        }

        // fresh superconductor samples: SCND is solid so the probes stay put
        void sim_create_scnd_cold(GameModel &gameModel, int x, int y)
        {
                auto *sim = gameModel.GetSimulation();
                int i = sim->create_part(-1, x, y, PT_SCND);
                if (i >= 0)
                {
                        sim->parts[i].temp = 4.0f; // deep below Tc
                }
        }
        void sim_create_scnd_warm(GameModel &gameModel, int x, int y)
        {
                auto *sim = gameModel.GetSimulation();
                int i = sim->create_part(-1, x, y, PT_SCND);
                if (i >= 0)
                {
                        sim->parts[i].temp = 295.0f; // room temperature, above Tc
                }
        }

        // spatial wavelength in px of the wave travelling along row y between
        // x0 and x1, from the zero crossings of az (2 crossings per period)
        double WaveLambdaPx(const EMField &emf, int y, int x0, int x1)
        {
                int crossings = 0;
                double firstX = -1, lastX = -1;
                for (int x = x0; x < x1; x++)
                {
                        double a = emf.cells[emf.CellIndex(x, y)].az;
                        double b = emf.cells[emf.CellIndex(x + 1, y)].az;
                        if ((a < 0 && b >= 0) || (a > 0 && b <= 0))
                        {
                                crossings++;
                                if (firstX < 0) firstX = x;
                                lastX = x;
                        }
                }
                if (crossings < 4)
                        return 0;
                return 2.0 * (lastX - firstX) / (crossings - 1);
        }

        double FieldEnergy(const EMField &emf)
        {
                double energy = 0;
                for (int y = 1; y < emf.gh - 1; y++)
                {
                        for (int x = 1; x < emf.gw - 1; x++)
                        {
                                int gi = x + y * emf.gw;
                                const auto &cell = emf.cells[gi];
                                energy += 0.5 * cell.dazdt * cell.dazdt + 0.5 * cell.az * cell.az;
                        }
                }
                return energy;
        }

        // inject a smooth planar WAVE-VELOCITY (dazdt) pulse centred on visible
        // cell column cx, confined to the material-free band y 40..170 so it never
        // crosses any conductor (the planar symmetry keeps it inside the band).
        // The profile sums to zero so the pulse is a pure pair of travelling waves
        // with no stationary (DC) component that would relax in place; returns the
        // field energy right after the injection
        double InjectPulse(EMField &emf, int cx)
        {
                // zero-sum derivative-shaped profile: smooth enough that its group
                // velocity stays near the lattice maximum (no slowly-crawling
                // high-k tail) and with no DC component (nothing relaxes in place)
                static const double profile[7] = { -0.125, -0.35, -0.75, 0, 0.75, 0.35, 0.125 };
                for (int y = 40; y < 170; ++y)
                {
                        for (int dx = 0; dx < 7; ++dx)
                        {
                                int px = (cx + dx - 3) * emf.cellSize;
                                emf.cells[emf.CellIndex(px, y)].dazdt = profile[dx];
                        }
                }
                return FieldEnergy(emf);
        }
}

void EMSelfTestTick(GameModel &gameModel, GameController &gameController)
{
        if (frame >= 4200)
        {
                // safety net: never hang forever
                std::cerr << "[EMSELFTEST] FAIL: timed out" << std::endl;
                std::exit(1);
        }
        // EMADJ continuous add/subtract stroke (frames 166..225): runs BEFORE
        // the sim tick of every held frame, exactly like the real
        // input-apply -> simulate order of the main loop
        if (frame >= 166 && frame <= 225)
        {
                EmAdjStrokeTick(gameModel.GetSimulation()->GetEMField());
        }
        switch (frame)
        {
        case 5:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                emf->enabled = true;
                gameModel.SetEMCellSize(4);
                gameModel.SetEMSourceMode(EMSRC_1S1F);
                gameModel.SetEMFrequency(10);
                gameModel.SetEMViewMode(EMVIEW_E_B_J);
                // the applet material modes ported as elements: conductors, magnets,
                // dielectric, resonant medium and the EMW source
                for (int x = 80; x < 130; ++x) sim->create_part(-1, x, 180, PT_EMPC);
                for (int x = 160; x < 210; ++x) sim->create_part(-1, x, 180, PT_EMEC);
                for (int x = 240; x < 290; ++x) sim->create_part(-1, x, 180, PT_EMFC);
                sim->create_part(-1, 320, 180, PT_EMW);
                for (int x = 350; x < 420; ++x) sim->create_part(-1, x, 250, PT_EMR);
                for (int x = 430; x < 480; ++x) sim->create_part(-1, x, 180, PT_EMFM);
                for (int x = 430; x < 480; ++x) sim->create_part(-1, x, 200, PT_EMDM);
                for (int x = 500; x < 550; ++x) sim->create_part(-1, x, 180, PT_EMJP);
                for (int x = 500; x < 550; ++x) sim->create_part(-1, x, 220, PT_EMJN);
                // a copper bar right below the positive current source picks up the
                // induced current and must heat up (Joule losses in real conductors)
                for (int x = 500; x < 550; ++x) sim->create_part(-1, x, 184, PT_CU);
                for (int x = 80; x < 140; ++x) sim->create_part(-1, x, 140, PT_EMMGD);
                for (int x = 80; x < 140; ++x) sim->create_part(-1, x, 240, PT_EMMGU);
                for (int x = 560; x < 610; ++x) sim->create_part(-1, x, 180, PT_EMDE);
                // real zone materials and powders
                for (int x = 80; x < 130; ++x) sim->create_part(-1, x, 280, PT_FE);
                for (int x = 160; x < 210; ++x) sim->create_part(-1, x, 280, PT_CU);
                for (int x = 240; x < 290; ++x) sim->create_part(-1, x, 280, PT_AG);
                for (int x = 320; x < 370; ++x) sim->create_part(-1, x, 280, PT_TI);
                for (int x = 400; x < 450; ++x) sim->create_part(-1, x, 280, PT_SCND);
                for (int x = 480; x < 530; ++x) sim->create_part(-1, x, 280, PT_PGRF);
                for (int x = 80; x < 130; ++x) sim->create_part(-1, x, 320, PT_FEPW);
                for (int x = 160; x < 210; ++x) sim->create_part(-1, x, 320, PT_PGPW);
                // superconductor must start below Tc to superconduct
                for (int x = 240; x < 290; ++x)
                {
                        int i = sim->create_part(-1, x, 320, PT_SCPW);
                        if (i >= 0)
                        {
                                sim->parts[i].temp = 4.0f;
                        }
                }
                std::cout << "[EMSELFTEST] setup complete" << std::endl;
                break;
        }
        case 25:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double maxdazdt = 0;
                int conductors = 0, mediums = 0, ferro = 0, diamag = 0, resonant = 0, sources = 0;
                for (auto &cell : emf->cells)
                {
                        maxdazdt = std::max(maxdazdt, std::abs(cell.dazdt));
                        if (cell.conductivity > 0) conductors++;
                        if (cell.medium > 0) mediums++;
                        if (cell.perm > 1) ferro++;
                        if (cell.perm < 1 && cell.perm != 1) diamag++;
                        if (cell.resonant) resonant++;
                        if (cell.jzext != 0) sources++;
                }
                Check(maxdazdt > 0.02, "EM wave energy present in the field");
                Check(conductors >= 150, "EM conductor elements + real conductors mapped to conductor cells");
                Check(mediums >= 10, "EMDE element mapped to dielectric cells");
                Check(ferro >= 10, "EMFM/FE/FEPW mapped to ferromagnet cells");
                Check(diamag >= 10, "EMDM/PGRF/PGPW/SCPW mapped to diamagnet cells");
                Check(resonant >= 10, "EMR element mapped to resonant cells");
                Check(sources >= 3, "global/EMW/EMJP/EMJN sources inject current");
                // the four magnet elements must map to magnetized cells
                int magnets = 0;
                for (auto &cell : emf->cells)
                {
                        if (cell.mx != 0 || cell.my != 0)
                                magnets++;
                }
                Check(magnets >= 20, "EMMGD/EMMGU elements mapped to magnet cells");
                gameModel.SetEMViewMode(EMVIEW_B_LINES); // exercise field line rendering
                break;
        }
        case 45:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                bool anyMag = false;
                bool anyAdj = false;
                for (int gi = 0; gi < emf->gw * emf->gh; ++gi)
                {
                        anyMag |= emf->ApplyMagDir(gi, 0.25f);
                        anyAdj |= emf->ApplyAdjust(gi, 0.5f);
                }
                Check(anyMag, "magnetization applies to ferromagnetic cells");
                Check(anyAdj, "combined parameter adjustment applies to material cells");
                // the unified EM adjust tool: all six properties and all three modes
                bool set[EMADJP_COUNT] = {};
                for (int prop = 0; prop < EMADJP_COUNT; ++prop)
                {
                        for (int gi = 0; gi < emf->gw * emf->gh; ++gi)
                        {
                                set[prop] |= emf->ApplyEMProperty(prop, EMADJA_SET, gi, 0.5f);
                        }
                }
                Check(set[EMADJP_CONDUCT], "EMADJ set mode: conductivity applies to conductors");
                Check(set[EMADJP_PERM], "EMADJ set mode: permeability applies to ferromagnets");
                Check(set[EMADJP_J], "EMADJ set mode: current applies to current sources");
                Check(set[EMADJP_MEDIUM], "EMADJ set mode: dielectric applies to medium cells");
                Check(set[EMADJP_MAG_DIR], "EMADJ set mode: mag direction applies to magnets");
                Check(set[EMADJP_MAG_STR], "EMADJ set mode: mag strength applies to magnets");
                // task 10: the J property also excites a current INSIDE a plain
                // EM-zone conductor (EMPC has conductivity 1 and no jzext)
                {
                        int condGi = emf->CellIndex(100, 180);
                        bool applied = emf->ApplyEMProperty(EMADJP_J, EMADJA_SET, condGi, 0.5f);
                        Check(applied, "EMADJ J mode accepts conductor cells (task 10)");
                }
                // add / subtract accumulate onto the effective value
                int gi = emf->CellIndex(100, 180); // an EMPC cell: conductivity 1
                emf->ApplyEMProperty(EMADJP_CONDUCT, EMADJA_SET, gi, 1.0f);
                emf->ApplyEMProperty(EMADJP_CONDUCT, EMADJA_SUB, gi, 0.25f);
                Check(std::abs(emf->cells[gi].ovConduct - 0.75f) < 1e-3, "EMADJ subtract mode accumulates onto the effective value");
                // task 10: after a SyncMaterials pass the conductor J override
                // must surface as an injected current in that cell
                {
                        // run one Update() manually so the override takes effect
                        // even though the EM field may be disabled in this stage
                        bool wasEnabled = emf->enabled;
                        emf->enabled = true;
                        emf->Update();
                        emf->enabled = wasEnabled;
                        float jz = float(emf->cells[emf->CellIndex(100, 180)].jzext);
                        Check(std::abs(jz - 0.5f) < 1e-3,
                                "conductor J override injects current inside the conductor (task 10)");
                }
                emf->ApplyEMProperty(EMADJP_CONDUCT, EMADJA_ADD, gi, 0.25f);
                Check(std::abs(emf->cells[gi].ovConduct - 1.0f) < 1e-3, "EMADJ add mode accumulates onto the effective value");
                // VacuumCell / ClearCellOverrides (EMCLR + erase integration)
                emf->VacuumCell(emf->CellIndex(100, 180));
                Check(emf->cells[emf->CellIndex(100, 180)].ovMask == 0, "EMCLR vacuum clears cell overrides");
                gameModel.SetEMViewMode(EMVIEW_FORCE); // exercise force flood fill
                break;
        }
        case 65:
        {
                gameModel.SetEMViewMode(EMVIEW_TYPE); // exercise material palette
                // fresh superconductor bars for the case-105 mapping checks
                // (created AFTER the case-45 adjust sweep so no tool override
                // pollutes them; SCND is solid and does not move)
                for (int x = 160; x < 200; ++x) sim_create_scnd_cold(gameModel, x, 250);
                for (int x = 220; x < 260; ++x) sim_create_scnd_warm(gameModel, x, 250);
                break;
        }
        case 85:
        {
                // THE CRASH REGRESSION TEST: opening the settings window is the exact
                // code path ("点击设置") that used to die with "Memory read/write error"
                // on handheld builds, where the language dropdown is never created.
                gameController.OpenOptions();
                Check(true, "options window constructed and survived NotifySettingsChanged");
                break;
        }
        case 104:
        {
                // pin the sample temperatures right before the mapping check:
                // TPT heat exchange with the 295 K ambient would otherwise drift
                // the 4 K bar towards Tc over the 40 frames in between
                auto *sim = gameModel.GetSimulation();
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        auto &p = sim->parts.data[i];
                        if (p.type != PT_SCND)
                        {
                                continue;
                        }
                        if (p.x >= 160 && p.x < 200 && std::abs(p.y - 250) < 3)
                        {
                                p.temp = 4.0f;
                        }
                        else if (p.x >= 220 && p.x < 260 && std::abs(p.y - 250) < 3)
                        {
                                p.temp = 295.0f;
                        }
                }
                break;
        }
        case 105:
        {
                std::cout << "[EMSELFTEST] options window survived NotifySettingsChanged" << std::endl;
                // superconductor temperature mapping, probed on the fresh bars
                // created at case 65 (temps re-pinned at case 104)
                {
                        auto *sim = gameModel.GetSimulation();
                        auto *emf = sim->GetEMField();
                        auto &cold = emf->cells[emf->CellIndex(180, 250)];
                        std::cout << "[EMSELFTEST] info: cold SCND cell cond=" << cold.conductivity
                                  << " perm=" << cold.perm << std::endl;
                        Check(cold.conductivity >= 1.0f - 1e-6f,
                                "superconductor below Tc maps to a perfect conductor");
                        Check(cold.perm <= EM_SC_PERM + 1e-6f && cold.perm > 0.0f,
                                "superconductor below Tc maps to the diamagnetic Meissner response");
                        auto &warm = emf->cells[emf->CellIndex(240, 250)];
                        std::cout << "[EMSELFTEST] info: warm SCND cell cond=" << warm.conductivity
                                  << " perm=" << warm.perm << std::endl;
                        Check(std::abs(warm.conductivity - EM_SC_QUENCH_CONDUCT) < 1e-6f,
                                "superconductor above Tc is quenched to a lossy conductor");
                        Check(std::abs(warm.perm - 1.0f) < 1e-6f,
                                "quenched superconductor carries no Meissner response");
                }
                break;
        }
        case 125:
        {
                // field -> particle interactions on the REAL system only: the vanilla
                // elements are reverted and must not be affected by the field
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double maxRealTemp = 0;
                auto &parts = sim->parts;
                for (int i = 0; i < parts.active; ++i)
                {
                        auto &p = parts.data[i];
                        if (p.type == PT_CU || p.type == PT_AG || p.type == PT_FE)
                        {
                                maxRealTemp = std::max(maxRealTemp, double(p.temp));
                        }
                }
                std::cout << "[EMSELFTEST] info: max real conductor temp " << maxRealTemp << "K" << std::endl;
                Check(maxRealTemp > 295.5, "Joule heating raises real conductor temperature");
                (void)emf;
                gameModel.SetEMViewMode(EMVIEW_E);
                break;
        }
        case 135:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                (void)sim;
                // energy boundedness: the field energy must stay finite over the run
                double energy = FieldEnergy(*emf);
                std::cout << "[EMSELFTEST] info: field energy " << energy << std::endl;
                Check(std::isfinite(energy) && energy < 1e6, "field energy stays bounded with sources on (no divergence)");
                break;
        }
        case 145:
        {
                // divergence soak test: a ferromagnet next to a driven source must
                // not blow up over a long run (the old self-excitation regression)
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                gameModel.SetEMSpeed(2); // fastest sub-step setting, worst CFL case
                double e0 = FieldEnergy(*emf);
                std::cout << "[EMSELFTEST] info: energy before soak " << e0 << std::endl;
                gameModel.SetEMViewMode(EMVIEW_E_B_J);
                break;
        }
        case 345:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double energy = FieldEnergy(*emf);
                std::cout << "[EMSELFTEST] info: energy after 200-frame soak at 2x speed " << energy << std::endl;
                Check(std::isfinite(energy) && energy < 1e6, "200-frame ferromagnet soak at 2x speed stays bounded");
                Check(emf->fieldClampHits == 0, "wave state never touched the safety clamp");
                gameModel.SetEMSpeed(1);
                break;
        }
        case 160:
        {
                // --- EMADJ continuous add/subtract test (frames 166..236) ---
                // one block of each material the adjust tool works on; the
                // stroke applies to all of them simultaneously with the tool's
                // 0.2s repeat cadence
                auto *sim = gameModel.GetSimulation();
                for (int y = 360; y < 372; ++y)
                {
                        for (int x = 120; x < 136; ++x) sim->create_part(-1, x, y, PT_EMFC); // fair conductor, cond .5
                        for (int x = 200; x < 216; ++x) sim->create_part(-1, x, y, PT_EMFM); // ferromagnetic conductor, perm 5 cond .5
                        for (int x = 280; x < 296; ++x) sim->create_part(-1, x, y, PT_EMPC); // perfect conductor, cond 1
                        for (int x = 360; x < 376; ++x) sim->create_part(-1, x, y, PT_EMDE); // dielectric, medium 191
                        for (int x = 440; x < 456; ++x) sim->create_part(-1, x, y, PT_EMJP); // +1 current source
                        for (int x = 520; x < 536; ++x) sim->create_part(-1, x, y, PT_EMMGR); // magnet pointing right
                }
                emAdjRepeatCounter = 0;
                emAdjStrokeActive = false;
                break;
        }
        case 226:
        {
                // 60 held frames with the 12-frame repeat -> 5 applications
                // (frames 166, 178, 190, 202, 214); every job must have moved
                // its property step by step, not just once ("加法模式无法连续
                // 增加" regression)
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                {
                        auto &c = emf->cells[emf->CellIndex(128, 366)];
                        Check(std::abs(c.conductivity - 0.85f) < 0.01f, "EMADJ add mode: 5x +0.07 on a fair conductor = 0.85");
                }
                {
                        auto &c = emf->cells[emf->CellIndex(208, 366)];
                        // the old else-if chain dropped the conduct override on
                        // every ferromagnet, so ADD could never raise it
                        Check(std::abs(c.conductivity - 0.85f) < 0.01f, "EMADJ add mode: conductivity rises continuously on a ferromagnetic conductor");
                        Check(std::abs(c.perm - 7.5f) < 0.05f, "EMADJ add mode: 5x +0.5 on a ferromagnet = 7.5");
                }
                {
                        auto &c = emf->cells[emf->CellIndex(288, 366)];
                        Check(std::abs(c.conductivity - 0.65f) < 0.01f, "EMADJ subtract mode: 5x -0.07 on a perfect conductor = 0.65");
                }
                {
                        auto &c = emf->cells[emf->CellIndex(368, 366)];
                        Check(c.medium == 141, "EMADJ subtract mode: 5x -10 on a dielectric = 141");
                }
                {
                        auto &c = emf->cells[emf->CellIndex(448, 366)];
                        Check(std::abs(float(c.jzext) - 1.5f) < 0.01f, "EMADJ add mode: 5x +0.1 on a current source = 1.5");
                }
                {
                        // 5x +90deg = +450deg = +90deg after wrap; EMMGR starts
                        // at mx=1 (angle 0), 90deg means mx=0, my=-|m|; the
                        // strength job raised |m| to 1.5 in the same stroke
                        auto &c = emf->cells[emf->CellIndex(528, 366)];
                        std::cout << "[EMSELFTEST] info: magnet cell mx=" << c.mx << " my=" << c.my
                                  << " ovDir=" << c.ovDir << " ovMask=" << int(c.ovMask)
                                  << " perm=" << c.perm << " cond=" << c.conductivity << std::endl;
                        Check(std::abs(c.mx) < 0.02f && std::abs(c.my + 1.5f) < 0.02f,
                                "EMADJ add mode: 5x +90deg rotates the magnet direction step by step");
                        Check(std::abs(std::sqrt(c.mx * c.mx + c.my * c.my) - 1.5f) < 0.02f,
                                "EMADJ add mode: 5x +0.1 magnet strength = 1.5, direction kept");
                }
                break;
        }
        case 236:
        {
                // remove the test materials again so they cannot pollute the
                // pulse / boundary / soak measurements below
                auto *sim = gameModel.GetSimulation();
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        int t = sim->parts.data[i].type;
                        if ((t == PT_EMFC || t == PT_EMFM || t == PT_EMPC ||
                                t == PT_EMDE || t == PT_EMJP || t == PT_EMMGR) &&
                                sim->parts.data[i].y >= 350)
                        {
                                sim->kill_part(i);
                        }
                }
                break;
        }
        case 2000:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double energy = FieldEnergy(*emf);
                std::cout << "[EMSELFTEST] info: energy after 1650-frame soak at 1x speed " << energy << std::endl;
                Check(std::isfinite(energy) && energy < 2e5, "long soak stays bounded (steady state, no divergence)");
                Check(emf->fieldClampHits == 0, "wave state never touched the safety clamp (long run)");
                break;
        }
        case 365:
        {
                // EMCLR must only clear the field, never delete particles
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                int before = 0;
                auto &parts = sim->parts;
                for (int i = 0; i < parts.active; ++i)
                {
                        if (parts.data[i].type)
                        {
                                before++;
                        }
                }
                int cx = emf->CellIndex(100, 180) % emf->gw;
                int cy = emf->CellIndex(100, 180) / emf->gw;
                int px = cx * emf->cellSize;
                int py = cy * emf->cellSize;
                // call VacuumCell exactly like the EMCLR tool does
                emf->VacuumCell(emf->CellIndex(px, py));
                int after = 0;
                for (int i = 0; i < parts.active; ++i)
                {
                        if (parts.data[i].type)
                        {
                                after++;
                        }
                }
                Check(before == after, "EMCLR clears only the field, no particles are deleted");
                break;
        }
        case 400:
        {
                // --- boundary condition tests (封闭/开放/循环) ---
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                // remove every field source AND every real charge so the pulse tests
                // are pure decay / propagation measurements: EMW/EMJP/EMJN inject
                // current, the permanent magnets establish a growing static field,
                // interop sparks and moving charges deposit currents - all of that
                // would drown the pulse
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        int t = sim->parts.data[i].type;
                        if (t == PT_EMW || t == PT_EMJP || t == PT_EMJN ||
                                t == PT_EMMGD || t == PT_EMMGU || t == PT_EMMGL ||
                                t == PT_EMMGR || t == PT_EMMG)
                        {
                                sim->kill_part(i);
                        }
                }
                gameModel.SetEMSpeed(2); // faster transit for the pulse tests
                gameModel.SetEMSourceMode(EMSRC_NONE); // no pumping, pure pulse
                gameModel.SetEMBoundaryMode(EMBND_CLOSED);
                Check(emf->padL == 0 && emf->visW == XRES / emf->cellSize,
                        "CLOSED boundary: grid == visible canvas");
                break;
        }
        case 650:
        {
                // inject after the startup filter window (filterCount > 200 by now)
                // so the smoothing cannot eat the pulse
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                closedE0 = InjectPulse(*emf, emf->visW / 2);
                break;
        }
        case 760:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double e = FieldEnergy(*emf);
                std::cout << "[EMSELFTEST] info: closed-boundary energy " << e << " (start " << closedE0 << ")" << std::endl;
                Check(e > 0.5 * closedE0, "CLOSED boundary conserves the pulse (full reflection)");
                break;
        }
        case 770:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                gameModel.SetEMBoundaryMode(EMBND_OPEN);
                Check(emf->padL > 0 && emf->gw == emf->visW + 2 * emf->padL,
                        "OPEN boundary: invisible absorber band outside the screen");
                break;
        }
        case 775:
        {
                // the freshly reallocated field only carries the 1e-10 vacuum
                // seed: the PML activity gate must rate the band quiet and skip
                // its passes (this is what makes OPEN as cheap as CLOSED)
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                Check(emf->pmlQuiet,
                        "PML activity gate: resting band is quiet (band passes skipped)");
                break;
        }
        case 780:
        {
                // debug benchmark: per-frame EM cost of the OPEN mode (quiet
                // gate, the common case) vs CLOSED vs OPEN with a wave forced
                // into the band, same grid resolution; information only
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                auto bench = [emf](const char *name) {
                        for (int f = 0; f < 30; f++) emf->Update(); // warmup
                        auto t0 = std::chrono::steady_clock::now();
                        for (int f = 0; f < 300; f++)
                        {
                                emf->Update();
                        }
                        auto t1 = std::chrono::steady_clock::now();
                        double us = std::chrono::duration<double, std::micro>(t1 - t0).count() / 300.0;
                        std::cout << "[EMSELFTEST] bench: " << name << " grid " << emf->gw << "x" << emf->gh
                                  << " pad " << emf->padL << " -> " << us << " us/frame" << std::endl;
                };
                bench("OPEN (quiet band, gate skipping)");
                gameModel.SetEMBoundaryMode(EMBND_CLOSED);
                bench("CLOSED (hard wall)             ");
                gameModel.SetEMBoundaryMode(EMBND_PERIODIC);
                bench("PERIODIC (ghost ring)          ");
                gameModel.SetEMBoundaryMode(EMBND_OPEN);
                // force the gate active so the PML passes actually run: put a
                // small wave in the guard strip
                for (int y = emf->padT + 2; y < emf->gh - emf->padT - 2; y += 3)
                {
                        emf->cells[emf->padL + 4 + y * emf->gw].az = 0.5;
                }
                emf->ScanPmlActivity();
                bench("OPEN (wave at the edge, PML on)");
                // the benchmark left wave energy in the field - reallocate the
                // grid (CLOSED -> OPEN switch) so the pulse tests below start
                // from a clean vacuum state
                gameModel.SetEMBoundaryMode(EMBND_CLOSED);
                gameModel.SetEMBoundaryMode(EMBND_OPEN);
                break;
        }
        case 985:
        {
                // gate timing probe: a tiny pulse at the exact canvas centre is
                // dozens of cells away from every guard strip, so the gate must
                // still rate the band quiet (the main pulse below overlaps the
                // top guard strip by design of InjectPulse's row band and trips
                // the gate immediately)
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                int gi = emf->CellIndex(XRES / 2, YRES / 2);
                emf->cells[gi].dazdt = 1.0;
                break;
        }
        case 995:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                Check(emf->pmlQuiet,
                        "PML activity gate: mid-canvas wave does not trip the band");
                break;
        }
        case 1000:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                openE0 = InjectPulse(*emf, emf->visW / 2);
                break;
        }
        case 1400:
        {
                // 400 frames after injection the pulse has crossed the visible
                // canvas (0.25 cells/frame at cs=4, speed 2) and entered the
                // band: the gate must have woken the PML up well before that
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                // TEMP DIAG
                double guardMax = 0, bandMax = 0, visMax = 0;
                int px=0, py=0; double pk=0;
                for (int j = 0; j < emf->gh; ++j)
                {
                        for (int i = 0; i < emf->gw; ++i)
                        {
                                int gi = i + j * emf->gw;
                                double v = std::max(std::abs(emf->cells[gi].az), std::abs(emf->cells[gi].dazdt));
                                bool inBand = i < emf->padL || i >= emf->gw - emf->padL || j < emf->padT || j >= emf->gh - emf->padT;
                                bool inGuard = (i >= emf->padL && i < emf->padL + 8) || (i >= emf->gw - emf->padL - 8 && i < emf->gw - emf->padL) ||
                                               (j >= emf->padT && j < emf->padT + 8) || (j >= emf->gh - emf->padT - 8 && j < emf->gh - emf->padT);
                                if (inBand) bandMax = std::max(bandMax, v);
                                else if (inGuard) guardMax = std::max(guardMax, v);
                                else visMax = std::max(visMax, v);
                                if (v > pk) { pk = v; px = i; py = j; }
                        }
                }
                std::cout << "[EMSELFTEST] info: gate=" << emf->pmlQuiet << " pad=" << emf->padL
                          << " gw=" << emf->gw << " guardMax=" << guardMax << " bandMax=" << bandMax
                          << " visMax=" << visMax << " peak=" << px << "," << py << " pk=" << pk << std::endl;
                Check(!emf->pmlQuiet,
                        "PML activity gate: wave reaching the edge wakes the band up");
                break;
        }
        case 1500:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double e = FieldEnergy(*emf);
                double vis = 0, pad = 0;
                double vmax = 0;
                int vmaxx = 0, vmaxy = 0;
                for (int y = 1; y < emf->gh - 1; ++y)
                {
                        for (int x = 1; x < emf->gw - 1; ++x)
                        {
                                const auto &c = emf->cells[x + y * emf->gw];
                                double e2 = 0.5 * c.dazdt * c.dazdt + 0.5 * c.az * c.az;
                                bool inView = x >= emf->padL && x < emf->padL + emf->visW &&
                                              y >= emf->padT && y < emf->padT + emf->visH;
                                if (inView) vis += e2; else pad += e2;
                                if (e2 > vmax) { vmax = e2; vmaxx = x; vmaxy = y; }
                        }
                }
                std::cout << "[EMSELFTEST] info: open-boundary energy " << e << " (start " << openE0 << ")"
                          << " vis=" << vis << " pad=" << pad
                          << " peak cell " << vmaxx << "," << vmaxy << " e2=" << vmax << std::endl;
                // The PML is a MATCHED layer: it absorbs gradually with depth, so
                // mid-band energy at this frame is by design (the pulse is still
                // inside the quartic profile). What must hold is:
                //  1. the wave LEFT the visible canvas (vis residue is the
                //     dispersed tail only) - a reflected wave would re-light it;
                //  2. the band is net-absorbing, not feeding (total decays).
                Check(vis < 0.1 * openE0 && e < 0.6 * openE0,
                        "OPEN boundary lets the pulse leave the screen without return (PML)");
                break;
        }
        case 1510:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                gameModel.SetEMBoundaryMode(EMBND_PERIODIC);
                Check(emf->padL == 1 && emf->gw == emf->visW + 2,
                        "PERIODIC boundary: one cell ghost ring");
                // inject near the LEFT edge: the left-moving half wraps almost
                // immediately and must reappear at the RIGHT edge, which in CLOSED
                // mode stays quiet until the wave crosses the whole screen
                InjectPulse(*emf, 6);
                break;
        }
        case 1545:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                // the pulse left through the left edge and must have wrapped to the
                // right edge (in CLOSED it would need a full screen traversal first)
                double wrapped = 0;
                for (int y = 1; y < emf->gh - 1; ++y)
                {
                        for (int x = emf->gw - 15; x < emf->gw - 1; ++x)
                        {
                                int gi = x + y * emf->gw;
                                wrapped = std::max(wrapped, std::abs(emf->cells[gi].az));
                                wrapped = std::max(wrapped, std::abs(emf->cells[gi].dazdt));
                        }
                }
                std::cout << "[EMSELFTEST] info: wrapped amplitude at left edge " << wrapped << std::endl;
                Check(wrapped > 0.05, "PERIODIC boundary wraps waves across the seam");
                gameModel.SetEMBoundaryMode(EMBND_OPEN); // restore the default
                gameModel.SetEMSpeed(1);
                break;
        }
        case 1610:
        {
                // region size setting (task 7): the simulated domain expands
                // beyond the visible canvas and stays centred on it
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                gameModel.SetEMRegionScale(2);
                Check(emf->visW == 2 * XRES / emf->cellSize && emf->visH == 2 * YRES / emf->cellSize,
                        "region scale 2 doubles the simulated domain (task 7)");
                Check(emf->renderOffX == emf->padL + (emf->visW - XRES / emf->cellSize) / 2 &&
                      emf->renderOffY == emf->padT + (emf->visH - YRES / emf->cellSize) / 2,
                        "visible canvas stays centred in the enlarged domain (task 7)");
                gameModel.SetEMRegionScale(1);
                Check(emf->visW == XRES / emf->cellSize,
                        "region scale 1 restores the visible canvas as the domain");
                // --- 0.5x: the domain is HALF the canvas per axis, rendered
                // magnified 2x so the field still covers the full screen
                // (0.5x fullscreen fix - the space is half, not the display)
                gameModel.SetEMRegionScale(0.5f);
                Check(emf->visW == XRES / 2 / emf->cellSize && emf->visH == YRES / 2 / emf->cellSize,
                        "region scale 0.5 shrinks the simulated space to half the canvas extent");
                Check(emf->renderScale == 2 &&
                      emf->renderW == XRES / (emf->cellSize * 2) &&
                      emf->renderH == YRES / (emf->cellSize * 2),
                        "region scale 0.5 renders the half domain magnified 2x (full-screen field)");
                Check(emf->renderOffX == emf->padL && emf->renderOffY == emf->padT,
                        "region scale 0.5 maps the whole domain onto the whole canvas");
                Check(emf->CellIndex(0, 0) == emf->padL + emf->padT * emf->gw,
                        "region scale 0.5: canvas corner maps onto the domain corner cell");
                Check(emf->CellIndex(XRES - 1, YRES - 1) ==
                      (emf->padL + emf->visW - 1) + (emf->padT + emf->visH - 1) * emf->gw,
                        "region scale 0.5: canvas maps onto every domain cell (zoomed 2x)");
                Check(emf->PixelInDomain(0, 0) && emf->PixelInDomain(XRES - 1, YRES - 1),
                        "region scale 0.5: the full screen carries the magnified field");
                gameModel.SetEMRegionScale(1);
                Check(emf->visW == XRES / emf->cellSize && emf->renderScale == 1,
                        "region scale 1 restores again after 0.5x");
                break;
        }
        case 1620:
        {
                // 1-pixel EM grid smoke test (设置 -> 电磁场网格大小 -> 1 像素):
                // one cell per particle, full alignment with matter
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                gameModel.SetEMCellSize(1);
                Check(emf->taddEff > EM_TADD_SUB + 1e-6f,
                        "1px EM grid: wave timestep shrunk (task 7 resolution decoupling)");
                Check(emf->visW == XRES && emf->visH == YRES &&
                      emf->renderOffX == emf->padL && emf->renderOffY == emf->padT,
                        "1px EM grid: one cell per particle, outflow band outside the canvas");
                emf->cells[emf->CellIndex(300, 100)].dazdt = 1.0; // small probe pulse
                break;
        }
        case 1680:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double e = FieldEnergy(*emf);
                Check(std::isfinite(e) && e < 1e5,
                        "1px EM grid stays bounded over 60 frames");
                gameModel.SetEMBoundaryMode(EMBND_OPEN);
                Check(emf->padL == EM_PAD_PX / emf->cellSize && emf->gw == XRES + 2 * emf->padL,
                        "1px EM grid + OPEN boundary: fixed 64px pixel-width pad");
                gameModel.SetEMBoundaryMode(EMBND_OPEN);
                gameModel.SetEMCellSize(EM_CELL_SIZE_DEFAULT);
                // every geometry switch above reallocated the field and wiped the
                // tool overrides; write one directly so the reset test below has
                // something to check for (ApplyAdjust would see a freshly cleared
                // field with no materials synced yet and do nothing)
                {
                        auto &cell = emf->cells[emf->CellIndex(100, 180)];
                        cell.ovMask |= EM_OV_CONDUCT;
                        cell.ovConduct = 0.5f;
                }
                break;
        }
        case 2100:
        {
                // new-save regression: a full simulation reset (new save, loaded
                // save, clear) must wipe every EM tool override and current, and
                // the erase tool must leave no invisible EM residue behind
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                bool anyOv = false;
                for (auto &cell : emf->cells)
                {
                        anyOv |= cell.ovMask != 0;
                }
                Check(anyOv, "EM overrides present before the reset");
                sim->clear_sim();
                anyOv = false;
                double leftoverCurrent = 0;
                for (auto &cell : emf->cells)
                {
                        anyOv |= cell.ovMask != 0;
                        leftoverCurrent += std::abs(cell.jz) + std::abs(cell.jzext);
                }
                Check(!anyOv, "clear_sim (new save) clears all EM overrides");
                Check(leftoverCurrent == 0, "clear_sim drops all currents");
                // --- superconductor rig: a cold SCPW grain falls past a vertical
                // magnet column; the diamagnetic (Meissner) pressure pushes it
                // AWAY from the stronger field, i.e. horizontally off its
                // straight-down path. A control grain falls the same distance
                // with no magnet nearby and drifts nowhere.
                for (int y = 140; y < 300; ++y)
                {
                        for (int x = 140; x < 147; ++x)
                        {
                                sim->create_part(-1, x, y, PT_EMMGD);
                        }
                }
                // catch floors: grains reaching the canvas bottom are deleted,
                // so land them on thick diamond shelves instead (a 1 px shelf can
                // be tunnelled through by a fast powder)
                for (int x = 140; x < 320; ++x)
                {
                        for (int y = 343; y < 348; ++y) sim->create_part(-1, x, y, PT_DMND);
                }
                for (int x = 330; x < 460; ++x)
                {
                        for (int y = 343; y < 348; ++y) sim->create_part(-1, x, y, PT_DMND);
                }
                {
                        int i = sim->create_part(-1, 152, 132, PT_SCPW);
                        if (i >= 0)
                        {
                                sim->parts[i].temp = 4.0f;
                                sim->parts[i].vx = 0;
                                sim->parts[i].vy = 0;
                        }
                }
                {
                        int i = sim->create_part(-1, 392, 132, PT_SCPW);
                        if (i >= 0)
                        {
                                sim->parts[i].temp = 4.0f;
                                sim->parts[i].vx = 0;
                                sim->parts[i].vy = 0;
                        }
                }
                scTestY0 = 152.0f;  // reuse: test start x
                scControlY0 = 392.0f; // control start x
                sim->emf->enabled = true;
                break;
        }
        case 2180:
        {
                // --- superconductor checks: the cold grain that fell past the
                // magnet column is deflected sideways; the control grain is not ---
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                (void)emf;
                float testX = 0, testN = 0, controlX = 0, controlN = 0;
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        auto &pt = sim->parts.data[i];
                        if (pt.type != PT_SCPW)
                        {
                                continue;
                        }
                        if (pt.x >= 147 && pt.x < 320)
                        {
                                testX += pt.x;
                                testN += 1;
                        }
                        else if (pt.x >= 330 && pt.x < 460)
                        {
                                controlX += pt.x;
                                controlN += 1;
                        }
                }
                if (testN > 0 && controlN > 0)
                {
                        testX /= testN;
                        controlX /= controlN;
                        float testD = testX - scTestY0;
                        float controlD = controlX - scControlY0;
                        std::cout << "[EMSELFTEST] info: SCPW past magnet x " << scTestY0 << " -> " << testX
                                  << " (drift " << testD << "), control x " << scControlY0 << " -> " << controlX
                                  << " (drift " << controlD << ")" << std::endl;
                        Check(testD > 3.0f && testD > controlD + 2.0f,
                                "cold superconductor powder is pushed away from a magnet (Meissner repulsion)");
                }
                else
                {
                        int shown = 0;
                        for (int i = 0; i < sim->parts.active && shown < 12; ++i)
                        {
                                auto &pt = sim->parts.data[i];
                                if (pt.type == PT_SCPW)
                                {
                                        std::cout << "[EMSELFTEST] info: SCPW at " << pt.x << "," << pt.y
                                                  << " temp=" << pt.temp << std::endl;
                                        shown++;
                                }
                        }
                        Check(false, "superconductor powders survived the rig");
                }
                // --- task 6 setup: vanilla -> EMTX transmitter rig ---
                // METL line sparked from its left end; EMTX watches its right end
                for (int x = 212; x <= 230; ++x) sim->create_part(-1, x, 80, PT_METL);
                int tx = sim->create_part(-1, 231, 80, PT_EMTX);
                if (tx >= 0)
                {
                        sim->parts[tx].ctype = 20; // carrier frequency 20 -> ~107 frame burst
                        emtxCellGi = emf->CellIndex(231, 80);
                }
                // spark the METL line (vanilla conversion of the conductor)
                {
                        unsigned r = sim->pmap[80][215];
                        if (r)
                        {
                                int ri = ID(r);
                                sim->parts[ri].ctype = PT_METL;
                                sim->part_change_type(ri, 215, 80, PT_SPRK);
                                sim->parts[ri].life = 4;
                        }
                }
                break;
        }
        case 2270:
        {
                // --- task 6 checks: burst started and reaches the field ---
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                // find the EMTX particle
                int found = -1;
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        if (sim->parts.data[i].type == PT_EMTX)
                        {
                                found = i;
                        }
                }
                Check(found >= 0, "EMTX transmitter exists");
                if (found >= 0)
                {
                        std::cout << "[EMSELFTEST] info: EMTX tmp=" << sim->parts.data[found].tmp << std::endl;
                        Check(sim->parts.data[found].tmp > 0,
                                "EMTX starts a burst when the adjacent conductor is sparked (task 6)");
                        Check(std::abs(emf->cells[emtxCellGi].jzext) > 1e-3,
                                "EMTX burst reaches the field as an injected current (task 6)");
                }
                break;
        }
        case 2400:
        {
                // --- task 6: the burst is finite and stops by itself ---
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                int found = -1;
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        if (sim->parts.data[i].type == PT_EMTX)
                                found = i;
                }
                if (found >= 0)
                {
                        Check(sim->parts.data[found].tmp == 0,
                                "EMTX burst ends by itself (finite wave packet, task 6)");
                }
                Check(std::abs(emf->cells[emtxCellGi].jzext) < 1e-9,
                        "EMTX injects nothing while idle (task 6)");
                // --- task 5 setup: EMJP source -> EMEC bar -> EMAN -> METL ---
                for (int x = 100; x <= 104; ++x) sim->create_part(-1, x, 300, PT_EMJP);
                for (int x = 108; x <= 141; ++x) sim->create_part(-1, x, 300, PT_EMEC);
                sim->create_part(-1, 143, 300, PT_EMAN);
                for (int x = 144; x <= 165; ++x) sim->create_part(-1, x, 300, PT_METL); // directly adjacent
                break;
        }
        case 2520:
        {
                // --- task 5 check: the antenna fired the vanilla conductor ---
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                (void)emf;
                int sparks = 0;
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        auto &pt = sim->parts.data[i];
                        if (pt.type == PT_SPRK && pt.ctype == PT_METL &&
                            pt.y > 295 && pt.y < 305 && pt.x > 140 && pt.x < 170)
                        {
                                sparks++;
                        }
                }
                std::cout << "[EMSELFTEST] info: sparks on the antenna-fed METL line: " << sparks << std::endl;
                Check(sparks > 0,
                        "EMAN sparks the adjacent vanilla conductor when the connected EM elements carry current (task 5)");
                // --- task 7 setup: wavelength probe at cellSize 2 ---
                sim->clear_sim();
                gameModel.SetEMBoundaryMode(EMBND_CLOSED);
                gameModel.SetEMCellSize(2);
                gameModel.SetEMSourceMode(EMSRC_NONE);
                sim->create_part(-1, 306, 192, PT_EMW, 20); // frequency 20
                sim->emf->enabled = true;
                break;
        }
        case 3020:
        {
                // 530 frames after the source went live: the wave has filled
                // x 306..~570 but has NOT yet reached the CLOSED side walls
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                lambdaPxSmall = WaveLambdaPx(*emf, 192, 400, 550);
                std::cout << "[EMSELFTEST] info: lambda at cellSize 2: " << lambdaPxSmall << "px" << std::endl;
                Check(lambdaPxSmall > 15 && lambdaPxSmall < 45,
                        "wave wavelength is in the expected pixel range at cellSize 2");
                // --- task 7: repeat the probe at cellSize 4 ---
                sim->clear_sim();
                gameModel.SetEMCellSize(4);
                Check(emf->taddEff < EM_TADD_SUB - 1e-6f,
                        "wave timestep shrinks for the coarser grid per the decoupling rule (task 7)");
                sim->create_part(-1, 306, 192, PT_EMW, 20);
                sim->emf->enabled = true;
                break;
        }
        case 3560:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double lambdaPxBig = WaveLambdaPx(*emf, 192, 400, 550);
                std::cout << "[EMSELFTEST] info: lambda at cellSize 4: " << lambdaPxBig << "px" << std::endl;
                Check(lambdaPxBig > 15 && lambdaPxBig < 45,
                        "wave wavelength is in the expected pixel range at cellSize 4");
                if (lambdaPxSmall > 0 && lambdaPxBig > 0)
                {
                        double rel = std::abs(lambdaPxBig - lambdaPxSmall) / lambdaPxSmall;
                        std::cout << "[EMSELFTEST] info: relative wavelength difference across resolutions: "
                                  << rel * 100 << "%" << std::endl;
                        Check(rel < 0.25,
                                "wave wavelength in px is (nearly) independent of the grid resolution (task 7)");
                }
                if (failCount)
                {
                        std::cerr << "[EMSELFTEST] RESULT: " << failCount << " FAILURES" << std::endl;
                        std::exit(1);
                }
                std::cout << "[EMSELFTEST] RESULT: ALL TESTS PASSED" << std::endl;
                ui::Engine::Ref().Exit();
                break;
        }
        }
        // cryostat for the superconductor levitation rig (frames 2100..2179):
        // an uncooled 4 K powder warms past Tc in about a second in the 295 K
        // ambient (correct physics - real superconductors need a coolant bath),
        // so the rig re-pins the temperature every frame to isolate the
        // Meissner force from the thermal quench
        if (frame >= 2100 && frame <= 2180)
        {
                auto *sim = gameModel.GetSimulation();
                for (int i = 0; i < sim->parts.active; ++i)
                {
                        auto &p = sim->parts.data[i];
                        if (p.type == PT_SCPW)
                        {
                                p.temp = 4.0f;
                        }
                }
        }
        frame++;
}

#endif // EMFIELD_DEBUG
