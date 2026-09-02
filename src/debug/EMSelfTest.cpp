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
        bool sawVanillaSpark = false;
        bool sawSparkFieldInject = false;
        float driftStartY = 0;

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
        if (frame >= 2200)
        {
                // safety net: never hang forever
                std::cerr << "[EMSELFTEST] FAIL: timed out" << std::endl;
                std::exit(1);
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
                // real charged particles
                sim->create_part(-1, 600, 100, PT_RELC);
                sim->create_part(-1, 600, 110, PT_RPRO);
                sim->create_part(-1, 600, 120, PT_RMON);
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
                // real particles registered
                Check(emf->realChargeCount == 3, "RPRO/RELC/RMON registered as real charges");
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
                // add / subtract accumulate onto the effective value
                int gi = emf->CellIndex(100, 180); // an EMPC cell: conductivity 1
                emf->ApplyEMProperty(EMADJP_CONDUCT, EMADJA_SET, gi, 1.0f);
                emf->ApplyEMProperty(EMADJP_CONDUCT, EMADJA_SUB, gi, 0.25f);
                Check(std::abs(emf->cells[gi].ovConduct - 0.75f) < 1e-3, "EMADJ subtract mode accumulates onto the effective value");
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
                auto *sim = gameModel.GetSimulation();
                // --- interop + conduction test setup ---
                // vanilla METL line with a real electron parked on it: the charge
                // must keep sparking the metal (EM -> vanilla direction)
                for (int x = 300; x < 340; ++x) sim->create_part(-1, x, 60, PT_METL);
                sim->create_part(-1, 320, 60, PT_RELC);
                // copper wire driven by an EMJP source at its top: a real electron
                // inside the wire must drift ALONG the wire (conduction)
                for (int y = 240; y < 300; ++y) sim->create_part(-1, 460, y, PT_CU);
                for (int y = 232; y < 240; ++y) sim->create_part(-1, 460, y, PT_EMJP);
                int driftEl = sim->create_part(-1, 460, 250, PT_RELC);
                if (driftEl >= 0)
                {
                        driftStartY = sim->parts.data[driftEl].y;
                }
                // annihilation: an electron/proton pair launched at each other must
                // neutralise when they meet
                int a1 = sim->create_part(-1, 200, 100, PT_RELC);
                int a2 = sim->create_part(-1, 210, 100, PT_RPRO);
                if (a1 >= 0) sim->parts[a1].vx = 0.3f;
                if (a2 >= 0) sim->parts[a2].vx = -0.3f;
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
        case 105:
        {
                std::cout << "[EMSELFTEST] options window survived NotifySettingsChanged" << std::endl;
                Check(sawVanillaSpark, "real charge on vanilla metal keeps it sparked (EM -> vanilla interop)");
                // Task 1: SPRK no longer excites the EM field - the test now verifies the
                // reverse direction is OFF, so the EMField is driven only by its own
                // sources and by real charges, not by vanilla sparks.
                Check(!sawSparkFieldInject, "powered vanilla spark does NOT inject current into the EM field (task 1)");
                // conduction drift: the electron on the driven copper wire moved along it
                {
                        auto *sim = gameModel.GetSimulation();
                        auto *emf = sim->GetEMField();
                        float moved = 0;
                        for (int i = 0; i < sim->parts.active; ++i)
                        {
                                auto &p = sim->parts.data[i];
                                if (p.type == PT_RELC && std::abs(p.x - 460.0f) < 3.0f)
                                {
                                        moved = std::max(moved, std::abs(p.y - driftStartY));
                                }
                        }
                        std::cout << "[EMSELFTEST] info: drift distance " << moved << "px" << std::endl;
                        Check(moved >= 2.0f, "real charge drifts along a driven conductor (current conducts)");
                        // annihilation: both carriers of the launched pair are gone
                        int left = 0;
                        for (int i = 0; i < sim->parts.active; ++i)
                        {
                                auto &p = sim->parts.data[i];
                                if ((p.type == PT_RELC || p.type == PT_RPRO) &&
                                        p.y > 90 && p.y < 110 && p.x > 180 && p.x < 230)
                                {
                                        left++;
                                }
                        }
                        Check(left == 0, "opposite carriers annihilate on contact");
                        (void)emf;
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
                // the velocity clamp: give an electron a huge kick and verify the
                // next frames clamp it back to the field propagation speed
                int el = -1;
                for (int i = 0; i < parts.active; ++i)
                {
                        if (parts.data[i].type == PT_RELC)
                        {
                                el = i;
                                break;
                        }
                }
                if (el >= 0)
                {
                        parts.data[el].vx = 50.0f; // way beyond the speed of light
                        parts.data[el].vy = 50.0f;
                        float vmax = emf->MaxParticleSpeed();
                        std::cout << "[EMSELFTEST] info: max particle speed " << vmax << " px/frame" << std::endl;
                }
                else
                {
                        Check(false, "real electron present for the FTL clamp test");
                }
                gameModel.SetEMViewMode(EMVIEW_E);
                break;
        }
        case 135:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                float vmax = emf->MaxParticleSpeed();
                auto &parts = sim->parts;
                bool clamped = true;
                for (int i = 0; i < parts.active; ++i)
                {
                        auto &p = parts.data[i];
                        if (p.type == PT_RPRO || p.type == PT_RELC || p.type == PT_RMON)
                        {
                                float sp = std::sqrt(p.vx * p.vx + p.vy * p.vy);
                                if (sp > vmax * 1.001f)
                                {
                                        clamped = false;
                                }
                        }
                }
                Check(clamped, "real particles never exceed the field propagation speed (FTL clamp)");
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
                                t == PT_RPRO || t == PT_RELC || t == PT_RMON ||
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
        case 1000:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                openE0 = InjectPulse(*emf, emf->visW / 2);
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
                Check(e < 0.15 * openE0, "OPEN boundary absorbs the pulse outside the screen (no return)");
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
                gameModel.SetEMBoundaryMode(EMBND_ABSORB); // restore the default
                gameModel.SetEMSpeed(1);
                break;
        }
        case 1620:
        {
                // 1-pixel EM grid smoke test (设置 -> 电磁场网格大小 -> 1 像素):
                // one cell per particle, full alignment with matter
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                gameModel.SetEMCellSize(1);
                Check(emf->visW == XRES && emf->visH == YRES && emf->padL == 0,
                        "1px EM grid: one cell per particle (absorb boundary)");
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
                Check(emf->padL == EM_OPEN_PAD_MAX && emf->gw == XRES + 2 * emf->padL,
                        "1px EM grid + OPEN boundary: capped invisible pad");
                gameModel.SetEMBoundaryMode(EMBND_ABSORB);
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
        // continuous interop observation between the setup (65) and the check (105):
        // a real charge parked on vanilla metal must keep it sparked, and a live
        // spark must inject current into the EM field
        if (frame >= 66 && frame < 105)
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                for (int x = 300; x < 340; ++x)
                {
                        unsigned r = sim->pmap[60][x];
                        if (r && TYP(r) == PT_SPRK)
                        {
                                sawVanillaSpark = true;
                                int gi = emf->CellIndex(x, 60);
                                if (emf->cells[gi].jzext > 0.01)
                                {
                                        sawSparkFieldInject = true;
                                }
                        }
                }
        }
        frame++;
}

#endif // EMFIELD_DEBUG
