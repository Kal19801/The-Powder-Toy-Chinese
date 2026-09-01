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
        frame++;
}
