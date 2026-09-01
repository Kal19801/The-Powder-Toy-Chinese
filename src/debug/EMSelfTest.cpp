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
}

void EMSelfTestTick(GameModel &gameModel, GameController &gameController)
{
        if (frame >= 200)
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
                // a bit of everything: metal bar, glass slab, iron bar, an EMW source,
                // a resonant slab and a small antenna
                for (int x = 80; x < 130; ++x) sim->create_part(-1, x, 180, PT_METL);
                for (int x = 160; x < 210; ++x) sim->create_part(-1, x, 180, PT_GLAS);
                for (int x = 240; x < 290; ++x) sim->create_part(-1, x, 180, PT_IRON);
                sim->create_part(-1, 320, 180, PT_EMW);
                for (int x = 350; x < 420; ++x) sim->create_part(-1, x, 250, PT_EMR);
                for (int x = 430; x < 480; ++x) sim->create_part(-1, x, 180, PT_WIRE);
                for (int x = 430; x < 480; ++x) sim->create_part(-1, x, 200, PT_METL);
                // the applet material modes ported as elements: current sources and magnets
                for (int x = 500; x < 550; ++x) sim->create_part(-1, x, 180, PT_EMJP);
                for (int x = 500; x < 550; ++x) sim->create_part(-1, x, 220, PT_EMJN);
                for (int x = 80; x < 140; ++x) sim->create_part(-1, x, 140, PT_EMMG);
                std::cout << "[EMSELFTEST] setup complete" << std::endl;
                break;
        }
        case 25:
        {
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double maxdazdt = 0;
                int conductors = 0, mediums = 0, ferro = 0, resonant = 0, sources = 0;
                for (auto &cell : emf->cells)
                {
                        maxdazdt = std::max(maxdazdt, std::abs(cell.dazdt));
                        if (cell.conductivity > 0) conductors++;
                        if (cell.medium > 0) mediums++;
                        if (cell.perm > 1) ferro++;
                        if (cell.resonant) resonant++;
                        if (cell.jzext != 0) sources++;
                }
                Check(maxdazdt > 0.02, "EM wave energy present in the field");
                Check(conductors >= 30, "TPT metal mapped to EM conductor cells");
                Check(mediums >= 10, "TPT glass mapped to EM dielectric cells");
                Check(ferro >= 10, "TPT iron mapped to EM ferromagnet cells");
                Check(resonant >= 10, "EMR element mapped to resonant cells");
                Check(sources >= 3, "global/EMW/EMJP/EMJN sources inject current");
                // EMMG magnets must map to magnetized cells (mx/my = +-1)
                int magnets = 0;
                for (auto &cell : emf->cells)
                {
                        if (cell.mx != 0 || cell.my != 0)
                                magnets++;
                }
                Check(magnets >= 10, "EMMG element mapped to magnet cells");
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
                Check(anyMag, "EMMD magnetization applies to ferromagnetic cells");
                Check(anyAdj, "EMADJ parameter adjustment applies to material cells");
                // the five separate applet adjust modes (EMAJ* tools)
                bool adjC = false, adjP = false, adjJ = false, adjD = false, adjS = false;
                for (int gi = 0; gi < emf->gw * emf->gh; ++gi)
                {
                        adjC |= emf->ApplyAdjustMode(EMADJM_CONDUCT, gi, 0.5f);
                        adjP |= emf->ApplyAdjustMode(EMADJM_PERM, gi, 0.5f);
                        adjJ |= emf->ApplyAdjustMode(EMADJM_J, gi, 0.5f);
                        adjD |= emf->ApplyAdjustMode(EMADJM_MEDIUM, gi, 0.5f);
                        adjS |= emf->ApplyAdjustMode(EMADJM_MAG_STR, gi, 0.5f);
                }
                Check(adjC, "EMAJC conductivity adjust applies to conductor cells");
                Check(adjP, "EMAJP permeability adjust applies to ferromagnet cells");
                Check(adjJ, "EMAJJ current adjust applies to current cells");
                Check(adjD, "EMAJD dielectric adjust applies to medium cells");
                Check(adjS, "EMAJS strength adjust applies to magnet cells");
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
                // The old crash happened synchronously inside OptionsController's
                // constructor (NotifySettingsChanged), so if OpenOptions() returns we
                // know the crash is gone.
                gameController.OpenOptions();
                Check(true, "options window constructed and survived NotifySettingsChanged");
                // note: while the modal window is open, GameController::Tick is not
                // called anymore (ShowWindow replaces the engine state); the window
                // closes itself again from OptionsView::OnTick's debug block below.
                break;
        }
        case 105:
        {
                std::cout << "[EMSELFTEST] options window survived NotifySettingsChanged" << std::endl;
                break;
        }
        case 125:
        {
                // field -> particle interactions
                auto *sim = gameModel.GetSimulation();
                auto *emf = sim->GetEMField();
                double maxMetalTemp = 0;
                int sparks = 0;
                auto &parts = sim->parts;
                for (int i = 0; i < parts.active; ++i)
                {
                        auto &p = parts.data[i];
                        if (p.type == PT_METL || p.type == PT_WIRE || p.type == PT_IRON)
                        {
                                maxMetalTemp = std::max(maxMetalTemp, double(p.temp));
                        }
                        if (p.type == PT_SPRK)
                        {
                                sparks++;
                        }
                }
                std::cout << "[EMSELFTEST] info: max metal temp " << maxMetalTemp
                          << "K, induced sparks so far: " << sparks << std::endl;
                Check(maxMetalTemp > 295.5, "Joule heating raises conductor temperature");
                gameModel.SetEMViewMode(EMVIEW_E);
                break;
        }
        case 145:
        {
                // new-save regression: a full simulation reset (new save, loaded
                // save, clear) must wipe every EM tool override, otherwise EMADJ/
                // EMMD painting leaks into the next save forever
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
                for (auto &cell : emf->cells)
                {
                        anyOv |= cell.ovMask != 0;
                }
                Check(!anyOv, "clear_sim (new save) clears all EM overrides");
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
