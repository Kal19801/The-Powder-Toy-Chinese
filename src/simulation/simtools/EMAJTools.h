#pragma once
#include "simulation/ToolCommon.h"
#include "simulation/EMField.h"
#include "simulation/Simulation.h"

// Shared implementation of the EMAJC / EMAJP / EMAJJ / EMAJD / EMAJS tools,
// the ports of the five non-directional adjust modes (MODE_ADJ_CONDUCT,
// MODE_ADJ_PERM, MODE_ADJ_J, MODE_ADJ_MEDIUM, MODE_ADJ_MAG_STR) of Paul
// Falstad's EMWave2 applet; MODE_ADJ_MAG_DIR is the EMMD tool. The tool
// strength slider plays the role of the applet's adjustBar.

namespace EMAJAdjust
{
        static inline bool perform(Simulation * sim, int mode, int x, int y, float strength)
        {
                auto *emf = sim->GetEMField();
                if (!emf || !emf->enabled)
                        return false;
                int gi = emf->CellIndex(x, y);
                if (emf->ApplyAdjustMode(mode, gi, strength))
                {
                        // the wave equation needs to know about the material change right away
                        // task 10: defer CalcBoundaries to the next frame (see EMADJ)
                        emf->NotifyCellChanged();
                        return true;
                }
                return false;
        }
}
