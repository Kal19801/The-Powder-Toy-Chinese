#pragma once
#include "simulation/ToolCommon.h"
#include "simulation/EMField.h"
#include "simulation/Simulation.h"

// Shared implementation of the EMMGD / EMMGU / EMMGL / EMMGR tools, the ports
// of the MODE_M_DOWN / MODE_M_UP / MODE_M_LEFT / MODE_M_RIGHT paint modes of
// Paul Falstad's EMWave2 applet; each draws EMMG particles with the matching
// .ctype magnetization direction (0 = down, 1 = up, 2 = left, 3 = right).

namespace EMMGPaint
{
        static inline int perform(Simulation * sim, Particle * cpart, int x, int y, int ctype)
        {
                if (cpart && cpart->type == PT_EMMG)
                {
                        if (cpart->ctype != ctype)
                        {
                                cpart->ctype = ctype;
                                return 1;
                        }
                        return 0;
                }
                if (cpart && cpart->type)
                {
                        // replace whatever is under the brush, just like drawing an element
                        sim->delete_part(x, y);
                }
                int i = sim->create_part(-1, x, y, PT_EMMG);
                if (i >= 0)
                {
                        sim->parts[i].ctype = ctype;
                        return 1;
                }
                return 0;
        }
}
