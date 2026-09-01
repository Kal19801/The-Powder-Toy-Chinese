#include "EMMGTools.h"

// EMMGR: paint permanent magnets pointing right, port of the MODE_M_RIGHT
// mode of Paul Falstad's EMWave2 applet.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMMGR()
{
        Identifier = "DEFAULT_TOOL_EMMGR";
        Name = "EMMGR";
        Colour = 0xC060E0_rgb;
        Description = ByteString("磁铁(向右):画出指向右方的永磁体,相当于applet的M RIGHT模式").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMMGPaint::perform(sim, cpart, x, y, 3); // right
}
