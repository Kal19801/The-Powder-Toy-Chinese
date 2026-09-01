#include "EMMGTools.h"

// EMMGU: paint permanent magnets pointing up, port of the MODE_M_UP mode of
// Paul Falstad's EMWave2 applet.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMMGU()
{
        Identifier = "DEFAULT_TOOL_EMMGU";
        Name = "EMMGU";
        Colour = 0xC060E0_rgb;
        Description = ByteString("磁铁(向上):画出指向上方的永磁体,相当于applet的M UP模式").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMMGPaint::perform(sim, cpart, x, y, 1); // up
}
