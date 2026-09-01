#include "EMMGTools.h"

// EMMGL: paint permanent magnets pointing left, port of the MODE_M_LEFT mode
// of Paul Falstad's EMWave2 applet.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMMGL()
{
        Identifier = "DEFAULT_TOOL_EMMGL";
        Name = "EMMGL";
        Colour = 0xC060E0_rgb;
        Description = ByteString("磁铁(向左):画出指向左方的永磁体,相当于applet的M LEFT模式").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMMGPaint::perform(sim, cpart, x, y, 2); // left
}
