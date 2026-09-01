#include "EMMGTools.h"

// EMMGD: paint permanent magnets pointing down, port of the MODE_M_DOWN mode
// of Paul Falstad's EMWave2 applet.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMMGD()
{
        Identifier = "DEFAULT_TOOL_EMMGD";
        Name = "EMMGD";
        Colour = 0xC060E0_rgb;
        Description = ByteString("磁铁(向下):画出指向下方的永磁体,相当于applet的M DOWN模式").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMMGPaint::perform(sim, cpart, x, y, 0); // down
}
