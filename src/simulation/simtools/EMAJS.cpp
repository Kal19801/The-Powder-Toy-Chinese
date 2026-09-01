#include "EMAJTools.h"

// EMAJS: adjust the strength (magnetization magnitude) of magnets under the
// brush, port of the MODE_ADJ_MAG_STR mode of Paul Falstad's EMWave2 applet
// (scale mx/my so that |m| = val). The applet only affects cells of type
// MAGNET; ferromagnets are accepted too, like in the combined EMADJ tool.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMAJS()
{
        Identifier = "DEFAULT_TOOL_EMAJS";
        Name = "EMAJS";
        Colour = 0x9060C0_rgb;
        Description = ByteString("调整磁强度:按工具强度设置永磁体(EMMG/磁化过的铁)的磁化强度").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMAJAdjust::perform(sim, EMADJM_MAG_STR, x, y, strength) ? 1 : 0;
}
