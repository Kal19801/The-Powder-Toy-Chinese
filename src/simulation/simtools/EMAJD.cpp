#include "EMAJTools.h"

// EMAJD: adjust the dielectric constant of dielectric media under the brush,
// port of the MODE_ADJ_MEDIUM mode of Paul Falstad's EMWave2 applet
// (medium = val * mediumMax). The applet only affects cells of type MEDIUM.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMAJD()
{
        Identifier = "DEFAULT_TOOL_EMAJD";
        Name = "EMAJD";
        Colour = 0xC0B060_rgb;
        Description = ByteString("调整介电常数:按工具强度设置电介质(玻璃等)的介电常数,仅对介质生效").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMAJAdjust::perform(sim, EMADJM_MEDIUM, x, y, strength) ? 1 : 0;
}
