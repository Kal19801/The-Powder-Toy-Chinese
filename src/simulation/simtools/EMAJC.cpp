#include "EMAJTools.h"

// EMAJC: adjust the conductivity of conductors under the brush, port of the
// MODE_ADJ_CONDUCT mode of Paul Falstad's EMWave2 applet. The tool strength
// slider picks the value; the applet only affects cells of type CONDUCTOR.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMAJC()
{
        Identifier = "DEFAULT_TOOL_EMAJC";
        Name = "EMAJC";
        Colour = 0x6090FF_rgb;
        Description = ByteString("调整电导率:按工具强度设置导体(金属等)的电导率,仅对导体生效").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMAJAdjust::perform(sim, EMADJM_CONDUCT, x, y, strength) ? 1 : 0;
}
