#include "EMAJTools.h"

// EMAJJ: adjust the current of current sources under the brush, port of the
// MODE_ADJ_J mode of Paul Falstad's EMWave2 applet (jz = sign * val). The
// applet only affects cells of type CURRENT; here it rescales the external
// current injected by sparks, electrons and the EMJP/EMJN/EMW sources.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMAJJ()
{
        Identifier = "DEFAULT_TOOL_EMAJJ";
        Name = "EMAJJ";
        Colour = 0xFF9060_rgb;
        Description = ByteString("调整电流:按工具强度设置电流源(火花/电子/EMJP/EMJN/EMW)的电流大小").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMAJAdjust::perform(sim, EMADJM_J, x, y, strength) ? 1 : 0;
}
