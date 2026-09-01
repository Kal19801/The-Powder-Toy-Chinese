#include "EMAJTools.h"

// EMAJP: adjust the permeability of ferromagnets under the brush, port of the
// MODE_ADJ_PERM mode of Paul Falstad's EMWave2 applet (perm = vali/2, where
// vali is the adjust slider 3..100, so perm 1.5..50). The applet only affects
// cells of type FERROMAGNET.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMAJP()
{
        Identifier = "DEFAULT_TOOL_EMAJP";
        Name = "EMAJP";
        Colour = 0x60C090_rgb;
        Description = ByteString("调整磁导率:按工具强度设置铁磁材料(铁等)的磁导率1.5~50,仅对铁磁体生效").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        return EMAJAdjust::perform(sim, EMADJM_PERM, x, y, strength) ? 1 : 0;
}
