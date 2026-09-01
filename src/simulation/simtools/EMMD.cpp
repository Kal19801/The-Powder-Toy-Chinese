#include "simulation/ToolCommon.h"
#include "simulation/EMField.h"
#include "simulation/Simulation.h"

// EMMD: magnetize ferromagnetic material under the brush (iron filings), a port
// of the "Adjust Mag Dir" mode of Paul Falstad's EMWave2 applet. The tool
// strength slider picks the magnetization direction as an angle from 0 to 2pi;
// the EMADJ tool picks the magnetization strength.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMMD()
{
	Identifier = "DEFAULT_TOOL_EMMD";
	Name = "EMMD";
	Colour = 0xCFA0FF_rgb;
	Description = ByteString("EM磁化:按工具强度设定的方向磁化铁等铁磁材料(形成永磁体)").FromUtf8();
	Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
	auto *emf = sim->GetEMField();
	if (!emf || !emf->enabled)
		return 0;
	int gi = emf->CellIndex(x, y);
	if (emf->ApplyMagDir(gi, strength))
	{
		emf->CalcBoundaries();
		return 1;
	}
	return 0;
}
