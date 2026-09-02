#include "simulation/ToolCommon.h"
#include "simulation/EMField.h"
#include "simulation/Simulation.h"

// EMADJ: adjust the EM parameter of the material under the brush, a port of the
// adjust modes (conductivity / permeability / current / dielectric / magnet
// strength) of Paul Falstad's EMWave2 applet. The tool strength slider picks the
// value; the override only applies while the underlying material still matches.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMADJ()
{
	Identifier = "DEFAULT_TOOL_EMADJ";
	Name = "EMADJ";
	Colour = 0x80BFFF_rgb;
	Description = ByteString("EM参数调整:按工具强度调整电磁材料参数(电导/磁导/介质/电流/磁强)").FromUtf8();
	Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
	auto *emf = sim->GetEMField();
	if (!emf || !emf->enabled)
		return 0;
	int gi = emf->CellIndex(x, y);
	if (emf->ApplyAdjust(gi, strength))
	{
		// the wave equation needs to know about the material change right away
		// task 10: defer the O(gw*gh) CalcBoundaries to the next
		// EMField::Update(); per-dab calls were the EMADJ lag
		emf->NotifyCellChanged();
		return 1;
	}
	return 0;
}
