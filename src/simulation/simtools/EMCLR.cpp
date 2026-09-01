#include "simulation/ToolCommon.h"
#include "simulation/EMField.h"
#include "simulation/Simulation.h"

// EMCLR: clear EM field state, port of the MODE_CLEAR mode of Paul Falstad's
// EMWave2 applet. Wipes the EM cell back to vacuum (material overrides and
// wave state) under the brush. Deliberately does NOT delete the underlying
// particles: this tool only clears the field, matter is removed with the
// normal eraser, exactly like the air eraser never touches particles either.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMCLR()
{
        Identifier = "DEFAULT_TOOL_EMCLR";
        Name = "EMCLR";
        Colour = 0xFFB060_rgb;
        Description = ByteString("EM清除:只清除刷子下的电磁场状态(恢复真空),不删除粒子;删除物质请用原版橡皮").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        auto *emf = sim->GetEMField();
        if (!emf)
        {
                return 0;
        }
        int gi = emf->CellIndex(x, y);
        auto &cell = emf->cells[gi];
        bool changed = cell.ovMask || cell.jz != 0 || cell.jzext != 0 || cell.jmext != 0 ||
                       cell.az != 1e-10 || cell.dazdt != 1e-10 || cell.epos != 0;
        emf->VacuumCell(gi);
        return changed ? 1 : 0;
}
