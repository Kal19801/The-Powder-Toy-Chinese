#include "simulation/ToolCommon.h"
#include "simulation/EMField.h"
#include "simulation/Simulation.h"

// EMCLR: clear EM material, port of the MODE_CLEAR mode of Paul Falstad's
// EMWave2 applet. Wipes the EM cell back to vacuum (material overrides and
// wave state) and deletes the particles under the brush, so anything the EM
// tools painted can be removed again from inside the game.

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength);

void SimTool::Tool_EMCLR()
{
        Identifier = "DEFAULT_TOOL_EMCLR";
        Name = "EMCLR";
        Colour = 0xFFB060_rgb;
        Description = ByteString("EM清除:删除刷子下的粒子并把电磁场该处恢复为真空(清除EMADJ/EMMD的效果)").FromUtf8();
        Perform = &perform;
}

static int perform(SimTool *tool, Simulation * sim, Particle * cpart, int x, int y, int brushX, int brushY, float strength)
{
        bool changed = false;
        if (cpart && cpart->type)
        {
                sim->delete_part(x, y);
                changed = true;
        }
        auto *emf = sim->GetEMField();
        if (emf)
        {
                int gi = emf->CellIndex(x, y);
                auto &cell = emf->cells[gi];
                if (!changed && (cell.ovMask || cell.jz != 0 || cell.jzext != 0))
                {
                        changed = true;
                }
                emf->VacuumCell(gi);
                return changed ? 1 : 0;
        }
        return changed ? 1 : 0;
}
