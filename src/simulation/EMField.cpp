#include "EMField.h"
#include "Simulation.h"
#include "ElementClasses.h"
#include "ElementDefs.h"
#include "common/tpt-rand.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numbers>

// NOTE: use std::numbers::pi instead of M_PI; M_PI is not defined by MSVC and
// not by MinGW under strict -std=c++20, which broke all Windows CI builds.

// Uncomment (or define through the build system) to trace the EM field internals
// while debugging; every block guarded by this macro is a debug block meant to be
// compiled out of release builds.
#ifndef EMFIELD_DEBUG
#define EMFIELD_DEBUG 0
#endif

#if EMFIELD_DEBUG
#define EMF_DBG(...) std::fprintf(stderr, __VA_ARGS__)
#else
#define EMF_DBG(...) do {} while (0)
#endif

namespace
{
        // Material properties of the EM field material elements and real zone
        // materials, ported from the painting modes of the applet:
        //   MODE_PERF_CONDUCTOR: conductivity 1
        //   MODE_GOOD_CONDUCTOR: conductivity .9
        //   MODE_FAIR_CONDUCTOR: conductivity .5
        //   MODE_FERROMAG: addPerm(5)  -> perm 5, conductivity .5
        //   MODE_DIAMAG:   addPerm(.5) -> perm .5, conductivity .5
        //   MODE_MEDIUM:   medium = mediumMax
        //   MODE_J_POS / MODE_J_NEG: jz = +-1
        // NOTE: original TPT elements are deliberately NOT mapped onto the EM
        // field anymore (their modifications were reverted); only the dedicated
        // elements of the EM and real zone categories participate.
        // SCND / SCPW are temperature dependent and handled separately.
        struct ElemMaterial
        {
                int type;
                float conductivity;
                float perm;      // relative permeability
                int medium;      // dielectric constant contribution
                bool source;     // injects current into the field
                float jz;        // current injected when source
        };

        const ElemMaterial elemMaterials[] = {
                // EM tab: the applet painting modes as elements
                { PT_EMPC, 1.0f, 1.0f, 0, false, 0 },
                { PT_EMEC, 0.9f, 1.0f, 0, false, 0 },
                { PT_EMFC, 0.5f, 1.0f, 0, false, 0 },
                { PT_EMFM, 0.5f, 5.0f, 0, false, 0 }, // applet addPerm(5)
                { PT_EMDM, 0.5f, 0.5f, 0, false, 0 }, // applet addPerm(.5)
                { PT_EMDE, 0.0f, 1.0f, EM_MEDIUM_MAX, false, 0 },
                // EM tab: current sources
                { PT_EMJP, 0.0f, 1.0f, 0, true,  1.0f },
                { PT_EMJN, 0.0f, 1.0f, 0, true, -1.0f },
                // real zone materials
                { PT_FE,   0.50f, 5.0f, 0, false, 0 },
                { PT_TI,   0.30f, 1.0f, 0, false, 0 },
                { PT_CU,   0.95f, 1.0f, 0, false, 0 },
                { PT_AG,   1.00f, 1.0f, 0, false, 0 },
                { PT_PGRF, 0.00f, 0.1f, 0, false, 0 },
                { PT_FEPW, 0.40f, 3.0f, 0, false, 0 },
                { PT_TIPW, 0.25f, 1.0f, 0, false, 0 },
                { PT_CUPW, 0.85f, 1.0f, 0, false, 0 },
                { PT_AGPW, 0.95f, 1.0f, 0, false, 0 },
                { PT_PGPW, 0.00f, 0.2f, 0, false, 0 },
        };

        // is this particle a real zone conductor that heats up from induced
        // currents? superconductors below their critical temperature are excluded
        // (zero resistance, no Joule losses)
        bool IsRealHeatable(int type, float temp)
        {
                switch (type)
                {
                case PT_EMPC: case PT_EMEC: case PT_EMFC:
                case PT_EMFM: case PT_EMDM:
                case PT_FE: case PT_TI: case PT_CU: case PT_AG:
                case PT_FEPW: case PT_TIPW: case PT_CUPW: case PT_AGPW:
                        return true;
                case PT_SCND: case PT_SCPW:
                        return temp >= EM_SC_TC; // quenched: resists again
                }
                return false;
        }

        // ferromagnetic / diamagnetic powders feel magnetic pressure;
        // superconducting powder only while it is still below Tc
        bool IsMagPowder(int type, float temp)
        {
                switch (type)
                {
                case PT_FEPW: case PT_PGPW:
                        return true;
                case PT_SCPW:
                        return temp < EM_SC_TC;
                }
                return false;
        }

        inline float ClampPerm(float perm)
        {
                // CFL stability bound of the leapfrog wave update, see
                // SimulationConfig.h; larger contrasts diverge (this is exactly
                // the old "ferromagnet self-excitation" divergence)
                return std::clamp(perm, EM_PERM_MIN, EM_PERM_MAX);
        }
}

EMField::EMField(Simulation & sim_) :
        sim(sim_)
{
        SetCellSize(EM_CELL_SIZE_DEFAULT);
        SetRegionScale(EM_REGION_SCALE_DEFAULT);
}

void EMField::SetCellSize(int newCellSize)
{
        bool valid = false;
        for (int i = 0; i < EM_CELL_SIZE_COUNT; ++i)
        {
                if (EM_CELL_SIZES[i] == newCellSize)
                {
                        valid = true;
                        break;
                }
        }
        if (!valid)
        {
                newCellSize = EM_CELL_SIZE_DEFAULT;
        }
        cellSize = newCellSize;
        ApplyGridGeometry();
}

void EMField::SetRegionScale(int newScale)
{
        bool valid = false;
        for (int i = 0; i < EM_REGION_SCALE_COUNT; ++i)
        {
                if (EM_REGION_SCALES[i] == newScale)
                {
                        valid = true;
                        break;
                }
        }
        if (!valid)
        {
                newScale = EM_REGION_SCALE_DEFAULT;
        }
        if (newScale == regionScale && cells.size() > 0)
        {
                return;
        }
        regionScale = newScale;
        ApplyGridGeometry();
}

void EMField::SetBoundaryMode(int newMode)
{
        if (newMode < 0 || newMode >= EMBND_COUNT)
        {
                newMode = EMBND_DEFAULT;
        }
        if (newMode == boundaryMode)
        {
                return;
        }
        boundaryMode = newMode;
        ApplyGridGeometry();
}

// Derive the full grid geometry from cellSize + regionScale + boundaryMode and
// reallocate the field. The simulation domain in pixels is regionScale * (XRES x
// YRES); the visible canvas is always XRES x YRES centred on the domain.
// - CLOSED: no padding (hard reflecting wall at the outermost simulated cell)
// - ABSORB / OPEN: pad by an invisible absorber band whose width is scaled so the
//   wave is attenuated below visibility before it can reflect off the hard outer
//   edge. ABSORB and OPEN now share the same geometry (the absorber sits OUTSIDE
//   the simulation domain in both), they differ only in ramp strength.
// - PERIODIC: pad by a one cell ghost ring which is refreshed from the opposite
//   edge every sub-step (true wrap-around without special-casing the wave update).
// Reallocating resets the wave state, the same trade-off as changing cellSize.
void EMField::ApplyGridGeometry()
{
        // total simulated cells (the domain); region > 1 means this extends past
        // the visible canvas in every direction
        visW = std::max(8, (XRES * regionScale) / cellSize);
        visH = std::max(8, (YRES * regionScale) / cellSize);
        padL = 0;
        padT = 0;
        margin = 0;
        switch (boundaryMode)
        {
        case EMBND_OPEN:
        case EMBND_ABSORB:
        {
                // absorber sits OUTSIDE the simulation domain (invisible padding).
                // Width scales with cellSize so the absorber is always ~the same
                // physical thickness in pixels; clamped to keep very fine grids
                // affordable. With regionScale > 1 the absorber is added on top of
                // the already-oversized domain, so the visible canvas stays bare.
                int pad = std::clamp(EM_MARGIN_AT_4 * EM_CELL_SIZE_DEFAULT / cellSize, 4, EM_OPEN_PAD_MAX);
                padL = pad;
                padT = pad;
                margin = pad;
                break;
        }
        case EMBND_PERIODIC:
                padL = 1;
                padT = 1;
                break;
        default: // EMBND_CLOSED: hard reflecting wall at the domain edge
                break;
        }
        gw = visW + 2 * padL;
        gh = visH + 2 * padT;
        cells.assign(gw * gh, Cell{});
        // visible canvas is always centred on the simulation domain
        int visCanvasW = XRES / cellSize; // cells of the visible canvas
        int visCanvasH = YRES / cellSize;
        renderOffX = padL + (visW - visCanvasW) / 2;
        renderOffY = padT + (visH - visCanvasH) / 2;
        if (renderOffX < padL) renderOffX = padL;
        if (renderOffY < padT) renderOffY = padT;
        forceBarValue = frequency;
        forceTimeZero = 0;
        t = 0;
        filterCount = 0;
        SetupSources();
        SetDamping();
        EMF_DBG("EMField: grid %dx%d cells (cell %dpx, region %dx, boundary %d, pad %d, vis %dx%d, margin %d, renderOff %d,%d)\n",
                gw, gh, cellSize, regionScale, boundaryMode, padL, visW, visH, margin, renderOffX, renderOffY);
}

void EMField::NotifyCellChanged()
{
        boundariesDirty = true;
}

void EMField::RunCalcBoundariesNow()
{
        CalcBoundaries();
        boundariesDirty = false;
}

void EMField::Clear()
{
        // The applet initialises the fields to 1e-10 instead of 0; do the same.
        for (auto &cell : cells)
        {
                cell.az = 1e-10;
                cell.dazdt = 1e-10;
                cell.epos = 0;
                cell.jz = 0;
                cell.jzext = 0;
                cell.jmext = 0;
                if (cell.resonant)
                {
                        cell.jz = 0;
                }
        }
        t = 0;
        forceTimeZero = 0;
        forceBarValue = frequency;
        filterCount = 0;
        filterCount %= 4;
}

void EMField::ClearAll()
{
        // applet doClearAll(): reset the wave state AND all material overrides; the
        // materials themselves are re-derived from the (now empty) particle list by
        // the next SyncMaterials(). Called whenever the whole simulation is reset
        // (new save, loaded save, clear), so no EM state leaks across saves.
        Clear();
        ClearOverrides();
        EMF_DBG("EMField: full reset (ClearAll)\n");
}

void EMField::ClearOverrides()
{
        for (auto &cell : cells)
        {
                cell.ovMask = 0;
                cell.ovConduct = 0;
                cell.ovPerm = 0;
                cell.ovJz = 0;
                cell.ovMedium = 0;
                cell.ovMag = 0;
                cell.ovDir = 0;
        }
        EMF_DBG("EMField: cleared all tool overrides\n");
}

void EMField::ClearCellOverrides(int gi)
{
        auto &cell = cells[gi];
        cell.ovMask = 0;
        cell.ovConduct = 0;
        cell.ovPerm = 0;
        cell.ovJz = 0;
        cell.ovMedium = 0;
        cell.ovMag = 0;
        cell.ovDir = 0;
}

void EMField::VacuumCell(int gi)
{
        // applet MODE_CLEAR / editFuncPoint preamble: wipe the cell back to vacuum;
        // used by the EMCLR tool, which only clears the field and never touches the
        // underlying particles (matter is removed with the normal eraser)
        ClearCellOverrides(gi);
        auto &cell = cells[gi];
        cell.jz = 0;
        cell.jzext = 0;
        cell.jmext = 0;
        cell.az = 1e-10;
        cell.dazdt = 1e-10;
        cell.epos = 0;
}

int EMField::CellIndex(int px, int py) const
{
        // map a pixel on the visible canvas to the corresponding cell of the
        // simulation domain. With regionScale > 1 the simulation extends beyond
        // the visible canvas, so the cell offset is renderOffX/Y (centred on
        // the domain); with regionScale = 1 renderOffX == padL so this reduces
        // to the previous behaviour.
        int cx = std::clamp(px / cellSize, 0, XRES / cellSize - 1) + renderOffX;
        int cy = std::clamp(py / cellSize, 0, YRES / cellSize - 1) + renderOffY;
        return cx + cy * gw;
}

// material type classification, faithful port of the applet's OscElement::getType();
// the order of the checks matches the applet exactly and must not be changed
int EMField::CellTypeOf(const Cell &oe)
{
        if (oe.perm < 1)
                return EMCT_DIAMAGNET;
        if (oe.perm > 1)
                return EMCT_FERROMAGNET;
        if (oe.mx != 0 || oe.my != 0)
                return EMCT_MAGNET;
        if (oe.medium > 0)
                return EMCT_MEDIUM;
        if (oe.conductivity > 0)
                return EMCT_CONDUCTOR;
        if (oe.jz != 0 || oe.jzext != 0)
                return EMCT_CURRENT;
        return EMCT_NONE;
}

// Derive the material properties of every EM cell from the TPT particles it covers.
void EMField::SyncMaterials()
{
        emwSources.clear();
        for (auto &cell : cells)
        {
                cell.perm = 1;
                cell.conductivity = 0;
                cell.mx = 0;
                cell.my = 0;
                cell.medium = 0;
                cell.jzext = 0;
                cell.jmext = 0;
                cell.resonant = false;
                cell.heatable = false;
                cell.magpowder = false;
                cell.magsolid = false;
        }

        auto &parts = sim.parts;
        for (int i = 0; i < parts.active; ++i)
        {
                auto &p = parts.data[i];
                if (!p.type)
                {
                        continue;
                }
                if (p.type == PT_EMW)
                {
                        // EMW: user-placeable EM source. .tmp overrides the global frequency,
                        // .ctype selects the waveform (0 = sine, 1 = packet).
                        float freq = frequency;
                        if (p.tmp >= 1 && p.tmp <= 40)
                        {
                                freq = float(p.tmp);
                        }
                        int waveform = (p.ctype == 1) ? EM_SWF_PACKET : EM_SWF_SIN;
                        emwSources.push_back({ CellIndex(int(p.x), int(p.y)), freq, waveform });
                        continue;
                }
                if (p.type == PT_EMR)
                {
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        cell.resonant = true;
                        cell.heatable = true; // resonant media dissipate absorbed energy
                        continue;
                }
                // --- vanilla SPRK no longer excites the EM field -------------------
                // Task 1: the vanilla -> EM current interop is disabled by design.
                // A vanilla spark on top of an EM-zone material used to inject
                // EM_SPRK_CURRENT into the cell's jzext, which coupled vanilla
                // circuits into the EM simulation. That coupling was unphysical
                // (a vanilla SPRK is a logical pulse, not a real current density)
                // and caused feedback loops with the EM-driven InteropParticles
                // path. The reverse direction (EM -> vanilla) is still active:
                // real charges touching a vanilla conductor spark it via
                // InteropParticles().
                //
                // If you want a vanilla circuit to drive the EM field, use the new
                // EMTX element (SPRK-triggered EM wave transmitter) or place an
                // EMJP/EMJN current source next to the conductor.
                // permanent magnets: the four dedicated direction elements, plus the
                // legacy EMMG element (.ctype 0..3) kept working for old saves
                switch (p.type)
                {
                case PT_EMMGD:
                        cells[CellIndex(int(p.x), int(p.y))].my =  1; // applet MODE_M_DOWN
                        continue;
                case PT_EMMGU:
                        cells[CellIndex(int(p.x), int(p.y))].my = -1; // applet MODE_M_UP
                        continue;
                case PT_EMMGL:
                        cells[CellIndex(int(p.x), int(p.y))].mx = -1; // applet MODE_M_LEFT
                        continue;
                case PT_EMMGR:
                        cells[CellIndex(int(p.x), int(p.y))].mx =  1; // applet MODE_M_RIGHT
                        continue;
                case PT_EMMG:
                {
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        switch (p.ctype)
                        {
                        case 1: cell.my = -1; break; // up
                        case 2: cell.mx = -1; break; // left
                        case 3: cell.mx =  1; break; // right
                        default: cell.my =  1; break; // down, like MODE_M_DOWN
                        }
                        continue;
                }
                default:
                        break;
                }
                // superconductors are temperature dependent: below the critical
                // temperature they are perfect conductors that expel the B field
                // (Meissner effect), above it they quench to a fair conductor
                if (p.type == PT_SCND || p.type == PT_SCPW)
                {
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        if (p.temp < EM_SC_TC)
                        {
                                cell.conductivity = std::max(cell.conductivity, 1.0f);
                                cell.perm = ClampPerm(0.01f); // Meissner: B expelled
                                cell.heatable = false;
                                cell.magpowder = (p.type == PT_SCPW); // levitates over magnets
                                cell.magsolid = (p.type == PT_SCND);  // solid Meissner body, mild push
                        }
                        else
                        {
                                cell.conductivity = std::max(cell.conductivity, 0.3f);
                                cell.perm = 1;
                                cell.heatable = true;
                        }
                        continue;
                }
                // --- Task 3: stationary magnetic monopole carries a static field ---
                // RMON deposits a constant magnetic current jmext in its cell so a
                // parked monopole radiates a 1/r B field of the same order as a
                // painted EMMG magnet (mx/my = +-1). The sign matches the pole
                // (ctype 0 = N = +1, ctype 1 = S = -1). Moving monopoles also
                // deposit jmext via DepositRealCharges(); both contributions sum.
                if (p.type == PT_RMON)
                {
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        float sign = (p.ctype == 1) ? -1.0f : 1.0f;
                        cell.jmext += double(sign * EM_MONO_STATIC);
                        continue;
                }
                // --- Task 5: EMAN antenna reads neighbour excitation and sparks ---
                // The antenna itself is non-material; it is processed in
                // InteropParticles() where it scans its 4-neighbourhood for any
                // EM-zone cell with |E|, |B|, |j| above a threshold and triggers
                // a vanilla SPRK on adjacent vanilla conductors. Nothing to do
                // here for the field state itself.
                if (p.type == PT_EMAN)
                {
                        continue;
                }
                // --- Task 6: EMTX transmitter injects a wave on SPRK ----------
                // The element is a passive field-material cell while not sparked;
                // while sparked (life > 0) its update() pushes a per-cell jzext
                // pulse into the field. The cell state is updated by the element's
                // own update function, not here.
                if (p.type == PT_EMTX)
                {
                        // The EMTX particle itself participates as a passive cell so
                        // neighbour EM cells see a slight impedance load; the
                        // actual wave injection happens in the element's update().
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        cell.perm = ClampPerm(1.0f); // vacuum-like, no disruption
                        continue;
                }
                // --- Task 10: EMJC current injector block --------------------
                // A powered (life > 0) EMJC injects a continuous current into its
                // cell, just like EMJP/EMJN but user-controlled via .tmp (signed
                // current -1..+1, sign sets direction) and .life (active gate).
                // Unpowered EMJC is inert. The element has PROP_CONDUCTS so SPRK
                // can drive it, mirroring the way BTRY drives a vanilla circuit.
                if (p.type == PT_EMJC)
                {
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        if (p.life > 0)
                        {
                                float j = std::clamp(p.tmp / 100.0f, -EM_JZEXT_MAX, EM_JZEXT_MAX);
                                if (j == 0.0f)
                                {
                                        j = EM_JZEXT_MAX; // default to full +current if .tmp unset
                                }
                                cell.jzext = double(j);
                        }
                        continue;
                }
                for (const auto &mat : elemMaterials)
                {
                        if (mat.type != p.type)
                        {
                                continue;
                        }
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        if (mat.source)
                        {
                                // applet semantics: painting a current source sets the
                                // cell current to +-1; stacked sources never add up,
                                // otherwise a full row would pump unbounded current
                                cell.jzext = mat.jz;
                        }
                        else
                        {
                                if (mat.conductivity > cell.conductivity)
                                {
                                        cell.conductivity = mat.conductivity;
                                }
                                if (IsRealHeatable(p.type, p.temp))
                                {
                                        cell.heatable = true;
                                }
                                if (IsMagPowder(p.type, p.temp))
                                {
                                        cell.magpowder = true;
                                }
                                // solid magnetic materials (FE/PGRF/EMFM/EMDM) get the weak
                                // version of the same magnetic pressure force
                                if (mat.perm != 1 && !IsMagPowder(p.type, p.temp) &&
                                        (p.type == PT_FE || p.type == PT_PGRF || p.type == PT_EMFM || p.type == PT_EMDM))
                                {
                                        cell.magsolid = true;
                                }
                        }
                        if (mat.perm != 1)
                        {
                                cell.perm = ClampPerm(mat.perm);
                        }
                        if (mat.medium > cell.medium)
                        {
                                cell.medium = mat.medium;
                        }
                        break;
                }
        }

        // Apply tool overrides, but only where the underlying material still matches,
        // exactly like the adjust modes of the applet which only affect matching cells.
        // Overrides on cells that lost all of their underlying material are wiped, so
        // erasing particles never leaves invisible EM leftovers behind.
        for (auto &cell : cells)
        {
                if (cell.ovMask & EM_OV_JZ && cell.jzext != 0)
                {
                        // the override carries the full signed current
                        cell.jzext = double(std::clamp(cell.ovJz, -EM_JZEXT_MAX, EM_JZEXT_MAX));
                }
                if (cell.perm > 1) // ferromagnet
                {
                        if (cell.ovMask & EM_OV_PERM)
                        {
                                cell.perm = ClampPerm(cell.ovPerm);
                        }
                }
                // magnetization overrides apply to ferromagnets AND to painted magnet
                // cells (applet TYPE_MAGNET), which carry perm == 1; the applet's
                // ADJ_MAG_STR sets the magnitude keeping the direction, while
                // ADJ_MAG_DIR sets the direction keeping the magnitude
                if (cell.perm > 1 || cell.mx != 0 || cell.my != 0)
                {
                        float mag = std::sqrt(cell.mx * cell.mx + cell.my * cell.my);
                        if (cell.ovMask & EM_OV_MAGSTR)
                        {
                                mag = std::clamp(cell.ovMag, 0.01f, 2.0f);
                        }
                        if (cell.ovMask & EM_OV_MAGDIR)
                        {
                                float angle = cell.ovDir * 2.0f * std::numbers::pi_v<float>;
                                cell.mx = mag * std::cos(angle);
                                cell.my = -mag * std::sin(angle);
                        }
                        else if ((cell.ovMask & EM_OV_MAGSTR) && (cell.mx != 0 || cell.my != 0))
                        {
                                // magnitude set by MAG_STR, direction preserved
                                float cur = std::sqrt(cell.mx * cell.mx + cell.my * cell.my);
                                if (cur > 0)
                                {
                                        cell.mx *= mag / cur;
                                        cell.my *= mag / cur;
                                }
                        }
                }
                else if (cell.ovMask & EM_OV_CONDUCT && cell.conductivity > 0)
                {
                        cell.conductivity = std::clamp(cell.ovConduct, 0.01f, 1.0f);
                }
                if (cell.medium > 0 && (cell.ovMask & EM_OV_MEDIUM))
                {
                        cell.medium = int(std::clamp(cell.ovMedium, 1.0f, float(EM_MEDIUM_MAX)));
                }
                // wipe overrides from cells whose material is gone; a painted
                // conductor whose particles were erased must not resurrect its
                // conductivity/adjustments when something else is drawn there later
                if (!(cell.perm != 1 || cell.conductivity > 0 || cell.medium > 0 ||
                      cell.mx != 0 || cell.my != 0 || cell.resonant || cell.jzext != 0))
                {
                        if (cell.ovMask)
                        {
                                ClearCellOverrides(int(&cell - cells.data()));
                        }
                }
                // induced current only exists inside conductors and resonant media;
                // when the conductor disappears (particles erased) its jz MUST be
                // dropped, otherwise the empty cell keeps radiating forever - that
                // was the "erased metal leaves a blob in the background" bug
                if (cell.conductivity <= 0 && !cell.resonant && cell.jz != 0)
                {
                        cell.jz = 0;
                }
        }

        CalcBoundaries();
}

void EMField::CalcBoundaries()
{
        // Mark all cells where the permeability, medium, magnetization or resonance
        // differs from one of the neighbours; the wave equation needs the hard path there.
        for (int y = 1; y < gh - 1; y++)
        {
                for (int x = 1; x < gw - 1; x++)
                {
                        int gi = x + y * gw;
                        auto &oe = cells[gi];
                        auto &e1 = cells[gi - 1];
                        auto &e2 = cells[gi + 1];
                        auto &e3 = cells[gi - gw];
                        auto &e4 = cells[gi + gw];
                        oe.gray = (oe.conductivity > 0 || oe.medium != 0 ||
                                   oe.perm != 1 || oe.mx != 0 || oe.my != 0 ||
                                   oe.resonant);
                        if (e1.perm != oe.perm || e2.perm != oe.perm ||
                            e3.perm != oe.perm || e4.perm != oe.perm ||
                            e1.medium != oe.medium || e2.medium != oe.medium ||
                            e3.medium != oe.medium || e4.medium != oe.medium ||
                            e1.mx != oe.mx || e2.mx != oe.mx || e3.mx != oe.mx || e4.mx != oe.mx ||
                            e1.my != oe.my || e2.my != oe.my || e3.my != oe.my || e4.my != oe.my ||
                            oe.resonant)
                        {
                                oe.boundary = true;
                        }
                        else
                        {
                                oe.boundary = false;
                        }
                }
        }
}

void EMField::SetDamping()
{
        for (auto &cell : cells)
        {
                cell.damp = 1;
                if (cell.medium > 0)
                {
                        cell.damp = .99; // need this to avoid reflections in dielectrics
                }
        }
        if (boundaryMode == EMBND_ABSORB || boundaryMode == EMBND_OPEN)
        {
                // exponential absorbing ramp in the invisible padding band; the
                // ramp starts WEAK at the inner edge (adjacent to the simulation
                // domain) and grows STRONG toward the outer edge, so the wave
                // sees a smooth impedance transition instead of a hard step that
                // would reflect it. The two boundary modes share the geometry
                // (padding outside the visible canvas); they differ only in ramp
                // strength - OPEN is the more aggressive absorber, ABSORB is the
                // softer one for users who want to see the wave gradually fade.
                float ramp = (boundaryMode == EMBND_OPEN) ? EM_OPEN_RAMP : EM_ABSORB_RAMP;
                for (int i = 0; i < margin; i++)
                {
                        // i=0 is the outermost cell of the absorber (next to the hard
                        // wall), i=margin-1 is the innermost (adjacent to the domain).
                        // The damping factor multiplies dazdt each sub-step, so a
                        // value close to 1 = weak damping (inner) and a smaller value
                        // = strong damping (outer). The factor is multiplicative
                        // across the band, giving total attenuation = product of all
                        // damp factors = exp(-sum).
                        double da = std::exp(-(margin - i) * ramp);
                        for (int x = 0; x < gw; x++)
                        {
                                cells[x + i * gw].damp = da;
                                cells[x + (gh - 1 - i) * gw].damp = da;
                        }
                        for (int y = 0; y < gh; y++)
                        {
                                cells[i + y * gw].damp = da;
                                cells[(gw - 1 - i) + y * gw].damp = da;
                        }
                }
        }
        // EMBND_CLOSED and EMBND_PERIODIC: no damping ramp at all; closed reflects
        // at the hard wall (never-updated outer ring), periodic wraps through the
        // ghost ring
}

void EMField::RefreshGhostRing()
{
        // PERIODIC boundary: copy the opposite edge into the ghost ring so the
        // interior update loops (1..gw-2) see true wrap-around neighbours. Both
        // az and dazdt are copied so the wave crosses the seam seamlessly.
        for (int j = 0; j < gh; j++)
        {
                int row = j * gw;
                cells[row] = cells[row + gw - 2];
                cells[row + gw - 1] = cells[row + 1];
        }
        for (int i = 0; i < gw; i++)
        {
                cells[i] = cells[i + (gh - 2) * gw];
                cells[i + (gh - 1) * gw] = cells[i + gw];
        }
}

void EMField::SetupSources()
{
        sourceMult = 1;
        sourceFreqCount = 1;
        sourcePlane = (sourceMode >= EMSRC_1S1F_PLANE);
        sourceWaveform = EM_SWF_SIN;
        sourceCount = 1;
        switch (sourceMode)
        {
        case EMSRC_NONE:          sourceCount = 0; break;
        case EMSRC_1S2F:          sourceFreqCount = 2; break;
        case EMSRC_2S1F:          sourceCount = 2; break;
        case EMSRC_2S2F:          sourceCount = 2; sourceFreqCount = 2; break;
        case EMSRC_3S1F:          sourceCount = 3; break;
        case EMSRC_4S1F:          sourceCount = 4; break;
        case EMSRC_1S1F_PACKET:   sourceWaveform = EM_SWF_PACKET; break;
        case EMSRC_1S2F_PLANE:    sourceFreqCount = 2; break;
        case EMSRC_2S1F_PLANE:    sourceCount = 2; break;
        case EMSRC_2S2F_PLANE:    sourceCount = 2; sourceFreqCount = 2; break;
        case EMSRC_1S1F_PLANE_PACKET: sourceWaveform = EM_SWF_PACKET; break;
        default: break;
        }
        if (sourceFreqCount >= 2)
        {
                auxFunction = EM_AUX_FREQ;
        }
        else if (sourceCount == 2 || sourceCount == 4)
        {
                auxFunction = EM_AUX_PHASE;
        }
        else
        {
                auxFunction = EM_AUX_NONE;
        }
        if (sourcePlane)
        {
                // plane sources are drawn between pairs of corner points of the
                // VISIBLE area, inset by one cell from the innermost usable ring
                sourceCount *= 2;
                int inset = (boundaryMode == EMBND_ABSORB) ? margin : 1;
                int x2 = padL + visW - inset - 1;
                int y2 = padT + visH - inset - 1;
                sources[0] = { padL + inset, padT + inset };
                sources[1] = { x2, padT + inset };
                sources[2] = { padL + inset, y2 };
                sources[3] = { x2, y2 };
        }
        else
        {
                // point sources sit around the centre of the visible canvas
                int inset = (boundaryMode == EMBND_ABSORB) ? margin : 1;
                sources[0] = { padL + visW / 2, padT + inset + 1 };
                sources[1] = { padL + visW / 2, padT + visH - inset - 2 };
                sources[2] = { padL + inset + 1, padT + visH / 2 };
                sources[3] = { padL + visW - inset - 2, padT + visH / 2 };
        }
        forceBarValue = frequency;
        forceTimeZero = 0;
        EMF_DBG("EMField: source mode %d, count %d plane %d waveform %d aux %d\n",
                sourceMode, sourceCount, int(sourcePlane), sourceWaveform, auxFunction);
}

void EMField::DoSources(double tadd, bool clear)
{
        double v = 0;
        double v2 = 0;
        if (sourceCount != 0)
        {
                double w = forceBarValue * (t - forceTimeZero) * EM_FREQ_MULT;
                double w2 = w;
                switch (auxFunction)
                {
                case EM_AUX_FREQ:
                        w2 = aux * t * EM_FREQ_MULT;
                        break;
                case EM_AUX_PHASE:
                {
                        // make sure the phase can be set all the way from 0 to pi
                        float au = aux - 1;
                        if (au > 38)
                        {
                                au = 38;
                        }
                        w2 = w + au * (std::numbers::pi / 38);
                        break;
                }
                }
                switch (sourceWaveform)
                {
                case EM_SWF_SIN:
                        v = std::sin(w);
                        if (sourceCount >= 2)
                        {
                                v2 = std::sin(w2);
                        }
                        else if (sourceFreqCount == 2)
                        {
                                v = (v + std::sin(w2)) * .5;
                        }
                        break;
                case EM_SWF_PACKET:
                {
                        double wp = std::fmod(w, std::numbers::pi * 2);
                        double adjw = wp / (EM_FREQ_MULT * forceBarValue);
                        adjw -= 10;
                        v = std::exp(-.01 * adjw * adjw) * std::sin(adjw * .2);
                        if (adjw < 0)
                        {
                                filterCount %= 4; // like the applet's doFilter()
                        }
                        break;
                }
                }
                if (clear)
                {
                        v = v2 = 0;
                }
                sources[0].v = sources[2].v = 2 * v * sourceMult;
                sources[1].v = sources[3].v = 2 * v2 * sourceMult;
                if (sourcePlane)
                {
                        for (int j = 0; j < sourceCount / 2; j++)
                        {
                                Source &src1 = sources[j * 2];
                                Source &src2 = sources[j * 2 + 1];
                                Source &src3 = sources[j];
                                // draw a line of current between the two corner points
                                int x1 = src1.x, y1 = src1.y;
                                int x2 = src2.x, y2 = src2.y;
                                if (y1 == y2)
                                {
                                        if (x1 == margin) x1 = 0;
                                        if (x2 == margin) x2 = 0;
                                        if (x1 == gw - margin - 1) x1 = gw - 1;
                                        if (x2 == gw - margin - 1) x2 = gw - 1;
                                }
                                if (x1 == x2)
                                {
                                        if (y1 == margin) y1 = 0;
                                        if (y2 == margin) y2 = 0;
                                        if (y1 == gh - margin - 1) y1 = gh - 1;
                                        if (y2 == gh - margin - 1) y2 = gh - 1;
                                }
                                double vv = src3.v * .1;
                                if (x1 == x2 && y1 == y2)
                                {
                                        cells[x1 + gw * y1].jzext = vv;
                                }
                                else if (std::abs(y2 - y1) > std::abs(x2 - x1))
                                {
                                        int sgn = (y2 > y1) ? 1 : -1;
                                        for (int y = y1; y != y2 + sgn; y += sgn)
                                        {
                                                int x = x1 + (x2 - x1) * (y - y1) / (y2 - y1);
                                                cells[x + gw * y].jzext = vv;
                                        }
                                }
                                else
                                {
                                        int sgn = (x2 > x1) ? 1 : -1;
                                        for (int x = x1; x != x2 + sgn; x += sgn)
                                        {
                                                int y = y1 + (y2 - y1) * (x - x1) / (x2 - x1);
                                                cells[x + gw * y].jzext = vv;
                                        }
                                }
                        }
                }
                else
                {
                        for (int i = 0; i < sourceCount; i++)
                        {
                                Source &src = sources[i];
                                cells[src.x + gw * src.y].jzext = src.v;
                        }
                }
        }

        // user-placed EMW sources
        for (const auto &es : emwSources)
        {
                double vew = 0;
                double w = es.freq * t * EM_FREQ_MULT;
                if (es.waveform == EM_SWF_PACKET)
                {
                        double wp = std::fmod(w, std::numbers::pi * 2);
                        double adjw = wp / (EM_FREQ_MULT * std::max(es.freq, 1.0f));
                        adjw -= 10;
                        vew = std::exp(-.01 * adjw * adjw) * std::sin(adjw * .2);
                }
                else
                {
                        vew = std::sin(w);
                }
                if (clear)
                {
                        vew = 0;
                }
                double amp = 2;
                auto &cell = cells[es.gi];
                if (cell.ovMask & EM_OV_JZ)
                {
                        amp = 2 * std::clamp(std::abs(cell.ovJz), 0.01f, 1.0f);
                }
                cell.jzext = vew * amp;
        }
}

void EMField::FilterGrid()
{
        // filter out high-frequency noise, but only right after the simulation has
        // been (re)started and only every 4th frame, exactly like the applet
        filterCount++;
        if ((filterCount & 3) != 0)
        {
                return;
        }
        if (filterCount > 200)
        {
                return;
        }
        // filter less aggressively if there is a source on the screen, to avoid
        // damping the waves
        double mult1 = (forceBarValue > 7 && sourceCount > 0 && sourceWaveform == EM_SWF_SIN) ? 40 : 8;
        double mult2 = 4 + mult1;
        for (int y = 1; y < gh - 1; y++)
        {
                for (int x = 1; x < gw - 1; x++)
                {
                        int gi = x + y * gw;
                        auto &oe = cells[gi];
                        if (oe.jz != 0 || oe.jzext != 0 || oe.jmext != 0 || oe.conductivity > 0)
                        {
                                continue;
                        }
                        if (oe.perm != cells[gi - 1].perm ||
                            oe.perm != cells[gi + 1].perm ||
                            oe.perm != cells[gi - gw].perm ||
                            oe.perm != cells[gi + gw].perm)
                        {
                                continue;
                        }
                        double jzm = cells[gi - 1].my - cells[gi + 1].my +
                                     cells[gi + gw].mx - cells[gi - gw].mx;
                        if (jzm != 0)
                        {
                                continue;
                        }
                        oe.az = (oe.az * mult1 + cells[gi - 1].az + cells[gi + 1].az +
                                 cells[gi - gw].az + cells[gi + gw].az) / mult2;
                }
        }
}

// --- rewritten current system ------------------------------------------------

// Gather all real charged particles into flat lists; called once per frame
// before the wave sub-steps.
void EMField::CollectRealCharges()
{
        realCharges.clear();
        realChargeCount = 0;
        auto &parts = sim.parts;
        for (int i = 0; i < parts.active; ++i)
        {
                auto &p = parts.data[i];
                float q = 0, g = 0, mass = 1;
                switch (p.type)
                {
                case PT_RPRO: q =  1;  mass = EM_PROTON_MASS; break;
                case PT_RELC: q = -1;  mass = 1;              break;
                case PT_RMON: g = (p.ctype == 1) ? -1 : 1;    break;
                default:      continue;
                }
                realChargeCount++;
                if (realCharges.size() < size_t(EM_PAIRWISE_LIMIT))
                {
                        // pairwise Coulomb only tracks a bounded number of charges
                        realCharges.push_back({ i, q, g, mass, 0, p.x, p.y, p.vx, p.vy });
                }
        }
        for (auto &rc : realCharges)
        {
                rc.gi = CellIndex(int(rc.x), int(rc.y));
        }
}

// Deposit the current of every real charge into the cells along the path it
// will travel this frame; a moving point charge IS a current density
// J = q v delta(x), and distributing it over the traversed cells keeps the
// deposit inside the light cone of the field (nothing is skipped).
void EMField::DepositRealCharges()
{
        float vmax = MaxParticleSpeed();
        for (const auto &rc : realCharges)
        {
                float vx = std::clamp(rc.vx, -vmax, vmax);
                float vy = std::clamp(rc.vy, -vmax, vmax);
                float speed = std::sqrt(vx * vx + vy * vy);
                if (speed < 1e-3f)
                {
                        continue; // a parked charge carries no current
                }
                // predicted pixel path of this frame
                float steps = std::max(1.0f, speed); // one sample per pixel
                int npath = 0;
                int pathGis[512];
                float px = rc.x, py = rc.y;
                float dx = vx / steps, dy = vy / steps;
                for (float s = 0; s < steps && npath < 512; s += 1.0f)
                {
                        int gi = CellIndex(int(px + dx * s), int(py + dy * s));
                        if (npath == 0 || pathGis[npath - 1] != gi)
                        {
                                pathGis[npath++] = gi;
                        }
                }
                float mag = speed / cellSize * EM_CHARGE_CURRENT / float(npath);
                for (int k = 0; k < npath; k++)
                {
                        auto &cell = cells[pathGis[k]];
                        if (rc.q != 0)
                        {
                                cell.jzext += rc.q * mag;
                        }
                        if (rc.g != 0)
                        {
                                // magnetic current drives dazdt directly, the dual of
                                // how the electric current jz drives the wave equation
                                cell.jmext += rc.g * speed / cellSize * EM_MONO_CURRENT / float(npath);
                        }
                }
        }
        // clamp the injected currents; the field can only carry a bounded amount
        // of energy per step and an over-injected cell is the fastest runaway
        for (auto &cell : cells)
        {
                if (cell.jzext > EM_JZEXT_MAX)
                {
                        cell.jzext = EM_JZEXT_MAX;
#if EMFIELD_DEBUG
                        jzClampHits++;
#endif
                }
                else if (cell.jzext < -EM_JZEXT_MAX)
                {
                        cell.jzext = -EM_JZEXT_MAX;
#if EMFIELD_DEBUG
                        jzClampHits++;
#endif
                }
        }
}

void EMField::InteractParticles(int substep, int substeps)
{
        auto &parts = sim.parts;
        float vmax = MaxParticleSpeed();
        // forces are integrated per sub-step so that a sub-step of the particle
        // motion always stays inside a sub-step of the field propagation;
        // every real particle in the simulation gets the field forces, while the
        // pairwise Coulomb acts between the tracked charges only (bounded work)
        for (int i = 0; i < parts.active; ++i)
        {
                auto &p = parts.data[i];
                float q = 0, g = 0, mass = 1;
                switch (p.type)
                {
                case PT_RPRO: q =  1; mass = EM_PROTON_MASS; break;
                case PT_RELC: q = -1; mass = 1;              break;
                case PT_RMON: g = (p.ctype == 1) ? -1 : 1;   break;
                default:      continue;
                }
                int gi = CellIndex(int(p.x), int(p.y));
                int cx = gi % gw;
                int cy = gi / gw;
                if (cx < 1 || cx >= gw - 1 || cy < 1 || cy >= gh - 1)
                {
                        continue;
                }
                auto b2 = [this](int gg) {
                        double dx = cells[gg - gw].az - cells[gg + gw].az;
                        double dy = cells[gg + 1].az - cells[gg - 1].az;
                        return dx * dx + dy * dy;
                };
                double fx = 0, fy = 0;
                // gradient force on electric charges: they are pulled toward
                // regions of strong field (polarity independent, like the
                // field-line pressure the applet shows in its force view)
                if (q != 0)
                {
                        double e2l = b2(gi - 1), e2r = b2(gi + 1);
                        double e2u = b2(gi - gw), e2d = b2(gi + gw);
                        fx += (e2r - e2l) * 0.5 * EM_GRAD_FORCE;
                        fy += (e2d - e2u) * 0.5 * EM_GRAD_FORCE;
                        // --- Lorentz-like force on individual moving charges ----
                        // In pure 2D TM mode the in-plane B field would give an
                        // out-of-plane v x B (unmodelled in 2D). We instead apply
                        // a fictitious in-plane force that treats |B| as if it
                        // were a perpendicular Bz field. The real 3D Lorentz force
                        // F = qv x B with v = (vx,vy,0) and B = (0,0,Bz) gives
                        // F = (q*vy*Bz, -q*vx*Bz, 0). We use |B| of the in-plane
                        // field as a proxy for Bz so a charge moving through a
                        // strong B region curves, the way a real charge would in
                        // a real Bz. Small gain keeps it a perturbation on top of
                        // the gradient force.
                        double bx = (cells[gi + gw].az - cells[gi - gw].az) * 0.5;
                        double by = (cells[gi - 1].az - cells[gi + 1].az) * 0.5;
                        double bmag = std::sqrt(bx * bx + by * by);
                        fx += q * double(p.vy) * bmag * EM_LORENTZ_FORCE;
                        fy -= q * double(p.vx) * bmag * EM_LORENTZ_FORCE;
                }
                // magnetic monopoles feel the in-plane B field directly: B = curl az
                if (g != 0)
                {
                        double bx = (cells[gi + gw].az - cells[gi - gw].az) * 0.5;
                        double by = (cells[gi - 1].az - cells[gi + 1].az) * 0.5;
                        fx += g * bx * EM_MONOPOLE_FORCE;
                        fy += g * by * EM_MONOPOLE_FORCE;
                }
                p.vx += float(fx / mass) * EM_TADD_SUB * 4.0f / float(substeps);
                p.vy += float(fy / mass) * EM_TADD_SUB * 4.0f / float(substeps);
                // strict speed of light: no particle may outrun the field, that is
                // the core anti-divergence invariant of the whole current system
                float sp = std::sqrt(p.vx * p.vx + p.vy * p.vy);
                if (sp > vmax)
                {
                        p.vx *= vmax / sp;
                        p.vy *= vmax / sp;
                }
                // conduction drift: inside a real-zone conductor the charge is
                // carried ALONG the wire by the local current, the EMWave2 jz
                // mechanism made directional. This is what makes our current
                // actually conduct: a charge entering one end of a copper wire
                // slides along it instead of rattling in place, and reaches the
                // other end where it can meet an opposite carrier or hand the
                // current over to a vanilla circuit. Monopoles do not ride wires.
                if (q != 0 && cells[gi].conductivity > 0)
                {
                        auto &cell = cells[gi];
                        // local wire tangent from conductor occupancy of the 4-neighbourhood
                        int cl = cells[gi - 1].conductivity > 0;
                        int cr = cells[gi + 1].conductivity > 0;
                        int cu = cells[gi - gw].conductivity > 0;
                        int cd = cells[gi + gw].conductivity > 0;
                        float tx = 0, ty = 0;
                        if (cl && cr) tx = 1;
                        else if (cu && cd) ty = 1;
                        else if (cl) tx = -1;
                        else if (cr) tx = 1;
                        else if (cu) ty = -1;
                        else if (cd) ty = 1;
                        if (tx != 0 || ty != 0)
                        {
                                // drive 1 (preferred): the field-intensity gradient along the
                                // wire - carriers are dragged away from the strong field side,
                                // i.e. away from whatever feeds the wire. Three intensities are
                                // summed because each carries the information in a different
                                // regime: |E|^2 (vacuum / weak conductors), |B|^2 (B penetrates
                                // conductors where E is screened) and the wire's own induced
                                // |jz|^2 (strongest right where the wire is fed). All three are
                                // sign-less so the direction stays steady for AC, DC and ringing.
                                double eL = cells[gi - 1].dazdt, eR = cells[gi + 1].dazdt;
                                double eU = cells[gi - gw].dazdt, eD = cells[gi + gw].dazdt;
                                auto b2f = [this](int gg) {
                                        double dx = cells[gg - gw].az - cells[gg + gw].az;
                                        double dy = cells[gg + 1].az - cells[gg - 1].az;
                                        return dx * dx + dy * dy;
                                };
                                double jL = cells[gi - 1].jz + cells[gi - 1].jzext;
                                double jR = cells[gi + 1].jz + cells[gi + 1].jzext;
                                double jU = cells[gi - gw].jz + cells[gi - gw].jzext;
                                double jD = cells[gi + gw].jz + cells[gi + gw].jzext;
                                double along = ((eR * eR - eL * eL) * tx + (eD * eD - eU * eU) * ty) * .5
                                             + ((b2f(gi + 1) - b2f(gi - 1)) * tx + (b2f(gi + gw) - b2f(gi - gw)) * ty) * .5
                                             + ((jR * jR - jL * jL) * tx + (jD * jD - jU * jU) * ty) * .5;
                                double drift = 0;
                                if (along > EM_DRIFT_NOISE || along < -EM_DRIFT_NOISE)
                                {
                                        drift = -std::clamp(along * 40.0, -1.0, 1.0);
                                }
                                else
                                {
                                        // drive 2 (fallback): the local conduction current itself
                                        // (external source or induced) sets the flow direction when
                                        // no usable gradient exists along the wire
                                        double drive = cell.jz + cell.jzext;
                                        if (std::abs(drive) > 1e-3)
                                        {
                                                drift = std::clamp(drive, -1.0, 1.0);
                                        }
                                }
                                if (drift != 0)
                                {
                                        float vdrift = float(drift) * EM_DRIFT_SPEED * cell.conductivity;
                                        // protons steer more sluggishly than electrons
                                        float blend = (q > 0) ? 0.35f : 1.0f;
                                        float nvx = p.vx * (1 - blend) + tx * vdrift * blend;
                                        float nvy = p.vy * (1 - blend) + ty * vdrift * blend;
                                        p.vx = std::clamp(nvx, -vmax, vmax);
                                        p.vy = std::clamp(nvy, -vmax, vmax);
                                }
                        }
                }
        }

        // pairwise Coulomb between the tracked charges (Newton's third law
        // honoured by pushing both partners); skipped entirely when there are
        // too many real particles to stay O(n^2)-cheap
        if (realCharges.size() > 1)
        {
                for (size_t a = 0; a < realCharges.size(); a++)
                {
                        auto &ra = realCharges[a];
                        if (ra.q == 0)
                        {
                                continue;
                        }
                        for (size_t b = a + 1; b < realCharges.size(); b++)
                        {
                                auto &rb = realCharges[b];
                                if (rb.q == 0)
                                {
                                        continue;
                                }
                                float ddx = ra.x - rb.x;
                                float ddy = ra.y - rb.y;
                                float r2 = ddx * ddx + ddy * ddy;
                                if (r2 > 64.0f)
                                {
                                        continue;
                                }
                                float r = std::sqrt(r2) + 1e-3f;
                                float f = EM_COULOMB * ra.q * rb.q / (r2 + 0.5f) / r; // f/r * (ddx,ddy)
                                auto &pa = parts.data[ra.i];
                                auto &pb = parts.data[rb.i];
                                pa.vx += f * ddx / ra.mass * EM_TADD_SUB * 4.0f / float(substeps);
                                pa.vy += f * ddy / ra.mass * EM_TADD_SUB * 4.0f / float(substeps);
                                pb.vx -= f * ddx / rb.mass * EM_TADD_SUB * 4.0f / float(substeps);
                                pb.vy -= f * ddy / rb.mass * EM_TADD_SUB * 4.0f / float(substeps);
                        }
                }
        }

        if (substep + 1 == substeps)
        {
                // once per frame: magnetic pressure on ferromagnetic / diamagnetic
                // powders (iron filings pull toward magnets, pyrolytic graphite and
                // superconductors are pushed away; solid magnetic materials respond
                // much more weakly), the motor effect (j x B force on current-
                // carrying conductors) and Joule heating of the real conductors
                // carrying induced current
                for (int cy = 1; cy < gh - 1; cy++)
                {
                        for (int cx = 1; cx < gw - 1; cx++)
                        {
                                int gi = cx + cy * gw;
                                auto &cell = cells[gi];
                                // pixel position of this cell on the visible canvas:
                                // the visible canvas is centred on the simulation
                                // domain, so cell (cx,cy) maps to pixel
                                // ((cx - renderOffX) * cs, (cy - renderOffY) * cs).
                                int px0 = (cx - renderOffX) * cellSize;
                                int py0 = (cy - renderOffY) * cellSize;
                                if (px0 < 0 || py0 < 0 || px0 >= XRES || py0 >= YRES)
                                {
                                        continue; // cell outside the visible canvas
                                }
                                if (cell.magpowder || cell.magsolid)
                                {
                                        auto b2f = [this](int gg) {
                                                double dx = cells[gg - gw].az - cells[gg + gw].az;
                                                double dy = cells[gg + 1].az - cells[gg - 1].az;
                                                return dx * dx + dy * dy;
                                        };
                                        double e2 = b2f(gi);
                                        double gx = (b2f(gi + 1) - b2f(gi - 1)) * 0.5;
                                        double gy = (b2f(gi + gw) - b2f(gi - gw)) * 0.5;
                                        // ferromagnets (perm > 1) move up the gradient of the
                                        // field energy, diamagnets (perm < 1) down it; the normaliser
                                        // keeps the force bounded while decaying like 1/L far away
                                        double sign = cell.perm > 1 ? 1.0 : -1.0;
                                        for (int sy = 0; sy < cellSize; sy++)
                                        {
                                                for (int sx = 0; sx < cellSize; sx++)
                                                {
                                                        int px = std::min(px0 + sx, XRES - 1);
                                                        int py = std::min(py0 + sy, YRES - 1);
                                                        int r = sim.pmap[py][px];
                                                        if (!r)
                                                        {
                                                                continue;
                                                        }
                                                        auto &p = parts.data[ID(r)];
                                                        bool powder = IsMagPowder(p.type, p.temp);
                                                        bool solid = p.type == PT_FE || p.type == PT_PGRF ||
                                                                     p.type == PT_EMFM || p.type == PT_EMDM ||
                                                                     (p.type == PT_SCND && p.temp < EM_SC_TC);
                                                        if (!powder && !solid)
                                                        {
                                                                continue;
                                                        }
                                                        float scale = powder ? EM_POWDER_FORCE : EM_SOLID_FORCE;
                                                        float fx = float(gx * sign * scale / (std::abs(e2) + EM_POWDER_NORM));
                                                        float fy = float(gy * sign * scale / (std::abs(e2) + EM_POWDER_NORM));
                                                        p.vx += std::clamp(fx, -0.8f, 0.8f);
                                                        p.vy += std::clamp(fy, -0.8f, 0.8f);
                                                }
                                        }
                                }
                                // --- motor effect: Lorentz force on a bulk current -------
                                // F = j x B. For our 2D TM field with j = (0,0,jz) and
                                // B = (Bx,By,0), F = (-jz*By, jz*Bx, 0), a real in-plane
                                // force on the current-carrying material. This is the
                                // force that makes a wire jump in a real magnetic field
                                // (the rail-gun / motor effect), and it was missing.
                                // Both the cell's induced jz and external jzext drive it;
                                // the result is applied to every particle of an EM-zone
                                // conductor in the cell, scaled by 1/substeps so the
                                // total per-frame impulse stays the same regardless of
                                // how many sub-steps the wave took.
                                double cellJz = cell.jz + cell.jzext;
                                if (cell.conductivity > 0 && std::abs(cellJz) > 0.005)
                                {
                                        double bx = (cells[gi + gw].az - cells[gi - gw].az) * 0.5;
                                        double by = (cells[gi - 1].az - cells[gi + 1].az) * 0.5;
                                        double mfx = -cellJz * by * EM_MOTOR_FORCE / double(substeps);
                                        double mfy =  cellJz * bx * EM_MOTOR_FORCE / double(substeps);
                                        float fxc = std::clamp(float(mfx), -0.5f, 0.5f);
                                        float fyc = std::clamp(float(mfy), -0.5f, 0.5f);
                                        if (std::abs(fxc) > 1e-4f || std::abs(fyc) > 1e-4f)
                                        {
                                                for (int sy = 0; sy < cellSize; sy++)
                                                {
                                                        for (int sx = 0; sx < cellSize; sx++)
                                                        {
                                                                int px = std::min(px0 + sx, XRES - 1);
                                                                int py = std::min(py0 + sy, YRES - 1);
                                                                int r = sim.pmap[py][px];
                                                                if (!r) continue;
                                                                auto &p = parts.data[ID(r)];
                                                                // apply to EM-zone conductors and real-zone
                                                                // conductors alike - the force is on the bulk
                                                                // material carrying the current, not on the
                                                                // individual charge
                                                                switch (p.type)
                                                                {
                                                                case PT_EMPC: case PT_EMEC: case PT_EMFC:
                                                                case PT_EMFM: case PT_EMDM:
                                                                case PT_FE: case PT_TI: case PT_CU: case PT_AG:
                                                                case PT_FEPW: case PT_TIPW: case PT_CUPW: case PT_AGPW:
                                                                case PT_SCND: case PT_SCPW:
                                                                        p.vx += fxc;
                                                                        p.vy += fyc;
                                                                        break;
                                                                default:
                                                                        break;
                                                                }
                                                        }
                                                }
                                        }
                                }
                                // Joule heating: induced current dissipates power in
                                // the real conductors; superconductors below Tc have
                                // no resistance and never heat up
                                if (cell.heatable && cell.conductivity > 0 && std::abs(cell.jz) > 0.01)
                                {
                                        float jzc = std::clamp(float(cell.jz), -2.0f, 2.0f);
                                        float heat = jzc * jzc * EM_JOULE_HEAT / float(substeps);
                                        if (heat >= 0.02f)
                                        {
                                                for (int sy = 0; sy < cellSize; sy++)
                                                {
                                                        for (int sx = 0; sx < cellSize; sx++)
                                                        {
                                                                int px = std::min(px0 + sx, XRES - 1);
                                                                int py = std::min(py0 + sy, YRES - 1);
                                                                int r = sim.pmap[py][px];
                                                                if (r && IsRealHeatable(parts.data[ID(r)].type, parts.data[ID(r)].temp))
                                                                {
                                                                        auto &p = parts.data[ID(r)];
                                                                        if (p.temp < MAX_TEMP)
                                                                        {
                                                                                p.temp = std::min(p.temp + heat, MAX_TEMP);
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}

// Once per frame, after the wave sub-steps: carrier neutralisation and the
// EM -> vanilla direction of the current interop.
//
// a) an opposite carrier pair meeting in the same spot annihilates (charge is
//    conserved: +1 and -1 vanish together);
// b) a real charge touching a vanilla conductor sparks it, handing the current
//    over to the vanilla conduction system - from there the spark propagates
//    and interacts exactly like any vanilla spark, so our current can do
//    everything a vanilla current can do. Particle indices are stable here
//    (kill_part uses a free list), so batch operations are safe.
void EMField::InteropParticles()
{
        auto &sd = SimulationData::CRef();
        auto &elements = sd.elements;
        auto &parts = sim.parts;

        // a) neutralisation over the tracked charges (bounded work); untracked
        //    charges simply never neutralise, like beyond the pairwise limit
        for (size_t a = 0; a + 1 < realCharges.size(); a++)
        {
                auto &ra = realCharges[a];
                if (parts.data[ra.i].type == PT_NONE)
                {
                        continue;
                }
                for (size_t b = a + 1; b < realCharges.size(); b++)
                {
                        auto &rb = realCharges[b];
                        if (parts.data[rb.i].type == PT_NONE)
                        {
                                continue;
                        }
                        bool opposite = (ra.q * rb.q < 0) || (ra.g * rb.g < 0);
                        if (!opposite)
                        {
                                continue;
                        }
                        float ddx = ra.x - rb.x;
                        float ddy = ra.y - rb.y;
                        if (ddx * ddx + ddy * ddy > 2.0f)
                        {
                                continue;
                        }
                        sim.kill_part(ra.i);
                        sim.kill_part(rb.i);
                        break;
                }
        }

        // b) EM -> vanilla: spark a vanilla conductor on contact
        for (int i = 0; i < parts.active; ++i)
        {
                auto &p = parts.data[i];
                if (p.type != PT_RPRO && p.type != PT_RELC)
                {
                        continue; // electric carriers only
                }
                int px = int(p.x + 0.5f);
                int py = int(p.y + 0.5f);
                if (px < 0 || py < 0 || px >= XRES || py >= YRES)
                {
                        continue;
                }
                unsigned r = sim.pmap[py][px];
                if (!r)
                {
                        continue;
                }
                int ri = ID(r);
                int tt = parts.data[ri].type;
                if (tt == PT_NONE || tt == PT_SPRK || !(elements[tt].Properties & PROP_CONDUCTS))
                {
                        continue;
                }
                // vanilla spark conversion, exactly like the photoelectric effect
                parts.data[ri].ctype = tt;
                sim.part_change_type(ri, px, py, PT_SPRK);
                parts.data[ri].life = 4;
        }

        // c) EMAN antenna: a coupler that bridges the EM field to vanilla SPRK.
        // For each EMAN particle, look at the EM cells in its 4-neighbourhood
        // (and its own cell). If any of them has |E|, |B| or |j| above an
        // "excitation" threshold (i.e. the EM zone is actively being driven),
        // the antenna fires a vanilla SPRK on every adjacent vanilla conductor
        // (PROP_CONDUCTS, not already sparked). The threshold scales with the
        // cellSize so a finer grid (smaller per-cell field magnitudes) still
        // triggers correctly.
        //
        // The antenna is the EM -> vanilla direction of "this EM zone is doing
        // something, push that signal out to the vanilla circuit". Pair it with
        // an EMTX transmitter (vanilla -> EM) for bidirectional coupling.
        double exciteThreshold = 0.02 * double(cellSize);
        for (int i = 0; i < parts.active; ++i)
        {
                auto &p = parts.data[i];
                if (p.type != PT_EMAN)
                {
                        continue;
                }
                int px = int(p.x + 0.5f);
                int py = int(p.y + 0.5f);
                if (px < 1 || py < 1 || px >= XRES - 1 || py >= YRES - 1)
                {
                        continue;
                }
                // gather excitation from the 3x3 neighbourhood of EM cells
                int gi = CellIndex(px, py);
                int cx = gi % gw;
                int cy = gi / gw;
                if (cx < 1 || cy < 1 || cx >= gw - 1 || cy >= gh - 1)
                {
                        continue;
                }
                double maxExcite = 0;
                for (int dy = -1; dy <= 1; dy++)
                {
                        for (int dx = -1; dx <= 1; dx++)
                        {
                                int ngi = gi + dx + dy * gw;
                                auto &nc = cells[ngi];
                                double e = std::abs(nc.dazdt); // |E|
                                double b = std::sqrt(std::pow(nc.az - cells[ngi - 1].az, 2) +
                                                     std::pow(nc.az - cells[ngi - gw].az, 2));
                                double j = std::abs(nc.jz + nc.jzext);
                                maxExcite = std::max({maxExcite, e, b, j});
                        }
                }
                if (maxExcite < exciteThreshold)
                {
                        continue;
                }
                // spark every adjacent vanilla conductor (4-neighbourhood)
                static const int dxs[4] = {-1, 1, 0, 0};
                static const int dys[4] = {0, 0, -1, 1};
                for (int k = 0; k < 4; k++)
                {
                        int nx = px + dxs[k];
                        int ny = py + dys[k];
                        unsigned r = sim.pmap[ny][nx];
                        if (!r) continue;
                        int ri = ID(r);
                        int tt = parts.data[ri].type;
                        if (tt == PT_NONE || tt == PT_SPRK) continue;
                        if (!(elements[tt].Properties & PROP_CONDUCTS)) continue;
                        // don't re-spark an active conductor (life > 0 means SPRK is ongoing)
                        if (parts.data[ri].life > 0) continue;
                        // fire SPRK on the vanilla conductor
                        parts.data[ri].ctype = tt;
                        sim.part_change_type(ri, nx, ny, PT_SPRK);
                        parts.data[ri].life = 4;
                }
        }
}

void EMField::Update()
{
        if (!enabled)
        {
                return;
        }
        // Consume any deferred boundary recompute the EMADJ tool scheduled. Doing
        // it once here, instead of after every brush dab, is the fix for the
        // EMADJ "lag" - CalcBoundaries is O(gw*gh) and used to run after every
        // single Apply(), making a brush stroke O(gw*gh * brush_area) per frame.
        if (boundariesDirty)
        {
                CalcBoundaries();
                boundariesDirty = false;
        }
        SyncMaterials();
        CollectRealCharges();
        DepositRealCharges();

        // the wave is integrated in fixed sub-steps of the exact applet timestep;
        // this keeps the CFL bound (perm contrast <= 2/tadd^2 = 32) satisfied for
        // every speed setting, matches the applet dynamics much more closely than
        // one large step per frame, and conserves energy in the linear regime
        int substeps = EM_SUBSTEPS[std::clamp(speed, 0, 4)];
        double tadd = EM_TADD_SUB;
        double tadd2 = tadd * tadd;

        for (int substep = 0; substep < substeps; substep++)
        {
                // PERIODIC boundary: refresh the ghost ring from the opposite edge
                // before every sub-step so the interior loops see wrap-around
                if (boundaryMode == EMBND_PERIODIC)
                {
                        RefreshGhostRing();
                }

                DoSources(tadd, false);

                // --- first pass: update dazdt from the neighbours (ported from EMWave2) ---
                double forcecoef = 1;
                int curMedium = 0;
                for (int j = 1; j < gh - 1; j++)
                {
                        for (int i = 1; i < gw - 1; i++)
                        {
                                int gi = i + j * gw;
                                auto &oe = cells[gi];
                                auto &oew = cells[gi - 1];
                                auto &oee = cells[gi + 1];
                                auto &oen = cells[gi - gw];
                                auto &oes = cells[gi + gw];

                                if (oe.conductivity > 0)
                                {
                                        oe.jz = 0;
                                }

                                double a;
                                if (oe.boundary)
                                {
                                        // we may be on the boundary between two media, so we have
                                        // to do some extra work
                                        if (oe.resonant)
                                        {
                                                oe.jz = oe.jz * .999 - oe.dazdt * .001 - oe.epos * .02;
                                                oe.epos += oe.jz * .2;
                                        }
                                        if (curMedium != oe.medium)
                                        {
                                                curMedium = oe.medium;
                                                forcecoef = (1 - (EM_MEDIUM_MAX_INDEX / EM_MEDIUM_MAX) * curMedium);
                                                forcecoef *= forcecoef;
                                        }
                                        double az = oe.az;
                                        double previ = (oew.az - az) / oew.perm;
                                        double nexti = (oee.az - az) / oee.perm;
                                        double prevj = (oen.az - az) / oen.perm;
                                        double nextj = (oes.az - az) / oes.perm;
                                        double basis = (nexti + previ + nextj + prevj) * .25;
                                        // calculate effective current
                                        double jz = oew.my - oee.my + oes.mx - oen.mx + oe.jz + oe.jzext;
                                        a = oe.perm * basis + jz + oe.jmext;
                                }
                                else
                                {
                                        // easy way
                                        double basis = (oew.az + oee.az + oen.az + oes.az) * .25;
                                        a = oe.jz + oe.jzext + oe.jmext - (oe.az - basis);
                                }
                                oe.dazdt = oe.dazdt * oe.damp + a * forcecoef;
                        }
                }

                // --- second pass: integrate az, and compute induced currents in conductors ---
                for (int j = 1; j < gh - 1; j++)
                {
                        for (int i = 1; i < gw - 1; i++)
                        {
                                int gi = i + j * gw;
                                auto &oe = cells[gi];
                                if (oe.conductivity > 0)
                                {
                                        double a = -oe.dazdt * oe.conductivity;
                                        oe.jz = a;
                                        oe.dazdt += a;
                                }
                                oe.az += oe.dazdt * tadd2;
                        }
                }
                t += tadd;

                InteractParticles(substep, substeps);
                FilterGrid();

#if EMFIELD_DEBUG
                // safety net bookkeeping: the clamp must never engage in normal
                // operation; if it does, the CFL bound somewhere above is broken
                for (auto &cell : cells)
                {
                        if (std::abs(cell.az) > EM_FIELD_CLAMP || std::abs(cell.dazdt) > EM_FIELD_CLAMP)
                        {
                                cell.az = std::clamp(cell.az, -double(EM_FIELD_CLAMP), double(EM_FIELD_CLAMP));
                                cell.dazdt = std::clamp(cell.dazdt, -double(EM_FIELD_CLAMP), double(EM_FIELD_CLAMP));
                                fieldClampHits++;
                        }
                }
                double maxdazdt = 0;
                for (const auto &cell : cells)
                {
                        maxdazdt = std::max(maxdazdt, std::abs(cell.dazdt));
                }
                if (maxdazdt > EM_FIELD_CLAMP * 0.5)
                {
                        EMF_DBG("EMField: WARNING dazdt close to clamp: %g\n", maxdazdt);
                }
#endif
        }

        // --- Task 9: anti-divergence safety net ----------------------------------
        // The leapfrog scheme conserves energy in the linear regime but feedback
        // between non-linear pieces (real charge -> field -> real charge, Joule
        // heating -> temperature-dependent conductivity, etc.) can pump energy
        // in. Catch the divergence early by detecting cells whose |az| or |dazdt|
        // has grown well beyond what the simulation normally carries (1e4 is the
        // hard safety clamp; here we use a softer 1e2 threshold and gently bleed
        // the excess off across a few frames instead of clamping hard). This is
        // applied unconditionally (not under EMFIELD_DEBUG) because it is the
        // user-facing stability net.
        {
                const double softLimit = 100.0;
                const double bleedRate = 0.95; // multiplicative decay per frame
                for (auto &cell : cells)
                {
                        if (std::abs(cell.az) > softLimit)
                        {
                                cell.az *= bleedRate;
                        }
                        if (std::abs(cell.dazdt) > softLimit)
                        {
                                cell.dazdt *= bleedRate;
                        }
                }
        }

        // once per frame: the EM -> vanilla direction of the current interop and
        // carrier neutralisation
        InteropParticles();

#if EMFIELD_DEBUG
        double energy = 0;
        for (int y = 1; y < gh - 1; y++)
        {
                for (int x = 1; x < gw - 1; x++)
                {
                        int gi = x + y * gw;
                        const auto &oe = cells[gi];
                        energy += 0.5 * oe.dazdt * oe.dazdt + 0.5 * oe.az * oe.az;
                }
        }
        EMF_DBG("EMField: t=%.2f energy=%.6g jzclamp=%lld fieldclamp=%lld\n",
                t, energy, jzClampHits, fieldClampHits);
#endif
}

double EMField::GetMagX(int gi) const
{
        auto &oe = cells[gi];
        double mm = 1 - 1 / oe.perm;
        return (cells[gi - gw].az - cells[gi + gw].az) * mm + oe.mx;
}

double EMField::GetMagY(int gi) const
{
        auto &oe = cells[gi];
        double mm = 1 - 1 / oe.perm;
        return (cells[gi + 1].az - cells[gi - 1].az) * mm + oe.my;
}

bool EMField::ApplyAdjust(int gi, float strength)
{
        // port of the applet's adjust modes; each cell is adjusted according to its
        // most specific material type, in the same priority order as getType()
        auto &cell = cells[gi];
        bool applied = false;
        if (cell.perm > 1) // ferromagnet
        {
                cell.ovMask |= EM_OV_PERM;
                cell.ovPerm = std::clamp(strength, 0.01f, 1.0f) * 32.0f;
                applied = true;
        }
        else if (cell.medium > 0)
        {
                cell.ovMask |= EM_OV_MEDIUM;
                cell.ovMedium = std::clamp(strength, 0.01f, 1.0f) * EM_MEDIUM_MAX;
                applied = true;
        }
        else if (cell.conductivity > 0)
        {
                cell.ovMask |= EM_OV_CONDUCT;
                cell.ovConduct = std::clamp(strength, 0.01f, 1.0f);
                applied = true;
        }
        else if (cell.mx != 0 || cell.my != 0)
        {
                cell.ovMask |= EM_OV_MAGSTR;
                cell.ovMag = std::clamp(strength, 0.01f, 1.0f);
                applied = true;
        }
        else if (cell.jzext != 0)
        {
                cell.ovMask |= EM_OV_JZ;
                cell.ovJz = (cell.jzext < 0 ? -1.0f : 1.0f) * std::clamp(strength, 0.01f, 1.0f);
                applied = true;
        }
        return applied;
}

bool EMField::ApplyMagDir(int gi, float strength)
{
        auto &cell = cells[gi];
        if (cell.perm > 1 || cell.mx != 0 || cell.my != 0)
        {
                cell.ovMask |= EM_OV_MAGDIR;
                cell.ovDir = std::clamp(strength, 0.0f, 1.0f);
                return true;
        }
        return false;
}

// One specific applet MODE_ADJ_* mode, faithful port of EMWave2's doAdjust():
// the value only lands on cells whose getType() matches the mode, the strength
// slider (0..1) plays the role of the applet's adjustBar (1..100 -> 0.01..1).
bool EMField::ApplyAdjustMode(int mode, int gi, float strength)
{
        auto &cell = cells[gi];
        float val = std::clamp(strength, 0.01f, 1.0f);
        switch (mode)
        {
        case EMADJM_CONDUCT:
                // applet: if (oe.getType() == TYPE_CONDUCTOR) oe.conductivity = val;
                if (CellTypeOf(cell) != EMCT_CONDUCTOR)
                        return false;
                cell.ovMask |= EM_OV_CONDUCT;
                cell.ovConduct = val;
                return true;
        case EMADJM_PERM:
        {
                // applet: vali clamped to >= 3 (0.03), perm = vali/2 (so 1.5 .. 50);
                // here the upper end is the CFL stability bound 32
                int vali = int(val * 100.0f);
                if (vali < 3)
                        vali = 3;
                if (CellTypeOf(cell) != EMCT_FERROMAGNET)
                        return false;
                cell.ovMask |= EM_OV_PERM;
                cell.ovPerm = ClampPerm(vali / 2.0f);
                return true;
        }
        case EMADJM_J:
                // applet: if (getType() == TYPE_CURRENT) oe.jz = sign * val; our
                // external currents live in jzext, the override replaces the value
                // keeping the sign of the underlying source
                if (CellTypeOf(cell) != EMCT_CURRENT)
                        return false;
                cell.ovMask |= EM_OV_JZ;
                cell.ovJz = (cell.jzext < 0 ? -1.0f : 1.0f) * val;
                return true;
        case EMADJM_MEDIUM:
                if (CellTypeOf(cell) != EMCT_MEDIUM)
                        return false;
                cell.ovMask |= EM_OV_MEDIUM;
                cell.ovMedium = val * EM_MEDIUM_MAX;
                return true;
        case EMADJM_MAG_DIR:
                // applet: if (getType() == TYPE_MAGNET) rotate; our version also
                // accepts ferromagnets, like the pre-existing EMMD tool
                return ApplyMagDir(gi, val);
        case EMADJM_MAG_STR:
                // applet: if (getType() == TYPE_MAGNET) scale |m| to val
                if (CellTypeOf(cell) != EMCT_MAGNET && CellTypeOf(cell) != EMCT_FERROMAGNET)
                        return false;
                cell.ovMask |= EM_OV_MAGSTR;
                cell.ovMag = val;
                return true;
        }
        return false;
}

// effective (override aware) value of one cell property
float EMField::EffectiveProperty(int property, int gi) const
{
        auto &cell = cells[gi];
        switch (property)
        {
        case EMADJP_CONDUCT:
                return cell.conductivity;
        case EMADJP_PERM:
                return cell.perm;
        case EMADJP_J:
                return float(cell.jzext);
        case EMADJP_MEDIUM:
                return float(cell.medium);
        case EMADJP_MAG_DIR:
        {
                // current direction as a fraction of 2pi; stored override wins
                if (cell.ovMask & EM_OV_MAGDIR)
                {
                        return cell.ovDir;
                }
                if (cell.mx == 0 && cell.my == 0)
                {
                        return 0;
                }
                float angle = std::atan2(-cell.my, cell.mx) / (2.0f * std::numbers::pi_v<float>);
                if (angle < 0)
                {
                        angle += 1;
                }
                return angle;
        }
        case EMADJP_MAG_STR:
                return std::sqrt(cell.mx * cell.mx + cell.my * cell.my);
        }
        return 0;
}

// Unified adjust tool: applies one property with one of the three modes
// (set / add / subtract) to one cell; the value only lands on cells whose
// material matches the property, exactly like the applet adjust modes.
bool EMField::ApplyEMProperty(int property, int applyMode, int gi, float value)
{
        auto &cell = cells[gi];
        float eff = EffectiveProperty(property, gi);
        float target = value;
        if (applyMode == EMADJA_ADD)
        {
                target = eff + value;
        }
        else if (applyMode == EMADJA_SUB)
        {
                target = eff - value;
        }
        switch (property)
        {
        case EMADJP_CONDUCT:
                if (cell.conductivity <= 0)
                        return false;
                cell.ovMask |= EM_OV_CONDUCT;
                cell.ovConduct = std::clamp(target, 0.0f, 1.0f);
                return true;
        case EMADJP_PERM:
                if (cell.perm <= 1)
                        return false;
                cell.ovMask |= EM_OV_PERM;
                cell.ovPerm = ClampPerm(target);
                return true;
        case EMADJP_J:
                if (cell.jzext == 0)
                        return false;
                cell.ovMask |= EM_OV_JZ;
                cell.ovJz = std::clamp(target, -EM_JZEXT_MAX, EM_JZEXT_MAX);
                return true;
        case EMADJP_MEDIUM:
                if (cell.medium <= 0)
                        return false;
                cell.ovMask |= EM_OV_MEDIUM;
                cell.ovMedium = std::clamp(target, 1.0f, float(EM_MEDIUM_MAX));
                return true;
        case EMADJP_MAG_DIR:
        {
                if (!(cell.perm > 1 || cell.mx != 0 || cell.my != 0))
                        return false;
                cell.ovMask |= EM_OV_MAGDIR;
                float frac = target / 360.0f; // the tool passes degrees
                frac -= std::floor(frac);
                cell.ovDir = frac;
                return true;
        }
        case EMADJP_MAG_STR:
        {
                if (!(cell.perm > 1 || cell.mx != 0 || cell.my != 0))
                        return false;
                cell.ovMask |= EM_OV_MAGSTR;
                cell.ovMag = std::clamp(target, 0.0f, 2.0f);
                return true;
        }
        }
        return false;
}
