#include "EMField.h"
#include "Simulation.h"
#include "ElementClasses.h"
#include "ElementDefs.h"
#include "common/tpt-rand.h"
#include <algorithm>
#include <cstring>
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

        inline float ClampPerm(float perm, float permMax)
        {
                // CFL stability bound of the leapfrog wave update, see
                // SimulationConfig.h; the bound is on the perm CONTRAST ratio
                // against the weakest neighbour, and it depends on the wave
                // timestep (which depends on the cell size). Larger contrasts
                // diverge (this is exactly the old "ferromagnet
                // self-excitation" divergence).
                return std::clamp(perm, EM_PERM_MIN, permMax);
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

void EMField::SetRegionScale(float newScale)
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
// - ABSORB / OPEN: pad by an invisible advective-outflow absorber band whose
//   width is fixed in PIXELS (EM_PAD_PX), so the band behaves identically at
//   every grid resolution; a wave crossing into the band keeps its impedance
//   and leaves the domain without reflecting (task 8)
// - PERIODIC: pad by a one cell ghost ring which is refreshed from the opposite
//   edge every sub-step (true wrap-around without special-casing the wave update).
// Reallocating resets the wave state, the same trade-off as changing cellSize.
void EMField::ApplyGridGeometry()
{
        // the wave timestep shrinks with the cell size so that the wave keeps
        // the same pixel speed and pixel wavelength at every resolution (task 7)
        taddEff = EmWaveTadd(cellSize);
        // total simulated cells (the domain); region > 1 means this extends past
        // the visible canvas in every direction, region < 1 means it covers only
        // the central part of the canvas
        visW = std::max(8, int(float(XRES) * regionScale) / cellSize);
        visH = std::max(8, int(float(YRES) * regionScale) / cellSize);
        padL = 0;
        padT = 0;
        margin = 0;
        switch (boundaryMode)
        {
        case EMBND_OPEN:
        {
                // absorber band OUTSIDE the simulation domain (invisible
                // padding), fixed thickness in pixels. With regionScale > 1 the
                // band is added on top of the already-oversized domain, so the
                // visible canvas stays bare. With regionScale < 1 the domain is
                // smaller than the canvas, so the band is capped to a quarter
                // of the domain to keep most of it usable (regionScale >= 1
                // keeps the exact original band width).
                int pad = std::clamp(EM_PAD_PX / cellSize, EM_PAD_MIN_CELLS, EM_PAD_MAX_CELLS);
                if (regionScale < 1.0f)
                {
                        pad = std::min(pad, std::min(visW, visH) / 4);
                        pad = std::max(pad, 2);
                }
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
        // guard against extreme regionScale+cellSize combinations exhausting memory
        while ((long long)(visW + 2 * padL) * (visH + 2 * padT) > EM_MAX_CELLS && regionScale > EM_REGION_SCALES[0])
        {
                // step down to the next smaller valid scale
                float next = EM_REGION_SCALES[0];
                for (int i = 0; i + 1 < EM_REGION_SCALE_COUNT; ++i)
                {
                        if (EM_REGION_SCALES[i] < regionScale && EM_REGION_SCALES[i] > next)
                        {
                                next = EM_REGION_SCALES[i];
                        }
                }
                regionScale = next;
                visW = std::max(8, int(float(XRES) * regionScale) / cellSize);
                visH = std::max(8, int(float(YRES) * regionScale) / cellSize);
                EMF_DBG("EMField: region scale reduced to %gx to fit the cell budget\n", double(regionScale));
        }
        gw = visW + 2 * padL;
        gh = visH + 2 * padT;
        cells.assign(gw * gh, Cell{});
        // CPML band state (OPEN only): face filter states start at zero, which
        // is consistent with the 1e-10 vacuum seed of az/dazdt (uniform field ->
        // zero face differences -> the recursion stays at zero).
        if (boundaryMode == EMBND_OPEN)
        {
                pmlPsiX.assign(gw * gh, 0.0);
                pmlPsiY.assign(gw * gh, 0.0);
        }
        else
        {
                pmlPsiX.clear();
                pmlPsiX.shrink_to_fit();
                pmlPsiY.clear();
                pmlPsiY.shrink_to_fit();
        }
        // visible canvas window in cells; at regionScale < 1 the domain is
        // smaller than the canvas and the renderer MAGNIFIES it by renderScale
        // (0.5x fix: the field still covers the full screen - the space is
        // half, not the display). At regionScale >= 1 the canvas is a 1:1
        // window into the (possibly larger) domain.
        renderScale = (regionScale < 1.0f) ? int(1.0f / regionScale + 0.5f) : 1;
        if (renderScale < 1)
        {
                renderScale = 1;
        }
        renderW = XRES / (cellSize * renderScale);
        renderH = YRES / (cellSize * renderScale);
        // visible window is always centred on the simulation domain
        renderOffX = padL + (visW - renderW) / 2;
        renderOffY = padT + (visH - renderH) / 2;
        forceBarValue = frequency;
        forceTimeZero = 0;
        t = 0;
        filterCount = 0;
        SetupSources();
        SetDamping();
        EMF_DBG("EMField: grid %dx%d cells (cell %dpx, region %g, boundary %d, pad %d, vis %dx%d, margin %d, renderScale %d, render window %dx%d, renderOff %d,%d)\n",
                gw, gh, cellSize, double(regionScale), boundaryMode, padL, visW, visH, margin, renderScale, renderW, renderH, renderOffX, renderOffY);
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
                if (cell.resonant)
                {
                        cell.jz = 0;
                }
        }
        // reset the CPML face filter states too (a uniform field has zero face
        // differences, so a zero recursion is the consistent rest state)
        if (!pmlPsiX.empty())
        {
                std::fill(pmlPsiX.begin(), pmlPsiX.end(), 0.0);
                std::fill(pmlPsiY.begin(), pmlPsiY.end(), 0.0);
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
        cell.az = 1e-10;
        cell.dazdt = 1e-10;
        cell.epos = 0;
}

int EMField::CellIndex(int px, int py) const
{
        // map a pixel on the visible canvas to the corresponding cell of the
        // simulation domain. Screen pixels per cell = cellSize * renderScale:
        // with regionScale > 1 the simulation extends beyond the visible
        // canvas, so the cell offset is renderOffX/Y (centred on the domain);
        // with regionScale = 1 renderOffX == padL so this reduces to the plain
        // pixel -> cell mapping; with regionScale = 0.5 the domain is half the
        // canvas and rendered magnified 2x, so a screen pixel maps to
        // px / (cellSize*2) in domain space; pixels whose domain position is
        // outside the domain clamp onto the nearest domain edge cell.
        int cs = cellSize * renderScale;
        int cx = std::clamp(px / cs, 0, renderW - 1) + renderOffX;
        int cy = std::clamp(py / cs, 0, renderH - 1) + renderOffY;
        cx = std::clamp(cx, padL, padL + visW - 1);
        cy = std::clamp(cy, padT, padT + visH - 1);
        return cx + cy * gw;
}

bool EMField::PixelInDomain(int px, int py) const
{
        // exact inverse of the render mapping: every pixel the renderer draws
        // from a cell counts as in-domain. The last partial cell column/row of
        // the canvas (XRES not divisible by the zoomed cell size, e.g. cs=8:
        // 612/8 = 76.5) clamps into the nearest cell, same as CellIndex - the
        // field visually covers the whole canvas and tools/particles at the
        // canvas edge agree on which cell they belong to.
        int cs = cellSize * renderScale;
        int cx = std::min(px / cs, renderW - 1) + renderOffX;
        int cy = std::min(py / cs, renderH - 1) + renderOffY;
        return cx >= padL && cx < padL + visW && cy >= padT && cy < padT + visH;
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
        // PERF: only the simulated domain rectangle can carry materials - the
        // OPEN absorber band and the regionScale > 1 padding are vacuum by
        // construction (every material writer checks PixelInDomain), and the
        // PERIODIC ghost ring is re-copied from the domain edge before every
        // sub-step. Iterating just the domain saves the full-grid sweep on the
        // padded geometries; behaviour is bit-identical.
        for (int y = padT; y < padT + visH; y++)
        {
                for (int x = padL; x < padL + visW; x++)
                {
                        auto &cell = cells[x + y * gw];
                        cell.perm = 1;
                        cell.conductivity = 0;
                        cell.mx = 0;
                        cell.my = 0;
                        cell.medium = 0;
                        cell.jzext = 0;
                        cell.resonant = false;
                        cell.heatable = false;
                        cell.magpowder = false;
                        cell.magsolid = false;
                }
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
                        if (!PixelInDomain(int(p.x), int(p.y)))
                        {
                                continue; // outside the simulated region (regionScale < 1)
                        }
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
                        if (PixelInDomain(int(p.x), int(p.y)))
                        {
                                auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                                cell.resonant = true;
                                cell.heatable = true; // resonant media dissipate absorbed energy
                        }
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
                if (!PixelInDomain(int(p.x), int(p.y)))
                {
                        // outside the simulated region (regionScale < 1): inert
                        continue;
                }
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
                // (Meissner effect), above it they quench to a fair conductor.
                // The mapping is ORDER-INDEPENDENT (deterministic even when several
                // particles share one cell, e.g. on coarse grids): conductivity
                // takes the maximum, permeability the minimum (the diamagnetic
                // Meissner response may only pull perm DOWN, never overwrite a
                // stronger ferromagnet), and the force flags accumulate.
                if (p.type == PT_SCND || p.type == PT_SCPW)
                {
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                        if (p.temp < EM_SC_TC)
                        {
                                cell.conductivity = std::max(cell.conductivity, 1.0f);
                                cell.perm = ClampPerm(std::min(cell.perm, EM_SC_PERM), PermMax());
                                cell.heatable = cell.heatable && false; // cold SC never heats
                                cell.magpowder = cell.magpowder || (p.type == PT_SCPW); // levitates over magnets
                                cell.magsolid = cell.magsolid || (p.type == PT_SCND);   // solid Meissner body, mild push
                        }
                        else
                        {
                                // quenched: resists again and dissipates; perm left
                                // as-is so a shared cell keeps its other materials
                                cell.conductivity = std::max(cell.conductivity, EM_SC_QUENCH_CONDUCT);
                                cell.heatable = true;
                        }
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
                // --- Task 6: EMTX transmitter radiates a burst on SPRK --------
                // The transmitter is NOT a vanilla conductor: it watches its
                // 4-neighbourhood from its own update() and starts a burst
                // (tmp = frames left) when an adjacent conductor is sparked.
                // Configuration: .ctype = carrier frequency 1..40 (0 = global
                // setting), .life = amplitude*4 (0 = unit amplitude). Keeping
                // the config in ctype/life is exactly what the task asks for;
                // because EMTX never becomes SPRK itself, a spark cycle cannot
                // clobber the user's values. The emission happens HERE (before
                // the wave sub-steps) - writing jzext from an element update
                // would be erased by the next SyncMaterials before the wave
                // ever sees it.
                if (p.type == PT_EMTX)
                {
                        if (p.tmp > 0 && PixelInDomain(int(p.x), int(p.y)))
                        {
                                // burst envelope: smooth Hann window over the
                                // remaining burst frames
                                float total = std::max(1, p.tmp2);
                                float env = std::sin(float(std::numbers::pi_v<float>) * (1.0f - float(p.tmp) / total));
                                float freq = frequency;
                                if (p.ctype >= 1 && p.ctype <= 40)
                                {
                                        freq = float(p.ctype);
                                }
                                float amp = (p.life > 0) ?
                                        std::clamp(float(p.life) / 4.0f, 0.05f, EM_JZEXT_MAX) : 1.0f;
                                auto &cell = cells[CellIndex(int(p.x), int(p.y))];
                                cell.jzext = double(amp * env) * std::sin(double(freq) * t * EM_FREQ_MULT);
                                // countdown runs once per frame (SyncMaterials is
                                // called once per frame from Update())
                                p.tmp -= 1;
                        }
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
                        if (p.life > 0 && PixelInDomain(int(p.x), int(p.y)))
                        {
                        auto &cell = cells[CellIndex(int(p.x), int(p.y))];
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
                        if (!PixelInDomain(int(p.x), int(p.y)))
                        {
                                break; // outside the simulated region: inert
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
                                cell.perm = ClampPerm(mat.perm, PermMax());
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
        // PERF: domain rectangle only, see the reset pass above - overrides can
        // never exist outside the domain.
        for (int oy = padT; oy < padT + visH; oy++)
        {
                for (int ox = padL; ox < padL + visW; ox++)
                {
                        auto &cell = cells[ox + oy * gw];
                        // task 10: the J override applies wherever a current can flow:
                        // dedicated current-source cells AND conductor cells (exciting a
                        // current inside an EM-zone conductor). Cells that lost their
                        // material get the override wiped below anyway.
                        if (cell.ovMask & EM_OV_JZ && (cell.jzext != 0 || cell.conductivity > 0))
                        {
                                // the override carries the full signed current
                                cell.jzext = double(std::clamp(cell.ovJz, -EM_JZEXT_MAX, EM_JZEXT_MAX));
                        }
                        if (cell.perm > 1) // ferromagnet
                        {
                                if (cell.ovMask & EM_OV_PERM)
                                {
                                        cell.perm = ClampPerm(cell.ovPerm, PermMax());
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
                        // EMADJ fix: the conductivity override used to sit in an `else if`
                        // behind the magnet branch, so it was silently dropped on every
                        // cell that ALSO carries magnetization/permeability (EMFM, FE,
                        // FEPW: perm 3..5 AND conductivity .4..5). On exactly those
                        // materials the EMADJ add mode could never raise the conductivity
                        // (the override was stored but never became effective), which is
                        // the "加法模式无法连续增加" report. Conductivity and magnetization
                        // are orthogonal cell properties, so apply it unconditionally.
                        if (cell.ovMask & EM_OV_CONDUCT && cell.conductivity > 0)
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
                                        ClearCellOverrides(ox + oy * gw);
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
        }

        CalcBoundaries();
}

void EMField::CalcBoundaries()
{
        // Mark all cells where the permeability, medium, magnetization or resonance
        // differs from one of the neighbours; the wave equation needs the hard path there.
        // PERF: the wave loops only ever read the flag inside the simulated domain
        // rectangle, so the sweep stays inside it too (the padded band is uniform
        // vacuum and its flag stays false from the Cell default).
        for (int y = std::max(1, padT); y < std::min(gh - 1, padT + visH); y++)
        {
                for (int x = std::max(1, padL); x < std::min(gw - 1, padL + visW); x++)
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
        // --- CPML face coefficient profiles (task 11) ---------------------------
        // The band cells do not use `damp` (their az/dazdt kinematics are the
        // undamped interior law); all absorption lives in the one-pole psi
        // recursion on the face differences. What is built here are the
        // per-face recursion coefficients from the standard CPML recipe:
        //   b = exp(-(sigma/kappa + alpha))
        //   c = sigma*(b-1) / (kappa*(kappa*sigma + alpha))
        // with a quartic sigma profile from the interface (sigma ~ 0, the
        // innermost face behaves like the interior, so the interface is
        // seamless) to the calibrated peak at the outer edge. alpha is a small
        // floor at the interface decaying to 0 outward: it is what makes the
        // stretched coordinate finite at DC, killing the low-frequency pole the
        // old split-field design had. kappa = 1 keeps the band wave speed
        // identical to the interior (raised above 1 only to attenuate grazing
        // evanescent tails if measurements demand it). The peak scales as
        // 1/cellSize (sigma_phys*cellSize = const) so the attenuation per PIXEL
        // is resolution-invariant, matching the task-7 speed/wavelength fix.
        pmlBPX.assign(gw, 1.0f);
        pmlCPX.assign(gw, 0.0f);
        pmlBPY.assign(gh, 1.0f);
        pmlCPY.assign(gh, 0.0f);
        if (boundaryMode == EMBND_OPEN)
        {
                float sigTau = EM_PML_SIGMA * taddEff / EM_TADD_SUB;
                float alTau = EM_PML_ALPHA * taddEff / EM_TADD_SUB;
                int D = std::max(2, padL);
                // face slot -> profile depth r (0 = interface face, 1 = outer
                // edge). psi slot i stores the +x face of cell i: the left band
                // occupies slots 0..padL-1 with the INTERFACE at slot padL-1,
                // the right band slots gw-padL-1..gw-2 with the INTERFACE at
                // slot gw-padL-1 (slot gw-2 faces the dead ring).
                auto buildFaces = [&](std::vector<float> &bp, std::vector<float> &cp,
                                      int nSlots, bool outwardPositive)
                {
                        for (int s = 0; s < nSlots; s++)
                        {
                                // s runs over the face slots; for the left/top
                                // bands slot 0 is the OUTERMOST face (against
                                // the dead ring) and slot nSlots-1 the interface
                                // face, for the right/bottom bands the other
                                // way round. r = profile depth: ~0 at the
                                // interface (seamless), ~1 at the outer edge.
                                double r = outwardPositive ? (s + 0.5) / double(nSlots)
                                                           : 1.0 - (s + 0.5) / double(nSlots);
                                double sig = sigTau * std::pow(r, double(EM_PML_POWER));
                                double al = alTau * (1.0 - r);
                                // bilinear-exact one-pole for the correction
                                // target -sigma/(alpha+sigma+i*omega):
                                //   b = (2-h)/(2+h), g = -2*sigma/(2+h)
                                // with h = sigma/kappa + alpha (per sub-step).
                                // DC check: psi/FD -> -sigma/h (the exact
                                // stretched-coordinate value).
                                double h = sig / EM_PML_KAPPA + al;
                                double b = (2.0 - h) / (2.0 + h);
                                double c = -2.0 * sig / (2.0 + h);
                                int slot = outwardPositive ? (gw - padL - 1 + s) : s;
                                bp[slot] = float(b);
                                cp[slot] = float(c);
                        }
                };
                D = padL;
                buildFaces(pmlBPX, pmlCPX, D, false);          // left band: slot 0 outer
                buildFaces(pmlBPX, pmlCPX, D, true);           // right band: slot n-1 outer
                D = std::max(2, padT);
                for (int s = 0; s < D; s++)                                    // top band
                {
                        double r = 1.0 - (s + 0.5) / double(D);
                        double sig = sigTau * std::pow(r, double(EM_PML_POWER));
                        double al = alTau * (1.0 - r);
                        double h = sig / EM_PML_KAPPA + al;
                        double b = (2.0 - h) / (2.0 + h);
                        double c = -2.0 * sig / (2.0 + h);
                        pmlBPY[s] = float(b);
                        pmlCPY[s] = float(c);
                }
                for (int s = 0; s < D; s++)                                    // bottom band
                {
                        int slot = gh - D - 1 + s;
                        double r = (s + 0.5) / double(D);
                        double sig = sigTau * std::pow(r, double(EM_PML_POWER));
                        double al = alTau * (1.0 - r);
                        double h = sig / EM_PML_KAPPA + al;
                        double b = (2.0 - h) / (2.0 + h);
                        double c = -2.0 * sig / (2.0 + h);
                        pmlBPY[slot] = float(b);
                        pmlCPY[slot] = float(c);
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
                // VISIBLE area, inset by one cell from its edge; the outflow
                // band lives entirely OUTSIDE the visible canvas, so no
                // margin-based inset is needed anymore
                sourceCount *= 2;
                int inset = 1;
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
                int inset = 1;
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
                        if (oe.jz != 0 || oe.jzext != 0 || oe.conductivity > 0)
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

// Once per frame, after the wave sub-steps: the force and heating couplings
// between the field and the matter inside it. All of these act on the CELL
// materials (never on individual charges - the real-zone charge carriers were
// removed from the mod), so a single pass over the domain is enough:
//  * magnetic pressure on ferromagnetic / diamagnetic powders (iron filings
//    pull toward magnets, pyrolytic graphite and superconductors are pushed
//    away) and the much milder version for solid magnetic materials;
//  * the motor effect (j x B force on current-carrying conductors);
//  * Joule heating of resistive real-zone conductors (superconductors below
//    Tc have no resistance and never heat up).
void EMField::ApplyFieldForces(int substeps)
{
        auto &parts = sim.parts;
        (void)parts;
        // PERF: iterate the visible canvas window only - cells outside it can
        // never touch particles (pmap is canvas-sized), the old full-grid loop
        // just burned a bounds check per band/padding cell per frame. Bit
        // identical: the skipped cells all hit the px0/py0 continue below.
        int cs = cellSize * renderScale; // screen px per cell (0.5x fix: zoomed)
        int cx0 = std::max(1, renderOffX);
        int cy0 = std::max(1, renderOffY);
        int cx1 = std::min(gw - 1, renderOffX + renderW);
        int cy1 = std::min(gh - 1, renderOffY + renderH);
        for (int cy = cy0; cy < cy1; cy++)
        {
                for (int cx = cx0; cx < cx1; cx++)
                {
                        int gi = cx + cy * gw;
                        auto &cell = cells[gi];
                        // pixel position of this cell on the visible canvas:
                        // the visible canvas is centred on the simulation
                        // domain, so cell (cx,cy) maps to pixel
                        // ((cx - renderOffX) * cs, (cy - renderOffY) * cs).
                        int px0 = (cx - renderOffX) * cs;
                        int py0 = (cy - renderOffY) * cs;
                        if (px0 < 0 || py0 < 0 || px0 >= XRES || py0 >= YRES)
                        {
                                continue; // cell outside the visible canvas
                        }
                        if (cell.magpowder || cell.magsolid)
                        {
                                // |B|^2 of the dynamic field: iron filings must respond to a
                                // magnet exactly like the applet diamagnet/ferromagnet pair
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
                                for (int sy = 0; sy < cs; sy++)
                                {
                                        for (int sx = 0; sx < cs; sx++)
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
                        // (the rail-gun / motor effect). Both the cell's induced
                        // jz and external jzext drive it.
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
                                        for (int sy = 0; sy < cs; sy++)
                                        {
                                                for (int sx = 0; sx < cs; sx++)
                                                {
                                                        int px = std::min(px0 + sx, XRES - 1);
                                                        int py = std::min(py0 + sy, YRES - 1);
                                                        int r = sim.pmap[py][px];
                                                        if (!r) continue;
                                                        auto &p = parts.data[ID(r)];
                                                        // apply to EM-zone conductors and real-zone
                                                        // conductors alike - the force is on the bulk
                                                        // material carrying the current
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
                                        for (int sy = 0; sy < cs; sy++)
                                        {
                                                for (int sx = 0; sx < cs; sx++)
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

// Once per frame, after the wave sub-steps: the EM -> vanilla direction of the
// current interop.
//
// EMAN antenna (task 5): a coupler that bridges the EM field to vanilla SPRK.
// The EM-zone elements CONNECTED to the antenna ARE its antenna: EMAN reads the
// excitation current |jz + jzext| (the induced current driven by the wave
// inside those elements) of the 3x3 EM-cell neighbourhood around itself and
// compares it against its threshold. The threshold is .ctype in raw current
// units (0 = default EMAN_THRESHOLD_DEFAULT, capped at EMAN_THRESHOLD_MAX); set
// it with the PROP tool or by drawing with a value. When the excitation
// exceeds the threshold, the antenna fires a vanilla SPRK on every adjacent
// vanilla conductor (PROP_CONDUCTS, not already sparked / not in cooldown).
//
// Pair with an EMTX transmitter (vanilla -> EM) for bidirectional coupling:
// EMTX radiates the wave, the wave drives currents inside the EM-zone
// conductors touching EMAN, EMAN sparks the next circuit.
void EMField::InteropParticles()
{
        auto &sd = SimulationData::CRef();
        auto &elements = sd.elements;
        auto &parts = sim.parts;

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
                // threshold from .ctype (raw current units)
                float thr = EMAN_THRESHOLD_DEFAULT;
                if (p.ctype > 0)
                {
                        thr = std::min(float(p.ctype), EMAN_THRESHOLD_MAX);
                }
                // gather the excitation current from the 3x3 neighbourhood of EM cells
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
                                double j = std::abs(nc.jz + nc.jzext);
                                maxExcite = std::max(maxExcite, j);
                        }
                }
                if (maxExcite < double(thr))
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

// One sub-step of the CPML absorbing band (task 11), VELOCITY pass.
// Two passes, both driven by the CURRENT state (az^n, dazdt = v^(n-1/2)):
//
//   1. psi recursion (one-pole CPML filter per band face, per axis),
//      trapezoidal in time for an exact bilinear mapping of the target
//      response -sigma/(alpha+sigma+i*omega):
//        psi' = b*psi + (g/2)*(FD + FD_prev)
//      FD_prev (the previous sub-step's face difference) is reconstructed
//      EXACTLY from the current state:
//        FD_prev = FD - tadd^2 * (dazdt[+]-dazdt)      (az integration law)
//      so no history arrays are needed. This pass must run BEFORE the drive
//      touches dazdt (the reconstruction needs the un-updated velocities).
//   2. band dazdt drive (the EXACT interior stencil with the face
//      differences replaced by their CPML-corrected versions):
//        dazdt += 0.25*( (FD+psi)_x+ - (FD+psi)_x- + (FD+psi)_y+ - (FD+psi)_y- )
//
// The band cells share the interior kinematics; the stretched-coordinate
// match keeps the band interface reflectionless at every frequency the
// one-pole can represent - including DC, where the old split-field design
// failed (|R|^2 up to 2e-2 at f=2).
void EMField::PmlStepA()
{
        // 1. psi recursion over every band face (in-place: each face update
        //    reads only its own psi state plus az/dazdt of its two cells).
        //    x faces of the left / right bands (full height, corners
        //    included), y faces of the top / bottom bands.
        double tadd2 = double(taddEff) * double(taddEff);
        for (int j = 1; j < gh - 1; j++)
        {
                for (int i = 0; i < padL; i++)
                {
                        int gi = i + j * gw;
                        double fd = cells[gi + 1].az - cells[gi].az;
                        double fd_old = fd - tadd2 * (cells[gi + 1].dazdt - cells[gi].dazdt);
                        pmlPsiX[gi] = float(double(pmlBPX[i]) * pmlPsiX[gi]
                                    + double(pmlCPX[i]) * 0.5 * (fd + fd_old));
                }
                for (int i = gw - padL - 1; i < gw - 1; i++)
                {
                        int gi = i + j * gw;
                        double fd = cells[gi + 1].az - cells[gi].az;
                        double fd_old = fd - tadd2 * (cells[gi + 1].dazdt - cells[gi].dazdt);
                        pmlPsiX[gi] = float(double(pmlBPX[i]) * pmlPsiX[gi]
                                    + double(pmlCPX[i]) * 0.5 * (fd + fd_old));
                }
        }
        for (int i = 1; i < gw - 1; i++)
        {
                for (int j = 0; j < padT; j++)
                {
                        int gi = i + j * gw;
                        double fd = cells[gi + gw].az - cells[gi].az;
                        double fd_old = fd - tadd2 * (cells[gi + gw].dazdt - cells[gi].dazdt);
                        pmlPsiY[gi] = float(double(pmlBPY[j]) * pmlPsiY[gi]
                                    + double(pmlCPY[j]) * 0.5 * (fd + fd_old));
                }
                for (int j = gh - padT - 1; j < gh - 1; j++)
                {
                        int gi = i + j * gw;
                        double fd = cells[gi + gw].az - cells[gi].az;
                        double fd_old = fd - tadd2 * (cells[gi + gw].dazdt - cells[gi].dazdt);
                        pmlPsiY[gi] = float(double(pmlBPY[j]) * pmlPsiY[gi]
                                    + double(pmlCPY[j]) * 0.5 * (fd + fd_old));
                }
        }
        // 2. band dazdt drive. Cells in a corner region are filtered on BOTH
        //    axes; cells in a single-axis band keep the plain interior
        //    difference on the other axis (their psi slots do not exist there,
        //    which is exactly the CPML convention: an unfiltered axis is a
        //    plain difference).
        auto stepA = [&](int i, int j, bool fx, bool fy)
        {
                int gi = i + j * gw;
                double lx;
                if (fx)
                {
                        lx = (cells[gi + 1].az - cells[gi].az + pmlPsiX[gi])
                           - (cells[gi].az - cells[gi - 1].az + pmlPsiX[gi - 1]);
                }
                else
                {
                        lx = cells[gi + 1].az - 2.0 * cells[gi].az + cells[gi - 1].az;
                }
                double ly;
                if (fy)
                {
                        ly = (cells[gi + gw].az - cells[gi].az + pmlPsiY[gi])
                           - (cells[gi].az - cells[gi - gw].az + pmlPsiY[gi - gw]);
                }
                else
                {
                        ly = cells[gi + gw].az - 2.0 * cells[gi].az + cells[gi - gw].az;
                }
                // accumulate, exactly like the interior leapfrog does (the
                // interior pass adds its acceleration to dazdt; the band cells
                // are skipped by the interior pass, so the band drive must add
                // the acceleration here - overwriting would zero the band
                // velocity every sub-step and the layer would reflect like a
                // hard wall)
                cells[gi].dazdt += 0.25 * (lx + ly);
        };
        // left / right bands, full height (they include the four corners);
        // y filtering applies only to rows inside the top / bottom bands
        for (int j = 1; j < gh - 1; j++)
        {
                bool fy = j < padT || j >= gh - padT;
                for (int i = 1; i < padL; i++)
                {
                        stepA(i, j, true, fy);
                }
                for (int i = gw - padL; i < gw - 1; i++)
                {
                        stepA(i, j, true, fy);
                }
        }
        // top / bottom bands, middle columns: x is unfiltered (interior columns)
        for (int i = padL; i < gw - padL; i++)
        {
                for (int j = 1; j < padT; j++)
                {
                        stepA(i, j, false, true);
                }
                for (int j = gh - padT; j < gh - 1; j++)
                {
                        stepA(i, j, false, true);
                }
        }
}

// One sub-step of the CPML absorbing band, POSITION pass: runs after the
// interior az sweep. Every band cell advances az with the SAME undamped
// kinematic law as the interior - the stretched coordinate never touches the
// kinematics, which is what keeps the layer matched (the absorption lives
// entirely in the psi memory of the velocity pass).
void EMField::PmlStepB()
{
        double tadd2 = double(taddEff) * double(taddEff);
        auto stepB = [&](int i, int j)
        {
                int gi = i + j * gw;
                cells[gi].az += cells[gi].dazdt * tadd2;
        };
        // left / right bands, full height (they include the four corners)
        for (int j = 1; j < gh - 1; j++)
        {
                for (int i = 1; i < padL; i++)
                {
                        stepB(i, j);
                }
                for (int i = gw - padL; i < gw - 1; i++)
                {
                        stepB(i, j);
                }
        }
        // top / bottom bands, middle columns
        for (int i = padL; i < gw - padL; i++)
        {
                for (int j = 1; j < padT; j++)
                {
                        stepB(i, j);
                }
                for (int j = gh - padT; j < gh - 1; j++)
                {
                        stepB(i, j);
                }
        }
}

// PERF (task 8 v3): once per frame, measure the wave activity in the PML band
// plus an interior guard strip. While everything is below EM_PML_QUIET the
// band passes are skipped for the whole frame - numerically exact, see the
// EMField.h comment. NaN/Inf states compare false against the threshold and
// therefore count as active, so the gate never hides a poisoned band from the
// clamp pass.
void EMField::ScanPmlActivity()
{
        if (boundaryMode != EMBND_OPEN || pmlPsiX.empty())
        {
                pmlQuiet = false;
                return;
        }
        // guard strip: interior cells within `guard` cells of the band
        // interface. 8 cells is deeper than the max per-frame wave travel
        // (16 sub-steps * taddEff/2 = 4 cells at cs=1, less at any other
        // resolution/speed), so a wave heading for the band always trips the
        // gate at least one frame before it can touch the interface.
        const int guard = 8;
        const int gx1 = std::min(padL + guard, gw);
        const int gx2 = std::max(gw - padL - guard, 0);
        const int gy1 = std::min(padT + guard, gh);
        const int gy2 = std::max(gh - padT - guard, 0);
        bool quiet = true;
        auto probe = [&](int gi)
        {
                const auto &oe = cells[gi];
                // NaN is never <= threshold, so poisoned cells read as active
                if (!(std::abs(oe.az) <= double(EM_PML_QUIET)) ||
                    !(std::abs(oe.dazdt) <= double(EM_PML_QUIET)) ||
                    !(std::abs(pmlPsiX[gi]) <= EM_PML_QUIET) ||
                    !(std::abs(pmlPsiY[gi]) <= EM_PML_QUIET))
                {
                        quiet = false;
                }
        };
        // left / right band + guard strips, full height
        for (int j = 0; j < gh && quiet; j++)
        {
                for (int i = 0; i < gx1; i++)
                {
                        probe(i + j * gw);
                        if (!quiet)
                        {
                                pmlQuiet = false;
                                return;
                        }
                }
                for (int i = gx2; i < gw; i++)
                {
                        probe(i + j * gw);
                        if (!quiet)
                        {
                                pmlQuiet = false;
                                return;
                        }
                }
        }
        // top / bottom band + guard strips, middle columns
        for (int i = gx1; i < gx2 && quiet; i++)
        {
                for (int j = 0; j < gy1; j++)
                {
                        probe(i + j * gw);
                        if (!quiet)
                        {
                                pmlQuiet = false;
                                return;
                        }
                }
                for (int j = gy2; j < gh; j++)
                {
                        probe(i + j * gw);
                        if (!quiet)
                        {
                                pmlQuiet = false;
                                return;
                        }
                }
        }
        pmlQuiet = quiet;
        EMF_DBG("EMField: PML activity gate %s\n", pmlQuiet ? "quiet (band passes skipped)" : "active");
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

        // The wave is integrated in fixed sub-steps. The per-sub-step timestep is
        // taddEff = EM_TADD_SUB*2/cellSize so the wave keeps the same pixel speed
        // and pixel wavelength at every grid resolution (task 7); the CFL bound
        // (perm ratio <= 2/taddEff^2) is enforced by the dynamic permeability
        // clamp. The clock t advances the CONSTANT EM_TADD_SUB per sub-step, so
        // source frequencies map to the same spatial wavelength at every
        // resolution as well. The traversal time of the simulated region depends
        // only on the region size, never on the grid resolution.
        int substeps = EM_SUBSTEPS[std::clamp(speed, 0, 4)];
        double tadd = taddEff;
        double tadd2 = double(taddEff) * double(taddEff);
        // band cells are integrated by the split-field PML (PmlStepA/B), not by
        // the interior wave update
        bool outflow = boundaryMode == EMBND_OPEN;
        // PERF (task 8 v3): skip the band passes entirely while the band and
        // its guard strip are quiet; see ScanPmlActivity()
        if (outflow)
        {
                ScanPmlActivity();
        }
        // perf: the wave loops visit the interior rectangle ONLY - the band is
        // handled by the dedicated PML passes, so no per-cell outflow branch is
        // needed inside the hot loop (bit-exact, the branch used to skip
        // exactly these cells)
        const int i0 = outflow ? padL : 1;
        const int i1 = outflow ? gw - padL : gw - 1;
        const int j0 = outflow ? padT : 1;
        const int j1 = outflow ? gh - padT : gh - 1;

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
                for (int j = j0; j < j1; j++)
                {
                        for (int i = i0; i < i1; i++)
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
                                        a = oe.perm * basis + jz;
                                }
                                else
                                {
                                        // easy way
                                        double basis = (oew.az + oee.az + oen.az + oes.az) * .25;
                                        a = oe.jz + oe.jzext - (oe.az - basis);
                                }
                                oe.dazdt = oe.dazdt * oe.damp + a * forcecoef;
                        }
                }

                // --- split-field PML band, velocity pass (task 8) --------------------
                // Must run BETWEEN the interior dazdt and az sweeps: the band
                // cells then read their neighbours' u^n values, exactly like the
                // interior stencil does, and the leapfrog stays centred.
                // Gated by the activity scan while the band is quiet (v3).
                if (outflow && !pmlQuiet)
                {
                        PmlStepA();
                }

                // --- second pass: integrate az, and compute induced currents in conductors ---
                for (int j = j0; j < j1; j++)
                {
                        for (int i = i0; i < i1; i++)
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
                // --- split-field PML band, position pass (task 8) --------------------
                // The split displacement is damped with the same profile as its
                // velocity, which makes the layer perfectly matched: a wave
                // crossing into the band keeps its impedance and leaves the
                // domain without reflecting (measured |R|^2 ~ 1e-11 at the
                // default frequency, ~5e-6 at f=5, vs ~0.01..0.19 for the old
                // advective outflow layer). No sweep-order constraints: the
                // velocity pass only reads az, this pass only touches own-cell
                // state. Gated by the activity scan while the band is quiet (v3).
                if (outflow && !pmlQuiet)
                {
                        PmlStepB();
                }
                t += EM_TADD_SUB;

#if EMFIELD_DEBUG
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

        // once per frame, like the applet: high-frequency noise filter
        FilterGrid();

        // once per frame: field -> matter forces and heating (magnetic pressure,
        // motor effect, Joule heat)
        ApplyFieldForces(substeps);

        // --- Task 9: energy conservation and the anti-divergence net -------------
        // In the linear regime the leapfrog scheme conserves energy exactly; the
        // historically divergent couplings are gone for good:
        //  * the permeability contrast is CFL-clamped (dynamic bound 2/taddEff^2),
        //    which removes the old ferromagnet self-excitation at its root;
        //  * all external currents are bounded by construction (every jzext
        //    writer assigns a value <= EM_JZEXT_MAX).
        // What remains is a two-level safety net, applied unconditionally:
        //  1. a hard clamp at EM_FIELD_CLAMP (a runaway must never reach inf),
        //  2. a soft bleed above EM_SOFT_LIMIT (1e3) that gently drains cells
        //     whose amplitude grew beyond anything legitimate feedback produces
        //     (legit f=2 standing waves peak near |az| ~ 900),
        //  3. a NaN/Inf guard: one poisoned cell used to be able to spread NaN
        //     through the whole grid permanently; now it is zeroed instead.
        // PERF: a cheap fused max-scan decides whether ANY cell needs attention;
        // in the normal (bounded) case the expensive per-cell pass is skipped
        // entirely. The pass itself is bit-for-bit the original one, so the
        // result is identical whether or not it runs.
        {
                // NOTE: double comparisons, exactly like the clamp pass below -
                // a float narrowing here could miss cells in the last ULP below
                // the soft limit and change the wave state
                // PERF: the band+guard cells were just probed by ScanPmlActivity
                // (quiet => finite and below 1e-9, nothing to clamp), so the
                // fast max-scan can stay inside the domain rectangle on a quiet
                // OPEN frame; the full grid is scanned whenever the band is active
                double maxAz = 0, maxDazdt = 0;
                bool poisoned = false;
                int scanN = (boundaryMode == EMBND_OPEN && pmlQuiet) ? padL + visW : gw * gh;
                if (boundaryMode == EMBND_OPEN && pmlQuiet)
                {
                        for (int y = padT; y < padT + visH; y++)
                        {
                                int gi = padL + y * gw;
                                for (int x = 0; x < visW; x++, gi++)
                                {
                                        const auto &cell = cells[gi];
                                        if (!std::isfinite(cell.az) || !std::isfinite(cell.dazdt))
                                        {
                                                poisoned = true;
                                                break;
                                        }
                                        maxAz = std::max(maxAz, std::abs(cell.az));
                                        maxDazdt = std::max(maxDazdt, std::abs(cell.dazdt));
                                }
                                if (poisoned)
                                {
                                        break;
                                }
                        }
                }
                else
                {
                        for (int gi = 0; gi < scanN; gi++)
                        {
                                const auto &cell = cells[gi];
                                if (!std::isfinite(cell.az) || !std::isfinite(cell.dazdt))
                                {
                                        poisoned = true;
                                        break;
                                }
                                maxAz = std::max(maxAz, std::abs(cell.az));
                                maxDazdt = std::max(maxDazdt, std::abs(cell.dazdt));
                        }
                }
                if (poisoned || maxAz > double(EM_SOFT_LIMIT) || maxDazdt > double(EM_SOFT_LIMIT))
                {
                        for (int gi = 0; gi < gw * gh; gi++)
                        {
                                auto &cell = cells[gi];
                                if (!std::isfinite(cell.az) || !std::isfinite(cell.dazdt))
                                {
                                        cell.az = 1e-10;
                                        cell.dazdt = 1e-10;
                                        cell.jz = 0;
                                        cell.jzext = 0;
                                        // keep the CPML band bookkeeping consistent
                                        if (!pmlPsiX.empty())
                                        {
                                                pmlPsiX[gi] = 0;
                                                pmlPsiY[gi] = 0;
                                        }
#if EMFIELD_DEBUG
                                        fieldClampHits++;
#endif
                                        continue;
                                }
                                if (std::abs(cell.az) > double(EM_FIELD_CLAMP))
                                {
                                        cell.az = std::clamp(cell.az, -double(EM_FIELD_CLAMP), double(EM_FIELD_CLAMP));
#if EMFIELD_DEBUG
                                        fieldClampHits++;
#endif
                                }
                                if (std::abs(cell.dazdt) > double(EM_FIELD_CLAMP))
                                {
                                        cell.dazdt = std::clamp(cell.dazdt, -double(EM_FIELD_CLAMP), double(EM_FIELD_CLAMP));
#if EMFIELD_DEBUG
                                        fieldClampHits++;
#endif
                                }
                                if (std::abs(cell.az) > double(EM_SOFT_LIMIT))
                                {
                                        cell.az *= EM_SOFT_BLEED;
                                }
                                if (std::abs(cell.dazdt) > double(EM_SOFT_LIMIT))
                                {
                                        cell.dazdt *= EM_SOFT_BLEED;
                                }
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
                cell.ovPerm = ClampPerm(vali / 2.0f, PermMax());
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
                // EMADJ fix: ALWAYS derive the current direction from the synced
                // magnetization (mx,my), never from the raw ovDir override. The
                // override is folded into mx/my by every SyncMaterials pass, so
                // the angle below already includes it - ADD/SUB accumulate once
                // per 0.2s stroke event like every other property. Reading the
                // raw override instead made the brush's per-pixel multi-hits
                // (one EM cell is hit cellSize^2 times per stroke) rotate the
                // magnet value*hits degrees within a SINGLE stroke - the
                // "加法模式无法连续增加" report for the direction property.
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
        if (property == EMADJP_MAG_DIR)
        {
                // EMADJ fix: EffectiveProperty returns the direction as a 0..1
                // fraction of 2pi while the tool passes DEGREES. Adding degrees
                // straight onto the fraction made ADD/SUB rotate by value/360
                // degrees per application (i.e. practically nothing) - convert
                // the current direction to degrees first so all three modes
                // operate in the same unit.
                eff *= 360.0f;
        }
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
                // lower clamp matches the SyncMaterials application clamp
                // (0.01, the applet's adjustBar minimum) so ADD/SUB cannot get
                // stuck between two different floors
                cell.ovConduct = std::clamp(target, 0.01f, 1.0f);
                return true;
        case EMADJP_PERM:
                if (cell.perm <= 1)
                        return false;
                cell.ovMask |= EM_OV_PERM;
                cell.ovPerm = ClampPerm(target, PermMax());
                return true;
        case EMADJP_J:
                // task 10: works on current sources AND on conductors (sets the
                // externally driven current flowing inside the conductor)
                if (cell.jzext == 0 && cell.conductivity <= 0)
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
                // lower clamp matches the SyncMaterials application clamp so
                // subtracting past the floor cannot oscillate between 0 and .01
                cell.ovMag = std::clamp(target, 0.01f, 2.0f);
                return true;
        }
        }
        return false;
}
