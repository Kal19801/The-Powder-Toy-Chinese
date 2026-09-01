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
#define EMF_DBG(...) do { if (false) { std::fprintf(stderr, __VA_ARGS__); } } while (0)
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
        gw = std::max(8, XRES / cellSize);
        gh = std::max(8, YRES / cellSize);
        cells.assign(gw * gh, Cell{});
        margin = std::clamp(EM_MARGIN_AT_4 * EM_CELL_SIZE_DEFAULT / cellSize, 2, std::min(gw, gh) / 4);
        forceBarValue = frequency;
        forceTimeZero = 0;
        t = 0;
        filterCount = 0;
        SetupSources();
        SetDamping();
        EMF_DBG("EMField: grid %dx%d cells (cell size %dpx, margin %d)\n", gw, gh, cellSize, margin);
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
        int cx = std::clamp(px / cellSize, 0, gw - 1);
        int cy = std::clamp(py / cellSize, 0, gh - 1);
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
                        }
                        else
                        {
                                cell.conductivity = std::max(cell.conductivity, 0.3f);
                                cell.perm = 1;
                                cell.heatable = true;
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
        // exponential absorbing ramp on all four edges, like the applet's hidden margin
        for (int i = 0; i < margin; i++)
        {
                double da = std::exp(-(margin - i) * .002);
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
                // plane sources are drawn between pairs of corner points
                sourceCount *= 2;
                int x2 = gw - margin - 1;
                int y2 = gh - margin - 1;
                sources[0] = { margin, margin };
                sources[1] = { x2, margin };
                sources[2] = { margin, y2 };
                sources[3] = { x2, y2 };
        }
        else
        {
                // point sources sit around the centre of the canvas
                sources[0] = { gw / 2, margin + 1 };
                sources[1] = { gw / 2, gh - margin - 2 };
                sources[2] = { margin + 1, gh / 2 };
                sources[3] = { gw - margin - 2, gh / 2 };
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
                // superconductors are pushed away) and Joule heating of the real
                // conductors carrying induced current
                for (int cy = 1; cy < gh - 1; cy++)
                {
                        for (int cx = 1; cx < gw - 1; cx++)
                        {
                                int gi = cx + cy * gw;
                                auto &cell = cells[gi];
                                if (cell.magpowder)
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
                                        // field energy, diamagnets (perm < 1) down it; the
                                        // soft normalisation keeps the force scale invariant
                                        double sign = cell.perm > 1 ? 1.0 : -1.0;
                                        float fx = float(gx * sign * EM_POWDER_FORCE / (std::abs(e2) + 0.05));
                                        float fy = float(gy * sign * EM_POWDER_FORCE / (std::abs(e2) + 0.05));
                                        int px0 = cx * cellSize;
                                        int py0 = cy * cellSize;
                                        for (int sy = 0; sy < cellSize; sy++)
                                        {
                                                for (int sx = 0; sx < cellSize; sx++)
                                                {
                                                        int px = std::min(px0 + sx, XRES - 1);
                                                        int py = std::min(py0 + sy, YRES - 1);
                                                        int r = sim.pmap[py][px];
                                                        if (r)
                                                        {
                                                                auto &p = parts.data[ID(r)];
                                                                if (IsMagPowder(p.type, p.temp))
                                                                {
                                                                        p.vx += std::clamp(fx, -0.3f, 0.3f);
                                                                        p.vy += std::clamp(fy, -0.3f, 0.3f);
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
                                                int px0 = cx * cellSize;
                                                int py0 = cy * cellSize;
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

void EMField::Update()
{
        if (!enabled)
        {
                return;
        }
        SyncMaterials();
        CollectRealCharges();
        DepositRealCharges();

        // the wave is integrated in fixed sub-steps of the exact applet timestep;
        // this keeps the CFL bound (perm contrast <= 2/tadd^2 = 32) satisfied for
        // every speed setting, matches the applet dynamics much more closely than
        // one large step per frame, and conserves energy in the linear regime
        int substeps = EM_SUBSTEPS[std::clamp(speed, 0, 2)];
        double tadd = EM_TADD_SUB;
        double tadd2 = tadd * tadd;

        for (int substep = 0; substep < substeps; substep++)
        {
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
