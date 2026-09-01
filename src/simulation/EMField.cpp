#include "EMField.h"
#include "Simulation.h"
#include "ElementClasses.h"
#include "ElementDefs.h"
#include "common/tpt-rand.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>

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
        // Conductivity of the TPT elements mapped onto the EM grid, ported from the
        // applet's conductor palette (1 = perfect conductor, 0 = vacuum).
        struct ElemMaterial
        {
                int type;
                float conductivity;
                float perm;      // relative permeability
                int medium;      // dielectric constant contribution
                bool source;     // injects current into the field
                float jz;        // current injected when source
        };

        // NOTE: keep in sync with the element table; EMW/EMR are handled separately.
        const ElemMaterial elemMaterials[] = {
                { PT_METL, 1.0f, 1.0f, 0, false, 0 }, // perfect conductor, like the applet's perf. conductor
                { PT_GOLD, 1.0f, 1.0f, 0, false, 0 },
                { PT_WIRE, 1.0f, 1.0f, 0, false, 0 },
                { PT_INWR, 0.9f, 1.0f, 0, false, 0 },
                { PT_IRON, 0.7f, 5.0f, 0, false, 0 }, // ferromagnet, like the applet's addPerm(5)
                { PT_TUNG, 0.6f, 1.0f, 0, false, 0 },
                { PT_BMTL, 0.5f, 1.0f, 0, false, 0 },
                { PT_PSCN, 0.3f, 1.0f, 0, false, 0 },
                { PT_NSCN, 0.3f, 1.0f, 0, false, 0 },
                { PT_BIZR, 0.5f, 0.5f, 0, false, 0 }, // bismuth is a diamagnet, like the applet's addPerm(.5)
                { PT_GLAS, 0.0f, 1.0f, EM_MEDIUM_MAX, false, 0 },
                { PT_BGLA, 0.0f, 1.0f, EM_MEDIUM_MAX, false, 0 },
                { PT_QRTZ, 0.0f, 1.0f, EM_MEDIUM_MAX, false, 0 },
                { PT_SPRK, 0.0f, 1.0f, 0, true, 1.0f },   // sparking cell = current source; the metal
                                                          // around it forms the antenna
                { PT_ELEC, 0.0f, 1.0f, 0, true, -1.0f },  // electrons are a negative current
        };

        bool IsConductorForSpark(int type)
        {
                switch (type)
                {
                case PT_METL: case PT_GOLD: case PT_WIRE: case PT_INWR:
                case PT_IRON: case PT_TUNG: case PT_BMTL:
                case PT_PSCN: case PT_NSCN:
                        return true;
                }
                return false;
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

int EMField::CellIndex(int px, int py) const
{
        int cx = std::clamp(px / cellSize, 0, gw - 1);
        int cy = std::clamp(py / cellSize, 0, gh - 1);
        return cx + cy * gw;
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
                cell.resonant = false;
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
                        cells[CellIndex(int(p.x), int(p.y))].resonant = true;
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
                                cell.jzext += mat.jz;
                        }
                        else
                        {
                                if (mat.conductivity > cell.conductivity)
                                {
                                        cell.conductivity = mat.conductivity;
                                }
                        }
                        if (mat.perm != 1)
                        {
                                cell.perm = mat.perm;
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
        for (auto &cell : cells)
        {
                if (cell.ovMask & EM_OV_JZ && cell.jzext != 0)
                {
                        cell.jzext = (cell.jzext < 0 ? -1.0 : 1.0) * double(std::clamp(cell.ovJz, 0.01f, 1.0f));
                }
                if (cell.perm > 1) // ferromagnet
                {
                        if (cell.ovMask & EM_OV_PERM)
                        {
                                cell.perm = std::clamp(cell.ovPerm, 0.5f, 50.0f);
                        }
                        if (cell.ovMask & (EM_OV_MAGDIR | EM_OV_MAGSTR))
                        {
                                float str = (cell.ovMask & EM_OV_MAGSTR) ? std::clamp(cell.ovMag, 0.01f, 2.0f) : 1.0f;
                                if (cell.ovMask & EM_OV_MAGDIR)
                                {
                                        float angle = cell.ovDir * 2.0f * float(M_PI);
                                        cell.mx = str * std::cos(angle);
                                        cell.my = -str * std::sin(angle);
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
                        w2 = w + au * (M_PI / 38);
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
                        double wp = std::fmod(w, M_PI * 2);
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
                        double wp = std::fmod(w, M_PI * 2);
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
                        amp = 2 * std::clamp(cell.ovJz, 0.01f, 1.0f);
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

void EMField::Update()
{
        if (!enabled)
        {
                return;
        }
        double tadd = EM_TADD[std::clamp(speed, 0, 2)];
        SyncMaterials();
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

        // --- second pass: integrate az, and compute induced currents in conductors ---
        double tadd2 = tadd * tadd;
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

        FilterGrid();
        InteractParticles();
}

void EMField::InteractParticles()
{
        auto &parts = sim.parts;
        // Joule heating and induced sparking in conductors
        for (int cy = 1; cy < gh - 1; cy++)
        {
                for (int cx = 1; cx < gw - 1; cx++)
                {
                        int gi = cx + cy * gw;
                        auto &cell = cells[gi];
                        // --- Joule heating: induced current dissipates power in conductors ---
                        if (cell.conductivity > 0 && std::abs(cell.jz) > 0.01)
                        {
                                float jzc = std::clamp(float(cell.jz), -2.0f, 2.0f);
                                float heat = jzc * jzc * EM_JOULE_HEAT;
                                if (heat >= 0.05f)
                                {
                                        int px0 = cx * cellSize;
                                        int py0 = cy * cellSize;
                                        for (int sy = 0; sy <= 1; sy++)
                                        {
                                                for (int sx = 0; sx <= 1; sx++)
                                                {
                                                        int px = std::min(px0 + sx * cellSize, XRES - 1);
                                                        int py = std::min(py0 + sy * cellSize, YRES - 1);
                                                        int r = sim.pmap[py][px];
                                                        if (r)
                                                        {
                                                                auto &p = parts.data[ID(r)];
                                                                if (IsConductorForSpark(p.type) && p.temp < MAX_TEMP)
                                                                {
                                                                        p.temp = std::min(p.temp + heat, MAX_TEMP);
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                        // --- strong E fields spark antennas (this is how EM waves are received) ---
                        if (std::abs(cell.dazdt) > EM_INDUCED_SPARK_THRESHOLD)
                        {
                                int px0 = cx * cellSize;
                                int py0 = cy * cellSize;
                                for (int sy = 0; sy <= 1; sy++)
                                {
                                        for (int sx = 0; sx <= 1; sx++)
                                        {
                                                int px = std::min(px0 + sx * cellSize, XRES - 1);
                                                int py = std::min(py0 + sy * cellSize, YRES - 1);
                                                int r = sim.pmap[py][px];
                                                if (r)
                                                {
                                                        auto &p = parts.data[ID(r)];
                                                        if (IsConductorForSpark(p.type) && p.life == 0 &&
                                                            sim.rng.chance(1, EM_INDUCED_SPARK_CHANCE))
                                                        {
                                                                p.ctype = p.type;
                                                                p.type = PT_SPRK;
                                                                p.life = 4;
                                                                EMF_DBG("EMField: induced spark at %d,%d\n", px, py);
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }

        // Magnetic pressure deflects moving charged particles (electron optics)
        for (int i = 0; i < parts.active; ++i)
        {
                auto &p = parts.data[i];
                if (p.type != PT_ELEC)
                {
                        continue;
                }
                int gi = CellIndex(int(p.x), int(p.y));
                int cx = gi % gw;
                int cy = gi / gw;
                if (cx < 1 || cx >= gw - 1 || cy < 1 || cy >= gh - 1)
                {
                        continue;
                }
                auto b2 = [this](int g) {
                        double dx = cells[g - gw].az - cells[g + gw].az;
                        double dy = cells[g + 1].az - cells[g - 1].az;
                        return dx * dx + dy * dy;
                };
                double b2l = b2(gi - 1);
                double b2r = b2(gi + 1);
                double b2u = b2(gi - gw);
                double b2d = b2(gi + gw);
                p.vx += std::clamp((b2r - b2l) * EM_MAG_FORCE, -0.5, 0.5);
                p.vy += std::clamp((b2d - b2u) * EM_MAG_FORCE, -0.5, 0.5);
        }
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
                cell.ovPerm = std::clamp(strength, 0.01f, 1.0f) * 50.0f;
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
                cell.ovJz = std::clamp(strength, 0.01f, 1.0f);
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
