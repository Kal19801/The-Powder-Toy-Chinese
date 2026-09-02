#pragma once
#include "SimulationConfig.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

class Simulation;
struct Particle;

// EM cell material override bits, set by the unified EM adjust tool; the
// override only applies while the underlying particles still provide a matching
// material, just like the adjust modes of the original applet only affect cells
// of the matching type. Overrides are wiped again as soon as a cell loses its
// underlying material, so nothing leaks after the particles are erased.
constexpr uint8_t EM_OV_CONDUCT = 1 << 0; // conductivity of conductors
constexpr uint8_t EM_OV_PERM    = 1 << 1; // permeability of ferromagnets
constexpr uint8_t EM_OV_JZ      = 1 << 2; // magnitude/sign of current sources
constexpr uint8_t EM_OV_MEDIUM  = 1 << 3; // dielectric constant of dielectrics
constexpr uint8_t EM_OV_MAGSTR  = 1 << 4; // strength of permanent magnets
constexpr uint8_t EM_OV_MAGDIR  = 1 << 5; // direction of permanent magnets (angle = value * 2pi)

// the five non-directional adjust modes of the applet (MODE_ADJ_*), kept for
// the applet-faithful ApplyAdjustMode() port
enum EmAdjustMode
{
        EMADJM_CONDUCT = 0, // MODE_ADJ_CONDUCT
        EMADJM_PERM,        // MODE_ADJ_PERM
        EMADJM_J,           // MODE_ADJ_J
        EMADJM_MEDIUM,      // MODE_ADJ_MEDIUM
        EMADJM_MAG_DIR,     // MODE_ADJ_MAG_DIR
        EMADJM_MAG_STR,     // MODE_ADJ_MAG_STR
};

// material type classification, a faithful port of the applet's OscElement::getType();
// the order of the checks matters and must not be changed
enum EmCellType
{
        EMCT_NONE = 0,
        EMCT_DIAMAGNET,
        EMCT_FERROMAGNET,
        EMCT_MAGNET,
        EMCT_MEDIUM,
        EMCT_CONDUCTOR,
        EMCT_CURRENT,
};

// the six adjustable EM properties of the unified adjust tool; the order is
// the order of the property dropdown of the tool (persisted in prefs)
enum EmAdjustProperty
{
        EMADJP_CONDUCT = 0, // conductivity 0..1 (conductors)
        EMADJP_PERM,        // permeability EM_PERM_MIN..EM_PERM_MAX (ferromagnets)
        EMADJP_J,           // current -EM_JZEXT_MAX..+EM_JZEXT_MAX (current sources)
        EMADJP_MEDIUM,      // dielectric constant 0..EM_MEDIUM_MAX (dielectrics)
        EMADJP_MAG_DIR,     // magnetization direction 0..1 (fraction of 2pi)
        EMADJP_MAG_STR,     // magnetization strength 0..2
        EMADJP_COUNT,
};

// application modes of the unified adjust tool
enum EmAdjustApply
{
        EMADJA_SET = 0, // set the value under the brush to the target
        EMADJA_ADD,     // add the target to the value under the brush, every 0.2s
        EMADJA_SUB,     // subtract the target from the value under the brush, every 0.2s
};

// TM-mode electrodynamics (EM wave) simulation, a port of Paul Falstad's
// EMWave2 applet (http://www.falstad.com/emwave/). The field lives on a coarse
// grid covering the whole TPT canvas; the material properties of each cell are
// derived every frame from the TPT particles it covers.
//
// The current system has been rewritten around the EMWave2 mechanism: the only
// currents are the ones the applet understands (external sources driving the
// wave equation and induced currents inside conductors), carried exclusively by
// the dedicated EM / real zone elements. Original TPT elements are NOT mapped
// onto the field anymore and keep their vanilla behaviour.
//
// The wave is integrated in fixed sub-steps of EM_TADD_SUB (the exact applet
// timestep), {1,2,4} of them per frame depending on the speed setting, which
// keeps the CFL bound 2/tadd^2 = 32 on the permeability contrast satisfied for
// every speed setting and conserves energy in the linear regime.
class EMField
{
public:
        Simulation & sim;

        bool enabled = false;
        int cellSize = EM_CELL_SIZE_DEFAULT; // EM cell edge in domain pixels
        // regionScale semantics (0.5x fix): the simulated domain measures
        // regionScale * (XRES x YRES) DOMAIN pixels. regionScale > 1 extends the
        // domain beyond the canvas and the canvas shows the central 1:1 window;
        // regionScale = 1 shows the domain 1:1; regionScale = 0.5 makes the
        // domain HALF the canvas in each dimension but the renderer MAGNIFIES it
        // by renderScale = 2 so the field still covers the FULL screen - the
        // space is half, not the display (it is never a half-screen box with
        // dead borders anymore). Waves cross the visible canvas renderScale
        // times faster and their on-screen wavelength is renderScale times
        // longer, which is exactly what "the simulated space is smaller" looks
        // like; the cell count still drops 4x, so it stays the cheap option.
        float regionScale = EM_REGION_SCALE_DEFAULT; // domain multiplier (0.5 = half the canvas extent per axis, rendered zoomed to full screen, 1 = 1:1, 2/4/8 = extends beyond)
        int boundaryMode = EMBND_DEFAULT;    // one of EmBoundaryMode
        int viewMode = EMVIEW_DEFAULT;       // one of EmViewMode, used by the renderer
        int sourceMode = EMSRC_DEFAULT;      // one of EmSourceMode
        float frequency = 10;                // source frequency parameter (applet forceBar, 1..40)
        float aux = 1;                       // phase difference / 2nd frequency (applet auxBar, 1..40)
        int brightness = 100;                // display gain (applet brightnessBar, 1..500)
        int lineDensity = 50;                // B field line density (applet lineDensityBar, 10..100)
        int speed = 1;                       // 0 = half, 1 = normal, 2 = double speed

        double t = 0;
        // effective wave timestep of the current grid: taddEff * cellSize stays
        // constant (EM_TADD_SUB * EM_CELL_SIZE_DEFAULT) so the wave keeps the same
        // pixel speed and pixel wavelength at every grid resolution; the clock t
        // always advances EM_TADD_SUB per sub-step (task 7)
        float taddEff = EM_TADD_SUB;

        struct Cell
        {
                // material properties, re-derived from TPT particles every frame
                float perm = 1;          // relative permeability (1 = vacuum, >1 ferromagnet, <1 diamagnet)
                float conductivity = 0;  // 0 = vacuum .. 1 = perfect conductor
                float mx = 0, my = 0;    // permanent magnetization
                int medium = 0;          // dielectric strength, 0 .. EM_MEDIUM_MAX
                float epos = 0;          // electron displacement of resonant media
                bool resonant = false;
                bool boundary = false;   // cell sits at a material boundary
                bool gray = false;       // cell carries any material
                // external current injected this frame (EM elements, tools)
                double jzext = 0;
                // induced current inside conductors (applet jz semantics); only
                // meaningful while the cell carries a conductor or a resonant medium
                double jz = 0;
                // true when a real zone conductor that actually resists (everything
                // but superconductors below Tc) sits in this cell; only those heat up
                bool heatable = false;
                // true when a ferromagnetic or superconducting (below Tc) powder sits
                // in this cell; those feel magnetic pressure from the field
                bool magpowder = false;
                // true when a solid magnetic material (FE/PGRF/EMFM/EMDM, or a
                // superconductor below Tc) sits in this cell; those feel the same
                // magnetic pressure much more weakly than powders do
                bool magsolid = false;
                // manual overrides applied by the unified adjust tool
                uint8_t ovMask = 0;
                float ovConduct = 0, ovPerm = 0, ovJz = 0, ovMedium = 0, ovMag = 0;
                float ovDir = 0; // magnetization direction, angle = value * 2pi
                // wave state
                double az = 1e-10, dazdt = 1e-10;
                double damp = 1;
                // last rendered colour, reused when drawing field lines
                uint32_t col = 0;
        };

        int gw = 0, gh = 0; // grid size in EM cells
        // geometry of the boundary padding: the visible canvas maps to the grid
        // rectangle [padL, padL+visW) x [padT, padT+visH). With regionScale = 1 the
        // visible area equals the simulation domain (apart from boundary padding).
        // With regionScale > 1 the simulation domain is regionScale * (XRES x YRES)
        // pixels, so visW > renderW: the cells outside the visible canvas are
        // simulated normally but never rendered, letting waves travel through the
        // invisible padding before reaching the actual boundary. With
        // regionScale = 0.5 the domain is HALF the canvas extent per axis and the
        // renderer magnifies it by renderScale = 2 to cover the full canvas
        // (0.5x fullscreen fix, see the regionScale comment above).
        // CLOSED uses no padding (hard wall at domain edge); OPEN pads
        // with an invisible split-field PML band; PERIODIC pads by a one cell
        // ghost ring that is refreshed from the opposite edge each sub-step.
        int padL = 0, padT = 0;
        int visW = 0, visH = 0; // total simulated grid extent in EM cells (may exceed the visible window)
        // visible canvas window in cells and its zoom factor:
        // renderScale = 1 for regionScale >= 1 (1:1 rendering), 2 for 0.5x.
        // One cell spans cellSize * renderScale SCREEN pixels.
        int renderScale = 1;
        int renderW = 0, renderH = 0; // visible canvas in cells (XRES / (cellSize*renderScale))
        int renderOffX = 0, renderOffY = 0; // cell-index offset: grid cell (renderOffX + i, renderOffY + j) renders at screen pixel (i*cellSize*renderScale, j*cellSize*renderScale)
        std::vector<Cell> cells;

        // dirty flag set by NotifyCellChanged, consumed at the start of Update()
        bool boundariesDirty = false;

        // one entry per EMW particle found during the last sync
        struct EmwSource
        {
                int gi;
                float freq;
                int waveform;
        };
        std::vector<EmwSource> emwSources;

        // applet source bookkeeping
        struct Source
        {
                int x = 0, y = 0;
                double v = 0;
        };
        Source sources[4];
        int sourceCount = 0;
        int sourceFreqCount = 1;
        bool sourcePlane = false;
        int sourceWaveform = EM_SWF_SIN;
        int auxFunction = EM_AUX_NONE;
        float sourceMult = 1;
        float forceBarValue = 10;
        double forceTimeZero = 0;

        int filterCount = 0;
        int margin = EM_MARGIN_AT_4; // absorbing edge width in cells

        // --- split-field PML state (task 8, OPEN only) ------------------------
        // The band replaces the wave equation with a damped split system:
        //   wx = wx*bX + Gx(az) ; wy = wy*bY + Gy(az)   (velocity split)
        //   ux = ux*bX + wx*tadd^2 ; uy = uy*bY + wy*tadd^2 (damped split
        //   displacement - damping ux with the SAME profile is what makes the
        //   layer perfectly matched)
        //   az = ux + uy ; dazdt = wx + wy
        // Only the y components are stored: ux = az - uy, wx = dazdt - wy.
        // pmlBX[i] / pmlBY[j] are the per-sub-step damping factors exp(-sigma),
        // graded quartically from the interface (1 at the inner edge) to the
        // profile peak at the outer edge; full-grid float arrays (8 B/cell).
        std::vector<float> pmlUy, pmlWy;
        std::vector<float> pmlBX, pmlBY;

        // PERF (task 8 v3): PML activity gate. ScanPmlActivity() runs once per
        // frame and measures the max wave state over the band plus an 8-cell
        // interior guard strip; while everything is below EM_PML_QUIET the two
        // band passes are skipped for the whole frame. Skipping is numerically
        // exact for a quiet band (damping a below-threshold state keeps it
        // below the threshold), and the guard strip is deeper than the max
        // per-frame wave travel (4 cells), so the layer always wakes up before
        // a wave can reach the interface. This is what makes the OPEN mode as
        // cheap as CLOSED in scenes without waves near the screen edges.
        bool pmlQuiet = false;
        void ScanPmlActivity();

        // debug counters, only meaningful with EMFIELD_DEBUG; the release build
        // keeps the memory (a few ints) but never reads them
        long long fieldClampHits = 0;
        long long jzClampHits = 0;

        explicit EMField(Simulation & sim);

        void SetCellSize(int newCellSize); // reallocates and clears the grid
        void SetRegionScale(float newScale); // reallocates and clears the grid
        void SetBoundaryMode(int newMode); // reallocates and clears the grid for the new boundary geometry
        void ApplyGridGeometry();          // derive gw/gh/pad/vis from cellSize + regionScale + boundaryMode and reallocate
        void Clear();                      // resets the wave state, keeps materials (applet doClear)
        void ClearAll();                   // resets the wave state and all tool overrides (applet doClearAll); used when the simulation itself is reset
        void ClearOverrides();             // removes all tool overrides
        void ClearCellOverrides(int gi);   // removes the tool overrides of one cell (used by the erase tools)
        void VacuumCell(int gi);           // applet MODE_CLEAR: remove all EM material and wave state from one cell
        void Update();                     // one simulation frame
        // NotifyCellChanged(): batched hint that one or more EM cells had their
        // material properties changed by an external tool (EMADJ); CalcBoundaries
        // is then run once at the start of the next Update() instead of after
        // every single brush dab, which is what made the EMADJ tool laggy.
        void NotifyCellChanged();
        // RunCalcBoundariesNow(): forces an immediate CalcBoundaries, used by the
        // EMADJ tool when the user lifts the brush so the next frame sees the
        // correct boundaries without waiting for an extra tick.
        void RunCalcBoundariesNow();

        void SyncMaterials();
        void CalcBoundaries();
        void SetDamping();
        void PmlStepA(); // split-field PML velocity update (band cells, task 8)
        void PmlStepB(); // split-field PML position update (band cells, task 8)
        void RefreshGhostRing(); // PERIODIC boundary: copy the opposite edge into the ghost ring
        void SetupSources();
        void DoSources(double tadd, bool clear);
        void FilterGrid();
        // once per frame, after the wave sub-steps: magnetic pressure on
        // ferromagnetic / diamagnetic materials (incl. superconductors below Tc),
        // the motor effect (j x B on current-carrying conductors) and Joule
        // heating of resistive real-zone conductors; the force is normalised by
        // the sub-step count so the per-frame impulse is speed-setting independent
        void ApplyFieldForces(int substeps);
        void InteropParticles(); // once per frame: EMAN antenna -> vanilla spark

        double GetMagX(int gi) const;
        double GetMagY(int gi) const;

        int CellIndex(int px, int py) const; // pixel coordinates -> cell index (clamped into the domain)
        // is this canvas pixel inside the simulated domain? with regionScale < 1
        // parts of the visible canvas simply have no EM field
        bool PixelInDomain(int px, int py) const;

        static int CellTypeOf(const Cell &oe); // applet OscElement::getType()

        // helpers used by the EM tools
        bool ApplyAdjust(int gi, float strength); // legacy combined adjust, picks the applet mode from the cell type
        bool ApplyMagDir(int gi, float strength); // applet MODE_ADJ_MAG_DIR equivalent
        bool ApplyAdjustMode(int mode, int gi, float strength); // one specific applet MODE_ADJ_* mode

        // unified adjust tool: one EMADJP_* property, one EMADJA_* apply mode and
        // the target value; ADD/SUB accumulate onto the effective value of the
        // cell, SET assigns the target; returns true when anything was changed
        bool ApplyEMProperty(int property, int applyMode, int gi, float value);

        // effective (override-aware) value of one cell property, used by ADD/SUB
        float EffectiveProperty(int property, int gi) const;

        // max real particle speed in px per frame, derived from the field
        // propagation speed; taddEff/2 cells per sub-step * cellSize px *
        // substeps = px/frame, which is cellSize-independent by design (task 7);
        // at regionScale < 1 the field is rendered magnified by renderScale, so
        // the on-screen speed scales with it (0.5x: waves cross the visible
        // canvas twice as fast - the space is half)
        float MaxParticleSpeed() const
        {
                int ss = EM_SUBSTEPS[speed < 0 || speed > 4 ? 1 : speed];
                return taddEff * 0.5f * cellSize * float(ss) * float(renderScale);
        }
        // CFL-safe permeability clamp for the current grid resolution:
        // the contrast ratio against the weakest perm must stay <= 2/taddEff^2
        float PermMax() const
        {
                return std::min(EM_PERM_MAX, EM_PERM_MIN * EmPermRatioMax(taddEff));
        }
};
