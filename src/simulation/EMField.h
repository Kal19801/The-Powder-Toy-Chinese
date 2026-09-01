#pragma once
#include "SimulationConfig.h"
#include <cstdint>
#include <cstdio>
#include <vector>

class Simulation;
struct Particle;

// EM cell material override bits, set by the EM adjustment tools; the override only
// applies while the underlying particles still provide a matching material, just like
// the adjust modes of the original applet only affect cells of the matching type.
constexpr uint8_t EM_OV_CONDUCT = 1 << 0; // conductivity of conductors
constexpr uint8_t EM_OV_PERM    = 1 << 1; // permeability of ferromagnets
constexpr uint8_t EM_OV_JZ      = 1 << 2; // magnitude of external current sources
constexpr uint8_t EM_OV_MEDIUM  = 1 << 3; // dielectric constant of dielectrics
constexpr uint8_t EM_OV_MAGSTR  = 1 << 4; // strength of permanent magnets
constexpr uint8_t EM_OV_MAGDIR  = 1 << 5; // direction of permanent magnets (angle = value * 2pi)

// the six separate adjust modes of the applet (MODE_ADJ_*); used by the EMAJ* tools
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

// TM-mode electrodynamics (EM wave) simulation, a faithful port of Paul Falstad's
// EMWave2 applet (http://www.falstad.com/emwave/). The field lives on a coarse grid
// covering the whole TPT canvas; the material properties of each cell are derived
// every frame from the TPT particles it covers, so original elements interact with
// the field: metals reflect and absorb waves, glass refracts them, sparks and EM
// sources radiate, waves heat conductors, exert magnetic pressure on ferromagnetic
// particles and can induce sparks in antennas.
class EMField
{
public:
        Simulation & sim;

        bool enabled = false;
        int cellSize = EM_CELL_SIZE_DEFAULT; // EM cell edge in pixels
        int viewMode = EMVIEW_DEFAULT;       // one of EmViewMode, used by the renderer
        int sourceMode = EMSRC_DEFAULT;      // one of EmSourceMode
        float frequency = 10;                // source frequency parameter (applet forceBar, 1..40)
        float aux = 1;                       // phase difference / 2nd frequency (applet auxBar, 1..40)
        int brightness = 100;                // display gain (applet brightnessBar, 1..500)
        int lineDensity = 50;                // B field line density (applet lineDensityBar, 10..100)
        int speed = 1;                       // 0 = half, 1 = normal, 2 = double speed

        double t = 0;

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
                // external current injected this frame (sparks, electrons, EM sources, tools)
                double jzext = 0;
                // induced current inside conductors (applet jz semantics)
                double jz = 0;
                // manual overrides applied by tools
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
        std::vector<Cell> cells;

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

        explicit EMField(Simulation & sim);

        void SetCellSize(int newCellSize); // reallocates and clears the grid
        void Clear();                      // resets the wave state, keeps materials (applet doClear)
        void ClearAll();                   // resets the wave state and all tool overrides (applet doClearAll); used when the simulation itself is reset
        void ClearOverrides();             // removes all tool overrides
        void ClearCellOverrides(int gi);   // removes the tool overrides of one cell (used by the erase tools)
        void VacuumCell(int gi);           // applet MODE_CLEAR: remove all EM material and wave state from one cell
        void Update();                     // one simulation frame

        void SyncMaterials();
        void CalcBoundaries();
        void SetDamping();
        void SetupSources();
        void DoSources(double tadd, bool clear);
        void FilterGrid();
        void InteractParticles();

        double GetMagX(int gi) const;
        double GetMagY(int gi) const;

        int CellIndex(int px, int py) const; // pixel coordinates -> cell index (clamped)

        static int CellTypeOf(const Cell &oe); // applet OscElement::getType()

        // helpers used by the EM tools
        bool ApplyAdjust(int gi, float strength); // combined adjust: picks the applet mode from the cell type
        bool ApplyMagDir(int gi, float strength); // applet MODE_ADJ_MAG_DIR equivalent
        bool ApplyAdjustMode(int mode, int gi, float strength); // one specific applet MODE_ADJ_* mode
};
