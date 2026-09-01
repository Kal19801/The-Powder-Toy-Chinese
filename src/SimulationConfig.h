#pragma once
#include <cstdint>
#include <common/Vec2.h>

constexpr int MENUSIZE = 40;
constexpr int BARSIZE  = 17;

constexpr float M_GRAV = 6.67300e-1f;

//CELL, the size of the pressure, gravity, and wall maps. Larger than 1 to prevent extreme lag
constexpr int CELL = 4;
constexpr Vec2<int> CELLS = Vec2(153, 96);
constexpr Vec2<int> RES = CELLS * CELL;

constexpr int XCELLS = CELLS.X;
constexpr int YCELLS = CELLS.Y;
constexpr int NCELL  = XCELLS * YCELLS;
constexpr int XRES   = RES.X;
constexpr int YRES   = RES.Y;
constexpr int NPART  = XRES * YRES;

constexpr int XCNTR = XRES / 2;
constexpr int YCNTR = YRES / 2;

constexpr Vec2<int> WINDOW = RES + Vec2(BARSIZE, MENUSIZE);

constexpr int WINDOWW = WINDOW.X;
constexpr int WINDOWH = WINDOW.Y;

constexpr int MAXSIGNS = 16;

constexpr int   ISTP            = CELL / 2;
constexpr float CFDS            = 4.0f / CELL;
constexpr float MAX_VELOCITY = 1e4f;

//Air constants
constexpr float AIR_TSTEPP = 0.3f;
constexpr float AIR_TSTEPV = 0.4f;
constexpr float AIR_VADV   = 0.3f;
constexpr float AIR_VLOSS  = 0.999f;
constexpr float AIR_PLOSS  = 0.9999f;

//EM field constants, ported from Paul Falstad's EMWave2 applet (TM electrodynamics)
//forceBar (source frequency) multiplier
constexpr float EM_FREQ_MULT   = .0233333f / 2;
//fudge factor used by the applet when drawing H/M fields
constexpr float EM_MHMULT      = 12.0f;
//maximum dielectric "medium" cell value and the refractive contribution of a fully saturated cell
constexpr int   EM_MEDIUM_MAX      = 191;
constexpr float EM_MEDIUM_MAX_INDEX = .5f;
//EM wave sub-step count per frame, indexed by the speed setting (0 = half, 1 = normal, 2 = double);
//see the rewritten current system constants below for the per-step timestep
//allowed EM cell edge sizes in pixels; the EM grid is then XRES/size x YRES/size cells
constexpr int EM_CELL_SIZES[] = { 2, 4, 8, 16 };
constexpr int EM_CELL_SIZE_COUNT = int(sizeof(EM_CELL_SIZES) / sizeof(EM_CELL_SIZES[0]));
constexpr int EM_CELL_SIZE_DEFAULT = 4;
//edge damping margin, in EM cells, scaled from the applet's fixed 20 cell margin at a 4px cell size
constexpr int EM_MARGIN_AT_4 = 20;
//scale of Joule heating applied to real conducting particles sitting in cells that carry induced current
constexpr float EM_JOULE_HEAT = 2.0f;

// --- rewritten current system (real zone), extending EMWave2 -------------------
// per-sub-step wave timestep, the exact value the applet integrates with
// (val = 5, tadd = val*.05 = 0.25); the port sub-steps {1,2,4}x this per frame
constexpr float EM_TADD_SUB = 0.25f;
// wave sub-steps per frame indexed by the speed setting (0 = 0.5x, 1 = 1x, 2 = 2x);
// 1x matches the applet, which integrates twice per frame at its default speed bar
constexpr int   EM_SUBSTEPS[3] = { 1, 2, 4 };
// stability bound of the leapfrog wave update: the per-cell acceleration scale
// c = perm_oe/perm_neighbour must satisfy c * |lambda_max| * tadd^2 <= 4 with
// |lambda_max| = 2 for the /4 scaled 5-point Laplacian, i.e. c <= 2/tadd^2 = 32
// at tadd = 0.25. perm above this diverges (the "ferromagnet self-excitation"
// of the old implementation; large perm amplified the Laplacian beyond CFL).
constexpr float EM_PERM_MAX  = 32.0f;
constexpr float EM_PERM_MIN  = 0.05f;
// external current clamp per cell; the field can only carry a bounded amount of
// energy per step, injecting more is unphysical and the fastest route to runaway
constexpr float EM_JZEXT_MAX = 2.0f;
// absolute clamp on the wave state, a pure safety net (the state should never
// get near this with the CFL bound respected); hit counter visible in debug builds
constexpr float EM_FIELD_CLAMP = 1e4f;
// field propagation speed: max group velocity of the leapfrog scheme is
// tadd/2 cells per sub-step (verified numerically over the whole Brillouin
// zone); the total per frame is EM_SUBSTEPS * EM_TADD_SUB / 2 cells, this is
// also the strict speed of light for real particles (nothing may outrun the
// field, otherwise currents appear outside their light cone and feedback
// diverges - the suspected cause of the old ferromagnet self-excitation)
constexpr float EM_CELLS_PER_SUBSTEP = EM_TADD_SUB * 0.5f;
// current injected per EM cell by a real charge crossing it, in units of
// charge * velocity; a moving point charge is a current density J = q v delta(x)
constexpr float EM_CHARGE_CURRENT = 0.06f;
// magnetic current injected by a moving monopole (drives dazdt directly,
// symmetric to how jz drives the electric wave equation)
constexpr float EM_MONO_CURRENT   = 0.06f;
// gradient-of-field-magnitude force on real charges (dielectrophoretic coupling
// of an in-plane charge to the TM field; polarity independent, like the
// field-line pressure the applet draws in its force view)
constexpr float EM_GRAD_FORCE     = 0.0006f;
// Coulomb-like pairwise force between real charges at short range (softened);
// electron/proton mass ratio = 1836 gives the proton its sluggish response
constexpr float EM_COULOMB        = 0.0045f;
constexpr float EM_PROTON_MASS    = 1836.0f;
// force of the in-plane B field (B = curl az) on a magnetic monopole, F = g B
constexpr float EM_MONOPOLE_FORCE = 0.05f;
// magnetic pressure on ferromagnetic / diamagnetic powders (iron filings vs
// pyrolytic graphite); sign chosen by permeability below/above 1
constexpr float EM_POWDER_FORCE   = 0.02f;
// pairwise Coulomb is O(n^2); above this many real particles only the field
// coupling is applied (still fully functional, just without charge-charge force)
constexpr int   EM_PAIRWISE_LIMIT = 400;
// superconductor critical temperature (YBCO-ish); above it the material
// quenches to a mere fair conductor and stops expelling the B field
constexpr float EM_SC_TC          = 93.0f;

//EM display modes, one entry of the applet's viewChooser each (order matters, persisted in prefs)
enum EmViewMode
{
        EMVIEW_OFF = 0,
        EMVIEW_E,
        EMVIEW_B,
        EMVIEW_B_LINES,
        EMVIEW_B_STRENGTH,
        EMVIEW_J,
        EMVIEW_E_B,
        EMVIEW_E_B_LINES,
        EMVIEW_E_B_J,
        EMVIEW_E_B_LINES_J,
        EMVIEW_H,
        EMVIEW_M,
        EMVIEW_TYPE,
        EMVIEW_A,
        EMVIEW_POYNTING,
        EMVIEW_ENERGY,
        EMVIEW_POYNTING_ENERGY,
        EMVIEW_FORCE,
        EMVIEW_EFF_CUR,
        EMVIEW_MAG_CHARGE,
        EMVIEW_CURL_E,
        EMVIEW_BX,
        EMVIEW_BY,
        EMVIEW_HX,
        EMVIEW_HY,
        EMVIEW_COUNT,
};

//global EM source layouts, port of the applet's sourceChooser (order matters, persisted in prefs)
enum EmSourceMode
{
        EMSRC_NONE = 0,
        EMSRC_1S1F,
        EMSRC_1S2F,
        EMSRC_2S1F,
        EMSRC_2S2F,
        EMSRC_3S1F,
        EMSRC_4S1F,
        EMSRC_1S1F_PACKET,
        EMSRC_1S1F_PLANE,
        EMSRC_1S2F_PLANE,
        EMSRC_2S1F_PLANE,
        EMSRC_2S2F_PLANE,
        EMSRC_1S1F_PLANE_PACKET,
        EMSRC_COUNT,
};

//source waveforms
constexpr int EM_SWF_SIN    = 0;
constexpr int EM_SWF_PACKET = 1;
//auxiliary slider functions
constexpr int EM_AUX_NONE  = 0;
constexpr int EM_AUX_PHASE = 1;
constexpr int EM_AUX_FREQ  = 2;

constexpr int EMVIEW_DEFAULT = EMVIEW_OFF;
constexpr int EMSRC_DEFAULT  = EMSRC_NONE;

constexpr int NGOL = 24;

enum DefaultBrushes
{
        BRUSH_CIRCLE,
        BRUSH_SQUARE,
        BRUSH_TRIANGLE,
        NUM_DEFAULTBRUSHES,
};

//Photon constants
constexpr int SURF_RANGE     = 10;
constexpr int NORMAL_MIN_EST =  3;
constexpr int NORMAL_INTERP  = 20;
constexpr int NORMAL_FRAC    = 16;

constexpr auto REFRACT = UINT32_C(0x80000000);

/* heavy flint glass, for awesome refraction/dispersion
   this way you can make roof prisms easily */
constexpr float GLASS_IOR  = 1.9f;
constexpr float GLASS_DISP = 0.07f;

constexpr float R_TEMP = 22;
