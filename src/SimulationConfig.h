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
//allowed EM cell edge sizes in pixels; the EM grid is then XRES/size x YRES/size cells;
//size 1 gives a field on the exact particle grid (1 cell = 1 particle, full alignment)
constexpr int EM_CELL_SIZES[] = { 1, 2, 4, 8, 16 };
constexpr int EM_CELL_SIZE_COUNT = int(sizeof(EM_CELL_SIZES) / sizeof(EM_CELL_SIZES[0]));
// Default cell size: 2 pixels per cell. The "1x" region (regionScale = 1) covers
// exactly the visible canvas at this granularity, matching the previous "2 pixel"
// setting; smaller cells give finer resolution at greater CPU cost, larger cells
// are cheaper but coarser. The region/grid decoupling (regionScale) is independent.
constexpr int EM_CELL_SIZE_DEFAULT = 2;

// Region size multiplier (regionScale): scales the simulated EM domain
// independently of the cell size. regionScale = 1 means the simulation domain
// equals the visible canvas (XRES x YRES) and renders 1:1; 2/4/8 expand the
// domain beyond the visible canvas, so waves can travel through (and partially
// into) the invisible padding before reaching the boundary while the canvas
// shows the central 1:1 window; 0.5 SHRINKS the domain to half the canvas
// extent per axis (XRES/2 x YRES/2, a quarter of the cells and CPU cost) but
// the renderer MAGNIFIES it 2x so the field still covers the FULL screen -
// the space is half, not the display. On screen the wave crosses the canvas
// twice as fast with twice the wavelength, which is exactly what a physically
// smaller world looks like. This decouples "grid resolution" from "physical
// domain size".
constexpr float EM_REGION_SCALES[] = { 0.5f, 1.0f, 2.0f, 4.0f, 8.0f };
constexpr int EM_REGION_SCALE_COUNT = int(sizeof(EM_REGION_SCALES) / sizeof(EM_REGION_SCALES[0]));
constexpr float EM_REGION_SCALE_DEFAULT = 1.0f;
// nearest valid region scale (used when loading prefs); -1 when not found
inline int EmRegionScaleIndex(float regionScale)
{
        for (int i = 0; i < EM_REGION_SCALE_COUNT; ++i)
        {
                if (EM_REGION_SCALES[i] == regionScale)
                {
                        return i;
                }
        }
        return -1;
}

//edge damping margin, in EM cells, scaled from the applet's fixed 20 cell margin at a 4px cell size
constexpr int EM_MARGIN_AT_4 = 20;
//scale of Joule heating applied to real conducting particles sitting in cells that carry induced current
constexpr float EM_JOULE_HEAT = 2.0f;

// --- wave timestep / grid decoupling (task 7) ---------------------------------
// The scheme integrates az with a timestep that must shrink proportionally to
// the cell size, so the wave keeps the SAME speed in pixels and the SAME
// wavelength in pixels at every grid resolution: taddEff * cellSize = const.
// The simulation clock t always advances EM_TADD_SUB per sub-step (the
// applet time unit), so source frequencies keep the same spatial wavelength
// at every resolution as well. Traversal time of the simulated region now
// depends only on the region size, never on the grid resolution.
constexpr float EM_TADD_SUB = 0.25f;
inline float EmWaveTadd(int cellSize)
{
        return EM_TADD_SUB * 2.0f / float(cellSize);
}
// wave sub-steps per frame indexed by the speed setting (0 = 0.5x, 1 = 1x, 2 = 2x);
// 1x matches the applet, which integrates twice per frame at its default speed bar
constexpr int   EM_SUBSTEPS[5] = { 1, 2, 4, 8, 16 };
// stability bound of the leapfrog wave update: the per-cell permeability
// CONTRAST ratio c = perm_oe/perm_neighbour must satisfy c <= 2/taddEff^2
// (verified numerically: at tadd 0.25 ratio 48 is stable, 640 diverges;
// at tadd 0.5 ratio 16 diverges, 8 is stable - both match the formula).
constexpr float EmPermRatioMax(float taddEff)
{
        return 2.0f / (taddEff * taddEff);
}
// usable permeability range: diamagnets need perm < 1 (applet MODE_DIAMAG = 0.5)
// and the contrast against the strongest ferromagnet must stay inside the CFL
// ratio bound, so [0.5, 16] (ratio 32, the bound at the reference tadd 0.25).
// Finer grids (cellSize 1, taddEff 0.5) further clamp the range dynamically.
constexpr float EM_PERM_MAX  = 16.0f;
constexpr float EM_PERM_MIN  = 0.5f;
// external current clamp per cell; the field can only carry a bounded amount of
// energy per step, injecting more is unphysical and the fastest route to runaway
constexpr float EM_JZEXT_MAX = 2.0f;
// absolute clamp on the wave state, the last-resort safety net (unconditional,
// unlike the debug counters); fields should never get near this because the
// CFL bound above is enforced on the permeability range
constexpr float EM_FIELD_CLAMP = 1e4f;
// soft anti-divergence bleed: cells whose |az|/|dazdt| exceed this lose 2% per
// frame; only true feedback runaways ever reach it (legit f=2 standing waves
// peak around |az| ~ 900)
constexpr float EM_SOFT_LIMIT  = 1e3f;
constexpr float EM_SOFT_BLEED  = 0.98f;
// field propagation speed: group velocity of the leapfrog scheme is
// taddEff/2 cells per sub-step; in pixels per frame that is
// taddEff/2 * cellSize * substeps, which is cellSize-independent now
constexpr float EM_CELLS_PER_SUBSTEP = EM_TADD_SUB * 0.5f; // reference (cellSize 2)
// hard cap on gw*gh so extreme regionScale+cellSize combos cannot OOM
constexpr long long EM_MAX_CELLS = 6000000LL;

// --- absorbing boundary (task 8; task 11: CPML matched band) ------------------
// OPEN pads the grid with an invisible band OUTSIDE the simulated region. The
// band is a convolutional PML (CFS-CPML) fitted to the exact leapfrog stencil
// of the interior wave update: the band cells keep the UNDAMPED interior
// kinematics (az += dazdt*tadd^2, dazdt from the same 5-point stencil) and the
// absorption lives in a one-pole recursive filter psi on the az face
// differences, per axis:
//   psi' = b*psi + (g/2)*(FD + FD_prev);  dazdt += 0.25*(FD + psi) differences
// The face correction target (1-sx)/sx = -sigma/(alpha+sigma+i*omega) is
// EXACTLY a one-pole, and the trapezoidal (bilinear) recursion
//   b = (2-h)/(2+h), g = -2*sigma/(2+h), h = sigma+alpha
// realises it with no mapping error at any frequency (the bilinear transform
// is exact for a one-pole; an impulse-invariant recursion leaves an O(omega)
// phase error measured at 10x worse reflection). FD_prev, the previous
// sub-step's face difference, is reconstructed exactly from the current state
// (FD_prev = FD - tadd^2*dazdt difference), so no extra history arrays exist.
// Two earlier designs of this band were REJECTED on hard evidence from the
// headless harness and must not come back:
//  - split-field Berenger (ux/uy + wx/wy): stable, but DC/low-frequency pole
//    AND a growing trapped mode (pmlUy measured growing ~1.2%/substep,
//    continuously re-radiating into the canvas);
//  - Berenger pair (V=Vx+Vy + face states Hx/Hy, impedance-matched pairs):
//    analytically matched at all frequencies, but the two-sided staggered
//    auxiliary leapfrog closes an amplifier loop at the staggered Nyquist
//    (measured |g| = 1.064/substep at kx=ky=pi) - a structural instability
//    that no per-axis damping profile removes (eigen scan: |g| > 1.06 for any
//    b when one axis is undamped).
// The CPML recursion is passive (|b| < 1, steady-state gain <= 1): no closed
// amplifier loop is possible, and being a low-pass memory on REAL az
// differences it cannot decouple from the wave.
//
// EMField::ScanPmlActivity() gates the two band passes per sub-step: while
// the band AND an 8-cell interior guard strip are entirely quiet (max wave
// state < EM_PML_QUIET), the passes are skipped - numerically exact, since
// filtering a zero state keeps it zero. A wave trips the gate at least 4
// cells (max per-frame travel) before it can reach the interface, so the
// layer always wakes up in time. Scenes with no wave near the screen edge
// run OPEN exactly as cheap as CLOSED.
constexpr int   EM_PAD_PX        = 64;   // PML band thickness in pixels
constexpr int   EM_PAD_MIN_CELLS = 16;
constexpr int   EM_PAD_MAX_CELLS = 96;
// CPML profile: sigma rises quartically from the interface face to the outer
// edge; the peak (per sub-step) scales with taddEff so the attenuation per
// PIXEL is resolution-invariant. The one-way optical depth of a quartic
// profile is ~sigma_peak*width/5 Np; 0.30 with D=32 gives ~2 Np one-way, and
// the round-trip residual is set by the (tiny) discrete stencil reflection,
// not by optical depth.
constexpr float EM_PML_SIGMA     = 0.60f;
constexpr float EM_PML_POWER     = 4.0f;
// alpha floor at the interface (decaying to 0 at the outer edge): makes the
// stretched coordinate finite at DC. Standard CPML practice; also damps the
// weakly-excited interface poles of the discrete recursion.
constexpr float EM_PML_ALPHA     = 0.03f;
// kappa = 1 keeps the band wave speed identical to the interior. Raise only
// if grazing evanescent tails are measured to leak (they attenuate the
// evanescent spectrum without touching propagating waves).
constexpr float EM_PML_KAPPA     = 1.0f;
// wave state below which a cell counts as quiet for the PML activity gate
// (the vacuum seed is 1e-10, so a resting band is quiet by construction)
constexpr float EM_PML_QUIET     = 1e-9f;

// --- EM boundary conditions (设置 -> 电磁场边界条件) ---------------------------
// NOTE: the numeric values are persisted in prefs, so they must not be renumbered.
// The old EMBND_ABSORB = 1 (identical split-field PML as OPEN) was removed from
// the UI; value 1 is migrated to OPEN when prefs are loaded (see GameModel).
enum EmBoundaryMode
{
        EMBND_CLOSED = 0,   // 封闭: perfectly conducting walls at the screen edge, full reflection
        EMBND_OPEN   = 2,   // 开放: split-field PML, waves leave like in infinite space
        EMBND_PERIODIC = 3, // 循环: waves leaving one edge re-enter from the opposite edge (ghost ring)
        EMBND_COUNT,
};
constexpr int EMBND_DEFAULT = EMBND_OPEN;
// legacy pref value of the removed 吸收 option (behaviourally identical to OPEN)
constexpr int EMBND_LEGACY_ABSORB = 1;

// Motor effect (Lorentz force on a current): the bulk force density on a
// current-carrying conductor is F = j x B. For our 2D TM field with j = (0,0,jz)
// and B = (Bx,By,0), this gives F = (-jz*By, jz*Bx, 0) - a real, in-plane force
// on the conductor. Gain scales the force so a visible motion results at the
// magnitudes of j and B the simulation typically produces.
constexpr float EM_MOTOR_FORCE     = 0.020f;
// magnetic pressure on ferromagnetic / diamagnetic powders (iron filings vs
// pyrolytic graphite); sign chosen by permeability below/above 1. Tuned against
// powder gravity (~0.4 px/frame^2) so the pull of a magnet clearly wins close by
// and decays like the field gradient further away.
constexpr float EM_POWDER_FORCE   = 0.55f;
constexpr float EM_POWDER_NORM    = 0.30f;
// same mechanism, much milder, for solid magnetic materials (FE/PGRF/EMFM/EMDM):
// blocks slowly slide instead of jumping, powders are the strong responders
constexpr float EM_SOLID_FORCE    = 0.08f;
// superconductor critical temperature (YBCO-ish); above it the material
// quenches to a mere fair conductor and stops expelling the B field
constexpr float EM_SC_TC          = 93.0f;
// permeability a superconductor maps to below Tc (the applet's diamagnet
// value); it makes the superconductor repel the field (magnetic pressure)
// the same way a painted diamagnet does, on top of its perfect conductivity
constexpr float EM_SC_PERM        = 0.5f;
// conductivity a quenched superconductor falls back to (between the applet's
// fair .5 and a lossy wire; high enough to screen, low enough to dissipate)
constexpr float EM_SC_QUENCH_CONDUCT = 0.3f;
// EMTX transmitter defaults (task 6): burst = Hann-windowed 2 carrier cycles
constexpr int   EMTX_BURST_MIN_FRAMES = 6;
constexpr int   EMTX_BURST_MAX_FRAMES = 600;
// EMAN antenna default excitation threshold (task 5), |jz+jzext| field units;
// .ctype overrides (raw units)
constexpr float EMAN_THRESHOLD_DEFAULT = 0.02f;
constexpr float EMAN_THRESHOLD_MAX     = 10.0f;

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

constexpr int EMVIEW_DEFAULT = EMVIEW_E_B_LINES_J;
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
