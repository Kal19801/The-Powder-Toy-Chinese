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
//EM wave time step per frame, indexed by the speed setting (0 = half speed, 1 = normal, 2 = double)
constexpr float EM_TADD[3] = { 0.125f, 0.25f, 0.5f };
//allowed EM cell edge sizes in pixels; the EM grid is then XRES/size x YRES/size cells
constexpr int EM_CELL_SIZES[] = { 2, 4, 8, 16 };
constexpr int EM_CELL_SIZE_COUNT = int(sizeof(EM_CELL_SIZES) / sizeof(EM_CELL_SIZES[0]));
constexpr int EM_CELL_SIZE_DEFAULT = 4;
//edge damping margin, in EM cells, scaled from the applet's fixed 20 cell margin at a 4px cell size
constexpr int EM_MARGIN_AT_4 = 20;
//scale of Joule heating applied to conducting particles sitting in cells that carry induced current
constexpr float EM_JOULE_HEAT = 2.0f;
//|dazdt| (E field amplitude) above which conductors may be induced into sparking, plus the chance per frame
constexpr float EM_INDUCED_SPARK_THRESHOLD = 1.8f;
constexpr int   EM_INDUCED_SPARK_CHANCE = 8;
//scale of the magnetic pressure force applied to ferromagnetic particles
constexpr float EM_MAG_FORCE = 0.03f;

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
