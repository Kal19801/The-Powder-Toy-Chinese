#pragma once
#include "common/String.h"

struct menu_section
{
        String::value_type icon;
        String name;
        int itemcount;
        int doshow;
};

constexpr int SC_WALL      =  0;
constexpr int SC_ELEC      =  1;
// TM-mode electromagnetic field elements (port of the EMWave2 applet painting modes)
constexpr int SC_EM        =  2;
// real-physics zone: real charged particles, real materials and their powders
constexpr int SC_REAL      =  3;
constexpr int SC_POWERED   =  4;
constexpr int SC_SENSOR    =  5;
constexpr int SC_FORCE     =  6;
constexpr int SC_EXPLOSIVE =  7;
constexpr int SC_GAS       =  8;
constexpr int SC_LIQUID    =  9;
constexpr int SC_POWDERS   = 10;
constexpr int SC_SOLIDS    = 11;
constexpr int SC_NUCLEAR   = 12;
constexpr int SC_SPECIAL   = 13;
constexpr int SC_LIFE      = 14;
constexpr int SC_TOOL      = 15;
constexpr int SC_FAVORITES = 16;
constexpr int SC_DECO      = 17;
