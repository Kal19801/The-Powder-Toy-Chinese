#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// EMAN: EM antenna. A coupler that bridges the EM field simulation to the
// vanilla TPT circuit system in the EM -> vanilla direction.
//
// Task 5: the EM-zone elements CONNECTED to the antenna are its antenna.
// Each frame (in EMField::InteropParticles) the antenna reads the excitation
// current |jz + jzext| - the current the wave drives inside those connected
// elements - of the 3x3 EM-cell neighbourhood around itself and compares it
// against its threshold. The threshold is .ctype in raw current units
// (0 = EMAN_THRESHOLD_DEFAULT); set it with the PROP tool or by drawing with
// a value. When the excitation exceeds the threshold, the antenna fires a
// vanilla SPRK on every adjacent vanilla conductor (PROP_CONDUCTS, not
// already sparked and not in spark cooldown).
//
// Pair with EMTX (the vanilla -> EM direction) for bidirectional coupling:
// a vanilla circuit drives EMTX, EMTX radiates an EM wave, the wave reaches
// EMAN, EMAN sparks the next vanilla circuit on the other side.
//
// The element itself is non-material - it doesn't add perm/conductivity/
// magnetisation to the cell. It is just a marker the EMField looks for.

static void create(ELEMENT_CREATE_FUNC_ARGS);

void Element::Element_EMAN()
{
        Identifier = "DEFAULT_PT_EMAN";
        Name = "EMAN";
        Colour = 0xFFC080_rgb;
        MenuVisible = 1;
        MenuSection = SC_EM;
        Enabled = 1;

        Advection = 0.0f;
        AirDrag = 0.00f * CFDS;
        AirLoss = 1.00f;
        Loss = 0.00f;
        Collision = 0.0f;
        Gravity = 0.0f;
        Diffusion = 0.00f;
        HotAir = 0.000f * CFDS;
        Falldown = 0;

        Flammable = 0;
        Explosive = 0;
        Meltable = 0;
        Hardness = 1;

        Weight = 100;

        HeatConduct = 251;
        Description = ByteString("电磁天线,把相连电磁区元素里的感应电流转换成SPRK,激发相邻原版导体;.ctype=电流阈值(0=默认0.05)").FromUtf8();

        Properties = TYPE_SOLID;

        LowPressure = IPL;
        LowPressureTransition = NT;
        HighPressure = IPH;
        HighPressureTransition = NT;
        LowTemperature = ITL;
        LowTemperatureTransition = NT;
        HighTemperature = ITH;
        HighTemperatureTransition = NT;

        Create = &create;
}

static void create(ELEMENT_CREATE_FUNC_ARGS)
{
        // placing any electromagnetic element enables the field simulation,
        // otherwise freshly drawn material would silently do nothing
        if (sim->emfOwner)
        {
                sim->emfOwner->enabled = true;
        }
        // the excitation threshold can be set while placing (draw value)
        if (v > 0)
        {
                sim->parts[i].ctype = std::min(v, int(EMAN_THRESHOLD_MAX));
        }
}
