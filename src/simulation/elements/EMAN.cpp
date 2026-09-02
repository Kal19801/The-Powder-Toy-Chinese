#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// EMAN: EM antenna. A coupler that bridges the EM field simulation to the
// vanilla TPT circuit system in the EM -> vanilla direction.
//
// While the particle exists, the EM field simulation scans a 3x3 neighbourhood
// of EM cells around the antenna each frame (in EMField::InteropParticles).
// If any of those cells carries |E|, |B| or |j| above a small excitation
// threshold (i.e. the EM zone next to the antenna is actively being driven,
// not just sitting at the 1e-10 vacuum floor), the antenna fires a vanilla
// SPRK on every adjacent vanilla conductor (PROP_CONDUCTS, not already
// sparked). The threshold scales with cellSize so a finer grid still triggers.
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
        Description = ByteString("电磁天线,把旁边电磁区的激发情况转换成SPRK脉冲,激发相邻的原版导体").FromUtf8();

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
}
