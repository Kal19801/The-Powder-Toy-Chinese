#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// Powder form of a real zone material; falls like a normal TPT powder and
// carries the (slightly weakened) EM material of its compact form. Ferro-
// magnetic and diamagnetic powders additionally feel magnetic pressure from
// the EM field, which is how iron filings line up and pyrolytic graphite
// levitates over magnets.
static void create(ELEMENT_CREATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS)
{
        if (sim->emfOwner)
        {
                sim->emfOwner->enabled = true;
        }
}

void Element::Element_FEPW()
{
        Identifier = "DEFAULT_PT_FEPW";
        Name = "FEPW";
        Colour = 0x8C8C8C_rgb;
        MenuVisible = 1;
        MenuSection = SC_REAL;
        Enabled = 1;

        Advection = 0.4f;
        AirDrag = 0.04f * CFDS;
        AirLoss = 0.94f;
        Loss = 0.95f;
        Collision = -0.1f;
        Gravity = 0.3f;
        Diffusion = 0.00f;
        HotAir = 0.000f * CFDS;
        Falldown = 1;

        Flammable = 0;
        Explosive = 0;
        Meltable = 10;
        Hardness = 15;

        Weight = 85;

        HeatConduct = 200;
        Description = ByteString("铁粉,铁磁性粉末,会被磁场吸引排列,高温熔化").FromUtf8();

        Properties = TYPE_PART | PROP_HOT_GLOW;

        LowPressure = IPL;
        LowPressureTransition = NT;
        HighPressure = IPH;
        HighPressureTransition = NT;
        LowTemperature = ITL;
        LowTemperatureTransition = NT;
        HighTemperature = ITH;
        HighTemperatureTransition = NT;

        HighTemperature = 1811.0f;
        HighTemperatureTransition = PT_LAVA; //@ FEPW -> LAVA(FEPW)

        Create = &create;
}
