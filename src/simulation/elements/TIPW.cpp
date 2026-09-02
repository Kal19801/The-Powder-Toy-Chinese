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

void Element::Element_TIPW()
{
        Identifier = "DEFAULT_PT_TIPW";
        Name = "TIPW";
        Colour = 0xA8A8BC_rgb;
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
        Meltable = 2;
        Hardness = 35;

        Weight = 55;

        HeatConduct = 90;
        Description = ByteString("钛粉,轻质耐腐蚀粉末,弱导电").FromUtf8();

        Properties = TYPE_PART | PROP_HOT_GLOW;

        LowPressure = IPL;
        LowPressureTransition = NT;
        HighPressure = IPH;
        HighPressureTransition = NT;
        LowTemperature = ITL;
        LowTemperatureTransition = NT;
        HighTemperature = ITH;
        HighTemperatureTransition = NT;

        HighTemperature = 1941.0f;
        HighTemperatureTransition = PT_LAVA; //@ TIPW -> LAVA(TIPW)

        Create = &create;
}
