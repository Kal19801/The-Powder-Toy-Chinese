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

void Element::Element_PGPW()
{
        Identifier = "DEFAULT_PT_PGPW";
        Name = "PGPW";
        Colour = 0x404040_rgb;
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
        Meltable = 0;
        Hardness = 10;

        Weight = 50;

        HeatConduct = 50;
        Description = ByteString("石墨粉,抗磁性粉末,会被磁场轻微排斥").FromUtf8();

        Properties = TYPE_PART;

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
