#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// Static field-material element; the TM-mode EM field simulation (EMField)
// derives its cell properties from this element every frame. The element
// itself has no particle physics - the field does the work, exactly like the
// corresponding painting mode of Paul Falstad's EMWave2 applet.

static void create(ELEMENT_CREATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS)
{
        // placing any electromagnetic element enables the field simulation,
        // otherwise freshly drawn material would silently do nothing
        if (sim->emfOwner)
        {
                sim->emfOwner->enabled = true;
        }
        sim->parts[i].ctype = 0;
}

void Element::Element_EMMGD()
{
        Identifier = "DEFAULT_PT_EMMGD";
        Name = "EMMGD";
        Colour = 0xB050D0_rgb;
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
        Description = ByteString("磁铁(下),产生向下的永磁化(EMWave2 Magnet Down)").FromUtf8();

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
