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
}

void Element::Element_EMEC()
{
        Identifier = "DEFAULT_PT_EMEC";
        Name = "EMEC";
        Colour = 0xB0B090_rgb;
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
        Description = ByteString("良导体,电磁场电导率=0.9,吸收大部分电磁波(EMWave2 良导体)").FromUtf8();

        Properties = TYPE_SOLID | PROP_HOT_GLOW;

        LowPressure = IPL;
        LowPressureTransition = NT;
        HighPressure = IPH;
        HighPressureTransition = NT;
        LowTemperature = ITL;
        LowTemperatureTransition = NT;
        // task 4: high melting point + red-hot glow approaching it, exactly
        // like vanilla FE (1811 K) / TTAN (1941 K); melts back to itself
        HighTemperature = 2843.0f;
        HighTemperatureTransition = PT_LAVA; //@ EMEC -> LAVA(EMEC)

        Create = &create;
}
