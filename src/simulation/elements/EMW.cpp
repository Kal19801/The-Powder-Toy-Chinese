#include "simulation/ElementCommon.h"

// EMW: EM wave source, ported from the sources of Paul Falstad's EMWave2 applet.
// Radiates a TM-mode wave into the EM field at the frequency set in the settings;
// .tmp overrides that frequency (1..40), .ctype selects the waveform (0 = sine,
// 1 = wave packet). The element itself is static; the field does the radiating.

static void create(ELEMENT_CREATE_FUNC_ARGS);

void Element::Element_EMW()
{
        Identifier = "DEFAULT_PT_EMW";
        Name = "EMW";
        Colour = 0xC0C0FF_rgb;
        MenuVisible = 1;
        MenuSection = SC_ELEC;
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
        Description = ByteString("电磁波源,向电磁场辐射正弦波(.tmp覆盖频率1~40,.ctype=1为波包)").FromUtf8();

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
        // default to the global frequency (tmp = 0) unless a specific one is drawn
        if (v > 0)
        {
                sim->parts[i].tmp = (v < 40) ? v : 40;
        }
}
