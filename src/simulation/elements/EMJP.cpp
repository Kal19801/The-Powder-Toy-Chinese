#include "simulation/ElementCommon.h"

// EMJP: positive DC current source, port of the MODE_J_POS paint mode of
// Paul Falstad's EMWave2 applet (oe.jz = +1). The element injects a steady
// positive current into the EM field, which creates the magnetic field around
// the source and radiates waves; conductors near it carry induced currents
// exactly like current sources do in the applet. The element itself is
// static; the field does the rest.

void Element::Element_EMJP()
{
        Identifier = "DEFAULT_PT_EMJP";
        Name = "EMJP";
        Colour = 0xFF5050_rgb;
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
        Description = ByteString("电流源(+),向电磁场持续注入正向电流,产生磁场并辐射电磁波").FromUtf8();

        Properties = TYPE_SOLID;

        LowPressure = IPL;
        LowPressureTransition = NT;
        HighPressure = IPH;
        HighPressureTransition = NT;
        LowTemperature = ITL;
        LowTemperatureTransition = NT;
        HighTemperature = ITH;
        HighTemperatureTransition = NT;
}
