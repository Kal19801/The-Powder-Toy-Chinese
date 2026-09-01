#include "simulation/ElementCommon.h"

// EMMG: permanent magnet, port of the MODE_M_DOWN / MODE_M_UP / MODE_M_LEFT /
// MODE_M_RIGHT paint modes of Paul Falstad's EMWave2 applet (mx/my = +-1).
// The magnetization direction is chosen with .ctype: 0 = down, 1 = up,
// 2 = left, 3 = right (the four applet tools are provided as the EMMGD /
// EMMGU / EMMGL / EMMGR drawing tools). The field exerts magnetic pressure on
// ferromagnetic particles around it and can be adjusted with the EMADJ and
// EMAJS tools. The element itself is static; the field does the physics.

void Element::Element_EMMG()
{
        Identifier = "DEFAULT_PT_EMMG";
        Name = "EMMG";
        Colour = 0xC060E0_rgb;
        MenuVisible = 0; // replaced by the four EMMGD/EMMGU/EMMGL/EMMGR elements; kept working for old saves
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
        Description = ByteString("永磁体,产生磁场(.ctype定方向:0下1上2左3右,EMMGD/EMMGU/EMMGL/EMMGR工具直接画)").FromUtf8();

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
