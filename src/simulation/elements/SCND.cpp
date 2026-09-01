#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// Real zone material; participates in the TM-mode EM field simulation through
// EMField::SyncMaterials(), which maps its particle type and temperature onto
// the EM cell material every frame.
static void create(ELEMENT_CREATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS)
{
        if (sim->emfOwner)
        {
                sim->emfOwner->enabled = true;
        }
}

void Element::Element_SCND()
{
        Identifier = "DEFAULT_PT_SCND";
        Name = "SCND";
        Colour = 0xE8F0FF_rgb;
        MenuVisible = 1;
        MenuSection = SC_REAL;
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
        Hardness = 10;

        Weight = 100;

        HeatConduct = 251;
        Description = ByteString("超导体,低于93K完全导电并排出磁场(迈斯纳效应),超温失超").FromUtf8();

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
