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

void Element::Element_TI()
{
        Identifier = "DEFAULT_PT_TI";
        Name = "TI";
        Colour = 0xC0C0D0_rgb;
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
        Meltable = 2;
        Hardness = 40;

        Weight = 100;

        HeatConduct = 100;
        Description = ByteString("钛,轻质耐腐蚀,弱导电,熔点极高").FromUtf8();

        Properties = TYPE_SOLID;

        LowPressure = IPL;
        LowPressureTransition = NT;
        HighPressure = IPH;
        HighPressureTransition = NT;
        LowTemperature = ITL;
        LowTemperatureTransition = NT;
        HighTemperature = ITH;
        HighTemperatureTransition = NT;

        HighTemperature = 1941.0f;
        HighTemperatureTransition = PT_LAVA; //@ TI -> LAVA(TI)

        Create = &create;
}
