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

void Element::Element_PGRF()
{
        Identifier = "DEFAULT_PT_PGRF";
        Name = "PGRF";
        Colour = 0x303030_rgb;
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
        Hardness = 15;

        Weight = 100;

        HeatConduct = 60;
        Description = ByteString("热解石墨,强抗磁性(磁导率0.1),可悬浮在强磁体上方").FromUtf8();

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
