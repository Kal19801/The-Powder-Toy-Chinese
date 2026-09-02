#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// EMJC: EM current injector block. A powered (life > 0) block that injects a
// continuous current into its EM cell, just like the EMJP/EMJN sources but
// user-tunable and gated by the vanilla SPRK system.
//
// .tmp sets the signed current magnitude as a percentage of EM_JZEXT_MAX
// (-100 = full negative, +100 = full positive, 0 = default full positive
// when sparked so a freshly placed block does something visible). The
// vanilla SPRK system drives .life, so a vanilla circuit controls the
// current just like it controls any other powered block (HSWC, PCLN, ...).
//
// This is the "current injector in EM zone conductors" requested by task 10.
// It complements the EMADJ tool: EMADJ adjusts the cell material properties
// (perm/conductivity/etc.) of the underlying particles, while EMJC injects
// a real external current into the field at a chosen point. Use EMJC to
// drive a circuit of EM-zone conductors from a vanilla switch.

static int update(UPDATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS);

void Element::Element_EMJC()
{
        Identifier = "DEFAULT_PT_EMJC";
        Name = "EMJC";
        Colour = 0xFFB060_rgb;
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
        Description = ByteString("电流注入块,受SPRK驱动(.tmp定电流大小及方向,-100~+100,0=默认正向全幅),向所在电磁元格注入持续电流").FromUtf8();

        // PROP_CONDUCTS so vanilla SPRK lands on it and lights up its .life
        // PROP_LIFE_DEC so the gate naturally turns off when the spark ends
        Properties = TYPE_SOLID | PROP_CONDUCTS | PROP_LIFE_DEC;

        LowPressure = IPL;
        LowPressureTransition = NT;
        HighPressure = IPH;
        HighPressureTransition = NT;
        LowTemperature = ITL;
        LowTemperatureTransition = NT;
        HighTemperature = ITH;
        HighTemperatureTransition = NT;

        Update = &update;
        Create = &create;
}

static void create(ELEMENT_CREATE_FUNC_ARGS)
{
        if (sim->emfOwner)
        {
                sim->emfOwner->enabled = true;
        }
}

static int update(UPDATE_FUNC_ARGS)
{
        // The actual injection happens in EMField::SyncMaterials (it reads
        // .life and .tmp directly). Here we only make sure the field is on
        // while an EMJC exists in the simulation, so a freshly placed block
        // doesn't sit idle until the user toggles the global enable.
        if (sim->emfOwner && !sim->emfOwner->enabled)
        {
                sim->emfOwner->enabled = true;
        }
        return 0;
}
