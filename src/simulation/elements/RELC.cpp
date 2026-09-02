#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// Real zone particle; the current carrier of the rewritten current system.
// Its motion is driven by the EM field inside EMField::InteractParticles()
// and it radiates by injecting current into the cells it crosses, the
// EMWave2 jz mechanism extended from static sources to moving charges.
//
// Task 2: this element used to be called "real electron"; it is now presented
// to the user as the negative charge (负电荷). The physics is unchanged
// (charge -1, unit mass, very nimble - it accelerates easily because F/m is
// 1836x larger than for the positive charge), but the motion parameters are
// identical to vanilla ELEC (AirLoss=1, Loss=1, Collision=-0.99) and the
// create() function gives the particle a random initial velocity like vanilla
// ELEC. The Update function is the explicit motion update requested by the
// task (see RPRO.cpp for the full rationale).
// Velocity is strictly clamped to the field propagation speed (CFL), so no
// current is ever deposited outside the light cone of the field.

static int update(UPDATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS);

void Element::Element_RELC()
{
        Identifier = "DEFAULT_PT_RELC";
        Name = "负电荷"; // task 2: real electron -> negative charge
        Colour = 0x80FFFF_rgb;
        MenuVisible = 1;
        MenuSection = SC_REAL;
        Enabled = 1;

        // motion parameters: identical to vanilla ELEC
        Advection = 0.0f;
        AirDrag = 0.00f * CFDS;
        AirLoss = 1.00f;
        Loss = 1.00f;
        Collision = -0.99f;
        Gravity = 0.0f;
        Diffusion = 0.00f;
        HotAir = 0.000f * CFDS;
        Falldown = 0;

        Flammable = 0;
        Explosive = 0;
        Meltable = 0;
        Hardness = 0;

        Weight = -1;

        HeatConduct = 251;
        Description = ByteString("负电荷,带负电的载流子(质量极小,极易加速),受电磁场驱动并辐射电磁波,放出来会动").FromUtf8();

        Properties = TYPE_ENERGY;

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
        // give the charge a random initial velocity like vanilla ELEC
        float a = sim->rng.between(0, 359) * std::numbers::pi_v<float> / 180.0f;
        sim->parts[i].vx = 2.0f * cosf(a);
        sim->parts[i].vy = 2.0f * sinf(a);
}

static int update(UPDATE_FUNC_ARGS)
{
        // Motion update function (task 2): see RPRO.cpp for the rationale.
        // A charge that came to rest needs a small nudge so the framework's
        // movement phase doesn't skip it permanently.
        if (parts[i].vx == 0.0f && parts[i].vy == 0.0f)
        {
                float a = sim->rng.between(0, 359) * std::numbers::pi_v<float> / 180.0f;
                parts[i].vx = 0.5f * cosf(a);
                parts[i].vy = 0.5f * sinf(a);
        }
        return 0;
}
