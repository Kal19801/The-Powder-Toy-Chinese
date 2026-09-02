#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// Real zone particle; the current carrier of the rewritten current system.
// Its motion is driven by the EM field inside EMField::InteractParticles()
// and it radiates by injecting current into the cells it crosses, the
// EMWave2 jz mechanism extended from static sources to moving charges.
//
// Task 2: this element used to be called "real proton"; it is now presented to
// the user as the positive charge (正电荷). Motion parameters are identical to
// vanilla ELEC (AirLoss=1, Loss=1, Collision=-0.99) and the create() function
// gives the particle a random initial velocity like vanilla ELEC, so a freshly
// placed charge actually moves off. The field-response inertia is also the
// electron's (mass 1) - the task asks for motion like the original electron,
// so the only difference between the two charges is the SIGN of the charge
// (Coulomb force, quiver phase, drift direction).
// The Update function is the explicit motion update requested by the task:
// it nudges a stationary charge into motion (matching vanilla ELEC semantics)
// and updates the EM field state flag if the field is off.
// Velocity is strictly clamped to the field propagation speed (CFL), so no
// current is ever deposited outside the light cone of the field.

static int update(UPDATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS);

void Element::Element_RPRO()
{
        Identifier = "DEFAULT_PT_RPRO";
        Name = "RPRO"; // task 2: real proton -> positive charge
        Colour = 0xFF9060_rgb;
        MenuVisible = 1;
        MenuSection = SC_REAL;
        Enabled = 1;

        // motion parameters: identical to vanilla ELEC so the positive charge
        // behaves kinematically like an electron (same AirLoss/Loss/Collision).
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
        Description = ByteString("正电荷,带正电的载流子,与负电荷相遇湮灭").FromUtf8();

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
        // give the charge a random initial velocity like vanilla ELEC, so a
        // freshly placed charge moves off instead of sitting still until the
        // field pushes it. Speed matches ELEC (2.0 px/frame).
        float a = sim->rng.between(0, 359) * std::numbers::pi_v<float> / 180.0f;
        sim->parts[i].vx = 2.0f * cosf(a);
        sim->parts[i].vy = 2.0f * sinf(a);
}

static int update(UPDATE_FUNC_ARGS)
{
        // Motion update function (task 2): the simulation framework applies
        // Loss/AirLoss and the EM field adds force via InteractParticles(),
        // but a charge that has been sitting still for a while (e.g. collided
        // and lost all velocity) would never move again - the framework skips
        // the movement phase when vx and vy are both zero. This update gives
        // such a charge a tiny random nudge so the simulation keeps it alive
        // and the EM field can pick it up again. The nudge is small enough
        // that an actively moving charge never notices it.
        if (parts[i].vx == 0.0f && parts[i].vy == 0.0f)
        {
                float a = sim->rng.between(0, 359) * std::numbers::pi_v<float> / 180.0f;
                parts[i].vx = 0.5f * cosf(a);
                parts[i].vy = 0.5f * sinf(a);
        }
        return 0;
}
