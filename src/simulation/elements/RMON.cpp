#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// Real zone particle; the current carrier of the rewritten current system.
// Its motion is driven by the EM field inside EMField::InteractParticles()
// and it radiates by injecting current into the cells it crosses, the
// EMWave2 jz mechanism extended from static sources to moving charges.
//
// Real magnetic monopole: .ctype selects the pole (0 = N, 1 = S).
// Task 3: a stationary monopole now also carries a static magnetic field of
// the same order as a painted EMMG magnet - the EMField deposits a constant
// jmext = +-EM_MONO_STATIC in the cell each frame, which is the dual of how
// EMJP deposits a constant jzext. Moving monopoles additionally deposit
// jmext along their path, so a moving monopole radiates waves.
// Velocity is strictly clamped to the field propagation speed (CFL), so no
// current is ever deposited outside the light cone of the field.

static int update(UPDATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS);

void Element::Element_RMON()
{
        Identifier = "DEFAULT_PT_RMON";
        Name = "RMON";
        Colour = 0x90FF90_rgb;
        MenuVisible = 1;
        MenuSection = SC_REAL;
        Enabled = 1;

        // motion parameters: identical to vanilla ELEC, matching the other
        // real charges - a freshly placed monopole actually moves off
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
        Description = ByteString("真实磁单极子,.ctype选极性(0=N极,1=S极),受磁场B直接驱动,静止时也带磁场(强度与磁铁相同)").FromUtf8();

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
        if (v > 0)
        {
                sim->parts[i].ctype = 1;
        }
        // give the monopole a random initial velocity so it moves off when
        // placed - matches the other real charges
        float a = sim->rng.between(0, 359) * std::numbers::pi_v<float> / 180.0f;
        sim->parts[i].vx = 2.0f * cosf(a);
        sim->parts[i].vy = 2.0f * sinf(a);
}

static int update(UPDATE_FUNC_ARGS)
{
        // Motion update function (task 2 / task 3): same nudge as RPRO/RELC
        // so a monopole that came to rest still gets picked up by the
        // framework's movement phase and can be pushed around by the B field.
        if (parts[i].vx == 0.0f && parts[i].vy == 0.0f)
        {
                float a = sim->rng.between(0, 359) * std::numbers::pi_v<float> / 180.0f;
                parts[i].vx = 0.5f * cosf(a);
                parts[i].vy = 0.5f * sinf(a);
        }
        return 0;
}
