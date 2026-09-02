#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// Real zone particle; the current carrier of the rewritten current system.
// Its motion is driven by the EM field inside EMField::InteractParticles()
// and it radiates by injecting current into the cells it crosses, the
// EMWave2 jz mechanism extended from static sources to moving charges.
//
// Real magnetic monopole: .ctype selects the pole (0 = N, 1 = S).
// Task 3: a monopole carries a true STATIC radial magnetic field
// B = g*EM_MONO_FIELD*r_hat/r (strength 1 at one cell distance, the same
// order as a painted EMMG magnet's near field), computed analytically in
// EMField::ComputeStaticB() and superposed on the dynamic field. It is kept
// out of the wave equation on purpose: Maxwell is linear, so a static field
// must not scatter waves, and driving it through the wave equation (the old
// jmext approach) would pump energy into the simulation forever. Moving
// monopoles additionally deposit jmext along their path, so a moving
// monopole radiates waves, exactly as a moving magnetic charge should.
// Velocity is strictly clamped to the field propagation speed (CFL), so no
// current is ever deposited outside the light cone of the field.

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
        Description = ByteString("真实磁单极子,.ctype选极性(0=N极,1=S极),受B场直接驱动,自带径向静态磁场(近场强度与磁铁相同),与磁铁一样吸引铁屑").FromUtf8();

        Properties = TYPE_ENERGY;

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

// No update() nudge: unlike the electric charges, a monopole at rest in a
// field-free region must STAY at rest, otherwise it would keep radiating
// (a moving magnetic charge radiates) and its static field would wander.
// A freshly placed monopole still moves off, because create() gives it a
// random initial velocity like the other real particles.
