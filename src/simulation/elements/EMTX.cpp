#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// EMTX: EM transmitter. The vanilla -> EM direction of the EM / vanilla
// circuit coupling - the inverse of EMAN.
//
// When the EMTX particle is sparked (life > 0), it injects a sinusoidal
// current pulse into its EM cell with frequency controlled by .ctype
// (1..40, cycles per applet time unit; values outside that range clamp to
// the global frequency) and amplitude controlled by .life (1..20; the
// amplitude is life/4 so a freshly sparked transmitter at life=4 emits at
// unit amplitude, and a long-lived spark driven by a BTRY emits stronger).
//
// The wave radiates out from the EMTX cell into the EM field exactly like
// an EMJP/EMJN source, except pulsed instead of continuous and user-tuned.
// Pair with an EMAN antenna on the other side of the field to relay the
// signal back into a vanilla circuit.
//
// While not sparked the EMTX is a passive vacuum-like cell - the field
// passes through it without disruption.

static int update(UPDATE_FUNC_ARGS);
static void create(ELEMENT_CREATE_FUNC_ARGS);

void Element::Element_EMTX()
{
        Identifier = "DEFAULT_PT_EMTX";
        Name = "EMTX";
        Colour = 0x80FFC0_rgb;
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
        Description = ByteString("电磁发射器,收到SPRK后向电磁场发射电磁波(.ctype=频率1~40,.life=增幅)").FromUtf8();

        // PROP_CONDUCTS so vanilla SPRK can land on EMTX and trigger it
        // PROP_LIFE_DEC so life counts down and the pulse is finite
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
        // default frequency = global (ctype = 0); user can override with the
        // PROP tool or by typing a value while placing
        if (v > 0)
        {
                sim->parts[i].ctype = std::min(v, 40);
        }
}

static int update(UPDATE_FUNC_ARGS)
{
        // Only fire while we are actively sparked (life > 0). The simulation
        // framework decrements life each frame (PROP_LIFE_DEC), so the pulse
        // naturally ends when the spark cycle does.
        if (parts[i].life <= 0)
        {
                return 0;
        }
        auto *emf = sim->emfOwner.get();
        if (!emf || !emf->enabled)
        {
                return 0;
        }
        int gi = emf->CellIndex(x, y);
        // Frequency: .ctype picks 1..40, 0 (default) means use the global.
        float freq = emf->frequency;
        if (parts[i].ctype >= 1 && parts[i].ctype <= 40)
        {
                freq = float(parts[i].ctype);
        }
        // Amplitude: .life drives the strength. life=4 (a freshly sparked
        // conductor) gives unit amplitude; higher life (e.g. from a long
        // BTRY-driven spark) gives proportionally more, clamped to the
        // field's per-cell current limit so we never inject runaway energy.
        float amp = std::clamp(float(parts[i].life) / 4.0f, 0.0f, EM_JZEXT_MAX);
        // Sinusoidal drive: same time base as the applet sources so the
        // emitted wave is in phase with the EMW sources and the EMJP/EMJN
        // current sources. The phase advances by freq * EM_FREQ_MULT per
        // simulation time unit; here we use emf->t directly so the pulse
        // stays phase-locked to the rest of the field.
        double phase = double(freq) * emf->t * EM_FREQ_MULT;
        double v = std::sin(phase) * double(amp);
        // Inject as an external current into our cell. This is exactly what
        // EMJP does (continuous +1) but pulsed and user-tuned. The clamp in
        // DepositRealCharges will keep the per-cell energy bounded.
        emf->cells[gi].jzext = v;
        return 0;
}
