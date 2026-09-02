#include "simulation/ElementCommon.h"
#include "simulation/EMField.h"

// EMTX: EM wave transmitter (task 6). The vanilla -> EM direction of the
// EM / vanilla coupling, the inverse of EMAN.
//
// "收到SPRK后能发射一道根据ctype设置频率，life设置振幅的电磁波":
//   .ctype = carrier frequency 1..40 (0 = use the global frequency setting)
//   .life  = amplitude * 4 (0 = unit amplitude, same strength as an EMJP
//            source; 8 = the field's maximum per-cell current)
//   .tmp   = burst frames remaining (runtime state, do not set by hand)
//   .tmp2  = burst frames total (runtime state, do not set by hand)
//
// When an ADJACENT conductor is sparked, the transmitter emits ONE wave
// packet: a Hann-windowed burst of two carrier cycles whose peak amplitude
// is life/4. The packet radiates into the EM field exactly like an EMJP/EMJN
// source, except it is a finite, smooth burst instead of a continuous drive,
// so it can carry an on/off signal.
//
// Design notes (why this shape):
//  * EMTX is deliberately NOT a vanilla conductor (no PROP_CONDUCTS). A
//    conductor gets its .ctype and .life overwritten by the spark cycle,
//    which would destroy the frequency and amplitude settings on the very
//    first spark. Instead EMTX never becomes SPRK itself: it watches its
//    neighbourhood and "receives" the spark of the conductor it touches.
//  * The actual field injection lives in EMField::SyncMaterials(), which
//    runs BEFORE the wave sub-steps. An element update() writing jzext
//    would have its value erased by the next SyncMaterials before the wave
//    ever sees it - that was the fatal flaw of the previous implementation.

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
        Description = ByteString("电磁发射器,相邻导体收到SPRK时发射一道电磁波(.ctype=频率1~40,0=全局;.life=振幅*4,0=单位振幅)").FromUtf8();

        // deliberately NOT PROP_CONDUCTS: the spark cycle would clobber the
        // user's frequency (.ctype) and amplitude (.life) settings
        Properties = TYPE_SOLID;

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
        // frequency can be set while placing (draw value 1..40)
        if (v > 0)
        {
                sim->parts[i].ctype = std::min(v, 40);
        }
        // draw value 2.. selects the amplitude (life = value*4, capped)
        if (v > 1)
        {
                sim->parts[i].life = std::min(v * 4, int(EM_JZEXT_MAX * 4));
        }
}

static int update(UPDATE_FUNC_ARGS)
{
        // The burst state machine:
        //   tmp == 0: idle. Watch for an adjacent SPRK; on sight, start a burst
        //             of 2 carrier cycles (frames computed from the speed and
        //             frequency so the packet carries a fixed number of cycles).
        //   tmp  > 0: bursting; the field injects in EMField::SyncMaterials and
        //             counts tmp down once per frame.
        // A continuous spark (e.g. BTRY feeding the neighbour) re-triggers the
        // burst as soon as the previous one ends, giving a continuous carrier.
        if (parts[i].tmp > 0)
        {
                return 0;
        }
        bool sparked = false;
        static constexpr int dxs[4] = { -1, 1, 0, 0 };
        static constexpr int dys[4] = { 0, 0, -1, 1 };
        for (int k = 0; k < 4 && !sparked; k++)
        {
                int nx = int(x + 0.5f) + dxs[k];
                int ny = int(y + 0.5f) + dys[k];
                if (nx < 0 || ny < 0 || nx >= XRES || ny >= YRES)
                {
                        continue;
                }
                unsigned r = pmap[ny][nx];
                if (r && TYP(r) == PT_SPRK && parts[ID(r)].life > 0)
                {
                        sparked = true;
                }
        }
        if (!sparked)
        {
                return 0;
        }
        auto *emf = sim->emfOwner.get();
        if (!emf || !emf->enabled)
        {
                return 0;
        }
        float freq = emf->frequency;
        if (parts[i].ctype >= 1 && parts[i].ctype <= 40)
        {
                freq = float(parts[i].ctype);
        }
        // frames per carrier cycle: 2pi / (freq * EM_FREQ_MULT * tau_per_frame),
        // tau_per_frame = EM_TADD_SUB * substeps; speed- and resolution-proof
        int substeps = EM_SUBSTEPS[std::clamp(emf->speed, 0, 4)];
        float tauPerFrame = EM_TADD_SUB * float(substeps);
        int framesPerCycle = int(2.0f * std::numbers::pi_v<float> / (freq * EM_FREQ_MULT * tauPerFrame));
        int burst = 2 * framesPerCycle;
        burst = std::clamp(burst, EMTX_BURST_MIN_FRAMES, EMTX_BURST_MAX_FRAMES);
        parts[i].tmp = burst;
        parts[i].tmp2 = burst;
        return 0;
}
