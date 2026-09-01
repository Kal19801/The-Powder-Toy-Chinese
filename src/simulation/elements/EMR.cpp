#include "simulation/ElementCommon.h"

// EMR: resonant EM medium, ported from the resonant medium of Paul Falstad's
// EMWave2 applet. A passive medium whose bound electrons behave like little
// springs: it stores and re-radiates energy, slowly absorbing waves that pass
// through it. The element itself is static; the field does the physics.

void Element::Element_EMR()
{
	Identifier = "DEFAULT_PT_EMR";
	Name = "EMR";
	Colour = 0xFFBF80_rgb;
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
	HotAir = 0.000f	* CFDS;
	Falldown = 0;

	Flammable = 0;
	Explosive = 0;
	Meltable = 0;
	Hardness = 1;

	Weight = 100;

	HeatConduct = 251;
	Description = ByteString("谐振介质,吸收并缓慢辐射电磁波,可用于制作电磁谐振腔").FromUtf8();

	Properties = TYPE_SOLID;

	LowPressure = IPL;
	LowPressureTransition = NT;
	HighPressure = IPH;
	HighPressureTransition = NT;
	LowTemperature = ITL;
	LowTemperatureTransition = NT;
	HighTemperature = ITH;
	HighTemperatureTransition = NT;
}
