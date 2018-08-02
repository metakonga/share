#include "materialCalculator.h"

/*using namespace MaterialCalculator;*/

VEC3D MaterialCalculator::CubeInertia(double m, double w, double h, double d)
{
	double x, y, z;
	if (d == 0)
	{
		d = 1.0;
	}
	
	x = (1.0 / 12.0) * m * (w * w + d * d);
	y = (1.0 / 12.0) * m * (h * h + d * d);
	z = (1.0 / 12.0) * m * (h * h + w * w);
	
	return VEC3D(x, y, z);
}

double MaterialCalculator::CircleInertia(double m, double r)
{
	return 0.5 * m * r * r;
}