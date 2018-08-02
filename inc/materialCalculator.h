#ifndef MATERIALCALCULATOR_H
#define MATERIALCALCULATOR_H

#include "algebraMath.h"

namespace MaterialCalculator
{
	//double rectangleMass(double )
	VEC3D CubeInertia(double m, double w, double h, double d);
	double CircleInertia(double m, double r);
}

#endif