#include "axis.h"

axis::axis(double min, double max)
	: QValueAxis(NULL)
{
	setMin(min);
	setMax(max);
}

axis::~axis()
{

}