#include "integrator.h"

integrator::integrator()
{

}

integrator::~integrator()
{

}

void integrator::calculatePositionFromVelocity(double p0, double dt, QVector<double> *vel, QList<double> *pos)
{
	int n = vel->size();
	pos->push_back(p0);
	for (unsigned int i = 1; i < n; i++)
	{
		double o = pos->at(i - 1) + 0.5*(vel->at(i - 1) + vel->at(i)) * dt;
		pos->push_back(o);
	}
}

void integrator::calculatePositionFromVelocity(double p0, double dt, QVector<QPointF>& points, QList<QPointF>& pos)
{
	max_p = DBL_MIN;
	min_p = DBL_MAX;
	int n = points.size();
	pos.push_back(QPointF(0.0, p0));
	for (unsigned int i = 1; i < n; i++)
	{
		double dx = points.at(i).x() - points.at(i - 1).x();
		double o = pos.at(i - 1).y() + 0.5*(points.at(i - 1).y() + points.at(i).y()) * dx;
		pos.push_back(QPointF(points.at(i).x(), o));
		if (max_p < o)
			max_p = o;
		if (min_p > o)
			min_p = o;
	}
}