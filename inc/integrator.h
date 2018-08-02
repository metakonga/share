#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <QVector>

class integrator
{
public:
	integrator();
	~integrator();

	void calculatePositionFromVelocity(double p0, double dt, QVector<double> *vel, QList<double> *pos);
	void calculatePositionFromVelocity(double p0, double dt, QVector<QPointF>& points, QList<QPointF>& pos);

	double max_p;
	double min_p;
};

#endif