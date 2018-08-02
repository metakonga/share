#ifndef AXIS_H
#define AXIS_H

#include <QValueAxis>

QT_CHARTS_USE_NAMESPACE

class axis : public QValueAxis
{
public:
	explicit axis(double min, double max);
	~axis();
};

#endif