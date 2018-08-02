#ifndef MCHSPLINE_H
#define MCHSPLINE_H

#include <QVector>
#include <QLineSeries>
#include <QScatterSeries>

QT_CHARTS_USE_NAMESPACE

class MCHSpline
{
public:
	MCHSpline();
	~MCHSpline();

	bool monotonic_cubic_Hermite_spline(double limit, QVector<double>& x_src, QVector<double>& y_src);
	bool monotonic_cubic_Hermite_spline(int idx, QPointF& new_p);
	bool monotonic_cubic_Hermite_spline(QPointF& new_p);
	void setData2LineSeries(QLineSeries* series);
	void calculate_curve_auto(QLineSeries* series, int sid, bool _base = false);
	void calculate_curve_auto(QLineSeries* series, int sid, int eid);
	void calculate_curve_auto(QVector<QPointF>& data);
	void calculate_curve(QScatterSeries* pseries, QList<QPointF>& rdata, double h);
	QVector<double>* xSource() { return &x_src; }
	QVector<double>* ySource() { return &y_src; }
	void update_curve(int idx, QPointF& new_p, QLineSeries* series);
	double limitation() { return limit; }
	double getInterpValue(unsigned int id, double rx);
	double getInterpValue(double rx);
	double calculate_derivative(double x);
	double calculate_dderivative(double x);

private:
	double h00(double t)
	{
		return 2 * t*t*t - 3 * t*t + 1;
	}
	double h10(double t)
	{
		return t*(1 - t)*(1 - t);
	}
	double h01(double t)
	{
		return t*t*(3 - 2 * t);
	}
	double h11(double t)
	{
		return t*t*(t - 1);
	}

	void monotonicity(int i, int j);

private:
	double limit;
	QVector<double> m;
	QVector<double> x_src;
	QVector<double> y_src;
	QVector<int> start_ids;
// 	QVector<double> x_data;
// 	QVector<double> y_data
};

#endif