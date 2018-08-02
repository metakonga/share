#ifndef CHARTVIEW_H
#define CHARTVIEW_H

//#include <QValueAxis>
#include <QChartView>
#include <QLineSeries>
#include <QScatterSeries>
#include <QGraphicsSimpleTextItem>
#include <QGraphicsPixmapItem>
#include <QGraphicsRectItem>
#include <QGraphicsLineItem>
#include <QTableWidget>
#include "MCHSpline.h"

//#include "chart.h"
#include "axis.h"

QT_BEGIN_NAMESPACE
class QGraphicsScene;
class QMouseEvent;
class QResizeEvent;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE
class QChart;
QT_CHARTS_END_NAMESPACE

class Callout;

QT_CHARTS_USE_NAMESPACE

class chartView : public QGraphicsView
{
	Q_OBJECT

public:
	enum ChartMode{ REALTIME_EDIT_CHART = 0, ONLY_DISPLAY_CHART };

	chartView(QWidget* parent = 0, ChartMode cm = ONLY_DISPLAY_CHART);
	~chartView();

	
	void addSeries(QLineSeries *_series);
	void setAxisRange(double x_min, double x_max, double y_min, double y_max, bool _mof = false);
	void setAxisXRange(double _min, double _max);
	void setAxisYRange(double _min, double _max);
	bool setChartData(QVector<double>* x, QVector<double>* y, int d = 1);
	void setTableWidget(QTableWidget* _table) { table = _table; }
	void setChartMode(ChartMode _cmode);

	void setTickCountX(int c);
	void setTickCountY(int c);
	void setAxisXLimit(double axl);
	void setAxisYLimit(double ayl);
	void setEnableClickEvent(bool b);
	void setControlLineSeries(int idx);
	void setControlScatterSeries(int idx);
	void setMCHS(int count, QPointF* points);
	void setMinMax(int ty, double minx, double maxx, double miny, double maxy);
	void setAxisBySeries(int ty);
	//void setAxisXY(int idx);

	void changeSeriesValue(unsigned int x, unsigned int y, double v);
	void exportSplineData(QString path);

	QMap<int, QLineSeries*>& LineSeries();
	QMap<int, QScatterSeries*>& ScatterSeries();
	QChart* Chart();
	MCHSpline* MCHS();
	QLineSeries* createLineSeries(int idx);
	QScatterSeries* createScatterSeries(int idx);

private:
	void mousePressEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseDoubleClickEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void keyPressEvent(QKeyEvent *event);
	void resizeEvent(QResizeEvent *event);

	int checkingNearPoint(double px, double py);
	void updateSeries(double newXmax, double newTmax);
	double calculationStrok();

	public slots:
	void keepCallout();
	void tooltip(QPointF, bool);


private:
	int ctr_line_idx;
	int ctr_scatter_idx;
	ChartMode cmode;
	bool onClickEvent;
	bool onMousePress;
	bool onMouseMiddleButton;
	QPointF movingPos;
	QPointF nonMovingPos;
	int checked_point_number;
	double previous_checking_distance;
	axis *ax;
	axis *ay;
	QGraphicsPixmapItem *m_checkItem;
	QGraphicsRectItem *m_rectHovered;
	QGraphicsLineItem *m_lineItemX;
	QGraphicsLineItem *m_lineItemY;
	QGraphicsSimpleTextItem *m_coordX;// = new QGraphicsSimpleTextItem(chart());
	QGraphicsSimpleTextItem *m_coordY;// = new QGraphicsSimpleTextItem(chart());
	QGraphicsSimpleTextItem *m_coordHoverX;
	QGraphicsSimpleTextItem *m_coordHoverY;


	QGraphicsScene *m_scene;
	QChart *m_chart;
	MCHSpline *mchs;
	QTableWidget *table;
	Callout *m_tooltip;
	QList<Callout *> m_callouts;
	QMap<int, QPointF> minMax_Axis_x;
	QMap<int, QPointF> minMax_Axis_y;
	QMap<int, QLineSeries*> lineSeries;
	QMap<int, QScatterSeries*> scatterSeries;
};

#endif