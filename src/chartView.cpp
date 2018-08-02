#include "chartView.h"
#include "callout.h"
#include "integrator.h"
#include "algebraMath.h"
#include <QString>
#include <QDebug>
#include <QPixmap>
#include <QImage>
#include <QTableWidgetItem>
#include <QtGui/QMouseEvent>
#include <QtWidgets/QGraphicsScene>

chartView::chartView(QWidget* parent, ChartMode cm)
	: QGraphicsView(new QGraphicsScene, parent)
	, ax(NULL)
	, ay(NULL)
	, m_coordX(NULL)
	, m_coordY(NULL)
	, m_coordHoverX(NULL)
	, m_coordHoverY(NULL)
	, m_rectHovered(NULL)
	, m_lineItemX(NULL)
	, m_lineItemY(NULL)
	, m_checkItem(NULL)
	, mchs(NULL)
	, table(NULL)
	, m_chart(NULL)
	, m_tooltip(NULL)
	, previous_checking_distance(0)
	, checked_point_number(-1)
	, onMousePress(false)
	, onMouseMiddleButton(false)
	, onClickEvent(true)
	, cmode(cm)
{
	setDragMode(QGraphicsView::NoDrag);
	setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	m_scene = scene();
	ax = new axis(0, 1);
	ay = new axis(0, 1);
	m_chart = new QChart;
	//m_chart->setMinimumSize(400, 400);
	m_chart->setMargins(QMargins(0, 0, 0, 0));
	m_coordX = new QGraphicsSimpleTextItem(m_chart);
	m_coordY = new QGraphicsSimpleTextItem(m_chart);
	QPen penBorder;
	penBorder.setColor(QColor(0, 0, 0));
	penBorder.setWidth(1);
	m_coordX->setPen(penBorder);
	m_coordY->setPen(penBorder);
	m_rectHovered = new QGraphicsRectItem(m_chart);
	m_rectHovered->setBrush(QBrush(Qt::yellow));
	m_coordHoverX = new QGraphicsSimpleTextItem(m_rectHovered);
	m_coordHoverY = new QGraphicsSimpleTextItem(m_rectHovered);
	penBorder.setColor(QColor(0, 0, 0));
	penBorder.setWidth(1);
	m_coordHoverX->setPen(penBorder);
	m_coordHoverY->setPen(penBorder);

	m_lineItemX = new QGraphicsLineItem(m_chart);
	m_lineItemY = new QGraphicsLineItem(m_chart);
	QPen penLine;
	penLine.setColor(QColor(0, 0, 0));
	penLine.setStyle(Qt::DotLine);
	m_lineItemX->setPen(penLine);
	m_lineItemY->setPen(penLine);

	QImage img;
	QPixmap buffer;
	img.load(":/Resources/check_icon.png");
	buffer = QPixmap::fromImage(img);
	buffer = buffer.scaled(20, 20);
	m_checkItem = new QGraphicsPixmapItem(buffer, m_chart);
	m_checkItem->setVisible(false);
	m_chart->setAcceptHoverEvents(true);
	setMouseTracking(true);
	m_chart->setZValue(50);
	m_coordHoverX->setZValue(20);
	setRenderHint(QPainter::Antialiasing);
	m_chart->setAxisX(ax);
	m_chart->setAxisY(ay);
// 	if (cmode == REALTIME_EDIT_CHART)
// 	{
// 		QLineSeries *series = new QLineSeries;
// 		QScatterSeries *p_series = new QScatterSeries;
// 		p_series->setMarkerSize(10);
// 
// 		m_chart->addSeries(series);
// 		m_chart->addSeries(p_series);
// 		m_chart->setAxisX(ax, series);
// 		m_chart->setAxisY(ay, series);
// 		p_series->attachAxis(ax);
// 		p_series->attachAxis(ay);
// 	}	
	m_chart->legend()->setVisible(false);
	m_scene->addItem(m_chart);
	mchs = new MCHSpline;
}

chartView::~chartView()
{
	if (ax) delete ax; ax = NULL;
	if (ay) delete ay; ay = NULL;
	if (m_coordX) delete m_coordX; m_coordX = NULL;
	if (m_coordY) delete m_coordY; m_coordY = NULL;
	if (m_coordHoverX) delete m_coordHoverX; m_coordHoverX = NULL;
	if (m_coordHoverY) delete m_coordHoverY; m_coordHoverY = NULL;
	if (m_rectHovered) delete m_rectHovered; m_rectHovered = NULL;
	if (m_checkItem) delete m_checkItem; m_checkItem = NULL;
	if (m_lineItemX) delete m_lineItemX; m_lineItemX = NULL;
	if (m_lineItemY) delete m_lineItemY; m_lineItemY = NULL;
	if (mchs) delete mchs; mchs = NULL;
	if (m_chart) m_chart->removeAllSeries();
	if (m_chart) delete m_chart; m_chart = NULL;
	if (m_scene) delete m_scene; m_scene = NULL;
}

void chartView::resizeEvent(QResizeEvent *event)
{
	if (scene())
	{
		scene()->setSceneRect(QRect(QPoint(0, 0), event->size()));
		m_chart->resize(event->size());
		foreach(Callout *callout, m_callouts)
			callout->updateGeometry();
	}
	QGraphicsView::resizeEvent(event);
}

void chartView::addSeries(QLineSeries *_series)
{
	m_chart->addSeries(_series);
	if (cmode == REALTIME_EDIT_CHART)
	{
		connect(_series, SIGNAL(clicked(QPointF)), this, SLOT(keepCallout()));
		connect(_series, SIGNAL(hovered(QPointF, bool)), this, SLOT(tooltip(QPointF, bool)));
	}
	m_chart->setAxisX(ax, _series);
	m_chart->setAxisY(ay, _series);
}

void chartView::keepCallout()
{
	m_callouts.append(m_tooltip);
	m_tooltip = new Callout(m_chart);
}

void chartView::tooltip(QPointF point, bool state)
{
	if (m_tooltip == 0)
		m_tooltip = new Callout(m_chart);
	
	if (state) {
		m_tooltip->setText(QString("X: %1 \nY: %2 ").arg(point.x(), 5, 'f', 3, '0').arg(point.y(), 5, 'f', 3, '0'));
		m_tooltip->setAnchor(point);
		m_tooltip->setZValue(1);
		m_tooltip->updateGeometry();
		m_tooltip->show();
	}
	else {
		m_tooltip->hide();
	}
}

void chartView::setChartMode(ChartMode _cmode)
{
	cmode = _cmode;
}

void chartView::setTickCountX(int c)
{
	ax->setTickCount(c);
// 	QFont font = ax->labelsFont();
// 	font.setPointSize(8);
// 	ax->setLabelsFont(font);
}

void chartView::setTickCountY(int c)
{
	ay->setTickCount(c);
}

bool chartView::setChartData(QVector<double>* x, QVector<double>* y, int d)
{
	QLineSeries *series = new QLineSeries;
	int n = x->size();
	//QVector<double>::iterator ix = x->begin();
	//QVector<double>::iterator iy = y->begin();
	double x_max = FLT_MIN;
	double x_min = FLT_MAX;
	double y_max = FLT_MIN;
	double y_min = FLT_MAX;
	int i = 0;
	for (i = 0; i < n; i += d)
	{
		double _x = x->at(i);
		double _y = y->at(i);
		series->append(_x, _y);
		x_min = _x < x_min ? _x : x_min;
		x_max = _x > x_max ? _x : x_max;
		y_min = _y < y_min ? _y : y_min;
		y_max = _y > y_max ? _y : y_max;
	}
 	if (i >= n)
 		series->append(x->at(n - 1), y->at(n -1));
	if (m_chart)
	{
		
	}
	else
		return false;
	return true;
}

void chartView::mouseReleaseEvent(QMouseEvent *event)
{
	if (cmode == REALTIME_EDIT_CHART)
	{
		if (onMousePress && checked_point_number >= 0)
		{
			QScatterSeries *p_series = scatterSeries[ctr_scatter_idx];// <QScatterSeries*>(m_chart->series().at(1));
			if (p_series->count() > 2)
			{
				QLineSeries *series = lineSeries[ctr_line_idx];// <QLineSeries*>(m_chart->series().at(0));
				QPointF new_p = m_chart->mapToValue(movingPos);
				QString nx = QString("%1").arg(new_p.x(), 5, 'f', 5, '0');
				QString ny = QString("%1").arg(new_p.y(), 5, 'f', 5, '0');
				new_p.setX(nx.toDouble());
				new_p.setY(ny.toDouble());
				p_series->replace(checked_point_number, new_p);
				mchs->update_curve(checked_point_number, new_p, series);
				mchs->calculate_curve_auto(series, checked_point_number, checked_point_number + 1);
				QTableWidgetItem* item = table->item(checked_point_number, 1);
				item->setText(QString("%1").arg(new_p.x()));
				item = table->item(checked_point_number, 2);
				item->setText(QString("%1").arg(new_p.y()));
			}
		}
	}
	else if (cmode == ONLY_DISPLAY_CHART)
	{

	}
	onMousePress = false;
	onMouseMiddleButton = false;
}

void chartView::updateSeries(double newXmax, double newYmax)
{
	double xMax = ax->max();
	double yMax = ay->max();
	QLineSeries *series = lineSeries[ctr_line_idx];
	QScatterSeries *p_series = scatterSeries[ctr_scatter_idx];
	QList<QPointF> series_list = series->points();
	QList<QPointF> new_series;
	QList<QPointF> p_series_list = p_series->points();
	QList<QPointF> new_p_series;
	foreach(QPointF value, series_list)
	{
		new_series.push_back(QPointF(value.x() * (xMax / newXmax), value.y() * (yMax / newYmax)));
	}
	foreach(QPointF value, p_series_list)
	{
		new_p_series.push_back(QPointF(value.x() * (xMax / newXmax), value.y() * (yMax / newYmax)));
	}
	series->clear();
	p_series->clear();
	series->replace(new_series);
	p_series->replace(new_p_series);
}

void chartView::setAxisRange(double x_min, double x_max, double y_min, double y_max, bool _mof)
{
	//if (m_chart->series().size() == 1)
	//{

	ax->setRange(x_min, x_max);
	ay->setRange(y_min, y_max);
	ax->setMax(x_max);
	ay->setMax(y_max);	
	//	return;
	//}
// 	double xmin = ax->min() > x_min ? x_min : ax->min();
// 	double xmax = ax->max() < x_max ? x_max : ax->max();
// 	double ymin = ay->min() > y_min ? y_min : ay->min();
// // 	double ymax = ay->max() < y_max ? y_max : ay->max();
// 	ax->setRange(xmin, xmax);
// 	ay->setRange(ymin, ymax);
// 	ax->setMax(xmax);
// 	ay->setMax(ymax);
}

void chartView::setAxisXRange(double _min, double _max)
{
	ax->setRange(_min, _max);
	ax->setMin(_min);
	ax->setMax(_max);
}

void chartView::setAxisYRange(double _min, double _max)
{
	ay->setRange(_min, _max);
	ay->setMin(_min);
	ay->setMax(_max);
}

void chartView::wheelEvent(QWheelEvent *event)
{
// 	QPoint p = event->angleDelta();
// 	double xmin = ax->min();
// 	double ymin = ay->min();
// 	double xmax = ax->max();
// 	double ymax = ay->max();
// 	if (p.y() > 0)
// 	{
// 		ax->setRange(0, xmax + 0.1);
// 		ay->setRange(ymin - 0.1, ymax + 0.1);
// 		ax->setMin(0);
// 		ay->setMin(ymin - 0.1);
// 		ax->setMax(xmax + 0.1);
// 		ay->setMax(ymax + 0.1);
// 	}
// 	else
// 	{
// 		ax->setRange(0, xmax - 0.1);
// 		ay->setRange(ymin + 0.1, ymax - 0.1);
// 		ax->setMin(0);
// 		ay->setMin(ymin + 0.1);
// 		ax->setMax(xmax - 0.1);
// 		ay->setMax(ymax - 0.1);
// 	}
}

void chartView::keyPressEvent(QKeyEvent *event)
{
	switch (event->key()){
	case Qt::Key_Up:

		break;
	case Qt::Key_Down:

		break;
	case Qt::Key_Left:

		break;
	case Qt::Key_Right:

		break;
	}
}

void chartView::mouseDoubleClickEvent(QMouseEvent *event)
{

}

void chartView::mousePressEvent(QMouseEvent *event)
{
	if (!onClickEvent)
		return;
	if (cmode == chartView::ONLY_DISPLAY_CHART)
		return;
	if (event->button() == Qt::MiddleButton)
	{
		onMouseMiddleButton = true;
	}
	else if (event->button() == Qt::LeftButton)
	{
		onMousePress = true;
		if (checked_point_number >= 0)
			return;
		QPointF xy = event->pos();
		QPointF p = m_chart->mapToValue(xy);
		QString spx = QString("%1").arg(p.x(), 5, 'f', 5, '0');
		QString spy = QString("%1").arg(p.y(), 5, 'f', 5, '0');
		QPointF p_r(spx.toDouble(), spy.toDouble());
		if (p_r.x() <= mchs->limitation())
			return;

		QLineSeries *series = NULL;
		QScatterSeries *p_series = NULL;
		if (m_chart->series().size() == 0)
		{

		}
		else if (m_chart->series().size() && cmode == REALTIME_EDIT_CHART)
		{
			series = lineSeries[ctr_line_idx];
			p_series = scatterSeries[ctr_scatter_idx];
			mchs->monotonic_cubic_Hermite_spline(p_r);
			mchs->setData2LineSeries(series);
			p_series->append(p_r.x(), p_r.y());
			if (table)
			{
				int nr = p_series->count() - 1;
				QString id = QString("%1").arg(nr);
				table->insertRow(nr);
				table->setItem(nr, 0, new QTableWidgetItem(id));
				table->setItem(nr, 1, new QTableWidgetItem(spx));
				table->setItem(nr, 2, new QTableWidgetItem(spy));
			}
			//calculationStrok();
		}
		
	}
	//qDebug() << "Area : " << calculationStrok();
	QGraphicsView::mousePressEvent(event);
}

double chartView::calculationStrok()
{
// 	integrator itor;
// 	QLineSeries* series = dynamic_cast<QLineSeries*>(m_chart->series().at(0));
// 	QVector<double> pos;
// 	itor.calculatePositionFromVelocity(0.0, 0.00001, series->pointsVector(), pos);
// 	return pos.at(pos.size() - 1);
	return 0;
}

QMap<int, QLineSeries*>& chartView::LineSeries()
{
	return lineSeries;
}

QMap<int, QScatterSeries*>& chartView::ScatterSeries()
{
	return scatterSeries;
}

QChart* chartView::Chart()
{
	return m_chart;
}

MCHSpline* chartView::MCHS()
{
	return mchs;
}

QLineSeries* chartView::createLineSeries(int idx)
{
	QLineSeries* ls = new QLineSeries;
	m_chart->addSeries(ls);
	ls->attachAxis(ax);
	ls->attachAxis(ay);
	lineSeries[idx] = ls;
	return ls;
}

QScatterSeries* chartView::createScatterSeries(int idx)
{
	QScatterSeries* ss = new QScatterSeries;
	m_chart->addSeries(ss);
	ss->attachAxis(ax);
	ss->attachAxis(ay);
	scatterSeries[idx] = ss;
	return ss;
}

void chartView::mouseMoveEvent(QMouseEvent *event)
{
	if (cmode == chartView::ONLY_DISPLAY_CHART)
		return;
//	QGraphicsView::mouseMoveEvent(event);
	if (!(m_chart->series().size()))
		return;
	QPointF xy = event->pos();
	QPointF p = m_chart->mapToValue(event->pos());
	
	if (p.x() <= ax->max() && p.x() >= ax->min() && p.y() <= ay->max() && p.y() >= ay->min())
	{
		m_coordHoverX->setVisible(true);
		m_coordHoverY->setVisible(true);
		m_rectHovered->setVisible(true);
		m_lineItemX->setVisible(true);
		m_lineItemY->setVisible(true);

		qreal x = m_chart->mapToPosition(p).x();
		qreal y = m_chart->mapToPosition(p).y();
		
		m_rectHovered->setRect(x, y - 31, 40, 30);

		qreal rectX = m_rectHovered->rect().x();
		qreal rectY = m_rectHovered->rect().y();
		qreal rectW = m_rectHovered->rect().width();
		qreal rectH = m_rectHovered->rect().height();

		/* We're setting the labels and nicely adjusting to chart axis labels (adjusting so the dot lines are centered on the label) */
		m_coordHoverX->setPos(rectX + rectW / 4 - 3, rectY + 1);
		m_coordHoverY->setPos(rectX + rectW / 4 - 3, rectY + rectH / 2 + 1);

		QPointF xp = m_chart->mapToPosition(QPointF(p.x(), 0));
		QPointF yp = m_chart->mapToPosition(QPointF(0, p.y()));
		m_lineItemX->setLine(xp.x(), xy.y(), xy.x(), xp.y()/* - 27*/);
		m_lineItemY->setLine(xp.x(), xy.y(), yp.x(), xy.y());

		/* Setting value to displayed with four digit max, float, 1 decimal */
		m_coordHoverX->setText(QString("%1").arg(p.x(), 5, 'f', 5, '0'));
		m_coordHoverY->setText(QString("%1").arg(p.y(), 5, 'f', 5, '0'));
		
		if (onMousePress && checked_point_number >= 0 && cmode == REALTIME_EDIT_CHART)
		{
			movingPos.setX(x);
			movingPos.setY(y);
			m_checkItem->setPos(QPointF(x, y));
			return;
		}
		else if (/*onMouseMiddleButton &&*/ cmode == REALTIME_EDIT_CHART)
		{
			QPointF dist = (p - nonMovingPos) * 0.2;
			double xmin = ax->min() + dist.x();
			double xmax = ax->max() + dist.x();
			double ymin = ay->min() + dist.y();
			double ymax = ay->max() + dist.y();
			if (xmin > 0 && ymin > 0)
			{
				ax->setRange(xmin, xmax);
				ay->setRange(ymin, ymax);
			}
			nonMovingPos = p;
			checkingNearPoint(x, y);
		}
		
	}
	else
	{
		m_coordHoverX->setVisible(false);
		m_coordHoverY->setVisible(false);
		m_rectHovered->setVisible(false);
		m_lineItemX->setVisible(false);
		m_lineItemY->setVisible(false);
	}
}

int chartView::checkingNearPoint(double px, double py)
{
	int cnumber = 0;
	//QLineSeries *series = lineSeries[ctr_line_idx];
	QScatterSeries *p_series = scatterSeries[ctr_scatter_idx];
	QVector<QPointF> points = p_series->pointsVector();
	//qDebug() << "checked_point_number : " << checked_point_number;
	if (checked_point_number >= 0)
	{
		double dx = px - m_chart->mapToPosition(points.at(checked_point_number)).x();
		double dy = py - m_chart->mapToPosition(points.at(checked_point_number)).y();
		double dist = sqrt(dx * dx + dy * dy);
		if (dist < 5)
		{
			m_checkItem->setPos(px, py);
			m_checkItem->setVisible(true);
		}
		else
		{
			m_checkItem->setVisible(false);
			checked_point_number = -1;
		}
		return checked_point_number;
	}
	previous_checking_distance = 0;
	checked_point_number = -1;
	if (!points.size())
		return -1;
	foreach(QPointF p, points)
	{
		double dx = px - m_chart->mapToPosition(p).x();
		double dy = py - m_chart->mapToPosition(p).y();
		double dist = sqrt(dx*dx + dy*dy);
		if (dist < 5)
		{
			checked_point_number = cnumber;
			//qDebug() << "checked_point_number : " << checked_point_number;
			return checked_point_number;
		}
		cnumber++;
	}
	return 0;
}

void chartView::setAxisXLimit(double axl)
{
	double xmin = ax->min();
	ax->setRange(xmin, axl);
	ax->setMax(axl);
}

void chartView::setAxisYLimit(double ayl)
{
	double ymin = ay->min();
	ay->setRange(ymin, ayl);
	ay->setMax(ayl);
}

void chartView::setEnableClickEvent(bool b)
{
	onClickEvent = b;
}

void chartView::setControlLineSeries(int idx)
{
	ctr_line_idx = idx;

}

void chartView::setControlScatterSeries(int idx)
{
	ctr_scatter_idx = idx;
}

void chartView::setMCHS(int count, QPointF* points)
{
	for (unsigned int i = 0; i < count; i++)
	{
		mchs->monotonic_cubic_Hermite_spline(points[i]);
	}
}

void chartView::setMinMax(int ty, double minx, double maxx, double miny, double maxy)
{
	minMax_Axis_x[ty] = QPointF(minx, maxx);
	minMax_Axis_y[ty] = QPointF(miny, maxy);
}

void chartView::setAxisBySeries(int ty)
{
	QPointF x_range = minMax_Axis_x[ty];
	QPointF y_range = minMax_Axis_y[ty];
	setAxisRange(x_range.x(), x_range.y(), y_range.x(), y_range.y());
}

// void chartView::setAxisXY()
// {
// 	m_chart->setAxisX(ax);
// 	m_chart->setAxisY(ay);
// }

void chartView::changeSeriesValue(unsigned int x, unsigned int y, double v)
{
	QLineSeries *series = lineSeries[ctr_line_idx];
	QScatterSeries *p_series = scatterSeries[ctr_scatter_idx];// <QScatterSeries*>(m_chart->series().at(1));
	if (p_series->count() > 2)
	{
		QVector<QPointF> ss = p_series->pointsVector();
		QPointF new_p;
		if(y == 1)
			new_p = QPointF(v, ss.at(x).y());
		else if(y == 2)
			new_p = QPointF(ss.at(x).x(), v);
		p_series->replace(x, new_p);
		mchs->update_curve(x, new_p, series);
		mchs->calculate_curve_auto(series, (int)x, (int)(x + 1));
	//	qDebug() << "Area : " << calculationArea();
	}
}

void chartView::exportSplineData(QString path)
{
	//integrator itor;
	//QString filePath = path + "splineData.txt";
	QLineSeries *series = lineSeries[ctr_line_idx];
	QScatterSeries *p_series = scatterSeries[ctr_scatter_idx];
	QList<QPointF> rdata;
	//QVector<double> pos;
	mchs->calculate_curve(p_series, rdata, 0.00001);
	QVector<QPointF> rd = series->pointsVector();
	//itor.calculatePositionFromVelocity(0.0, 0.00001, rdata, pos);
	QFile qf(path);
	qf.open(QIODevice::WriteOnly);
	QTextStream qts(&qf);
 	//unsigned int sz = series->pointsVector().size();
// 	qts 
	int pVal = 1;
	for (unsigned int i = 0; i < rdata.size(); i++)
	{
		QPointF v = rdata.at(i);
		if (pVal == static_cast<int>(100000 * (v.x() + 1e-9)))
			continue;
		//double p = pos.at(i);
		qts << v.x() << " " << v.y() << endl;
		pVal = static_cast<int>(100000 * (v.x() + 1e-9));
	}
	qts << 0.035 << " " << 0 << endl;
	qf.close();
	
}