#ifndef PLOTCONTROL_H
#define PLOTCONTROL_H

#include <QDockWidget>
#include <QLineEdit>
#include <QLabel>
#include <QGridLayout>
#include <QTabWidget>
#include <QList>

class chartView;

typedef struct
{
	QString xt;
	QString yt;
	int xtick;
	int ytick;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
}plotControlData;

class plotControl : public QDockWidget
{
	Q_OBJECT
	enum lineEditType{ LEXT = 0, LEXMIN, LEXMAX, LEXTICK, LEYT, LEYMIN, LEYMAX, LEYTICK };
public:
	plotControl(QWidget* parent);
	~plotControl();

	void bindChartView(chartView* pw);
	void setPlotProperty(QString nm, QString xt, QString yt, double xmin, double xmax, double ymin, double ymax);
	void setPlotProperty(QString pn);
	plotControlData getPlotProperty(QString pn);
	QList<QLineEdit*> LEs;

	private slots:
	void changeProperty();

private:
	QWidget *dockWidgetContents;
	QGridLayout *gridLayout;
	QTabWidget *TW_CONTROL;
	QWidget *W_AXIS_TAB;
	QGridLayout *axisTabGrid;
	QLabel *LTITLE;
	QLabel *LMINR;
	QLabel *LMAXR;
	QLabel *LTICK;
	QLabel *L_X;
	QLabel *L_Y;
	QWidget *W_VOID_TAB;
	chartView* vcht;
	QString curPlot;
	QMap<QString, plotControlData> pcds;
};

#endif