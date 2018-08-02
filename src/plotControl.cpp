#include "plotControl.h"
#include "chartView.h"
#include "messageBox.h"
#include <QtWidgets>

plotControl::plotControl(QWidget* parent)
	: QDockWidget(parent)
{
	setObjectName("Plot Property");
	dockWidgetContents = new QWidget();
	dockWidgetContents->setObjectName("Dock Contents");
	gridLayout = new QGridLayout(dockWidgetContents);
	TW_CONTROL = new QTabWidget(dockWidgetContents);
	TW_CONTROL->setObjectName("Property Tab");
	TW_CONTROL->setTabShape(QTabWidget::Triangular);
	W_AXIS_TAB = new QWidget();
	W_AXIS_TAB->setObjectName("Axis Tab");
	axisTabGrid = new QGridLayout(W_AXIS_TAB);
	
	LTITLE = new QLabel(W_AXIS_TAB);
	LTITLE->setObjectName(QStringLiteral("LTITLE"));
	LTITLE->setAlignment(Qt::AlignCenter);
	LTITLE->setText("Title");

	axisTabGrid->addWidget(LTITLE, 0, 1, 1, 1);

	LMINR = new QLabel(W_AXIS_TAB);
	LMINR->setObjectName(QStringLiteral("LMINR"));
	LMINR->setAlignment(Qt::AlignCenter);
	LMINR->setText("Min. range");

	axisTabGrid->addWidget(LMINR, 0, 2, 1, 1);

	LMAXR = new QLabel(W_AXIS_TAB);
	LMAXR->setObjectName(QStringLiteral("LMAXR"));
	LMAXR->setAlignment(Qt::AlignCenter);
	LMAXR->setText("Max. range");

	axisTabGrid->addWidget(LMAXR, 0, 3, 1, 1);

	LTICK = new QLabel(W_AXIS_TAB);
	LTICK->setObjectName(QStringLiteral("LTICK"));
	LTICK->setAlignment(Qt::AlignCenter);
	LTICK->setText("Tick count");

	axisTabGrid->addWidget(LTICK, 0, 4, 1, 1);

	L_X = new QLabel(W_AXIS_TAB);
	L_X->setObjectName(QStringLiteral("L_X"));
	L_X->setText("X");

	axisTabGrid->addWidget(L_X, 1, 0, 1, 1);
	QLineEdit* le;
	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_XTITLE"));
	le->setText("None");
	axisTabGrid->addWidget(le, 1, 1, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEXT, le);

	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_XMINR"));
	le->setText("0");
	axisTabGrid->addWidget(le, 1, 2, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEXMIN, le);

	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_XMAXR"));
	le->setText("1.0");
	axisTabGrid->addWidget(le, 1, 3, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEXMAX, le);

	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_XTICKCOUNT"));
	le->setText("5");
	axisTabGrid->addWidget(le, 1, 4, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEXTICK, le);

	L_Y = new QLabel(W_AXIS_TAB);
	L_Y->setObjectName(QStringLiteral("L_Y"));
	L_Y->setText("Y");
	axisTabGrid->addWidget(L_Y, 2, 0, 1, 1);

	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_YTITLE"));
	le->setText("None");
	axisTabGrid->addWidget(le, 2, 1, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEYT, le);

	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_YMINR"));
	le->setText("0");
	axisTabGrid->addWidget(le, 2, 2, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEYMIN, le);

	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_YMAXR"));
	le->setText("1.0");
	axisTabGrid->addWidget(le, 2, 3, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEYMAX, le);

	le = new QLineEdit(W_AXIS_TAB);
	le->setObjectName(QStringLiteral("LE_YTICKCOUNT"));
	le->setText("5");
	axisTabGrid->addWidget(le, 2, 4, 1, 1);
	connect(le, SIGNAL(editingFinished()), this, SLOT(changeProperty()));
	LEs.insert(LEYTICK, le);

	TW_CONTROL->addTab(W_AXIS_TAB, QString());
	W_VOID_TAB = new QWidget();
	W_VOID_TAB->setObjectName(QStringLiteral("Void Tab"));
	TW_CONTROL->addTab(W_VOID_TAB, QString());

	gridLayout->addWidget(TW_CONTROL, 0, 1, 1, 1);

	setWidget(dockWidgetContents);
}

plotControl::~plotControl()
{

}

void plotControl::bindChartView(chartView* pw)
{
	vcht = pw;
	vcht->Chart()->axisX()->setTitleText(LEs.at(LEXT)->text());
	vcht->Chart()->axisY()->setTitleText(LEs.at(LEYT)->text());
}


void plotControl::setPlotProperty(QString nm, QString xt, QString yt, double xmin, double xmax, double ymin, double ymax)
{
	LEs.at(LEXT)->setText(xt);
	LEs.at(LEYT)->setText(yt);
	LEs.at(LEXMIN)->setText(QString("%1").arg(xmin));
	LEs.at(LEXMAX)->setText(QString("%1").arg(xmax));
	LEs.at(LEYMIN)->setText(QString("%1").arg(ymin));
	LEs.at(LEYMAX)->setText(QString("%1").arg(ymax));
	vcht->Chart()->axisX()->setTitleText(xt);
	vcht->Chart()->axisY()->setTitleText(yt);
	vcht->setAxisXRange(xmin, xmax);
	vcht->setAxisYRange(ymin, ymax);
	plotControlData pcd = { xt, yt, 5, 5, xmin, xmax, ymin, ymax };
	pcds[nm] = pcd;
	curPlot = nm;
}

void plotControl::setPlotProperty(QString pn)
{
	plotControlData pcd = pcds[pn];
	LEs.at(LEXT)->setText(pcd.xt);
	LEs.at(LEYT)->setText(pcd.yt);
	LEs.at(LEXMIN)->setText(QString("%1").arg(pcd.xmin));
	LEs.at(LEXMAX)->setText(QString("%1").arg(pcd.xmax));
	LEs.at(LEYMIN)->setText(QString("%1").arg(pcd.ymin));
	LEs.at(LEYMAX)->setText(QString("%1").arg(pcd.ymax));
	vcht->Chart()->axisX()->setTitleText(pcd.xt);
	vcht->Chart()->axisY()->setTitleText(pcd.yt);
	vcht->setAxisXRange(pcd.xmin, pcd.xmax);
	vcht->setAxisYRange(pcd.ymin, pcd.ymax);
	vcht->setTickCountX(pcd.xtick);
	vcht->setTickCountY(pcd.ytick);
	curPlot = pn;
}

plotControlData plotControl::getPlotProperty(QString pn)
{
	return pcds[pn];
}

void plotControl::changeProperty()
{
	QLineEdit* le = (QLineEdit*)sender();
	QString str = le->text();
	if (le == LEs.at(LEXT))
	{
		vcht->Chart()->axisX()->setTitleText(str);
		pcds[curPlot].xt = str;
	}
	else if (le == LEs.at(LEYT))
	{
		vcht->Chart()->axisY()->setTitleText(str);
		pcds[curPlot].yt = str;
	}
	else if (le == LEs.at(LEXMIN) || le == LEs.at(LEXMAX) || le == LEs.at(LEXTICK))
	{
		double xmin = LEs.at(LEXMIN)->text().toDouble();
		double xmax = LEs.at(LEXMAX)->text().toDouble();
		int xtick = LEs.at(LEXTICK)->text().toInt();
		if (xmin >= xmax)
		{
			messageBox::run("Maximum range must be larger than minimum range.");
			return;
		}
			
		vcht->setAxisXRange(xmin, xmax);
		vcht->setTickCountX(xtick);
		pcds[curPlot].xmin = xmin;
		pcds[curPlot].xmax = xmax;
		pcds[curPlot].xtick = xtick;
	}
	else if (le == LEs.at(LEYMIN) || le == LEs.at(LEYMAX) || le == LEs.at(LEYTICK))
	{
		double ymin = LEs.at(LEYMIN)->text().toDouble();
		double ymax = LEs.at(LEYMAX)->text().toDouble();
		int ytick = LEs.at(LEYTICK)->text().toInt();
		if (ymin >= ymax)
		{
			messageBox::run("Maximum range must be larger than minimum range.");
			return;
		}
		vcht->setAxisYRange(ymin, ymax);
		vcht->setTickCountY(ytick);
		pcds[curPlot].ymin = ymin;
		pcds[curPlot].ymax = ymax;
		pcds[curPlot].ytick = ytick;
	}
}
