#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

#include <QMainWindow>
#include "chartView.h"
//#include <QColumnView>
#include "resultDatabase.h"
#include "plotControl.h"
#include "resultStorage.h"

//class model;
class QStandardItem;
class QStandardItemModel;
class Callout;
class QToolBar;
class QLineEdit;
class waveHeightSensor;

class plotWindow : public QMainWindow
{
	Q_OBJECT

	enum ACTION_ID { WAVE_HEIGHT = 0 };

	struct 
	{ 
		int begin; 
		int end; 
		double location;
	}waveHeightInputData;

public:
	plotWindow(QWidget* parent = NULL);
	virtual ~plotWindow();

	void setResultStorage(resultStorage* _rs);
	void uploadingResults();
	void sensorItemPlot();
	void pointMassItemPlot();
	void closeEvent(QCloseEvent *event);

	static bool IsActivate() { return isActivate; }
	static bool isActivate;

	public slots:
	void updatePlot();

	private slots:
	void click_waveHeight();
	void changeComboBoxItem(int);
	void editingCommand();
	
private:
	QLineSeries* createLineSeries(QString n);
	//void columnViewAdjustSameSize(int ncol);
	//void setBodyRoot(QStandardItem *p);
	/*void setJointForceRoot(QStandardItem *p);*/
	//void resizeEvent(QResizeEvent *event);

private:
	bool isAutoUpdateProperties;
	bool isEditingCommand;
	int commandStatus;
	int wWidth;
	int wHeight;
	int xSize;
	int ySize;
	int openColumnCount;
	QString curPlotName;
	chartView *vcht;
	resultDatabase *tree;
	plotControl *prop;
	QStandardItemModel *cmodel;
	resultStorage* rs;
	Callout *m_tooltip;
	QList<Callout *> m_callouts;
	QToolBar *mainToolBar;
	QMap<ACTION_ID, QAction*> actions;
	//QComboBox* plot_item;
	QMap<QString, QLineSeries*> seriesMap;
	QDockWidget *commDock;
	QLineEdit *comm;
	waveHeightSensor *whs;



};
#endif