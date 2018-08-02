#ifndef RESULT_DATABASE_H
#define RESULT_DATABASE_H

#include <QDockWidget>
#include <QTreeWidget>
#include <QComboBox>

#include "sph_model.h"
#include "sensor.h"

class resultDatabase : public QDockWidget
{
public:
	enum tRoot { PART_ROOT = 0, SENSOR_ROOT, PMASS_ROOT, REACTION_ROOT };
	resultDatabase(QWidget* parent);
	~resultDatabase();

	void addChild(tRoot, QString& _nm);
	void bindItemComboBox(QComboBox* t);
	QStringList selectedLists();
	tRoot selectedType();
	QComboBox* plotItemComboBox();
	QString plotTarget();

	private slots:
	void contextMenu(const QPoint&);
	void clickItem(QTreeWidgetItem*, int);

private:
	QTreeWidget *tree;
	QMap<tRoot, QTreeWidgetItem*> roots;
	QStringList sLists;
	tRoot tSelected;
	QComboBox *plot_item;
	QString target;
};

#endif