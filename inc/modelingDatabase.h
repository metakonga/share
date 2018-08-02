#ifndef MODELINGDATABASE_H
#define MODELINGDATABASE_H

#include <QDockWidget>
#include <QTreeWidget>
#include <QComboBox>

#include "sph_model.h"
#include "sensor.h"

class modelingDatabase : public QDockWidget
{
public:
	enum tRoot { EXPRESSION_ROOT = 0, SENSOR_ROOT, RIGIDBODY_ROOT };
	modelingDatabase(QWidget* parent);
	~modelingDatabase();

	void addChild(tRoot, QString& _nm);
	QStringList selectedLists();
	tRoot selectedType();

	private slots:
	void contextMenu(const QPoint&);
	void clickItem(QTreeWidgetItem*, int);

private:
	QTreeWidget *tree;
	QMap<tRoot, QTreeWidgetItem*> roots;
	QStringList sLists;
	tRoot tSelected;
};

#endif