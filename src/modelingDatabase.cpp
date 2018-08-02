#include "modelingDatabase.h"

#include "resultDatabase.h"
#include <QAction>
#include <QMenu>
#include <QDebug>

modelingDatabase::modelingDatabase(QWidget* parent)
	: QDockWidget(parent)
{
	setWindowTitle("Result database");
	tree = new QTreeWidget;
	tree->setColumnCount(1);
	tree->setContextMenuPolicy(Qt::CustomContextMenu);
	roots[EXPRESSION_ROOT] = new QTreeWidgetItem(tree);
	roots[SENSOR_ROOT] = new QTreeWidgetItem(tree);
	roots[RIGIDBODY_ROOT] = new QTreeWidgetItem(tree);
	//roots[CUBE_ROOT] = new QTreeWidgetItem(vtree);
	//roots[CYLINDER_ROOT] = new QTreeWidgetItem(vtree);
	//roots[POLYGON_ROOT] = new QTreeWidgetItem(vtree);
	//roots[MASS_ROOT] = new QTreeWidgetItem(vtree);

	roots[EXPRESSION_ROOT]->setText(0, "Expression");
	roots[SENSOR_ROOT]->setText(0, "Sensor");
	roots[RIGIDBODY_ROOT]->setText(0, "Rigid body");
	/*roots[CUBE_ROOT]->setText(0, "Cube");
	roots[CYLINDER_ROOT]->setText(0, "Cylinder");
	roots[POLYGON_ROOT]->setText(0, "Polygon");
	roots[MASS_ROOT]->setText(0, "Mass");
	*/
	roots[EXPRESSION_ROOT]->setIcon(0, QIcon(":/Resources/expression.png"));
	roots[SENSOR_ROOT]->setIcon(0, QIcon(":/Resources/Sensor_icon.png"));
	roots[RIGIDBODY_ROOT]->setIcon(0, QIcon(":/Resources/mass.png"));
	/*roots[CUBE_ROOT]->setIcon(0, QIcon(":/Resources/pRec.png"));
	roots[CYLINDER_ROOT]->setIcon(0, QIcon(":/Resources/cylinder.png"));
	roots[POLYGON_ROOT]->setIcon(0, QIcon(":/Resources/icPolygon.png"));
	roots[MASS_ROOT]->setIcon(0, QIcon(":/Resources/mass.png"));*/
	connect(tree, &QTreeWidget::customContextMenuRequested, this, &modelingDatabase::contextMenu);
	connect(tree, &QTreeWidget::itemClicked, this, &modelingDatabase::clickItem);

	tree->setSelectionMode(QAbstractItemView::SelectionMode::ContiguousSelection);
	setWidget(tree);
}

void modelingDatabase::contextMenu(const QPoint& pos)
{
	QTreeWidgetItem* item = tree->itemAt(pos);
	if (!item->parent())
		return;
	QAction *act = new QAction(tr("Property"), this);
	act->setStatusTip(tr("property menu"));
	//connect(act, SIGNAL(triggered()), this, SLOT(actProperty()));

	QMenu menu(this);
	menu.addAction(act);

	QPoint pt(pos);
	menu.exec(tree->mapToGlobal(pos));
}

void modelingDatabase::clickItem(QTreeWidgetItem* item, int col)
{
// 	if (!item)
// 		return;
// 	QList<QTreeWidgetItem*> items = tree->selectedItems();
// 	if (items.size() > 1)
// 	{
// 		foreach(QTreeWidgetItem* it, items)
// 		{
// 			sLists.push_back(it->text(col));
// 		}
// 		return;
// 	}
// 	QString target = item->text(col);
// 	QTreeWidgetItem* parent = item->parent();
// 	if (!parent)
// 	{
// 		return;
// 	}
// 	QString sp = parent->text(0);
// 	if (sp == "Part")
// 	{
// 
// 	}
// 	else if (sp == "Sensor")
// 	{
// 		tSelected = SENSOR_ROOT;
// 		sLists.push_back(target);
// 		plot_item->clear();
// 		sensor* s = sph_model::SPHModel()->Sensors()[target];
// 		switch (s->sensorType())
// 		{
// 		case sensor::WAVE_HEIGHT_SENSOR:
// 			plot_item->insertItem(0, "Wave height");
// 			break;
// 		}
// 
// 	}

}

modelingDatabase::~modelingDatabase()
{
	if (tree) delete tree; tree = NULL;
}

void modelingDatabase::addChild(tRoot tr, QString& _nm)
{
	QTreeWidgetItem* child = new QTreeWidgetItem();
	child->setText(0, _nm);
	roots[tr]->addChild(child);
}

QStringList modelingDatabase::selectedLists()
{
	return sLists;
}

modelingDatabase::tRoot modelingDatabase::selectedType()
{
	return tSelected;
}
