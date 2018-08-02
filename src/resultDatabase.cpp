#include "resultDatabase.h"
#include "cmdWindow.h"
#include <QAction>
#include <QMenu>
#include <QDebug>

resultDatabase::resultDatabase(QWidget* parent)
	: QDockWidget(parent)
{
	setWindowTitle("Result database");
	tree = new QTreeWidget;
	tree->setColumnCount(1);
	tree->setContextMenuPolicy(Qt::CustomContextMenu);
	roots[PART_ROOT] = new QTreeWidgetItem(tree);
	roots[SENSOR_ROOT] = new QTreeWidgetItem(tree);
	roots[PMASS_ROOT] = new QTreeWidgetItem(tree);
	roots[REACTION_ROOT] = new QTreeWidgetItem(tree);
	//roots[CUBE_ROOT] = new QTreeWidgetItem(vtree);
	//roots[CYLINDER_ROOT] = new QTreeWidgetItem(vtree);
	//roots[POLYGON_ROOT] = new QTreeWidgetItem(vtree);
	//roots[MASS_ROOT] = new QTreeWidgetItem(vtree);

	roots[PART_ROOT]->setText(0, "Part");
	roots[SENSOR_ROOT]->setText(0, "Sensor");
	roots[PMASS_ROOT]->setText(0, "Mass");
	roots[REACTION_ROOT]->setText(0, "Reaction");
	/*roots[CUBE_ROOT]->setText(0, "Cube");
	roots[CYLINDER_ROOT]->setText(0, "Cylinder");
	roots[POLYGON_ROOT]->setText(0, "Polygon");
	roots[MASS_ROOT]->setText(0, "Mass");
*/
	roots[PART_ROOT]->setIcon(0, QIcon(":/Resources/parts-256.png"));
	roots[SENSOR_ROOT]->setIcon(0, QIcon(":/Resources/Sensor_icon.png"));
	roots[PMASS_ROOT]->setIcon(0, QIcon(":/Resources/mass.png"));
	roots[REACTION_ROOT]->setIcon(0, QIcon(":/Resources/revolute.png"));
	/*roots[CUBE_ROOT]->setIcon(0, QIcon(":/Resources/pRec.png"));
	roots[CYLINDER_ROOT]->setIcon(0, QIcon(":/Resources/cylinder.png"));
	roots[POLYGON_ROOT]->setIcon(0, QIcon(":/Resources/icPolygon.png"));
	roots[MASS_ROOT]->setIcon(0, QIcon(":/Resources/mass.png"));*/
	plot_item = new QComboBox(this);
	connect(tree, &QTreeWidget::customContextMenuRequested, this, &resultDatabase::contextMenu);
	connect(tree, &QTreeWidget::itemClicked, this, &resultDatabase::clickItem);
	
	tree->setSelectionMode(QAbstractItemView::SelectionMode::ContiguousSelection);
	setWidget(tree);
}

// resultDatabase::~resultDatabase()
// {
// 	
// }

void resultDatabase::contextMenu(const QPoint& pos)
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

void resultDatabase::clickItem(QTreeWidgetItem* item, int col)
{
	if (!item)
		return;
	QList<QTreeWidgetItem*> items = tree->selectedItems();
	cmdWindow::write(CMD_INFO, "The number of selected items : " + QString("%1").arg(items.size()));
	if (items.size() > 1)
	{
		foreach(QTreeWidgetItem* it, items)
		{
			sLists.push_back(it->text(col));
		}
		return;
	}
	QString tg = item->text(col);
	if (target == tg)
		return;
	target = tg;
	QTreeWidgetItem* parent = item->parent();
	if (!parent)
	{
		return;
	}
	QString sp = parent->text(0);
	if (sp == "Part")
	{

	}
	else if (sp == "Sensor")
	{
		tSelected = SENSOR_ROOT;
		sLists.push_back(target);
		plot_item->clear();
		sensor* s = sph_model::SPHModel()->Sensors()[target];
		switch (s->sensorType())
		{
		case sensor::WAVE_HEIGHT_SENSOR:
			plot_item->addItem("Wave height");
			break;
		}
	}
	else if (sp == "Mass")
	{
		tSelected = PMASS_ROOT;
		plot_item->clear();
		//plot_item->addItem("PX");
		plot_item->addItems(getPMResultString());
	}
	else if (sp == "Reaction")
	{
		tSelected = REACTION_ROOT;
		plot_item->clear();
		plot_item->addItems(getRFResultString());
	}

}

resultDatabase::~resultDatabase()
{
	if (tree) delete tree; tree = NULL;
}

void resultDatabase::addChild(tRoot tr, QString& _nm)
{
	QTreeWidgetItem* child = new QTreeWidgetItem();
	child->setText(0, _nm);
	roots[tr]->addChild(child);
	//cmdWindow::write(CMD_INFO, "Added root item : " + _nm);
}

void resultDatabase::bindItemComboBox(QComboBox* t)
{
	plot_item = t;
}

QStringList resultDatabase::selectedLists()
{
	return sLists;
}

resultDatabase::tRoot resultDatabase::selectedType()
{
	return tSelected;
}

QComboBox* resultDatabase::plotItemComboBox()
{
	return plot_item;
}

QString resultDatabase::plotTarget()
{
	return target;
}
