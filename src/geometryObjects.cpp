#include "geometryObjects.h"
#include "glwidget.h"
#include "Objects.h"
#include "database.h"

unsigned int geometryObjects::nCube = 0;
unsigned int geometryObjects::nPlane = 0;
unsigned int geometryObjects::nCylinder = 0;
unsigned int geometryObjects::nPolyObject = 0;

geometryObjects::geometryObjects()
{

}

geometryObjects::~geometryObjects()
{
	qDeleteAll(objs);
}

void geometryObjects::Save(QTextStream& qts)
{
	qts << endl
		<< "GEOMETRY_OBJECTS_DATA" << endl;
	foreach(QString log, logs)
	{
		qts << log;
	}
	qts << "END_DATA" << endl;	
}

void geometryObjects::Open(QTextStream& qts)
{
	QString ch;
	while (ch != "END_DATA")
	{
		qts >> ch;
		if (ch == "OBJECT")
			qts >> ch;
		if (ch == "cube")
		{
			QString _name;
			int mt, gu;
			double e, p, d, s;
			VEC3D start, end;
			qts >> ch >> _name
				>> ch >> mt
				>> ch >> gu
				>> ch >> start.x >> start.y >> start.z
				>> ch >> end.x >> end.y >> end.z
				>> ch >> e >> p >> d >> s;
			makeCube(
				_name, (material_type)mt, (geometry_use)gu,
				start, end, e, p, d, s);
		}
		else if (ch == "plane")
		{
			QString _name;
			int mt, gu;
			double e, p, d, s;
			VEC3D pa, pb, pc, pd;
			qts >> ch >> _name
				>> ch >> mt
				>> ch >> gu
				>> ch >> pa.x >> pa.y >> pa.z
				>> ch >> pb.x >> pb.y >> pb.z
				>> ch >> pc.x >> pc.y >> pc.z
				>> ch >> pd.x >> pd.y >> pd.z
				>> ch >> e >> p >> d >> s;
			makePlane(
				_name, (material_type)mt, (geometry_use)gu,
				pa, pb, pc, pd, e, p, d, s);
		}
		else if (ch == "polygon")
		{
			QString _name, file;
			int mt, gu, ist;
			double e, p, d, s;
			qts >> ch >> _name
				>> ch >> file
				>> ch >> mt
				>> ch >> gu
				>> ch >> ist
				>> ch >> e >> p >> d >> s;
			makePolygonObject(
				(import_shape_type)ist, file, (geometry_use)gu,
				(material_type)mt, e, p, d, s);
		}
	}
}

void geometryObjects::insertObject(object* _obj)
{
	objs[_obj->Name()] = _obj;
}

object* geometryObjects::Object(QString n)
{
	QStringList l = objs.keys();
	QStringList::const_iterator it = qFind(l, n);
	if (it == l.end())
		return NULL;
	return objs[n];
}

cube* geometryObjects::makeCube(
	QString _name, material_type mt, geometry_use gu, 
	VEC3D& start, VEC3D& end, double e, double p, double d, double s)
{
	cube* c = new cube(_name, gu);
	c->define(start, end);
	c->setMaterial(mt, e, d, p, s);
	objs[_name] = c;
	GLWidget::GLObject()->makeCube(c);
	database::DB()->addChild(database::CUBE_ROOT, c->Name());

	QString log;
	QTextStream qts(&log);
	qts << "OBJECT " << "cube" << endl
		<< "NAME " << _name << endl
		<< "MATERIAL_TYPE " << (int)mt << endl
		<< "GEOMETRY_USE " << (int)gu << endl
		<< "START_POINT " << start.x << " " << start.y << " " << start.z << endl
		<< "END_POINT " << end.x << " " << end.y << " " << end.z << endl
		<< "MATERIAL " << e << " " << p << " " << d << " " << s << endl;
	logs[c->Name()] = log;

	return c;
}

plane* geometryObjects::makePlane(
	QString _name, material_type mt, geometry_use gu, 
	VEC3D& pa, VEC3D& pb, VEC3D& pc, VEC3D& pd, 
	double e, double p, double d, double s)
{
	plane* pl = new plane(_name, gu);
	pl->define(pa, pb, pc, pd);
	pl->setMaterial(mt, e, d, p, s);
	objs[_name] = pl;
	GLWidget::GLObject()->makePlane(pl);
	database::DB()->addChild(database::PLANE_ROOT, pl->Name());

	QString log;
	QTextStream qts(&log);
	qts << "OBJECT " << "plane" << endl
		<< "NAME " << _name << endl
		<< "MATERIAL_TYPE " << (int)mt << endl
		<< "GEOMETRY_USE " << (int)gu << endl
		<< "POINT_A " << pa.x << " " << pa.y << " " << pa.z << endl
		<< "POINT_B " << pb.x << " " << pb.y << " " << pb.z << endl
		<< "POINT_C " << pc.x << " " << pc.y << " " << pc.z << endl
		<< "POINT_D " << pd.x << " " << pd.y << " " << pd.z << endl
		<< "MATERIAL " << e << " " << p << " " << d << " " << s << endl;
	logs[pl->Name()] = log;
	return pl;
}

polygonObject* geometryObjects::makePolygonObject
	(import_shape_type t, QString file, geometry_use gu
	,material_type mt, double e, double p, double d, double s)
{
	int begin = file.lastIndexOf("/") + 1;
	int end = file.lastIndexOf(".");
	QString _nm = file.mid(begin, end - begin);
	polygonObject* po = new polygonObject(_nm, gu);
	po->define(t, file);
	po->setMaterial(mt, e, d, p, s);
	objs[_nm] = po;
	GLWidget::GLObject()->makePolygonObject(po);
	database::DB()->addChild(database::POLYGON_ROOT, po->Name());

	QString log;
	QTextStream qts(&log);
	qts << "OBJECT " << "polygon" << endl
		<< "NAME " << _nm << endl
		<< "FILE " << file << endl
		<< "MATERIAL_TYPE " << (int)mt << endl
		<< "GEOMETRY_USE " << (int)gu << endl
		<< "IMPORT_SHAPE_TYPE " << (int)t << endl
		<< "MATERIAL " << e << " " << p << " " << d << " " << s << endl;
	logs[po->Name()] = log;
	return po;
}