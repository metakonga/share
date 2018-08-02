#ifndef OBJECT_H
#define OBJECT_H

#include "types.h"
#include "algebraMath.h"
#include <QString>
#include <QTextStream>

class pointMass;

class object
{
public:
	object();
	object(QString _name, geometry_type gt, geometry_use gu);
	object(const object& obj);
	virtual ~object();

	virtual void update(pointMass* pm){};
	virtual void Save(QTextStream& ts){};

	unsigned int ID() const { return id; }
	QString Name() const { return name; }
	geometry_type ObjectType() const { return obj_type; }
	material_type MaterialType() const { return mat_type; }
	geometry_use RollType() const { return roll_type; }
	double Density() const { return d; }
	double Youngs() const { return y; }
	double Poisson() const { return p; }
	double Shear() const { return sm; }
	unsigned int particleCount() const { return count; }

	void setRoll(geometry_use tr)  { roll_type = tr; }
	void setID(unsigned int _id) { id = _id; }
	void setMaterial(material_type _tm, double _y = 0, double _d = 0, double _p = 0, double _s = 0);

protected:
	static unsigned int count;
	unsigned int id;	
	QString name;
	geometry_use roll_type;
	geometry_type obj_type;
	material_type mat_type;
	double d;		// density
	double y;		// young's modulus
	double p;		// poisson ratio
	double sm;		// shear modulus

};

#endif