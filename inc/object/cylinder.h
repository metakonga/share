#ifndef CYLINDER_H
#define CYLINDER_H

#include "object.h"

struct device_cylinder_info;

class cylinder : public object
{
public:
	cylinder();
	cylinder(QString& _name, geometry_use _roll);
	cylinder(const cylinder& _cube);
	virtual ~cylinder();

	virtual void update(pointMass* pm);
	virtual void Save(QTextStream& ts);

	bool define(double _br, double _tr, VEC3D _bpos, VEC3D _tpos);
	double baseRadius() const { return br; }
	double topRadius() const { return tr; }
	double length() const { return len; }
	VEC3D origin() const { return org; }
	VEC3D basePos() const { return bpos; }
	VEC3D topPos() const { return tpos; }
	void setOrientation(double e1, double e2, double e3);
	EPD orientation() const { return ep; }

private:
	double br;				// base radius
	double tr;				// top radius
	double len;				// length
	VEC3D org;		// origin
	EPD ep;
	VEC3D bpos;
	VEC3D tpos;
	VEC3D loc_bpos;
	VEC3D loc_tpos;

	device_cylinder_info* dci;
	//MAT33F A;
	//VEC3F spos1;
	//VEC3F spos2;
};

// inline
// std::ostream& operator<<(std::ostream& oss, const cube& my){
// 	oss << "OBJECT CUBE " << my.objectName() << " " << my.rolltype() << " " << my.materialType() << std::endl;
// 	oss << my.origin() << std::endl
// 		<< my.min_point() << std::endl
// 		<< my.max_point() << std::endl
// 		<< my.cube_size() << std::endl;
// 	return oss;
// }

#endif