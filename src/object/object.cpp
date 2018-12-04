#include "object.h"
#include "vobject.h"
#include <cuda_runtime.h>

//unsigned int object::sid = 0;
unsigned int object::count = 0;
//std::ofstream object::io_object;

object::object()
	: name("")
	, obj_type(NO_GEOMETRY_TYPE)
	, roll_type(NO_USE_TYPE)
	, mat_type(NO_MATERIAL)
	, dim(DIM3)
	, d(0)
	, y(0)
	, p(0)
	, sm(0)
	, id(0)
	, vol(0)
	, vobj(NULL)
	, marker(NULL)
{
	gmc.enable = false;
}

object::object(QString _name, geometry_type _tobj, geometry_use _roll)
	: name(_name)
	, obj_type(_tobj)
	, roll_type(_roll)
	, mat_type(NO_MATERIAL)
	, dim(DIM3)
	, d(0)
	, y(0)
	, p(0)
	, sm(0)
	, id(0)
	, vol(0)
	, vobj(NULL)
	, marker(NULL)
{
	id = count++;
	gmc.enable = false;
}

object::object(const object& obj)
	: name(obj.Name())
	, id(obj.ID())
	, obj_type(obj.ObjectType())
	, mat_type(obj.MaterialType())
	, roll_type(obj.RollType())
	, d(obj.Density())
	, y(obj.Youngs())
	, p(obj.Poisson())
	, sm(obj.Shear())
	, gmc(obj.MotionCondition())
{

}

object::~object()
{

}

void object::setMotionCondition(double st, double et, double cv, VEC3D u)
{
	gmc.enable = true;
	gmc.st = st;
	gmc.et = et;
	gmc.cv = cv;
	gmc.unit = u;
}

void object::updateView(VEC3D p, VEC3D r)
{
	r.x = RAD2DEG(r.x);
	r.y = RAD2DEG(r.y);
	r.z = RAD2DEG(r.z);
	vobj->updateView(p, r);
	marker->updateView(p, r);
// 	vobj->setInitialPosition(p);
// 	marker->setInitialPosition(p);

// 	vobj->setInitialAngle(r);
// 	marker->setInitialAngle(r);
}

// void object::setMaterial(tMaterial _tm)
// {
// 	mat_type = _tm;
// 	d = material::getDensity(mat_type);
// 	y = material::getYoungs(mat_type);
// 	p = material::getPoisson(mat_type);
// }

void object::setMaterial(material_type _tm, double _y /*= 0*/, double _d /*= 0*/, double _p /*= 0*/, double _s)
{
	mat_type = _tm;
	if (mat_type == USER_INPUT)
	{
		y = _y;
		d = _d;
		p = _p;
		sm = _s;
	}
	else
	{
		cmaterialType cmt = getMaterialConstant(mat_type);
		y = cmt.youngs;
		d = cmt.density;
		p = cmt.poisson;
		sm = cmt.shear;
		//ms = vol * d;
	}

}
