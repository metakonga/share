#ifndef CUBE_H
#define CUBE_H

#include "object.h"
#include "plane.h"

class cube : public object
{
public:
	cube(){}
	cube(QString& _name, geometry_use _roll);
	cube(const cube& _cube);
	virtual ~cube();

	virtual void update(pointMass* pm){};
	virtual void save_object_data(QTextStream& ts);

	plane* Planes() const { return planes; }
	device_plane_info* deviceCubeInfo() { return dpi; }

	bool define(VEC3D& min, VEC3D& max);
	VEC3D origin() { return ori; }
	VEC3D origin() const { return ori; }
	VEC3D min_point() { return min_p; }
	VEC3D min_point() const { return min_p; }
	VEC3D max_point() { return max_p; }
	VEC3D max_point() const { return max_p; }
	VEC3D cube_size() { return size; }
	VEC3D cube_size() const { return size; }
	plane planes_data(int i) const { return planes[i]; }

private:
	VEC3D ori;
	VEC3D min_p;
	VEC3D max_p;
	VEC3D size;
	plane *planes;

	device_plane_info *dpi;
};

#endif