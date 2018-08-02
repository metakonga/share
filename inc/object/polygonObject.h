#ifndef POLYGONOBJECT_H
#define POLYGONOBJECT_H

#include "object.h"

class polygonObject : public object
{
public:
	polygonObject();
	polygonObject(QString file, geometry_use _roll);
	//polygonObject(QString file);
	polygonObject(const polygonObject& _poly);
	virtual ~polygonObject();

	//virtual unsigned int makeParticles(double rad, VEC3UI &_size, VEC3D& spacing, unsigned int nstack, bool isOnlyCount, VEC4D_PTR pos = NULL, unsigned int sid = 0);
	//virtual void cuAllocData(unsigned int _np);
	virtual void update(pointMass* pm);
	virtual void Save(QTextStream& ts);
	bool define(import_shape_type t, QString file);
	void updateDeviceFromHost();

	QString meshDataFile() const { return filePath; }
	unsigned int numVertex() const { return nvertex; }
	unsigned int numIndex() const { return nindex; }
	double maxRadius() const { return maxRadii; }
	VEC3D* vertexSet() const { return vertice; }
	VEC3D* normalSet() const { return normals; }
	VEC3UI* indexSet() const { return indice; }
	host_polygon_mass_info* hostMassInfo() const { return h_mass; }
	//device_polygon_mass_info* deviceMassInfo() const { return d_mass; }
	VEC4D* hostSphereSet() const { return h_sph; }
	double* deviceSphereSet() const { return d_sph; }
	host_polygon_info* hostPolygonInfo() const { return h_poly; }
	//device_polygon_info* devicePolygonInfo() const { return d_poly; }
	VEC3D getOrigin() const { return org; }

private:
	void fromMS3DASCII(QTextStream& file);

private:
	unsigned int nvertex;
	unsigned int nindex;
	double maxRadii;
	unsigned int *id;
	VEC3D org;
	VEC3D *vertice;
	VEC3UI *indice;
	VEC3D *normals;
	VEC4D *h_sph;
	double *d_sph;
	host_polygon_info* h_poly;
//	device_polygon_info* d_poly;
	host_polygon_mass_info* h_mass;
//	device_polygon_mass_info* d_mass;

	QString filePath;
};

#endif