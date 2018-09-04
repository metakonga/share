#include "polygonObject.h"
#include "mphysics_cuda_dec.cuh"
#include "mass.h"

#include <QFile>
#include <QTextStream>

unsigned int polygonObject::nPolygonObject = 0;

polygonObject::polygonObject()
	: pointMass()
	, vertexList(NULL)
	, indexList(NULL)
//	, local_vertice(NULL)
// 	, h_poly(NULL)
// 	, h_sph(NULL)
// 	, d_sph(NULL)
// //	, h_mass(NULL)
// 	//, d_mass(NULL)
// 	, vertice(NULL)
// 	, indice(NULL)
// 	, normals(NULL)
// 	, id(NULL)
// 	, nvertex(0)
// 	, nindex(0)
	, ntriangle(0)
	, maxRadii(0)
	, filePath("")
{

}

polygonObject::polygonObject(QString _name, geometry_use _roll)
	: pointMass(_name, POLYGON_SHAPE, _roll)
	, vertexList(NULL)
	, indexList(NULL)
//	, local_vertice(NULL)
// 	, h_poly(NULL)
// 	, h_sph(NULL)
// // 	, d_sph(NULL)
// 	, vertice(NULL)
// 	, normals(NULL)
// 	, id(NULL)
// //	, h_mass(NULL)
// 	//, d_mass(NULL)
// 	, indice(NULL)
// 	, nvertex(0)
/*	, nindex(0)*/
	, maxRadii(0)
	, ntriangle(0)
	, filePath("")
{
// 	int begin = file.lastIndexOf("/") + 1;
// 	int end = file.lastIndexOf(".");
// 	QString _nm = file.mid(begin, end - begin);
// 	name = _nm;
}

polygonObject::polygonObject(const polygonObject& _poly)
	: pointMass(_poly)
	, vertexList(NULL)
	, indexList(NULL)
//	, local_vertice(NULL)
// 	, h_poly(NULL)
// 	, h_sph(NULL)
//	, d_sph(NULL)
//	, vertice(NULL)
//	, indice(NULL)
//	, normals(NULL)
//	, id(NULL)
//	, h_mass(NULL)
//	//, d_mass(NULL)
// 	, nvertex(_poly.numVertex())
// 	, nindex(_poly.numIndex())
	, ntriangle(_poly.NumTriangle())
	, maxRadii(_poly.maxRadius())
	, filePath(_poly.meshDataFile())
	//, org(_poly.getOrigin())
{
// 	if (_poly.hostMassInfo())
// 	{
// 		h_mass = new host_polygon_mass_info;
// 		memcpy(h_mass, _poly.hostMassInfo(), sizeof(host_polygon_mass_info));
// 	}
// 	if (_poly.hostPolygonInfo())
// 	{
// 		h_poly = new host_polygon_info[nindex];
// 		memcpy(h_poly, _poly.hostPolygonInfo(), sizeof(host_polygon_info) * nindex);
// 	}
// 	if (_poly.hostSphereSet())
// 	{
// 		h_sph = new VEC4D[ntriangle];
// 		memcpy(h_sph, _poly.hostSphereSet(), sizeof(VEC4D) * ntriangle);
// 	}
// 	if (_poly.vertexSet())
// 	{
// 		vertice = new VEC3D[nvertex];
// 		memcpy(vertice, _poly.vertexSet(), sizeof(VEC3D) * nvertex);
// 	}
// 	if (_poly.indexSet())
// 	{
// 		indice = new VEC3UI[nindex];
// 		memcpy(indice, _poly.indexSet(), sizeof(VEC3UI) * nindex);
// 	}
// 	if (_poly.devicePolygonInfo())
// 	{
// 		checkCudaErrors(cudaMalloc((void**)&d_poly, sizeof(device_polygon_info) * nindex));
// 		checkCudaErrors(cudaMemcpy(d_poly,  _poly.devicePolygonInfo(), sizeof(device_polygon_info) * nindex, cudaMemcpyDeviceToDevice));
// 	}
// 	if (_poly.deviceSphereSet())
// 	{
// 		checkCudaErrors(cudaMalloc((void**)&d_sph, sizeof(VEC4D) * nindex));
// 		checkCudaErrors(cudaMemcpy(d_sph, _poly.deviceSphereSet(), sizeof(VEC4D) * nindex, cudaMemcpyDeviceToDevice));
// 	}
}

polygonObject::~polygonObject()
{
// 	if (vertice) delete[] vertice; vertice = NULL;
// 	if (indice) delete[] indice; indice = NULL;
// 	if (h_poly) delete[] h_poly; h_poly = NULL;
// 	if (h_sph) delete[] h_sph; h_sph = NULL;
// 	if (local_vertice) delete[] local_vertice; local_vertice = NULL;
// 	if (normals) delete[] normals; normals = NULL;
// 	if (id) delete[] id; id = NULL;
	//if (h_mass) delete h_mass; h_mass = NULL;
	//if (d_mass) checkCudaErrors(cudaFree(h_mass)); h_mass = NULL;
// 	if (d_sph) checkCudaErrors(cudaFree(d_sph)); d_sph = NULL;
	//if (d_poly) checkCudaErrors(cudaFree(d_poly)); d_poly = NULL;
	nPolygonObject--;
}

bool polygonObject::define(import_shape_type t, int _ntriangle, double* vList, unsigned int *iList)
{
	switch (t)
	{
	case MILKSHAPE_3D_ASCII: _fromMS3DASCII(_ntriangle, vList, iList); break;
	case STL_ASCII: _fromSTLASCII(_ntriangle, vList); break;
	}
	return true;
}

void polygonObject::_fromMS3DASCII(int _ntriangle, double* vList, unsigned int *iList)
{
	ntriangle = _ntriangle;
	//	unsigned int a, b, c;
	vertexList = vList;
	indexList = iList;
	//	h_poly = new host_polygon_info[ntriangle];
	//	h_sph = new VEC4D[ntriangle];
	/// 	VEC3D com;
	unsigned int a, b, c;
	VEC3D P, Q, R, V, W, N;
	VEC3D com;
	for (unsigned int i = 0; i < ntriangle; i++)
	{
		a = iList[i * 3 + 0];
		b = iList[i * 3 + 1];
		c = iList[i * 3 + 2];
		P = VEC3D(vList[a * 3 + 0], vList[a * 3 + 1], vList[a * 3 + 2]);
		Q = VEC3D(vList[b * 3 + 0], vList[b * 3 + 1], vList[b * 3 + 2]);
		R = VEC3D(vList[c * 3 + 0], vList[c * 3 + 1], vList[c * 3 + 2]);
		V = Q - P;
		W = R - P;
		N = V.cross(W);
		N = N / N.length();
		VEC3D M1 = (Q + P) / 2;
		VEC3D M2 = (R + P) / 2;
		VEC3D D1 = N.cross(V);
		VEC3D D2 = N.cross(W);
		// 		if (N.z < 0.0)
		// 		{
		// 			bool p = true;
		// 		}
		double t = 0;
		if (abs(D1.x*D2.y - D1.y*D2.x) > 1E-13)
		{
			t = (D2.x*(M1.y - M2.y)) / (D1.x*D2.y - D1.y*D2.x) - (D2.y*(M1.x - M2.x)) / (D1.x*D2.y - D1.y*D2.x);
		}
		else if (abs(D1.x*D2.z - D1.z*D2.x) > 1E-13)
		{
			t = (D2.x*(M1.z - M2.z)) / (D1.x*D2.z - D1.z*D2.x) - (D2.z*(M1.x - M2.x)) / (D1.x*D2.z - D1.z*D2.x);
		}
		else if (abs(D1.y*D2.z - D1.z*D2.y) > 1E-13)
		{
			t = (D2.y*(M1.z - M2.z)) / (D1.y*D2.z - D1.z*D2.y) - (D2.z*(M1.y - M2.y)) / (D1.y*D2.z - D1.z*D2.y);
		}
		VEC3D Ctri = M1 + t * D1;
		VEC4D sph;
		double rad = (Ctri - P).length();
		com += Ctri;
		// 		while (abs(fc - ft) > 0.00001)
		// 		{
		// 			d = ft * sph.w;
		// 			double p = d / po.N.length();
		// 			VEC3D _c = Ctri - p * po.N;
		// 			sph.x = _c.x; sph.y = _c.y; sph.z = _c.z;
		// 			sph.w = (_c - po.P).length();
		// 			fc = d / sph.w;
		// 		}
		if (rad > maxRadii)
			maxRadii = rad;
	}
	pointMass::pos = com / ntriangle;
	// 	QFile qf(file);
	// 	qf.open(QIODevice::ReadOnly);
	// 	QTextStream qts(&qf);
	// 	switch (tm)
	// 	{
	// 	case MILKSHAPE_3D_ASCII:
	// 		fromMS3DASCII(qts);
	// 		break;
	// 	}
	nPolygonObject++;
// 	QString ch;
// 	//unsigned int nvertex = 0;
// 	//unsigned int npoly = 0;
// 	qs >> ch >> ch >> ch >> ch >> ch >> ch >> ch >> ch >> ch >> ch;
// 	unsigned int ui;
// 	qs >> ch;
// 	int begin = ch.indexOf("\"");
// 	int end = ch.lastIndexOf("\"");
// 	QString _name = ch.mid(begin + 1, end - 1);
// 	qs >> ch >> ch >> nvertex;
// 	vertice = new VEC3D[nvertex];
// 	local_vertice = new VEC3D[nvertex];
// 	id = new unsigned int[nvertex];
// 	double x_max = FLT_MIN; double x_min = FLT_MAX;
// 	double y_max = FLT_MIN; double y_min = FLT_MAX;
// 	double z_max = FLT_MIN; double z_min = FLT_MAX;
// 	
// 	unsigned int real_nvert = 1;
// 	id[0] = 0;
// 	qs >> ch >> vertice[0].x >> vertice[0].y >> vertice[0].z >> ch >> ch >> ch;
// 	for (unsigned int i = 1; i < nvertex; i++){
// 		bool isBreak = false;
// 		double x, y, z;
// 		qs >> ch >> x >> y >> z >> ch >> ch >> ch;
// 		for (unsigned int j = 0; j < real_nvert; j++)
// 		{
// 			if (vertice[j].x == x && vertice[j].y == y && vertice[j].z == z)
// 			{
// 				id[i] = j;
// 				isBreak = true;
// 			}
// 				
// 		}
// 		if (isBreak)
// 			continue;
// 		id[i] = real_nvert;
// 		vertice[real_nvert] = VEC3D(x, y, z);
// 		unsigned int k = real_nvert;
// 		if (x_min > vertice[k].x) x_min = vertice[k].x;
// 		if (x_max < vertice[k].x) x_max = vertice[k].x;
// 		if (y_min > vertice[k].y) y_min = vertice[k].y;
// 		if (y_max < vertice[k].y) y_max = vertice[k].y;
// 		if (z_min > vertice[k].z) z_min = vertice[k].z;
// 		if (z_max < vertice[k].z) z_max = vertice[k].z;
// 		real_nvert++;
// 	}
// 	com = VEC3D(0.5 * (x_min + x_max), 0.5 * (y_min + y_max), 0.5 * (z_min + z_max));
// 	com = 0.001 * com;
// 	for (unsigned int i = 0; i < real_nvert; i++)
// 	{
// 		vertice[i].x *= 0.001;
// 		vertice[i].y *= 0.001;
// 		vertice[i].z *= 0.001;
// 		local_vertice[i] = vertice[i] - com;
// 	}
// 	
// 	qs >> ui;
// 	VEC3D* nid = new VEC3D[ui];
// 	for (unsigned int i = 0; i < ui; i++)
// 		qs >> nid[i].x >> nid[i].y >> nid[i].z;
// 	qs >> nindex;
// 	indice = new VEC3UI[nindex];
// 	normals = new VEC3D[nindex];
// 	//unsigned int *num_id = new unsigned int[real_nvert];
// //	memset(num_id, 0, sizeof(unsigned int) * real_nvert);
// 	memset(normals, 0, sizeof(double) * real_nvert * 3);
// 	unsigned int x, y, z;
// 	VEC3D P, Q, R, V, W, N;
// 	QMap<unsigned int, QMap<unsigned int, unsigned int>> clist;
// 	for (unsigned int i = 0; i < nindex; i++){
// 		qs >> ch >> indice[i].x >> indice[i].y >> indice[i].z >> x >> y >> z >> ch;
// 		unsigned int ix, iy, iz;
// 		ix = id[indice[i].x];
// 		iy = id[indice[i].y];
// 		iz = id[indice[i].z];
// 		indice[i].x = ix;
// 		indice[i].y = iy;
// 		indice[i].z = iz;
// // 		P = vertice[ix];
// // 		Q = vertice[iy];
// // 		R = vertice[iz];
// // 		V = Q - P;
// // 		W = R - P;
// // 		N = V.cross(W);
// // 		N = N / N.length();
// 		QMap<unsigned int, unsigned int> *cx = &clist[ix];
// 		QMap<unsigned int, unsigned int> *cy = &clist[iy];
// 		QMap<unsigned int, unsigned int> *cz = &clist[iz];
// 		(*cx)[x] += 1;
// 		(*cy)[y] += 1;
// 		(*cz)[z] += 1;
// 		
// // 		normals[ix/*indice[i].x*/] += nid[x];
// // 		normals[iy/*indice[i].y*/] += nid[y];
// // 		normals[iz/*indice[i].z*/] += nid[z];
// // 		num_id[ix]++;
// // 		num_id[iy]++;
// // 		num_id[iz]++;
// // 		normals[indice[i].y] += nid[y];
// // 		normals[indice[i].z] += nid[z];
// 	}
// 	QMapIterator<unsigned int, QMap<unsigned int, unsigned int>> it(clist);
// 	while (it.hasNext())
// 	{
// 		it.next();
// 		unsigned int k = it.key();
// 		QMap<unsigned int, unsigned int> c = it.value();
// 		QList<unsigned int> v = c.keys();
// 		foreach(unsigned n, v)
// 		{
// 			normals[k] += nid[n];
// 		}
// 		normals[k] = normals[k] / normals[k].length();
// 	}
// // 	for (unsigned int i = 0; i < real_nvert; i++)
// // 	{
// // 		normals[i] = normals[i] / normals[i].length();
// // 	}
// 	delete[] nid;
// //	h_poly = new host_polygon_info[nindex];
// 	//h_loc_poly = new host_polygon_info[nindex];
// 	h_sph = new VEC4D[nindex];
// 	//h_loc_sph = new VEC4D[nindex];
// 	double fc = 0;
// 	double ft = 0.7;
// 	for (unsigned int i = 0; i < nindex; i++){
// 		fc = 0;
// 		host_polygon_info po;
// 		po.P = vertice[ indice[i].x ];
// 		po.Q = vertice[ indice[i].y ];
// 		po.R = vertice[ indice[i].z ];
// 		po.V = po.Q - po.P;
// 		po.W = po.R - po.P;
// 		po.N = po.V.cross(po.W);
// 		po.N = po.N / po.N.length();
// //		h_poly[i] = po;
// 		VEC3D M1 = (po.Q + po.P) / 2;
// 		VEC3D M2 = (po.R + po.P) / 2;
// 		VEC3D D1 = po.N.cross(po.V);
// 		VEC3D D2 = po.N.cross(po.W);
// 		double t = 0;//t = (D2.x*(M1.y - M2.y)) / (D1.x*D2.y - D1.y*D2.x) - (D2.y*(M1.x - M2.x)) / (D1.x*D2.y - D1.y*D2.x);
// 		if (abs(D1.x*D2.y - D1.y*D2.x) > 1E-13)
// 		{
// 			t = (D2.x*(M1.y - M2.y)) / (D1.x*D2.y - D1.y*D2.x) - (D2.y*(M1.x - M2.x)) / (D1.x*D2.y - D1.y*D2.x);
// 		}
// 		else if (abs(D1.x*D2.z - D1.z*D2.x) > 1E-13)
// 		{
// 			t = (D2.x*(M1.z - M2.z)) / (D1.x*D2.z - D1.z*D2.x) - (D2.z*(M1.x - M2.x)) / (D1.x*D2.z - D1.z*D2.x);
// 		}
// 		else if (abs(D1.y*D2.z - D1.z*D2.y) > 1E-13)
// 		{
// 			t = (D2.y*(M1.z - M2.z)) / (D1.y*D2.z - D1.z*D2.y) - (D2.z*(M1.y - M2.y)) / (D1.y*D2.z - D1.z*D2.y);
// 		}
// 		VEC3D Ctri = M1 + t * D1;
// 		VEC4D sph;
// 		sph.w = (Ctri - po.P).length();
// 		sph.x = Ctri.x; sph.y = Ctri.y; sph.z = Ctri.z;
// // 		while (abs(fc - ft) > 0.00001)
// // 		{
// // 			d = ft * sph.w;
// // 			double p = d / po.N.length();
// // 			VEC3D _c = Ctri - p * po.N;
// // 			sph.x = _c.x; sph.y = _c.y; sph.z = _c.z;
// // 			sph.w = (_c - po.P).length();
// // 			fc = d / sph.w;
// // 		}
// 		if (sph.w > maxRadii)
// 			maxRadii = sph.w;
// 		h_sph[i] = sph;
// 	}
// // 	for (unsigned int i = 0; i < real_nvert; i++){
// // 		vertice[i] = vertice[i] - org;
// // 	}
// 	nvertex = real_nvert;
// // 	h_mass = new host_polygon_mass_info;
// // 	h_mass->origin = org;
// // 	h_mass->omega = 0;
// // 	h_mass->vel = 0;
// // 	h_mass->ep = EPD(1.0, 0, 0, 0);
}

void polygonObject::_fromSTLASCII(int _ntriangle, double* vList)
{
	ntriangle = _ntriangle;
	vertexList = vList;
	unsigned int a, b, c;
	VEC3D P, Q, R, V, W, N;
	VEC3D com;
	for (unsigned int i = 0; i < ntriangle; i++)
	{
		P = VEC3D(vList[i * 9 + 0], vList[i * 9 + 1], vList[i * 9 + 2]);
		Q = VEC3D(vList[i * 9 + 3], vList[i * 9 + 4], vList[i * 9 + 5]);
		R = VEC3D(vList[i * 9 + 6], vList[i * 9 + 7], vList[i * 9 + 8]);
		V = Q - P;
		W = R - P;
		N = V.cross(W);
		N = N / N.length();
		VEC3D M1 = (Q + P) / 2;
		VEC3D M2 = (R + P) / 2;
		VEC3D D1 = N.cross(V);
		VEC3D D2 = N.cross(W);
		double t = 0;
		if (abs(D1.x*D2.y - D1.y*D2.x) > 1E-13)
		{
			t = (D2.x*(M1.y - M2.y)) / (D1.x*D2.y - D1.y*D2.x) - (D2.y*(M1.x - M2.x)) / (D1.x*D2.y - D1.y*D2.x);
		}
		else if (abs(D1.x*D2.z - D1.z*D2.x) > 1E-13)
		{
			t = (D2.x*(M1.z - M2.z)) / (D1.x*D2.z - D1.z*D2.x) - (D2.z*(M1.x - M2.x)) / (D1.x*D2.z - D1.z*D2.x);
		}
		else if (abs(D1.y*D2.z - D1.z*D2.y) > 1E-13)
		{
			t = (D2.y*(M1.z - M2.z)) / (D1.y*D2.z - D1.z*D2.y) - (D2.z*(M1.y - M2.y)) / (D1.y*D2.z - D1.z*D2.y);
		}
		VEC3D Ctri = M1 + t * D1;
		VEC4D sph;
		double rad = (Ctri - P).length();
		com += Ctri;
		if (rad > maxRadii)
			maxRadii = rad;
	}
	pointMass::pos = com / ntriangle;
	nPolygonObject++;
}

//void polygonObject::cuAllocData(unsigned int _np)
//{
// 	if (!d_poly)
// 		checkCudaErrors(cudaMalloc((void**)&d_poly, sizeof(device_polygon_info) * nindex));
// 	if (!d_sph)
// 		checkCudaErrors(cudaMalloc((void**)&d_sph, sizeof(double4) * nindex));
// 	if (!d_mass)
// 		checkCudaErrors(cudaMalloc((void**)&d_mass, sizeof(device_polygon_mass_info)));
// 
// 	updateDeviceFromHost();
//}

void polygonObject::updateDeviceFromHost()
{
	//checkCudaErrors(cudaMemcpy(d_poly, h_poly, sizeof(device_polygon_info) * nindex, cudaMemcpyHostToDevice));
	//checkCudaErrors(cudaMemcpy(d_sph, h_sph, sizeof(double) * 4 * nindex, cudaMemcpyHostToDevice));
	//checkCudaErrors(cudaMemcpy(d_mass, h_mass, sizeof(device_polygon_mass_info), cudaMemcpyHostToDevice));
}

//void polygonObject::update(pointMass* pm)
//{
//// 	mass* m = object::pointMass();
//// 	VEC3D vel = m->getVelocity();
//// 	EPD ep = m->getEP();
//// 	EPD ev = m->getEV();
//// 	org = m->getPosition();
//// 	VEC3D *g_vertice = new VEC3D[nvertex];
//// 	for(unsigned int i = 0; i < nvertex; i++){
//// 		g_vertice[i] = org + m->toGlobal(vertice[i]);
//// 	}
//// 	for (unsigned int i = 0; i < nindex; i++){
//// 		//fc = 0;
//// 		//VEC3D psph = org + m->toGlobal(VEC3D(h_sph[i].x, h_sph[i].y, h_sph[i].z));
//// 		host_polygon_info po;
//// 		po.P = g_vertice[indice[i].x];
//// 		po.Q = g_vertice[indice[i].y];
//// 		po.R = g_vertice[indice[i].z];
//// 		po.V = po.Q - po.P;
//// 		po.W = po.R - po.P;
//// 		po.N = po.V.cross(po.W);
//// 		po.N = po.N / po.N.length();
//// 		h_poly[i] = po;
//// 		VEC3D M1 = (po.Q + po.P) / 2;
//// 		VEC3D M2 = (po.R + po.P) / 2;
//// 		VEC3D D1 = po.N.cross(po.V);
//// 		VEC3D D2 = po.N.cross(po.W);
//// 		double t = 0;//t = (D2.x*(M1.y - M2.y)) / (D1.x*D2.y - D1.y*D2.x) - (D2.y*(M1.x - M2.x)) / (D1.x*D2.y - D1.y*D2.x);
//// 		if (abs(D1.x*D2.y - D1.y*D2.x) > 1E-13)
//// 		{
//// 			t = (D2.x*(M1.y - M2.y)) / (D1.x*D2.y - D1.y*D2.x) - (D2.y*(M1.x - M2.x)) / (D1.x*D2.y - D1.y*D2.x);
//// 		}
//// 		else if (abs(D1.x*D2.z - D1.z*D2.x) > 1E-13)
//// 		{
//// 			t = (D2.x*(M1.z - M2.z)) / (D1.x*D2.z - D1.z*D2.x) - (D2.z*(M1.x - M2.x)) / (D1.x*D2.z - D1.z*D2.x);
//// 		}
//// 		else if (abs(D1.y*D2.z - D1.z*D2.y) > 1E-13)
//// 		{
//// 			t = (D2.y*(M1.z - M2.z)) / (D1.y*D2.z - D1.z*D2.y) - (D2.z*(M1.y - M2.y)) / (D1.y*D2.z - D1.z*D2.y);
//// 		}
//// 		VEC3D Ctri = M1 + t * D1;
//// 		VEC4D sph;
//// 		sph.w = (Ctri - po.P).length();
//// 		sph.x = Ctri.x; sph.y = Ctri.y; sph.z = Ctri.z;
//// 		h_sph[i] = sph;
//// 		//h_sph[i] = VEC4D(psph.x, psph.y, psph.z, h_sph[i].w);
//// 	}
//// 	h_mass->origin = org;
//// 	h_mass->vel = vel;
//// 	h_mass->omega = VEC3D(
//// 		2.0*(-ep.e1 * ev.e0 + ep.e0 * ev.e1 - ep.e3 * ev.e2 + ep.e2 * ev.e3),
//// 		2.0*(-ep.e2 * ev.e0 + ep.e3 * ev.e1 + ep.e0 * ev.e2 - ep.e1 * ev.e3),
//// 		2.0*(-ep.e3 * ev.e0 - ep.e2 * ev.e1 + ep.e1 * ev.e2 + ep.e0 * ev.e3));
//// 	h_mass->ep = ep;
//// 	updateDeviceFromHost();
//// 	delete[] g_vertice; g_vertice = NULL;
//}