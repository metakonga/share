#include "plane.h"
#include "mphysics_cuda_dec.cuh"
#include "pointMass.h"

plane::plane()
	:object()
	, l1(0)
	, l2(0)
{

}

plane::plane(QString& _name,  geometry_use _roll)
	: object(_name, PLANE, _roll)
	, l1(0)
	, l2(0)
{

}

plane::plane(const plane& _plane)
	: object(_plane)
	, l1(_plane.L1())
	, l2(_plane.L2())
	, xw(_plane.XW())
	, uw(_plane.UW())
	, u1(_plane.U1())
	, u2(_plane.U2())
	, pa(_plane.PA())
	, pb(_plane.PB())
	, w2(_plane.W2())
	, w3(_plane.W3())
	, w4(_plane.W4())
{

}

plane::~plane()
{
//	if (dpi) checkCudaErrors(cudaFree(dpi));
}

bool plane::define(VEC3D& _xw, VEC3D& _pa, VEC3D& _pc, VEC3D& _pb)
{
	w2 = _pa; 
	minp.x = w2.x < minp.x ? w2.x : minp.x; minp.y = w2.y < minp.y ? w2.y : minp.y; minp.z = w2.z < minp.z ? w2.z : minp.z;
	maxp.x = w2.x < maxp.x ? w2.x : maxp.x; maxp.y = w2.y < maxp.y ? w2.y : maxp.y; maxp.z = w2.z < maxp.z ? w2.z : maxp.z;
	w3 = _pc;
	minp.x = w3.x < minp.x ? w3.x : minp.x; minp.y = w3.y < minp.y ? w3.y : minp.y; minp.z = w3.z < minp.z ? w3.z : minp.z;
	maxp.x = w3.x < maxp.x ? w3.x : maxp.x; maxp.y = w3.y < maxp.y ? w3.y : maxp.y; maxp.z = w3.z < maxp.z ? w3.z : maxp.z;
	w4 = _pb;
	minp.x = w4.x < minp.x ? w4.x : minp.x; minp.y = w4.y < minp.y ? w4.y : minp.y; minp.z = w4.z < minp.z ? w4.z : minp.z;
	maxp.x = w4.x < maxp.x ? w4.x : maxp.x; maxp.y = w4.y < maxp.y ? w4.y : maxp.y; maxp.z = w4.z < maxp.z ? w4.z : maxp.z;

	xw = _xw;
	minp.x = xw.x < minp.x ? xw.x : minp.x; minp.y = xw.y < minp.y ? xw.y : minp.y; minp.z = xw.z < minp.z ? xw.z : minp.z;
	maxp.x = xw.x < maxp.x ? xw.x : maxp.x; maxp.y = xw.y < maxp.y ? xw.y : maxp.y; maxp.z = xw.z < maxp.z ? xw.z : maxp.z;
	pa = _pa;
	pb = _pb;

	pa -= xw;
	pb -= xw;
	l1 = pa.length();
	l2 = pb.length();
	u1 = pa / l1;
	u2 = pb / l2;
	uw = u1.cross(u2);

	return true;
}

bool plane::define(VEC3D& _xw, VEC3D& _pa, VEC3D& _pb)
{
	xw = _xw;
	pa = _pa;
	pb = _pb;

	pa -= xw;
	pb -= xw;
	l1 = pa.length();
	l2 = pb.length();
	u1 = pa / l1;
	u2 = pb / l2;
	uw = u1.cross(u2);

	return true;
}

void plane::update(pointMass* pm)
{
// 	if (tsd == CPU){
// 		xw += dt * VEC3D(0.0, 0.0, 0.0);	//plane motion setting, m/s, CPU code
// 		w2 += dt * VEC3D(0.0, 0.0, 0.0);
// 		w3 += dt * VEC3D(0.0, 0.0, 0.0);
// 		w4 += dt * VEC3D(0.0, 0.0, 0.0);
// 	}
// 	else if(tsd == GPU){
// 		device_plane_info *h_dpi = new device_plane_info;
// 		checkCudaErrors(cudaMemcpy(h_dpi, dpi, sizeof(device_plane_info), cudaMemcpyDeviceToHost));
// 		// plane motion setting, m/s, GPU code, mde파일에서 plane 설정값중 마지막값 '1'로 바꾸기
// 		h_dpi->xw = make_double3(h_dpi->xw.x + dt * 0.0, h_dpi->xw.y + dt *  0.2, h_dpi->xw.z + dt *  0.0);
// 		h_dpi->w2 = make_double3(h_dpi->w2.x + dt * 0.0, h_dpi->w2.y + dt *  0.2, h_dpi->w2.z + dt *  0.0);
// 		h_dpi->w3 = make_double3(h_dpi->w3.x + dt * 0.0, h_dpi->w3.y + dt *  0.2, h_dpi->w3.z + dt *  0.0);
// 		h_dpi->w4 = make_double3(h_dpi->w4.x + dt * 0.0, h_dpi->w4.y + dt *  0.2, h_dpi->w4.z + dt *  0.0);
// 		checkCudaErrors(cudaMemcpy(dpi, h_dpi, sizeof(device_plane_info), cudaMemcpyHostToDevice));
// 		delete h_dpi;
// 	}
	
}

void plane::Save(QTextStream& ts)
{

	//QTextStream ts(&(md->modelStream()));
//	bool isExistMass = ms ? true : false;

// 	ts << "OBJECT PLANE " << id << " " << name << " " << roll_type << " " << mat_type << " " << (int)_expression << " " << isExistMass << endl
// 		<< xw.x << " " << xw.y << " " << xw.z << endl;
// 	VEC3D ap = xw + pa;
// 	ts << ap.x << " " << ap.y << " " << ap.z << endl;
// 	VEC3D cp = xw + pb + l1 * u1;
// 	ts << cp.x << " " << cp.y << " " << cp.z << endl;
// 	VEC3D bp = xw + pb;
// 	ts << bp.x << " " << bp.y << " " << bp.z << endl;
// 
// 	if (isExistMass)
// 	{
// 		save_mass_data(ts);
// 	}
}