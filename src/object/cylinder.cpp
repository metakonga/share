#include "cylinder.h"

cylinder::cylinder()
	: pointMass()
	, dci(NULL)
	, br(0)
	, tr(0)
{

}

cylinder::cylinder(QString& _name, geometry_use _roll)
	: pointMass(_name, CYLINDER, _roll)
	, dci(NULL)
	, br(0)
	, tr(0)
{

}

cylinder::cylinder(const cylinder& _cy)
	: org(_cy.origin())
	, br(_cy.baseRadius())
	, tr(_cy.topRadius())
	, len(_cy.length())
	, bpos(_cy.basePos())
	, tpos(_cy.topPos())
	, pointMass(_cy)
	, ep(_cy.orientation())
	, dci(NULL)
{
	
}

cylinder::~cylinder()
{
	
}
//
//void cylinder::update(pointMass* pm)
//{
//// 	mass* m = object::pointMass();
//// 	ep = m->getEP();
//// 	EPD ev = m->getEV(); 
//// 	org = m->Position();
//// 	bpos = org + m->toGlobal(loc_bpos);
//// 	tpos = org + m->toGlobal(loc_tpos);
//// 	if (dci){
//// 		device_cylinder_info* h_dci = new device_cylinder_info;
//// 		checkCudaErrors(cudaMemcpy(h_dci, dci, sizeof(device_cylinder_info), cudaMemcpyDeviceToHost));
//// 		h_dci->origin = make_double3(org.x, org.y, org.z);
//// 		h_dci->vel = make_double3(m->getVelocity().x, m->getVelocity().y, m->getVelocity().z);
//// 		h_dci->omega = make_double3(
//// 			2.0*(-ep.e1 * ev.e0 + ep.e0 * ev.e1 - ep.e3 * ev.e2 + ep.e2 * ev.e3), 
//// 			2.0*(-ep.e2 * ev.e0 + ep.e3 * ev.e1 + ep.e0 * ev.e2 - ep.e1 * ev.e3), 
//// 			2.0*(-ep.e3 * ev.e0 - ep.e2 * ev.e1 + ep.e1 * ev.e2 + ep.e0 * ev.e3));
//// 		h_dci->ep = make_double4(ep.e0, ep.e1, ep.e2, ep.e3);
////  		h_dci->pbase = make_double3(bpos.x, bpos.y, bpos.z);
////  		h_dci->ptop = make_double3(tpos.x, tpos.y, tpos.z);
//// 		checkCudaErrors(cudaMemcpy(dci, h_dci, sizeof(device_cylinder_info), cudaMemcpyHostToDevice));
//// 		delete h_dci; h_dci = NULL;
//// 	}
//}

// void cylinder::cuAllocData(unsigned int _np)
// { 
// 	device_cylinder_info *_dci = new device_cylinder_info;
// 	_dci->len = len;
// 	_dci->rbase = br;
// 	_dci->rtop = tr;
// 	_dci->pbase = make_double3(bpos.x, bpos.y, bpos.z);
// 	_dci->ptop = make_double3(tpos.x, tpos.y, tpos.z);
// 	//_dci->obase = make_float4(b_ep.e0, b_ep.e1, b_ep.e2, b_ep.e3);
// 	//_dci->otop = make_float4(t_ep.e0, t_ep.e1, t_ep.e2, t_ep.e3);
// 	_dci->ep = make_double4(ep.e0, ep.e1, ep.e2, ep.e3);
// 	_dci->origin = make_double3(org.x, org.y, org.z);
// 	checkCudaErrors(cudaMalloc((void**)&dci, sizeof(device_cylinder_info)));
// 	checkCudaErrors(cudaMemcpy(dci, _dci, sizeof(device_cylinder_info), cudaMemcpyHostToDevice));
// 	delete _dci;
// }

bool cylinder::define(double _br, double _tr, VEC3D _bpos, VEC3D _tpos)
{
	br = _br;
	tr = _tr;
	//len = _len;
	//org = _org;
	bpos = _bpos;
	tpos = _tpos;
	ep = EPD(1.0, 0.0, 0.0, 0.0);
	len = (tpos - bpos).length();
	org = 0.5 * (tpos + bpos);
	loc_tpos = ep.toLocal(tpos - org);/* VEC3D(0, -len * 0.5, 0);*/
	loc_bpos = ep.toLocal(bpos - org);// VEC3D(0, len * 0.5, 0);
	//p1 = org + 0.5f * len * dpos;
	//p2 = org - 0.5f * len * dpos;

// 	float h_len = len * 0.5f;
// 	float th = M_PI * 0.5f;
// 	float ap = acos((tpos.z - org.z) / h_len);
// 	float xi = asin((org.y - tpos.y) / (sin(ap) * h_len));
// 	euler.x = xi / M_PI;
// 	euler.y = th / M_PI;
	// 	euler.z = ap / M_PI;
// 	VEC3D bo = bpos - org;
// 	VEC3D u = bo / bo.length();
// 	VEC3D bv = bo.length() * 1.1 * u;

	//loc_bpos = len * 0.5 * u;
	//loc_tpos = -len * 0.5 * u;
// 	double h_len = len * 0.5;
// 	double dth = M_PI * 0.5;
// 	double th = M_PI * 0.5;
// 	double ap = acos(u.z);
// 	double xi = asin(-u.y);
// 	b_ep.setFromEuler(xi, th, ap);
// 	double _ap = ap + M_PI;
//  	if (_ap > M_PI)
//  		ap = ap - M_PI;
// 	t_ep.setFromEuler(xi, th, ap);
// 	ep = t_ep;
	
	return true;
}

//unsigned int cylinder::makeParticles(double rad, VEC3UI &_size, VEC3D& _spacing, unsigned int nstack, bool isOnlyCount, VEC4D_PTR pos, unsigned int sid)
// {
// 	unsigned int np = 0;
// 	return np;
// }

void cylinder::setOrientation(double e1, double e2, double e3)
{
	ep.setFromEuler(e1, e2, e3);
}
// std::fstream& cube::operator<<(std::fstream& oss)
//{
//
//}