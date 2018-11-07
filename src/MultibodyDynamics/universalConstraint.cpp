#include "universalConstraint.h"

universalConstraint::universalConstraint()
	: kinematicConstraint()
{

}

universalConstraint::universalConstraint(mbd_model* _md, QString& _nm, kinematicConstraint::Type kt,
	pointMass* ip, VEC3D& _spi, VEC3D& _fi, VEC3D& _gi,
	pointMass* jp, VEC3D& _spj, VEC3D& _fj, VEC3D& _gj)
	: kinematicConstraint(_md, _nm, kt, ip, _spi, _fi, _gi, jp, _spj, _fj, _gj)
{

}

universalConstraint::~universalConstraint()
{

}

void universalConstraint::constraintEquation(double m, double* rhs)
{
	if (ib->NumDOF() == DIM2)
	{
		return;
	}
	else
	{
		VEC3D v3 = jb->Position() + jb->toGlobal(spj) - ib->Position() - ib->toGlobal(spi);
		rhs[srow + 0] = m * v3.x;
		rhs[srow + 1] = m * v3.y;
		rhs[srow + 2] = m * v3.z;
		v3 = jb->toGlobal(hj);
		//v = ib->toGlobal(kconst->g_i());
		rhs[srow + 3] = m * v3.dot(ib->toGlobal(hi));
	}
}

void universalConstraint::constraintEquation(double *rhs, double* q, double* dq)
{
	int i = ib->ID() * 7;
	int j = jb->ID() * 7;
	VEC3D ri = i < 0 ? ib->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
	VEC3D rj = j < 0 ? jb->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
	EPD epi = i < 0 ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD epj = j < 0 ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	VEC3D v3 = (rj + epj.toGlobal(spj)) - (ri + epi.toGlobal(spi));
	//VEC3D v3 = jb->Position() + jb->toGlobal(spj) - ib->Position() - ib->toGlobal(spi);
	rhs[srow + 0] = v3.x;
	rhs[srow + 1] = v3.y;
	rhs[srow + 2] = v3.z;
	v3 = epj.toGlobal(hj);
	//v = ib->toGlobal(kconst->g_i());
	rhs[srow + 3] = v3.dot(epi.toGlobal(hi));
}

void universalConstraint::constraintJacobian(SMATD& cjaco)
{
	if (ib->NumDOF() == DIM2)
	{
		return;
	}
	else
	{
		if (ib->MassType() != pointMass::GROUND)
		{
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, icol + i) = -1;
			EPD ep = ib->getEP();// m->getParameterv(ib);
			cjaco.extraction(srow + 0, icol + 3, POINTER(B(ep, -spi)), MAT3X4);
			cjaco.extraction(srow + 3, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ep, hi))), VEC4);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			//if (!ib->ID()) ib = &ground;
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, jcol + i) = 1;
			EPD ep = jb->getEP();
			cjaco.extraction(srow + 0, jcol + 3, POINTER(B(ep, spj)), MAT3X4);
			cjaco.extraction(srow + 3, jcol + 3, POINTER(transpose(ib->toGlobal(hi), B(ep, hj))), VEC4);
		}
	}
}

void universalConstraint::constraintJacobian(SMATD& cjaco, double* q, double* dq)
{
	int i = ib->ID() * 7;
	int j = jb->ID() * 7;
	EPD epi = i < 0 ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD epj = j < 0 ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	if (ib->MassType() != pointMass::GROUND)
	{
		//unsigned int id = ib->ID() * 7;
		for (unsigned i(0); i < 3; i++) cjaco(srow + i, icol + i) = -1;
		//EPD ep(q[id + 3], q[id + 4], q[id + 5], q[id + 6]); //EPD ep = ib->getEP();// m->getParameterv(ib);
		cjaco.extraction(srow + 0, icol + 3, POINTER(B(epi, -spi)), MAT3X4);
		cjaco.extraction(srow + 3, icol + 3, POINTER(transpose(epj.toGlobal(hj), B(epi, hi))), VEC4);
	}
	if (jb->MassType() != pointMass::GROUND)
	{
		//if (!ib->ID()) ib = &ground;
		//unsigned int id = jb->ID() * 7;
		for (unsigned i(0); i < 3; i++) cjaco(srow + i, jcol + i) = 1;
		//EPD ep(q[id + 3], q[id + 4], q[id + 5], q[id + 6]); //EPD ep = jb->getEP();
		cjaco.extraction(srow + 0, jcol + 3, POINTER(B(epj, spj)), MAT3X4);
		cjaco.extraction(srow + 3, jcol + 3, POINTER(transpose(epi.toGlobal(hi), B(epj, hj))), VEC4);
	}
}

void universalConstraint::derivate(MATD& lhs, double mul)
{

}

void universalConstraint::differentialEquation(double *rhs, double* q, double *dq, double t)
{
	bool ig = ib->MassType() == pointMass::GROUND;
	bool jg = jb->MassType() == pointMass::GROUND;
	int i = ig ? -1 : ib->ID() * 7;
	int j = jg ? -1 : jb->ID() * 7;
	EPD pi = ig ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD pj = jg ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	EPD dpi = ig ? ib->getEV() : EPD(dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	EPD dpj = jg ? jb->getEV() : EPD(dq[j + 3], dq[j + 4], dq[j + 5], dq[j + 6]);
	VEC3D gamma = spherical_differential(dpi, dpj);
	rhs[srow + 0] = gamma.x;
	rhs[srow + 1] = gamma.y;
	rhs[srow + 2] = gamma.z;
	double d1_0 = dot_1_differential(hi, hj, pi, pj, dpi, dpj);
	rhs[srow + 3] = d1_0;
}
