#include "sphericalConstraint.h"

sphericalConstraint::sphericalConstraint()
	: kinematicConstraint()
{

}

sphericalConstraint::sphericalConstraint(mbd_model* _md, QString& _nm, kinematicConstraint::Type kt,
	pointMass* ip, VEC3D& _spi, VEC3D& _fi, VEC3D& _gi,
	pointMass* jp, VEC3D& _spj, VEC3D& _fj, VEC3D& _gj)
	: kinematicConstraint(_md, _nm, kt, ip, _spi, _fi, _gi, jp, _spj, _fj, _gj)
{

}

sphericalConstraint::~sphericalConstraint()
{

}

void sphericalConstraint::constraintEquation(double m, double* rhs)
{
	if (ib->NumDOF() == DIM2)
		return;
	else
	{
		VEC3D v3 = jb->Position() + jb->toGlobal(spj) - ib->Position() - ib->toGlobal(spi);
		rhs[srow + 0] = m * v3.x;
		rhs[srow + 1] = m * v3.y;
		rhs[srow + 2] = m * v3.z;
	}
}

void sphericalConstraint::constraintEquation(double *rhs, double* q, double* dq)
{
	int i = ib->ID() * 7;
	int j = jb->ID() * 7;
	VEC3D ri = i < 0 ? ib->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
	VEC3D rj = j < 0 ? jb->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
	EPD epi = i < 0 ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD epj = j < 0 ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	VEC3D v3 = (rj + epj.toGlobal(spj)) - (ri + epi.toGlobal(spi));
	rhs[srow + 0] = v3.x;
	rhs[srow + 1] = v3.y;
	rhs[srow + 2] = v3.z;
}

void sphericalConstraint::constraintJacobian(SMATD& cjaco)
{
	if (ib->NumDOF() == DIM2)
		return;
	else
	{
		if (ib->MassType() != pointMass::GROUND)
		{
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, icol + i) = -1;
			EPD ep = ib->getEP();// m->getParameterv(ib);
			cjaco.extraction(srow + 0, icol + 3, POINTER(B(ep, -spi)), MAT3X4);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			//if (!ib->ID()) ib = &ground;
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, jcol + i) = 1;
			EPD ep = jb->getEP();
			cjaco.extraction(srow + 0, jcol + 3, POINTER(B(ep, spj)), MAT3X4);
		}
	}
}

void sphericalConstraint::constraintJacobian(SMATD& cjaco, double* q, double* dq)
{
	if (ib->MassType() != pointMass::GROUND)
	{
		unsigned int id = ib->ID() * 7;
		for (unsigned i(0); i < 3; i++) cjaco(srow + i, icol + i) = -1;
		EPD ep(q[id + 3], q[id + 4], q[id + 5], q[id + 6]); //EPD ep = ib->getEP();// m->getParameterv(ib);
		cjaco.extraction(srow + 0, icol + 3, POINTER(B(ep, -spi)), MAT3X4);
	}
	if (jb->MassType() != pointMass::GROUND)
	{
		//if (!ib->ID()) ib = &ground;
		unsigned int jd = jb->ID() * 7;
		for (unsigned i(0); i < 3; i++) cjaco(srow + i, jcol + i) = 1;
		EPD ep(q[jd + 3], q[jd + 4], q[jd + 5], q[jd + 6]); //EPD ep = jb->getEP();
		cjaco.extraction(srow + 0, jcol + 3, POINTER(B(ep, spj)), MAT3X4);
	}
}

void sphericalConstraint::derivate(MATD& lhs, double mul)
{
	VEC3D L(lm[0], lm[1], lm[2]);
	MAT44D Di = -D(spi, L);
	MAT44D Dj = D(spj, L);

	lhs.plus(icol + 3, icol + 3, POINTER(Di), MAT4x4, mul);
	lhs.plus(jcol + 3, jcol + 3, POINTER(Dj), MAT4x4, mul);
}

void sphericalConstraint::differentialEquation(double *rhs, double* q, double *dq, double t)
{
	bool ig = ib->MassType() == pointMass::GROUND;
	bool jg = jb->MassType() == pointMass::GROUND;
	int i = ig ? -1 : ib->ID() * 7;
	int j = jg ? -1 : jb->ID() * 7;
	EPD dpi = ig ? ib->getEV() : EPD(dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	EPD dpj = jg ? jb->getEV() : EPD(dq[j + 3], dq[j + 4], dq[j + 5], dq[j + 6]);
	VEC3D v3 = spherical_differential(dpi, dpj);
	rhs[srow + 0] = v3.x;
	rhs[srow + 1] = v3.y;
	rhs[srow + 2] = v3.z;
}
