#include "revoluteConstraint.h"

revoluteConstraint::revoluteConstraint()
	: kinematicConstraint()
{

}

revoluteConstraint::revoluteConstraint(mbd_model* _md, QString& _nm, kinematicConstraint::Type kt,
	pointMass* ip, VEC3D& _spi, VEC3D& _fi, VEC3D& _gi, 
	pointMass* jp, VEC3D& _spj, VEC3D& _fj, VEC3D& _gj)
	: kinematicConstraint(_md, _nm, kt, ip, _spi, _fi, _gi, jp, _spj, _fj, _gj)
{

}

revoluteConstraint::~revoluteConstraint()
{

}

void revoluteConstraint::constraintEquation(double m, double* rhs)
{
	if (ib->NumDOF() == DIM2)
	{
		VEC3D v3 = jb->getPosition() + jb->toGlobal2D(spj) - ib->getPosition() - ib->toGlobal2D(spi);
		rhs[srow + 0] = m * v3.x;
		rhs[srow + 1] = m * v3.y;
	}
	else
	{
		VEC3D v3;
		v3 = jb->getPosition() + jb->toGlobal(spj) - ib->getPosition() - ib->toGlobal(spi);
		rhs[srow + 0] = m * v3.x;
		rhs[srow + 1] = m * v3.y;
		rhs[srow + 2] = m * v3.z;
		v3 = jb->toGlobal(hj);
		//v = ib->toGlobal(kconst->g_i());
		rhs[srow + 3] = m * v3.dot(ib->toGlobal(fi));
		rhs[srow + 4] = m * v3.dot(ib->toGlobal(gi));
	}
}

void revoluteConstraint::constraintJacobian(SMATD& cjaco)
{
	if (ib->NumDOF() == DIM2)
	{
		if (ib->MassType() != pointMass::GROUND)
		{
			int ic = ib->ID() * 3;
			for (unsigned t(0); t < 2; t++)
				cjaco(srow + t, ic + t) = -1;
			double th = ib->getPosition().z;
			cjaco(srow + 0, ic + 2) = -spi.x * sin(th) - spi.y * cos(th);
			cjaco(srow + 1, ic + 2) = spi.x * cos(th) - spi.y * sin(th);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			int jc = jb->ID() * 3;
			for (unsigned t(0); t < 2; t++)
				cjaco(srow + t, jc + t) = 1;
			double th = jb->getPosition().z;
			cjaco(srow + 0, jc + 2) = -spj.x * sin(th) - spj.y * cos(th);
			cjaco(srow + 1, jc + 2) = spj.x * cos(th) - spj.y * sin(th);
		}
	}
	else
	{
		if (ib->MassType() != pointMass::GROUND)
		{
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, icol + i) = -1;
			EPD ep = ib->getEP();// m->getParameterv(ib);
			cjaco.extraction(srow + 0, icol + 3, POINTER(B(ep, -spi)), MAT3X4);
			cjaco.extraction(srow + 3, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ep, fi))), VEC4);
			cjaco.extraction(srow + 4, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ep, gi))), VEC4);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, jcol + i) = 1;
			EPD ep = jb->getEP();
		//	VEC3D spj = jb->toGlobal(spj);
			cjaco.extraction(srow + 0, jcol + 3, POINTER(B(ep, spj)), MAT3X4);
			cjaco.extraction(srow + 3, jcol + 3, POINTER(transpose(ib->toGlobal(fi), B(ep, hj))), VEC4);
			cjaco.extraction(srow + 4, jcol + 3, POINTER(transpose(ib->toGlobal(gi), B(ep, hj))), VEC4);
		}
	}
}