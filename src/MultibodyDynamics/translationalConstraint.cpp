#include "translationalConstraint.h"

translationalConstraint::translationalConstraint()
	: kinematicConstraint()
{

}

translationalConstraint::translationalConstraint(
	mbd_model* _md, QString& _nm, kinematicConstraint::Type kt,
	pointMass* ip, VEC3D& _spi, VEC3D& _fi, VEC3D& _gi, 
	pointMass* jp, VEC3D& _spj, VEC3D& _fj, VEC3D& _gj)
	: kinematicConstraint(_md, _nm, kt, ip, _spi, _fi, _gi, jp, _spj, _fj, _gj)
{

}

translationalConstraint::~translationalConstraint()
{

}

void translationalConstraint::constraintEquation(double m, double* rhs)
{
	if (ib->NumDOF() == DIM2)
	{

	}
	else
	{
		VEC3D v3 = jb->toGlobal(hj);
		VEC3D v3g = ib->toGlobal(gi);
		VEC3D v3f = ib->toGlobal(fi);
		rhs[srow + 0] = m * v3.dot(v3f);
		rhs[srow + 1] = m * v3.dot(v3g);
		v3 = jb->Position() + jb->toGlobal(spj) - ib->Position();// -ib->toGlobal(kconst->sp_i());
		rhs[srow + 2] = m * (v3.dot(v3f) - spi.dot(fi));
		rhs[srow + 3] = m * (v3.dot(v3g) - spi.dot(gi));
		rhs[srow + 4] = m * v3f.dot(jb->toGlobal(fj));
	}
}

void translationalConstraint::constraintJacobian(SMATD& cjaco)
{
	if (ib->NumDOF() == DIM2)
	{

	}
	else
	{
		VEC3D dij = (jb->Position() + jb->toGlobal(spj)) - (ib->Position() + ib->toGlobal(spi));
		if (ib->MassType() != pointMass::GROUND)
		{
			cjaco.extraction(srow + 0, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ib->getEP(), fi))), VEC4);
			cjaco.extraction(srow + 1, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ib->getEP(), gi))), VEC4);
			cjaco.extraction(srow + 2, icol + 0, POINTER((-fi)), POINTER(transpose(dij + ib->toGlobal(spi), B(ib->getEP(), fi))), VEC3_4);
			cjaco.extraction(srow + 3, icol + 0, POINTER((-gi)), POINTER(transpose(dij + ib->toGlobal(spi), B(ib->getEP(), gi))), VEC3_4);
			cjaco.extraction(srow + 4, icol + 3, POINTER(transpose(jb->toGlobal(fj), B(ib->getEP(), fi))), VEC4);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			cjaco.extraction(srow + 0, jcol + 3, POINTER(transpose(ib->toGlobal(fi), B(jb->getEP(), hj))), VEC4);
			cjaco.extraction(srow + 1, jcol + 3, POINTER(transpose(ib->toGlobal(gi), B(jb->getEP(), hj))), VEC4);
			cjaco.extraction(srow + 2, jcol + 0, POINTER(ib->toGlobal(fi)), POINTER(transpose(ib->toGlobal(fi), B(jb->getEP(), spj))), VEC3_4);
			cjaco.extraction(srow + 3, jcol + 0, POINTER(ib->toGlobal(gi)), POINTER(transpose(ib->toGlobal(gi), B(jb->getEP(), spj))), VEC3_4);
			cjaco.extraction(srow + 4, jcol + 3, POINTER(transpose(ib->toGlobal(fi), B(jb->getEP(), fj))), VEC4);
		}
	}
}
