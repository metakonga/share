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

void revoluteConstraint::constraintEquation2D(double m, double* rhs)
{
	VEC3D v3 = j->getPosition() + j->toGlobal2D(spj) - i->getPosition() - i->toGlobal2D(spi);
	rhs[srow + 0] = m * v3.x;
	rhs[srow + 1] = m * v3.y;
}

void revoluteConstraint::constraintJacobian2D(SMATD& cjaco)
{
	if (i->MassType() != pointMass::GROUND)
	{
		int ic = i->ID() * 3;
		for (unsigned t(0); t < 2; t++) 
			cjaco(srow + t, ic + t) = -1;
		double th = i->getPosition().z;
		cjaco(srow + 0, ic + 2) = -spi.x * sin(th) - spi.y * cos(th);
		cjaco(srow + 1, ic + 2) = spi.x * cos(th) - spi.y * sin(th);
	}
	if (j->MassType() != pointMass::GROUND)
	{
		int jc = j->ID() * 3;
		for (unsigned t(0); t < 2; t++) 
			cjaco(12 + srow + t, jc + t) = 1;
		double th = j->getPosition().z;
		cjaco(12 + srow + 0, jc + 2) = -spj.x * sin(th) - spj.y * cos(th);
		cjaco(12 + srow + 1, jc + 2) = spj.x * cos(th) - spj.y * sin(th);
	}
}