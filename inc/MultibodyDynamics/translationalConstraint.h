#ifndef TRANSLATIONALCONSTRAINT_H
#define TRANSLATIONALCONSTRAINT_H

#include "kinematicConstraint.h"

class translationalConstraint : public kinematicConstraint
{
public:
	translationalConstraint();
	translationalConstraint(mbd_model* _md, QString& _nm, kinematicConstraint::Type kt,
		pointMass* ip, VEC3D& spi, VEC3D& _fi, VEC3D& _gi,
		pointMass* jp, VEC3D& spj, VEC3D& _fj, VEC3D& _gj);
	virtual ~translationalConstraint();

	virtual void calculation_reaction_force(double ct);
	virtual void constraintEquation(double m, double* rhs);
	virtual void constraintEquation(double *ce, double* q, double* dq);
	virtual void constraintJacobian(SMATD& cjaco);
	virtual void constraintJacobian(SMATD& cjaco, double* q, double* dq);
	virtual void derivate(MATD& lhs, double mul);
	virtual void differentialEquation(double *rhs, double* q, double *dq, double t);
};

#endif