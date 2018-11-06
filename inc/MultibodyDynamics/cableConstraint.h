#ifndef CABLECONSTRAINT_H
#define CABLECONSTRAINT_H

#include "kinematicConstraint.h"

class cableConstraint : public kinematicConstraint
{
public:
	cableConstraint(
		QString _name,
		pointMass* _fi, VEC3D& _fspi, pointMass* _fj, VEC3D& _fspj,
		pointMass* _si, VEC3D& _sspi, pointMass* _sj, VEC3D& _sspj);
	virtual ~cableConstraint();

	void setApplyCableForceMass(pointMass* acm);
	//void constraintEquation(int sr, double* rhs, double mul);
	//void constraintJacobian(int sr, SMATD& cjaco);
	void saveCableConstraintData(QTextStream& qts);
	void updateCableInitLength();

	void calculation_reaction_force2body();
	virtual void calculation_reaction_force(double ct);
	virtual void constraintEquation(double m, double* rhs);
	virtual void constraintJacobian(SMATD& cjaco);
	virtual void derivate(MATD& lhs, double mul);
	virtual void differentialEquation(double *rhs, double* q, double *dq, double t){};

private:
	double f_c;
	double s_c;
	bool ispulling;
	double pre_f_c;
	pointMass* fi;
	pointMass* fj;
	pointMass* si;
	pointMass* sj;

	pointMass* a_cf_mass;

	VEC3D fspi, fspj, sspi, sspj;
};

#endif