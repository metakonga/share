#ifndef GEARCONSTRAINT_H
#define GEARCONSTRAINT_H

#include "kinematicConstraint.h"

class artificialCoordinate;

class gearConstraint : public kinematicConstraint
{
	struct rotationInfo
	{
		double angi;
		double angj;
		double pangi;
		double pangj;
	};
public:
	gearConstraint();
	gearConstraint(QString _nm, pointMass* i, kinematicConstraint* ik, pointMass* j, kinematicConstraint* jk, double r);
	gearConstraint(QString _nm);
	virtual ~gearConstraint();

	void setFirstBody_i(pointMass* _fi) { fmi = _fi; }
	void setFirstBody_j(pointMass* _fj) { fmj = _fj; }
	void setSecondBody_i(pointMass* _si) { smi = _si; }
	void setSecondBody_j(pointMass* _sj) { smj = _sj; }
	void setGearRatio(double gr) { ratio = gr; }
	void setContactPoint(VEC3D& p) { cpoint = p; }
	void setFirstLocalPosition(VEC3D& p) { fs = p; }
	void setSecondLocalPosition(VEC3D& p) { ss = p; }
	void setFirstArtificialCoordinate(artificialCoordinate* ac) { fac = ac; }
	void setSecondArtificialCoordinate(artificialCoordinate* ac) { sac = ac; }

	//void constraintEquation(int sr, double* rhs, double mul);
	//void constraintJacobian(int sr, SMATD& cjaco);
	void saveGearConstraintData(QTextStream& qts);
	void setGroundPointer(pointMass* _g) { g = _g; }
	void decisionPreviousData();

	virtual void calculation_reaction_force(double ct) {};
	virtual void constraintEquation(double m, double* rhs);
	virtual void constraintJacobian(SMATD& cjaco);


private:
	pointMass *g;
	pointMass *fmi;
	pointMass *fmj;
	pointMass *smi;
	pointMass *smj;
	artificialCoordinate* fac;
	artificialCoordinate* sac;
	double ratio;
	VEC3D cpoint;
	VEC3D init_f_i;
	VEC3D init_f_j;

	VEC3D init_g_i;
	VEC3D init_g_j;
	int ki;
	int kj;
	double dthi;
	double dthj;
 	double ci;
 	double cj;
 	VEC3D fs;
 	VEC3D ss;
	int nri;
	int nrj;

	rotationInfo roti;
	rotationInfo rotj;
};

#endif