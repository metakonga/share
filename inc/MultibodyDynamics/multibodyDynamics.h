#ifndef MULTIBODYDYNAMICS_H
#define MULTIBODYDYNAMICS_H

//#include "msimulation.h"
#include "simulation.h"
#include "kinematicConstraint.h"
#include "mbd_model.h"
#include "pointMass.h"
#include "utility/parabola_predictor.h"
#include <QFile>

class waveGeneratorModel;

class multibodyDynamics : public simulation
{
	typedef QMap<QString, kinematicConstraint*>::iterator kConstIterator;
	//typedef QMap<QString, drivingConstraint*>::iterator dConstIterator;
	typedef QMap<QString, pointMass*>::iterator massIterator;

public:
	multibodyDynamics();
	multibodyDynamics(mbd_model *_md);
	virtual ~multibodyDynamics();
	bool initialize(startingModel* stm = NULL);

	void prediction(unsigned int cs);
	int correction(unsigned int cs);
	double oneStepCorrection();

	bool saveResult(double ct);
	void saveFinalResult(QFile& qf);
	void setIntegratorType(mbd_integrator_type itype) { itor_type = itype; }
	unsigned int getOutCount() { return outCount; }
	unsigned int N_NR_Iteration();

	mbd_model* MbdModel() { return md; }

	int oneStepAnalysis(double ct, unsigned int cstep);
	void calcReactionForce(double ct);
	//bool cpuRun();

private:
	unsigned int mdim;
	unsigned int tdim;
	unsigned int nnz;
	unsigned int n_NR_iteration;
	int dof;				// degree of freedom

	// 	MKL_INT lapack_one;
	// 	MKL_INT lapack_info;
	// 	MKL_INT *ptDof;
	// 	MKL_INT *permutation;

	double dt2accp;
	double dt2accv;
	double dt2acc;
	double divalpha;
	double divbeta;
	double alpha;
	double beta;
	double gamma;
	double eps;
	double *lagMul;

	VECD rhs;
	MATD lhs;
	VECD pre;
	VECD ipp;
	VECD ipv;
	VECD ee;
	VECD cEQ;
	SMATD cjaco;
	//pointMass ground;

	void calculateRhs();
	void solveInitialValue(VECD& y, double t);
	void FULL_LEOM();
	void FULL_LEOM(VECD& y, double t);
	bool initialize_implicit_hht();
	bool initialize_explicit_rk4();
	void calcMassMatrix(double mul = 1.0);
	void calcMassMatrix(VECD& y);
	void calcForceVector(VECD *vec = NULL);
	void calcForceVector(VECD& y, double t);
	VEC4D calcInertiaForce(EPD& ev, MAT33D& J, EPD& ep);
	VEC4D calcMoment(EPD& ep, VEC3D& m);
	void sparseConstraintJacobian();
	void sparseConstraintJacobian(VECD& y, double t);
	int runke_kutta_method(unsigned int cs);

	void calcMassSystemJacobian(double mul);
	void calcForceSystemJacobian(double gamma, double beta);
	void calcConstraintSystemJacobian(double mul);

	MAT44D D(VEC3D& a, VEC3D& lag);
	void constraintEquation();
	void constraintEquation(VECD& y, double *ce, double t);

	mbd_integrator_type itor_type;
	QFile qf_out;
	unsigned int outCount;
	mbd_model* md;
	utility::parabola_predictor paraPredictor;
};

#endif