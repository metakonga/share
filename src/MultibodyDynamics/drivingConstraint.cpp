#include "drivingConstraint.h"
#include "kinematicConstraint.h"
#include "simulation.h"
#include "numeric_utility.h"
#include <QStringList>
#include <QDebug>
drivingConstraint::drivingConstraint()
	: maxnnz(0)
	, theta(0)
	, start_time(0.0)
	, plus_time(0.0)
{

}

drivingConstraint::drivingConstraint(QString _name)
	: name(_name)
	, maxnnz(0)
	, theta(0)
	, start_time(0.0)
	, plus_time(0.0)
{

}

drivingConstraint::~drivingConstraint()
{

}

void drivingConstraint::define(
	kinematicConstraint* kc, drivingConstraint::Type td, double init, double cont)
{
//	m = kc->jMass();
//	pos0 = m->Position();
	kconst = kc;
// 	qDebug() << kconst->name();
// 	qDebug() << (int)td;
// 	qDebug() << init;
// 	qDebug() << cont;
	type = td;
	init_v = init;
	cons_v = cont;
	if (td == DRIVING_ROTATION)
	{
		pointMass* im = kconst->iMass();
		pointMass* jm = kconst->jMass();
		VEC3D fi = im->toGlobal(kconst->f_i());
		VEC3D fj = jm->toGlobal(kconst->f_j());
		init_v = acos(fi.dot(fj));;
	}
// 	else
// 	{
// 		VEC3D dist = kconst->CurrentDistance();
// 		init_v = dist.length();
// 	}
//	cons_v = cont;
}

void drivingConstraint::updateInitialCondition()
{
	if (type == DRIVING_ROTATION)
	{
		pointMass* im = kconst->iMass();
		pointMass* jm = kconst->jMass();
		VEC3D fi = im->toGlobal(kconst->f_i());
		VEC3D fj = jm->toGlobal(kconst->f_j());
		init_v = acos(fi.dot(fj));
	}
	else
	{
		VEC3D dist = kconst->CurrentDistance();
		init_v = dist.length();
	}
}

void drivingConstraint::constraintEquation(double mul, double* rhs)
{
	double v = 0.0;
	if (type == DRIVING_TRANSLATION)
	{
		VEC3D dist = kconst->CurrentDistance();
		VEC3D hi = kconst->iMass()->toGlobal(kconst->h_i());
		if (start_time > simulation::ctime + plus_time)
			v = hi.dot(dist) - (init_v + 0.0 * simulation::ctime);
		else
			v = hi.dot(dist) - (init_v + cons_v * (simulation::ctime - start_time + plus_time));
		rhs[srow] = mul * v;
	}
	else if (type == DRIVING_ROTATION)
	{
		if (start_time > simulation::ctime + plus_time)
			v = theta - (init_v + 0.0 * simulation::ctime);
		else
			v = theta - (init_v + cons_v * (simulation::ctime - start_time + plus_time));
		rhs[srow] = mul * v;
	}
}

void drivingConstraint::constraintEquation(double *rhs, double* q, double* dq, double ct)
{
	pointMass* im = kconst->iMass();
	pointMass* jm = kconst->jMass();
	int i = im->ID() * 7;
	int j = jm->ID() * 7;
	double v = 0.0;
	if (type == DRIVING_TRANSLATION)
	{
		VEC3D ri = i < 0 ? im->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
		VEC3D rj = j < 0 ? jm->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
		EPD epi = i < 0 ? im->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
		EPD epj = j < 0 ? jm->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
		VEC3D dist = (rj + epj.toGlobal(kconst->sp_j())) - (ri + epi.toGlobal(kconst->sp_i())); //kconst->CurrentDistance();
		//VEC3D dist = kconst->CurrentDistance();
		VEC3D hi = epi.toGlobal(kconst->h_i());
		if (start_time > simulation::ctime + plus_time)
			v = hi.dot(dist) - (init_v + 0.0 * simulation::ctime);
		else
			v = hi.dot(dist) - (init_v + cons_v * (simulation::ctime - start_time + plus_time));
		rhs[srow] = v;
	}
	else if (type == DRIVING_ROTATION)
	{
		if (start_time > simulation::ctime + plus_time)
			v = theta - (init_v + 0.0 * simulation::ctime);
		else
			v = theta - (init_v + cons_v * (simulation::ctime - start_time + plus_time));
		rhs[srow] = v;
	}
}

void drivingConstraint::constraintJacobian(SMATD& cjaco, double ct)
{
	if (type == DRIVING_TRANSLATION)
	{
		VEC3D fdij = kconst->CurrentDistance();
		pointMass* im = kconst->iMass();
		pointMass* jm = kconst->jMass();
		VEC3D hi = im->toGlobal(kconst->h_i());
		VEC3D D1;
		VEC4D D2;
		int ic = kconst->iColumn();
		int jc = kconst->jColumn();
		if (im->MassType() != pointMass::GROUND)
		{
			D1 = -hi;
			D2 = transpose(fdij, B(im->getEP(), kconst->h_i())) - transpose(hi, B(im->getEP(), kconst->sp_i()));
			cjaco.extraction(srow, ic, POINTER(D1), POINTER(D2), VEC3_4);
		}
		if (jm->MassType() != pointMass::GROUND)
		{
			D1 = hi;
			transpose(hi, B(jm->getEP(), kconst->sp_j()));
			cjaco.extraction(srow, jc, POINTER(D1), POINTER(D2), VEC3_4);
		}
	}
	else if (type == DRIVING_ROTATION)
	{
		pointMass* im = kconst->iMass();
		pointMass* jm = kconst->jMass();
		int ic = kconst->iColumn();
		int jc = kconst->jColumn();
		VEC3D fi = im->toGlobal(kconst->f_i());
		VEC3D fj = jm->toGlobal(kconst->f_j());
		VEC3D gi = im->toGlobal(kconst->g_i());
		double stheta = acos(fi.dot(fj));
		double prad = 0.0;
// 		qDebug() << "start_time : " << start_time;
// 		qDebug() << "simulation_time : " << simulation::ctime;
		if (start_time > ct + plus_time)
			prad = init_v + 0.0 * ct;
		else
			prad = init_v + cons_v * (ct - start_time + plus_time);
		if (prad > M_PI)
			stheta = stheta;
	//	qDebug() << "prad : " << prad;
		theta = numeric::utility::angle_coefficient(prad, stheta);
		VEC3D zv;
		VEC4D D1 = transpose(fj, cos(theta) * B(im->getEP(), kconst->g_i()));
		VEC4D D2 = transpose((cos(theta) * gi - sin(theta) * fi), B(jm->getEP(), kconst->f_j()));
		if (im->MassType() != pointMass::GROUND)
			cjaco.extraction(srow, ic, POINTER(zv), POINTER(D1), VEC3_4);
		if (jm->MassType() != pointMass::GROUND)
			cjaco.extraction(srow, jc, POINTER(zv), POINTER(D2), VEC3_4);
	}
	
}

void drivingConstraint::constraintJacobian(SMATD& cjaco, double *q, double *dq, double ct)
{
	pointMass* im = kconst->iMass();
	pointMass* jm = kconst->jMass();
	int i = im->ID() * 7;
	int j = jm->ID() * 7;
	if (type == DRIVING_TRANSLATION)
	{
		VEC3D ri = i < 0 ? im->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
		VEC3D rj = j < 0 ? jm->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
		EPD epi = i < 0 ? im->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
		EPD epj = j < 0 ? jm->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
		VEC3D fdij = (rj + epj.toGlobal(kconst->sp_j())) - (ri + epi.toGlobal(kconst->sp_i())); //kconst->CurrentDistance();
		VEC3D hi = epi.toGlobal(kconst->h_i());
		VEC3D D1;
		VEC4D D2;
		int ic = kconst->iColumn();
		int jc = kconst->jColumn();
		if (im->MassType() != pointMass::GROUND)
		{
			D1 = -hi;
			D2 = transpose(fdij, B(epi, kconst->h_i())) - transpose(hi, B(epi, kconst->sp_i()));
			cjaco.extraction(srow, ic, POINTER(D1), POINTER(D2), VEC3_4);
		}
		if (jm->MassType() != pointMass::GROUND)
		{
			D1 = hi;
			transpose(hi, B(epj, kconst->sp_j()));
			cjaco.extraction(srow, jc, POINTER(D1), POINTER(D2), VEC3_4);
		}
	}
	else if (type == DRIVING_ROTATION)
	{
// 		pointMass* im = kconst->iMass();
// 		pointMass* jm = kconst->jMass();
		EPD epi = i < 0 ? im->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
		EPD epj = j < 0 ? jm->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
		int ic = kconst->iColumn();
		int jc = kconst->jColumn();
		VEC3D fi = epi.toGlobal(kconst->f_i());
		VEC3D fj = epj.toGlobal(kconst->f_j());
		VEC3D gi = epi.toGlobal(kconst->g_i());
		double stheta = acos(fi.dot(fj));
		double prad = 0.0;
		// 		qDebug() << "start_time : " << start_time;
		// 		qDebug() << "simulation_time : " << simulation::ctime;
		if (start_time > ct + plus_time)
			prad = init_v + 0.0 * ct;
		else
			prad = init_v + cons_v * (ct - start_time + plus_time);
		if (prad > M_PI)
			stheta = stheta;
		//	qDebug() << "prad : " << prad;
		theta = numeric::utility::angle_coefficient(prad, stheta);
		VEC3D zv;
		VEC4D D1 = transpose(fj, cos(theta) * B(epi, kconst->g_i()));
		VEC4D D2 = transpose((cos(theta) * gi - sin(theta) * fi), B(epj, kconst->f_j()));
		if (im->MassType() != pointMass::GROUND)
			cjaco.extraction(srow, ic, POINTER(zv), POINTER(D1), VEC3_4);
		if (jm->MassType() != pointMass::GROUND)
			cjaco.extraction(srow, jc, POINTER(zv), POINTER(D2), VEC3_4);
	}
}

void drivingConstraint::differentialEquation(double* q, double *dq, double *rhs)
{
	pointMass* im = kconst->iMass();
	pointMass* jm = kconst->jMass();
	int i = im->ID() * 7;
	int j = jm->ID() * 7;
	EPD pi = i < 0 ? im->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD pj = j < 0 ? jm->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	EPD dpi = i < 0 ? im->getEV() : EPD(dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	EPD dpj = j < 0 ? jm->getEV() : EPD(dq[j + 3], dq[j + 4], dq[j + 5], dq[j + 6]);
	if (type == DRIVING_TRANSLATION)
	{
		VEC3D ri = i < 0 ? im->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
		VEC3D rj = j < 0 ? jm->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
		VEC3D dri = i < 0 ? im->getVelocity() : VEC3D(dq[i + 0], dq[i + 1], dq[i + 2]);
		VEC3D drj = j < 0 ? jm->getVelocity() : VEC3D(dq[j + 0], dq[j + 1], dq[j + 2]);
		VEC3D dij = rj + pj.toGlobal(kconst->sp_j()) - ri - pi.toGlobal(kconst->sp_i());
		VEC3D spi = kconst->sp_i();
		VEC3D spj = kconst->sp_j();
		VEC3D _ai = kconst->h_i();
		VEC3D ai = pi.toGlobal(_ai);
		double v = dri.dot(B(pi, _ai) * dpi) - dij.dot(B(dpi, _ai) * dpi) + ai.dot(B(dpi, spi) * dpi) + dpi.dot(transpose(B(pi, spi), B(pi, _ai) * dpi))
			- drj.dot(B(pi, _ai) * dpi) - ai.dot(B(dpj, spj) * dpj) - dpj.dot(transpose(B(pj, spj), B(pi, _ai) * dpi))
			- dpi.dot(transpose(B(pi, _ai), drj + B(pj, spj) * dpj - dri - B(pi, _ai) * dpi));
		rhs[srow] = v;
	}
	else if (type = DRIVING_ROTATION)
	{
		VEC3D _gi = kconst->g_i();
		VEC3D _fi = kconst->f_i();
		VEC3D _fj = kconst->f_j();
		VEC3D fi = pi.toGlobal(kconst->f_i());
		VEC3D fj = pj.toGlobal(kconst->f_j());
		VEC3D gi = pi.toGlobal(kconst->g_i());
		double v = -cos(theta) * (fj.dot(B(dpi, _gi) * dpi) + gi.dot(B(dpj, _fj) * dpj) + 2.0 * dpi.dot(transpose(B(pi, _gi), B(pj, _fj) * dpj)))
			+ sin(theta) * (fj.dot(B(dpi, _fi) * dpi) + fi.dot(B(dpj, _fj) * dpj) + 2.0 * dpi.dot(transpose(B(pi, _fi), B(pj, _fj) * dpj)));
		rhs[srow] = v;
	}
}

void drivingConstraint::saveData(QTextStream& qts)
{
	qts << "ELEMENT " << "drive_constraint" << endl
		<< "NAME " << name << endl
		<< "TYPE " << type << endl
		<< "TARGET " << kconst->name() << endl
		<< "START_TIME " << start_time << endl
		<< "PARAMETES " << init_v << " " << cons_v << endl;
}
