#include "springDamperModel.h"

springDamperModel::springDamperModel()
	: forceElement()
	, init_l(0)
	, k(0)
	, c(0)
	, f(0)
	, l(0)
	, dl(0)
{
	
}

springDamperModel::springDamperModel(
	QString _n, mbd_model* _md, 
	pointMass* _b, VEC3D _baseLoc,
	pointMass* _a, VEC3D _actionLoc,
	double _k, double _c)
	: forceElement(_n, _md, forceElement::TSDA, _b, _a)
	, init_l(0)
	, k(_k)
	, c(_c)
	, baseLoc(_baseLoc)
	, actionLoc(_actionLoc)
	, f(0)
	, l(0)
	, dl(0)
{
	init_l = (actionLoc - baseLoc).length();
	spi = base->toLocal(baseLoc - base->Position());
	spj = action->toLocal(actionLoc - action->Position());
}

springDamperModel::~springDamperModel()
{

}

void springDamperModel::calcForce(VECD* rhs)
{
	VEC3D Qi;
	VEC4D QRi;
	VEC3D Qj;
	VEC4D QRj;
	L = action->Position() + action->toGlobal(spj) - base->Position() - base->toGlobal(spi);
	l = L.length();
	VEC3D dL = action->getVelocity() + B(action->getEV(), spj) * action->getEP() - base->getVelocity() - B(base->getEV(), spi) * base->getEP();
	dl = L.dot(dL) / l;
	f = k * (l - init_l) + c * dl;
	Qi = (f / l) * L;
	Qj = -Qi;
	QRi = (f / l) * transpose(B(base->getEP(), spi), L);
	QRj = -(f / l) * transpose(B(action->getEP(), spj), L);
		
	if (base->MassType() != pointMass::GROUND)
	{
		if (base->NumDOF() == DIM2){}
		else
		{
			int irc = base->ID() * 7;
			rhs->plus(irc, POINTER3(Qi), POINTER(QRi), 3, 4);
		}		
	}
	if (action->MassType() != pointMass::GROUND)
	{
		if (action->NumDOF() == DIM2){}
		else
		{
			int jrc = action->ID() * 7;
			rhs->plus(jrc, POINTER3(Qj), POINTER(QRj), 3, 4);
		}
	}
	//rhs->plus(irc, POINTER3)
}

void springDamperModel::calcForce(double* rhs, double* q, double *dq, double t)
{
	VEC3D Qi;
	VEC4D QRi;
	VEC3D Qj;
	VEC4D QRj;
	VEC3D ri, rj, vi, vj;
	EPD epi, epj, evi, evj;
	if (base->MassType() == pointMass::GROUND)
	{
		ri = base->Position();
		epi = base->getEP();// (q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
		vi = base->getVelocity();// (dq[i + 0], dq[i + 1], dq[i + 2]);
		evi = base->getEV();// (dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	}
	else
	{
		unsigned int i = base->ID() * 7;
		ri = VEC3D(q[i + 0], q[i + 1], q[i + 2]);
		epi = EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
		vi = VEC3D(dq[i + 0], dq[i + 1], dq[i + 2]);
		evi = EPD(dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	}
	if (action->MassType() == pointMass::GROUND)
	{
		rj = action->Position();
		epj = action->getEP();// (q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
		vj = action->getVelocity();// (dq[i + 0], dq[i + 1], dq[i + 2]);
		evj = action->getEV();// (dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	}
	else
	{
		unsigned int j = action->ID() * 7;
		rj = VEC3D(q[j + 0], q[j + 1], q[j + 2]);
		epj = EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
		vj = VEC3D(dq[j + 0], dq[j + 1], dq[j + 2]);
		evj = EPD(dq[j + 3], dq[j + 4], dq[j + 5], dq[j + 6]);
	}
	L = rj + epj.toGlobal(spj) - ri - epi.toGlobal(spi);// action->Position() + action->toGlobal(spj) - base->Position() - base->toGlobal(spi);
	l = L.length();
	VEC3D dL = vj + B(evj, spj) * epj - vi - B(evi, spi) * epi;
	dl = L.dot(dL) / l;
	f = k * (l - init_l) + c * dl;
	Qi = (f / l) * L;
	Qj = -Qi;
	QRi = (f / l) * transpose(B(epi, spi), L);
	QRj = -(f / l) * transpose(B(epj, spj), L);

	if (base->MassType() != pointMass::GROUND)
	{
		if (base->NumDOF() == DIM2){}
		else
		{
			int irc = base->ID() * 7;
			rhs[irc + 0] = Qi.x;
			rhs[irc + 1] = Qi.y;
			rhs[irc + 2] = Qi.z;

			rhs[irc + 3] = QRi.x;
			rhs[irc + 4] = QRi.y;
			rhs[irc + 5] = QRi.z;
			rhs[irc + 6] = QRi.w;
			//rhs->plus(jrc, POINTER3(Qj), POINTER(QRj), 3, 4);
		}
	}
	if (action->MassType() != pointMass::GROUND)
	{
		if (action->NumDOF() == DIM2){}
		else
		{
			int jrc = action->ID() * 7;
			rhs[jrc + 0] = Qj.x;
			rhs[jrc + 1] = Qj.y;
			rhs[jrc + 2] = Qj.z;

			rhs[jrc + 3] = QRj.x;
			rhs[jrc + 4] = QRj.y;
			rhs[jrc + 5] = QRj.z;
			rhs[jrc + 6] = QRj.w;
			//rhs->plus(jrc, POINTER3(Qj), POINTER(QRj), 3, 4);
		}
	}
}

void springDamperModel::derivate(MATD& lhs, double mul)
{
	if (base->NumDOF() != action->NumDOF())
		return;
	else if (base->NumDOF() == DIM2)
		return;

	MAT33D N1;
	MAT34D NP1;
	MAT33D N2;
	MAT34D NP2;
	double c1 = (k - (f / l) - c*(dl / l));
	MAT33D v1 = transpose(L, L) / pow(l, 2.0);
	VEC3D v = action->getVelocity() + action->toDGlobal(spj) - base->getVelocity() - base->toDGlobal(spi);
	MAT33D v2 = transpose(L, v) / pow(l, 2.0);
	MAT33D diag = MAT33D(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1);
	N1 = -(c1 * v1 + c * v2 + (f / l) * diag);
	N2 = c * v1;
	NP1 = -N1 * B(base->getEP(), spi) - N2 * B(base->getEV(), spi);
	NP2 = N1 * B(action->getEP(), spj) + N2 * B(action->getEV(), spj);
	int irc = base->ID() * 7;
	int jrc = action->ID() * 7;
	lhs.plus(irc, irc, POINTER(N1), MAT3X3, -mul);
	//N1 = -N1;
	lhs.plus(irc, jrc, POINTER((-N1)), MAT3X3, -mul);
	lhs.plus(irc, irc + 3, POINTER(NP1), MAT3X4, -mul);
	lhs.plus(irc, jrc + 3, POINTER(NP2), MAT3X4, -mul);

	lhs.plus(jrc, irc, POINTER((-N1)), MAT3X3, -mul);
	//N1 = -N1;
	lhs.plus(jrc, jrc, POINTER(N1), MAT3X3, -mul);
	//NP1 = -NP1;
	lhs.plus(jrc, irc + 3, POINTER((-NP1)), MAT3X4, -mul);
	//NP2 = -NP2;
	lhs.plus(jrc, jrc + 3, POINTER((-NP2)), MAT3X4, -mul);
	MAT43D QqN1 = -transpose(B(base->getEP(), spi), N1);
	MAT44D QqN2 = -transpose(B(base->getEP(), spi), (N1 * B(base->getEP(), spi) + N2 * B(base->getEV(), spi))) + (f / l) * D(spi, L);
	MAT43D QqN3 = -QqN1;
	MAT44D QqN4 = transpose(B(base->getEP(), spi), N1 * B(action->getEP(), spj) + N2 * B(action->getEV(), spj));
	lhs.plus(irc + 3, irc, POINTER((QqN1)), MAT4x3, -mul);
	lhs.plus(irc + 3, irc + 3, POINTER(QqN2), MAT4x4, -mul);
	lhs.plus(irc + 3, jrc, POINTER(QqN3), MAT4x3, -mul);
	lhs.plus(irc + 3, jrc + 3, POINTER(QqN4), MAT4x4, -mul);

	QqN1 = -transpose(B(action->getEP(), spj), N1);
	QqN2 = -transpose(B(base->getEP(), spi), (N1 * B(base->getEP(), spi) + N2 * B(base->getEV(), spi)));
	QqN3 = -QqN1;
	QqN4 = transpose(B(base->getEP(), spi), N1 * B(action->getEP(), spj) + N2 * B(action->getEV(), spj)) + (f / l) * D(spj, L);
	lhs.plus(jrc + 3, irc, POINTER(QqN1), MAT4x3, -mul);
	lhs.plus(jrc + 3, irc + 3, POINTER(QqN2), MAT4x4, -mul);
	lhs.plus(jrc + 3, jrc, POINTER(QqN3), MAT4x3, -mul);
	lhs.plus(jrc + 3, jrc + 3, POINTER(QqN4), MAT4x4, -mul);
}

void springDamperModel::derivate_velocity(MATD& lhs, double mul)
{
	MAT33D coe = c * transpose(L, L) / pow(l, 2.0);
	MAT33D QdqF1 = -coe;
	MAT34D QdqF2 = -coe * B(base->getEP(), spi);
	MAT33D QdqF3 = coe;
	MAT34D QdqF4 = coe * B(action->getEP(), spj);

	MAT43D coe1 = c * transpose(B(base->getEP(), spi), transpose(L, L) / pow(l, 2.0));
	MAT43D QdqN1 = -coe1;
	MAT44D QdqN2 = -(coe1 * B(base->getEP(), spi));
	MAT43D QdqN3 = coe1;
	MAT44D QdqN4 = coe1 * B(action->getEP(), spj);

	int irc = base->ID() * 7;
	int jrc = action->ID() * 7;
	lhs.plus(irc, irc, POINTER(QdqF1), MAT3X3, mul);
	lhs.plus(irc, irc + 3, POINTER(QdqF2), MAT3X4, mul);
	lhs.plus(irc, jrc, POINTER(QdqF3), MAT3X3, mul);
	lhs.plus(irc, jrc + 3, POINTER(QdqF4), MAT3X4, mul);


	lhs.plus(irc + 3, irc, POINTER(QdqN1), MAT4x3, mul);
	lhs.plus(irc + 3, irc + 3, POINTER(QdqN2), MAT4x4, mul);
	lhs.plus(irc + 3, jrc, POINTER(QdqN3), MAT4x3, mul);
	lhs.plus(irc + 3, jrc + 3, POINTER(QdqN4), MAT4x4, mul);

	QdqF1 = -QdqF1;
	QdqF2 = -QdqF2;
	QdqF3 = -QdqF3;
	QdqF4 = -QdqF4;

	coe1 = c * transpose(B(action->getEP(), spj), transpose(L, L) / pow(l, 2.0));
	QdqN1 = coe1;
	QdqN2 = coe1 * B(base->getEP(), spi);
	QdqN3 = -coe1;
	QdqN4 = -(coe1 * B(action->getEP(), spj));

	lhs.plus(jrc, irc, POINTER(QdqF1), MAT3X3, mul);
	lhs.plus(jrc, irc + 3, POINTER(QdqF2), MAT3X4, mul);
	lhs.plus(jrc, jrc, POINTER(QdqF3), MAT3X3, mul);
	lhs.plus(jrc, jrc + 3, POINTER(QdqF4), MAT3X4, mul);


	lhs.plus(jrc + 3, irc, POINTER(QdqN1), MAT4x3, mul);
	lhs.plus(jrc + 3, irc + 3, POINTER(QdqN2), MAT4x4, mul);
	lhs.plus(jrc + 3, jrc, POINTER(QdqN3), MAT4x3, mul);
	lhs.plus(jrc + 3, jrc + 3, POINTER(QdqN4), MAT4x4, mul);
}

void springDamperModel::saveData(QTextStream& qts)
{
	qts << "ELEMENT " << "tsda" << endl
		<< "NAME " << name << endl
		<< "FIRST_BODY " << base->Name() << endl
		<< "SECOND_BODY " << action->Name() << endl
		<< "FIRST_LOCATION " << baseLoc.x << " " << baseLoc.y << " " << baseLoc.z << endl
		<< "SECOND_LOCATION " << actionLoc.x << " " << actionLoc.y << " " << actionLoc.z << endl
		<< "COEFF_SPRING " << k << endl << "COEFF_DAMPING " << c << endl;
// 	qts << endl
// 		<< "FORCE_ELEMENT " << type << endl
// 		<< "NAME " << name << endl
// 		<< "BASE " << base->Name() << endl
// 		<< "ACTION " << action->Name() << endl
// 		<< "BASELOC " << baseLoc.x << " " << baseLoc.y << " " << baseLoc.z << endl
// 		<< "ACTIONLOC " << actionLoc.x << " " << actionLoc.y << " " << actionLoc.z << endl
// 		<< "BASELOCAL " << spi.x << " " << spi.y << " " << spi.z << endl
// 		<< "ACTIONLOCAL " << spj.x << " " << spj.y << " " << spj.z << endl
// 		<< "PARAMETER " << init_l << " " << k << " " << c << endl;
}

void springDamperModel::saveResult(double t)
{

}
