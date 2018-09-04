#include "multibodyDynamics.h"
#include "kinematicConstraint.h"
#include "forceElement.h"
#include "linearSolver.hpp"
#include "errors.h"
#include <QDebug>
#include <QTime>
#include <QTextStream>
#include <QFile>
#include <cusolverDn.h>

multibodyDynamics::multibodyDynamics()
	: md(NULL)
	, mdim(0)
	, outCount(0)
	//, permutation(NULL)
{

}

multibodyDynamics::multibodyDynamics(mbd_model *_md)
	: md(_md)
	, mdim(0)
	, outCount(0)
	//, permutation(NULL)
{
	//ground.setID(-1);
	//ground.makeTransformationMatrix();
}

multibodyDynamics::~multibodyDynamics()
{
	//qf_out.close();
	if (qf_out.isOpen())
		qf_out.close();
	//if (permutation) delete[] permutation; permutation = NULL;
}

bool multibodyDynamics::initialize()
{
	outCount = 0;
	int sdim = 0;
	multibodyDynamics::nnz = 0;
	alpha = -0.3;
	beta = (1 - alpha) * (1 - alpha) / 4;
	gamma = 0.5 - alpha;
	eps = 1E-6;
	unsigned int nm = md->pointMasses().size();
// 	int od = md->mode2D() ? 3 : 7;
// 	int minus2D = md->mode2D() ? -3 : 0;
// 	if (nm)
// 	{
// 		mdim = nm * od;
// 	}
	unsigned int nr = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		pm->setID(nr++);
		pm->setInertia();
		mdim += pm->NumDOF();
		if (pm->NumDOF() == DIM3)
			sdim++;
	}
	unsigned int srow = 0;
	foreach(kinematicConstraint* kconst, md->kinConstraint())
	{
		kconst->setStartRow(srow);
		nnz += kconst->maxNNZ();
		switch (kconst->constType())
		{
		case kinematicConstraint::FIXED:		srow += 6; break;
		case kinematicConstraint::REVOLUTE:		srow += 5; break;
		case kinematicConstraint::TRANSLATIONAL:	srow += 5; break;
		case kinematicConstraint::SPHERICAL:		srow += 3; break;
		case kinematicConstraint::UNIVERSAL:		srow += 4; break;
		case kinematicConstraint::CABLE:			srow += 1; continue;
		case kinematicConstraint::GEAR:				srow += 1; break;
		case kinematicConstraint::COINCIDE:			srow += 1; break;
		default:
			break;
		}
		kconst->setFirstColumn(kconst->iMass()->ID() * kconst->iMass()->NumDOF());
		kconst->setSecondColumn(kconst->jMass()->ID() * kconst->jMass()->NumDOF());
	}
	sdim += srow;
	dof = mdim - sdim;
// 	if (!(md->mode2D()))
// 		sdim += md->pointMasses().size();
	tdim = mdim + sdim;

	lhs.alloc(tdim, tdim); lhs.zeros();
	rhs.alloc(tdim); rhs.zeros();
	pre.alloc(mdim); pre.zeros();
	ipp.alloc(mdim); ipp.zeros();
	ipv.alloc(mdim); ipv.zeros();
	ee.alloc(tdim); ee.zeros();
	cEQ.alloc(tdim - mdim); cEQ.zeros();
	dt2accp = simulation::dt*simulation::dt*(1 - 2 * beta)*0.5;
	dt2accv = simulation::dt*(1 - gamma);
	dt2acc = simulation::dt*simulation::dt*beta;
	divalpha = 1 / (1 + alpha);
	divbeta = -1 / (beta*simulation::dt*simulation::dt);
	cjaco.alloc(nnz + (nm)* 4, tdim - mdim, mdim);
	FULL_LEOM();
	 	QFile qf(model::path + "/" + "mat.txt");
	 	qf.open(QIODevice::WriteOnly);
	 	QTextStream qts(&qf);
	 	for (unsigned int i = 0; i < lhs.rows(); i++)
	 	{
	 		for (unsigned int j = 0; j < lhs.cols(); j++)
	 		{
	 			qts << lhs(i, j) << " ";
	 		}
	 		qts << endl;
	 	}
	 	qf.close();
	linearSolver ls(LAPACK_COL_MAJOR, tdim, 1, tdim, tdim);
	int info = ls.solve(lhs.getDataPointer(), rhs.get_ptr());

	unsigned int idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		pm->setAcceleration(VEC3D(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2)));
		if (pm->NumDOF() == DIM3)
			pm->setEA(EPD(rhs(idx + 3), rhs(idx + 4), rhs(idx + 5), rhs(idx + 6)));
		idx += pm->NumDOF();
	}
// 	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++)
// 	{
// 		pointMass* m = it.value();
// 		int idx = m->ID() * od;
// 		m->setAcceleration(VEC3D(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2)));
// 		if (md->mode2D() == false)
// 			m->setEA(EPD(rhs(idx + 3), rhs(idx + 4), rhs(idx + 5), rhs(idx + 6)));
// 		i++;
// 	}

	lagMul = rhs.get_ptr() + mdim;
	double* m_lag = rhs.get_ptr() + mdim;
	foreach(kinematicConstraint* kconst, md->kinConstraint())
	{
		kconst->setLagrangeMultiplierPointer(m_lag);
		//m_lag += kconst->numConst();
	}
	//constraintEquation();
	qf_out.setFileName(model::path + "/" + model::name + ".mrf");
	return true;
}

void multibodyDynamics::calcMassMatrix(double mul /* = 1.0 */)
{
	int cnt = 0;
//	lhs.resize(tdim, tdim);
	lhs.zeros();
	unsigned int idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		if (pm->NumDOF() == DIM2)
		{
			lhs(idx, idx) = lhs(idx + 1, idx + 1) = mul * pm->getMass();
			lhs(idx + 2, idx + 2) = mul * pm->getPriInertia().z;
		}
		else
		{
			lhs(idx, idx) = lhs(idx + 1, idx + 1) = lhs(idx + 2, idx + 2) = mul * pm->getMass();
		//	idx += 3;
			MAT44D LTJL = 4.0*mul*transpose(pm->getEP().G(), pm->getInertia() * pm->getEP().G());
			for (int j(0); j < 4; j++){
				for (int k(0); k < 4; k++){
					if (LTJL(j, k))
						lhs(idx+3 + j, idx+3 + k) = LTJL(j, k);
				}
			}
		}
		idx += pm->NumDOF();
	}
// 	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
// 		pointMass* m = it.value();
// 
// 		if (md->mode2D())
// 		{
// 			int idx = m->ID() * 3;
// 			lhs(idx, idx) = lhs(idx + 1, idx + 1) = mul * m->getMass();
// 			lhs(idx + 2, idx + 2) = mul * m->getPriInertia().z;
// 			//	cnt += 3;
// 		}
// 		else
// 		{
// 			int idx = m->ID() * 7;
// 			lhs(idx, idx) = lhs(idx + 1, idx + 1) = lhs(idx + 2, idx + 2) = mul * m->getMass();
// 			idx += 3;
// 			MAT44D LTJL = 4.0*mul*transpose(m->getEP().G(), m->getInertia() * m->getEP().G());
// 			for (int j(0); j < 4; j++){
// 				for (int k(0); k < 4; k++){
// 					if (LTJL(j, k))
// 						lhs(idx + j, idx + k) = LTJL(j, k);
// 				}
// 			}
// 			//cnt += 4;
// 		}
// 	}
}

VEC4D multibodyDynamics::calcInertiaForce(EPD& ev, MAT33D& J, EPD& ep)
{
	double GvP0 = -ev.e1*ep.e0 + ev.e0*ep.e1 + ev.e3*ep.e2 - ev.e2*ep.e3;
	double GvP1 = -ev.e2*ep.e0 - ev.e3*ep.e1 + ev.e0*ep.e2 + ev.e1*ep.e3;
	double GvP2 = -ev.e3*ep.e0 + ev.e2*ep.e1 - ev.e1*ep.e2 + ev.e0*ep.e3;
	return VEC4D(
		8 * (-ev.e1*J.a00*GvP0 - ev.e2*J.a11*GvP1 - ev.e3*J.a22*GvP2),
		8 * (ev.e0*J.a00*GvP0 - ev.e3*J.a11*GvP1 + ev.e2*J.a22*GvP2),
		8 * (ev.e3*J.a00*GvP0 + ev.e0*J.a11*GvP1 - ev.e1*J.a22*GvP2),
		8 * (-ev.e2*J.a00*GvP0 + ev.e1*J.a11*GvP1 + ev.e0*J.a22*GvP2));
}

VEC4D multibodyDynamics::calcMoment(EPD& ep, VEC3D& m)
{
	VEC4D out;
	out = transpose(2 * m, ep.G());
	///*VEC4D*/ out = transpose(2 * m, ep.G());
	return out;
}

void multibodyDynamics::calcForceVector(VECD* vec)
{
	VECD *out = vec ? vec : &rhs;
	out->zeros();
	int nmass = md->pointMasses().size();
	int i = 1;
	int cnt = 0;
	VEC3D nf;
	VEC4D mm, im, rf;
	VEC3D g = model::gravity;
	VEC3D cm;
	unsigned int idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		nf = pm->getMass() * g + pm->getCollisionForce() + pm->getExternalForce() + pm->getHydroForce();
		mm = calcMoment(pm->getEP(), pm->getCollisionMoment() + pm->getExternalMoment() + pm->getHydroMoment());
		rf = calcInertiaForce(pm->getEV(), pm->getInertia(), pm->getEP());
		//rf.x += mm.x; rf.y += mm.y; rf.z += mm.z; rf.w += mm.w;
		rf += mm;
		if (pm->NumDOF() == DIM2)
		{
			double m3 = 
				pm->getCollisionMoment().z + 
				pm->getExternalMoment().z + 
				pm->getHydroMoment().z;
			(*out)(idx + 0) = nf.x;
			(*out)(idx + 1) = nf.y;
			(*out)(idx + 2) = m3;
		}
		else
		{
			out->insert(idx, POINTER3(nf), POINTER4(rf), 3, 4);
		}
		idx += pm->NumDOF();
	}
	foreach(forceElement* fe, md->forceElements())
	{
		fe->calcForce(out);
	}
}

void multibodyDynamics::sparseConstraintJacobian()
{
	cjaco.zeroCount();
	unsigned int sr = 0;
	foreach(kinematicConstraint* kc, md->kinConstraint())
	{
		kc->constraintJacobian(cjaco);
		sr += kc->numConst();
	}
	unsigned int idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		if (pm->NumDOF() == DIM2)
			continue;
		EPD ep2 = 2 * pm->getEP();
		cjaco.extraction(sr++, idx + 3, ep2.Pointer(), 4);
		idx += pm->NumDOF();
	}
}

void multibodyDynamics::FULL_LEOM()
{
	sparseConstraintJacobian();
	calcMassMatrix();
	calcForceVector();

	for (int i(0); i < cjaco.nnz(); i++){
		lhs(cjaco.ridx[i] + mdim, cjaco.cidx[i]) = lhs(cjaco.cidx[i], cjaco.ridx[i] + mdim) = cjaco.value[i];
	}
	//lhs.display();
}

bool multibodyDynamics::saveResult(double ct)
{
	unsigned int cnt = 0;
	//	QTextStream qts(&qf_out);
	//qf_out.write((char*)&p, sizeof(unsigned int));
	foreach(pointMass* m, md->pointMasses())
	{
		VEC3D aa = 2.0 * m->getEP().G() * m->getEA();
		VEC3D av = 2.0 * m->getEP().G() * m->getEV();
		resultStorage::pointMassResultData pmr =
		{
			ct,
			m->getPosition(),
			m->getEP(),
			m->getVelocity(),
			av,
			m->getAcceleration(),
			aa,
			m->getEA(),
			m->getMass() * model::gravity + m->getCollisionForce() + m->getExternalForce() + m->getHydroForce()
		};
		model::rs->insertPointMassResult(m->Name(), pmr);
		//m->appendResultData(pmr);
	}
	foreach(kinematicConstraint* k, md->kinConstraint())
	{
		k->calculation_reaction_force(ct);
	}
	// 	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++)
	// 	{
	// 		
	// 		pointMass *m = it.value();
	// 		if (m->ID() == 0)
	// 			continue;
	// 		unsigned int name_size = m->name().size();
	// 		qf_out.write((char*)&name_size, sizeof(unsigned int));
	// 		qf_out.write((char*)m->name().toStdString().c_str(), sizeof(char) * m->name().size());
	// 		qf_out.write((char*)&cnt, sizeof(unsigned int));
	// 		qf_out.write((char*)&ct, sizeof(double));
	// 		qf_out.write((char*)&(m->getPosition()), sizeof(VEC3D));
	// 		qf_out.write((char*)&(m->getEP()), sizeof(EPD));
	// 		qf_out.write((char*)&(m->getVelocity()), sizeof(VEC3D));
	// 		qf_out.write((char*)&(m->getEV()), sizeof(EPD));
	// 		qf_out.write((char*)&(m->getAcceleration()), sizeof(VEC3D));
	// 		qf_out.write((char*)&(m->getEA()), sizeof(EPD));
	// 	//	m->getGeometryObject()->insertResultData(outCount, m->getPosition(), m->getEP());
	// 		cnt++;
	// 	}
	//outCount++;
	return true;
}

bool multibodyDynamics::oneStepAnalysis(double ct, unsigned int cstep)
{
	//qDebug() << "current time : " << ct;
	prediction(cstep);
	return correction(cstep);
}

void multibodyDynamics::calcReactionForce(double ct)
{
	foreach(kinematicConstraint* kconst, md->kinConstraint())
	{
		kconst->calculation_reaction_force(ct);
	}
}

void multibodyDynamics::prediction(unsigned int cs)
{
	//	VEC3D aa = 2.0 * md->pointMasses()["drum0"]->getEP().L() * md->pointMasses()["drum0"]->getEA();
	
	calcForceVector(&pre);
	sparseConstraintJacobian();
	for (unsigned int i = 0; i < cjaco.nnz(); i++){
		pre(cjaco.cidx[i]) -= cjaco.value[i] * lagMul[cjaco.ridx[i]];
	}
	pre *= alpha / (1 + alpha);
	unsigned int idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		for (int j = 0; j < 3; j++){
			ipp(idx + j) = pm->getPosition()(j)+pm->getVelocity()(j)* simulation::dt + rhs(idx + j) * dt2accp;
			ipv(idx + j) = pm->getVelocity()(j)+rhs(idx + j) * dt2accv;
		}
		if (!(md->mode2D()))
		{
			unsigned int sr = idx+3;
			for (int j = 0; j < 4; j++){
				ipp(sr + j) = pm->getEP()(j)+pm->getEV()(j)* simulation::dt + rhs(sr + j) * dt2accp;
				ipv(sr + j) = pm->getEV()(j)+rhs(sr + j) * dt2accv;
			}
		}
		idx += pm->NumDOF();
	}
	idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		VEC3D p, v;
		EPD ep, ev;
		for (int j = 0; j < 3; j++)
		{
			p(j) = ipp(idx + j) + rhs(idx + j) * dt2acc;
			v(j) = ipv(idx + j) + rhs(idx + j) * dt * gamma;
		}
		pm->setPosition(p);
		pm->setVelocity(v);
		if (pm->NumDOF() == DIM3)
		{
			unsigned int sr = idx + 3;
			for (int j = 0; j < 4; j++){
				ep(j) = ipp(sr + j) + rhs(sr + j) * dt2acc;
				ev(j) = ipv(sr + j) + rhs(sr + j) * dt * gamma;
			}
			pm->setEP(ep);
			pm->setEV(ev);
			pm->makeTransformationMatrix();
		}
		else
		{
			pm->makeTransformationMatrix2D();
		}
		idx += pm->NumDOF();
	}
}

void multibodyDynamics::calcMassSystemJacobian(double mul)
{
	EPD e;
	EPD edd;
	MAT33D inertia;
	MAT44D data;
	unsigned int idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		e = pm->getEP();
		edd = pm->getEA();
		inertia = pm->getInertia();
		data = 4.0 * mul * (-transpose(e.G(), inertia * edd.G()) + opMiner(inertia * (e.G() * edd)));
		lhs.plus(idx + 3, idx + 3, POINTER(data), MAT4x4);
		idx += pm->NumDOF();
	}
}

void multibodyDynamics::calcForceSystemJacobian()
{
	EPD e;
	EPD ed;
	MAT33D inertia;
	MAT44D data; 
	unsigned int idx = 0;
	foreach(pointMass* pm, md->pointMasses())
	{
		e = pm->getEP();
		ed = pm->getEV();
		inertia = pm->getInertia();
		data = -8.0 * gamma * dt * transpose(ed.G(), inertia * ed.G());
		lhs.plus(idx + 3, idx + 3, POINTER(data), MAT4x4);
		data = -8.0 * beta * dt * dt * (transpose(ed.G(), inertia * e.G()) + opMiner(inertia * (ed.G() * e)));
		lhs.plus(idx + 3, idx + 3, POINTER(data), MAT4x4);
		idx += pm->NumDOF();
	}
}

MAT44D multibodyDynamics::D(VEC3D& a, VEC3D& lag)
{
	double one = 2 * (a.x*lag.x + a.y*lag.y + a.z*lag.z);
	VEC3D upper = transpose(a, tilde(lag));
	VEC3D left = a.cross(lag);
	MAT33D mid(
		-a.z*lag.z - a.y*lag.y - lag.z*a.z - lag.y*a.y + one, a.y*lag.y + lag.y*a.x, a.z*lag.x + lag.z*a.x,
		a.x*lag.y + lag.x*a.y, -a.z*lag.z - a.x*lag.x - lag.z*a.z - lag.x*a.x + one, a.z*lag.y + lag.z*a.y,
		a.x*lag.z + lag.x*a.z, a.y*lag.z + lag.y*a.z, -a.y*lag.y - a.x*lag.x - lag.y*a.y - lag.x*a.x + one);
	return MAT44D(
		2 * one, upper.x, upper.y, upper.z,
		left.x, mid.a00, mid.a01, mid.a02,
		left.y, mid.a10, mid.a11, mid.a12,
		left.z, mid.a20, mid.a21, mid.a22);
}

void multibodyDynamics::calcConstraintSystemJacobian(double mul)
{
	int sr = 0;
	int ic = 0;
	int jc = 0;
	pointMass* ib = NULL;
	pointMass* jb = NULL;
	VEC3D dij;
	MAT44D Dv;
	MAT34D Bv;
	MATD m_lhs(lhs.rows(), lhs.cols()); m_lhs.zeros();
	for (kConstIterator it = md->kinConstraint().begin(); it != md->kinConstraint().end(); it++){
		kinematicConstraint *kconst = it.value();
		ib = kconst->iMass();
		jb = kconst->jMass();
		ic = kconst->iColumn();
		jc = kconst->jColumn();
		switch (kconst->constType())
		{
		case kinematicConstraint::REVOLUTE:
			if (ib->MassType() != pointMass::GROUND){
				//if (!jb->ID()) jb = &ground;
				Dv = -D(kconst->sp_i(), VEC3D(lagMul[sr + 0], lagMul[sr + 1], lagMul[sr + 2]));
				m_lhs.plus(ic + 3, ic + 3, POINTER(Dv), MAT4x4);
				Dv = lagMul[sr + 3] * D(kconst->g_i(), jb->toGlobal(kconst->h_j())) + lagMul[sr + 4] * D(kconst->f_i(), jb->toGlobal(kconst->h_j()));
				m_lhs.plus(ic + 3, ic + 3, POINTER(Dv), MAT4x4);
				if (jb->MassType() != pointMass::GROUND)
				{
					Dv = lagMul[sr + 3] * transpose(B(jb->getEP(), kconst->h_j()), B(ib->getEP(), kconst->g_i())) + lagMul[sr + 4] * transpose(B(jb->getEP(), kconst->h_j()), B(ib->getEP(), kconst->f_i()));
					m_lhs.plus(jc + 3, ic + 3, POINTER(Dv), MAT4x4);
				}
			}
			if (jb->MassType() != pointMass::GROUND){
				//if (!ib->ID()) ib = &ground;
				Dv = D(kconst->sp_j(), VEC3D(lagMul[sr + 0], lagMul[sr + 1], lagMul[sr + 2]));
				m_lhs.plus(jc + 3, jc + 3, POINTER(Dv), MAT4x4);
				Dv = lagMul[sr + 3] * D(kconst->h_j(), ib->toGlobal(kconst->g_i())) + lagMul[sr + 4] * D(kconst->h_j(), ib->toGlobal(kconst->f_i()));
				m_lhs.plus(jc + 3, jc + 3, POINTER(Dv), MAT4x4);
				if (ib->MassType() != pointMass::GROUND)
				{
					Dv = lagMul[sr + 3] * transpose(B(ib->getEP(), kconst->g_i()), B(jb->getEP(), kconst->h_j())) + lagMul[sr + 4] * transpose(B(ib->getEP(), kconst->f_i()), B(jb->getEP(), kconst->h_j()));
					m_lhs.plus(ic + 3, jc + 3, POINTER(Dv), MAT4x4);
				}
			}
			sr += 5;
			break;
		case kinematicConstraint::TRANSLATIONAL:
			if (ib->MassType() != pointMass::GROUND)
			{
				//if (!jb->ID()) jb = &ground;
				Dv = lagMul[sr + 0] * D(kconst->g_i(), jb->toGlobal(kconst->h_j())) + lagMul[sr + 1] * D(kconst->f_i(), jb->toGlobal(kconst->h_j())) + lagMul[sr + 4] * D(kconst->f_i(), jb->toGlobal(kconst->f_j()));
				m_lhs.plus(ic + 3, ic + 3, POINTER(Dv), MAT4x4);
				if (jb->MassType() != pointMass::GROUND)
				{
					Dv = lagMul[sr + 0] * transpose(B(ib->getEP(), kconst->g_i()), B(jb->getEP(), kconst->h_j())) + lagMul[sr + 1] * transpose(B(ib->getEP(), kconst->f_i()), B(jb->getEP(), kconst->h_j())) + lagMul[sr + 4] * transpose(B(ib->getEP(), kconst->f_i()), B(jb->getEP(), kconst->f_j()));
					m_lhs.plus(ic + 3, jc + 3, POINTER(Dv), MAT4x4);
				}
				dij = (jb->getPosition() + jb->toGlobal(kconst->sp_j())) - (ib->getPosition() + ib->toGlobal(kconst->sp_i()));
				Bv = -lagMul[sr + 2] * B(ib->getEP(), kconst->f_i()) - lagMul[sr + 3] * B(ib->getEP(), kconst->g_i());
				m_lhs.plus(ic + 0, ic + 3, POINTER(Bv), MAT3X4);
				m_lhs.plus(ic + 3, ic + 0, POINTER(Bv), MAT4x3);
				Dv = lagMul[sr + 2] * D(kconst->f_i(), dij + ib->toGlobal(kconst->sp_i())) + lagMul[sr + 3] * D(kconst->g_i(), dij + ib->toGlobal(kconst->sp_i()));
				m_lhs.plus(ic + 3, ic + 3, POINTER(Dv), MAT4x4);
				if (jb->MassType() != pointMass::GROUND)
				{
					Bv = -Bv;
					m_lhs.plus(ic + 3, jc + 0, POINTER(Bv), MAT4x3);
					Dv = lagMul[sr + 2] * transpose(B(ib->getEP(), kconst->f_i()), B(jb->getEP(), kconst->sp_j())) + lagMul[sr + 3] * transpose(B(ib->getEP(), kconst->g_i()), B(jb->getEP(), kconst->sp_j()));
					m_lhs.plus(ic + 3, jc + 3, POINTER(Dv), MAT4x4);
				}
			}
			if (jb->MassType() != pointMass::GROUND)
			{
				//if (!ib->ID()) ib = &ground;
				Dv = lagMul[sr + 0] * D(kconst->h_j(), jb->toGlobal(kconst->g_i())) + lagMul[sr + 1] * D(kconst->h_j(), jb->toGlobal(kconst->f_i())) + lagMul[sr + 4] * D(kconst->f_j(), jb->toGlobal(kconst->f_i()));
				m_lhs.plus(jc + 3, jc + 3, POINTER(Dv), MAT4x4);
				if (ib->MassType() != pointMass::GROUND)
				{
					Dv = lagMul[sr + 0] * transpose(B(jb->getEP(), kconst->h_j()), B(ib->getEP(), kconst->g_i())) + lagMul[sr + 1] * transpose(B(jb->getEP(), kconst->h_j()), B(ib->getEP(), kconst->f_i())) + lagMul[sr + 4] * transpose(B(jb->getEP(), kconst->f_j()), B(ib->getEP(), kconst->f_i()));
					m_lhs.plus(jc + 3, ic + 3, POINTER(Dv), MAT4x4);
				}
				Dv = lagMul[sr + 2] * D(kconst->sp_j(), ib->toGlobal(kconst->f_i())) + lagMul[sr + 3] * D(kconst->sp_j(), ib->toGlobal(kconst->g_i()));
				m_lhs.plus(jc + 3, jc + 3, POINTER(Dv), MAT4x4);
				if (ib->MassType() != pointMass::GROUND)
				{
					Bv = lagMul[sr + 2] * B(ib->getEP(), kconst->f_i()) + lagMul[sr + 3] * B(ib->getEP(), kconst->g_i());
					m_lhs.plus(jc + 0, ic + 3, POINTER(Bv), MAT3X4);
					Dv = lagMul[sr + 2] * transpose(B(jb->getEP(), kconst->sp_j()), B(ib->getEP(), kconst->f_i())) + lagMul[sr + 3] * transpose(B(jb->getEP(), kconst->sp_j()), B(ib->getEP(), kconst->g_i()));
					m_lhs.plus(jc + 3, ic + 3, POINTER(Dv), MAT4x4);
				}
			}
			sr += 5;
		}
	}
	for (massIterator m = md->pointMasses().begin(); m != md->pointMasses().end(); m++){
		pointMass* ms = m.value();
		if (ms->MassType() == pointMass::GROUND)
			continue;
		size_t id = (ms->ID() - 1) * 7;
		MAT44D Dv;
		size_t sr = tdim - mdim;
		Dv = lagMul[sr - ms->ID()] * 2.0;
		m_lhs.plus(id + 3, id + 3, POINTER(Dv), MAT4x4);
	}
	for (unsigned int i = 0; i < m_lhs.rows(); i++){
		for (unsigned int j = 0; j < m_lhs.cols(); j++){
			lhs(i, j) += mul * m_lhs(i, j);
		}
	}
}

void multibodyDynamics::constraintEquation()
{
	double* rhs = ee.get_ptr() + mdim;
	unsigned int sr = 0;
	foreach(kinematicConstraint *kc, md->kinConstraint())
	{
		kc->constraintEquation(divbeta, rhs);
		sr += kc->numConst();
	}
	foreach(pointMass* pm, md->pointMasses())
	{
		double d = pm->getEP().dot();
		rhs[sr++] = divbeta * (d - 1.0);
	}
}

bool multibodyDynamics::correction(unsigned int cs)
{
	int niter = 0;
	double e_norm = 1;
	linearSolver ls(LAPACK_COL_MAJOR, tdim, 1, tdim, tdim);
	while (1){
		niter++;
		if (niter > 100)
		{
			errors::setError(errors::MBD_EXCEED_NR_ITERATION);
			n_NR_iteration = niter;
			return false;
		}
		calcForceVector(&ee);
		sparseConstraintJacobian();
		for (int i = 0; i < cjaco.nnz(); i++){
			ee(cjaco.cidx[i]) -= cjaco.value[i] * lagMul[cjaco.ridx[i]];
		}

		calcMassMatrix(divalpha);
		unsigned int k = 0;
		foreach(pointMass* pm, md->pointMasses())
		{
			ee(k + 0) += -lhs(k + 0, k + 0) * rhs(k + 0);
			ee(k + 1) += -lhs(k + 1, k + 1) * rhs(k + 1);
			ee(k + 2) += -lhs(k + 2, k + 2) * rhs(k + 2);
			if (pm->NumDOF() == DIM3)
			{
				for (unsigned int j = 3; j < 7; j++){
					ee(k + j) += -(lhs(k + j, k + 3) * rhs(k + 3) + lhs(k + j, k + 4)*rhs(k + 4) + lhs(k + j, k + 5)*rhs(k + 5) + lhs(k + j, k + 6)*rhs(k + 6));
				}
			}
			k += pm->NumDOF();
		}
		calcMassSystemJacobian(divalpha * beta * dt * dt);
		calcForceSystemJacobian();
		//calcConstraintSystemJacobian(beta * dt * dt);
		for (int i(0); i < cjaco.nnz(); i++){
			lhs(cjaco.ridx[i] + mdim, cjaco.cidx[i]) = lhs(cjaco.cidx[i], cjaco.ridx[i] + mdim) = cjaco.value[i];
		}
		for (int i(0); i < mdim; i++) ee(i) -= pre(i);
		constraintEquation();
		e_norm = ee.norm();
		int info = ls.solve(lhs.getDataPointer(), ee.get_ptr());
		rhs += ee;
		unsigned int idx = 0;
		foreach(pointMass* pm, md->pointMasses())
		{
			VEC3D p, v;
			EPD ep, ev;
			for (int i = 0; i < 3; i++)
			{
				p(i) = ipp(idx + i) + dt2acc * rhs(idx + i);
				v(i) = ipv(idx + i) + dt * gamma * rhs(idx + i);
			}
			pm->setAcceleration(VEC3D(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2)));
			pm->setPosition(p);
			pm->setVelocity(v);
			idx += 3;
			if (pm->NumDOF() == DIM3)
			{
				for (int i = 0; i < 4; i++)
				{
					ep(i) = ipp(idx + i) + dt2acc * rhs(idx + i);
					ev(i) = ipv(idx + i) + dt * gamma * rhs(idx + i);
				}
				pm->setEA(EPD(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2), rhs(idx + 3)));
				pm->setEP(ep);
				pm->setEV(ev);
				idx += 4;
				pm->makeTransformationMatrix();
			}
			else
			{
				pm->makeTransformationMatrix2D();
			}
		}
		if (e_norm <= 1e-5)
		{
			n_NR_iteration += niter;
			break;
		}
	}
	return true;
}

double multibodyDynamics::oneStepCorrection()
{
	linearSolver ls(LAPACK_COL_MAJOR, tdim, 1, tdim, tdim);
	double e_norm = 1;
	calcForceVector(&ee);
	sparseConstraintJacobian();
	for (int i = 0; i < cjaco.nnz(); i++){
		ee(cjaco.cidx[i]) -= cjaco.value[i] * lagMul[cjaco.ridx[i] - mdim];
	}

	calcMassMatrix(divalpha);
	for (unsigned int i = 0, k = 0; i < md->pointMasses().size() - 1; i++, k += 7){
		ee(k + 0) += -lhs(k + 0, k + 0) * rhs(k + 0);
		ee(k + 1) += -lhs(k + 1, k + 1) * rhs(k + 1);
		ee(k + 2) += -lhs(k + 2, k + 2) * rhs(k + 2);
		for (unsigned int j = 3; j < 7; j++){
			ee(k + j) += -(lhs(k + j, k + 3) * rhs(k + 3) + lhs(k + j, k + 4)*rhs(k + 4) + lhs(k + j, k + 5)*rhs(k + 5) + lhs(k + j, k + 6)*rhs(k + 6));
		}
	}
	calcMassSystemJacobian(divalpha * beta * dt * dt);
	calcConstraintSystemJacobian(beta * dt * dt);
	for (int i(0); i < cjaco.nnz(); i++){
		lhs(cjaco.ridx[i], cjaco.cidx[i]) = lhs(cjaco.cidx[i], cjaco.ridx[i]) = cjaco.value[i];
	}
	for (int i(0); i < mdim; i++) ee(i) -= pre(i);
	constraintEquation();
	e_norm = ee.norm();
	//lhs.display();
	int info = ls.solve(lhs.getDataPointer(), ee.get_ptr());
	//dgesv_(ptDof, &lapack_one, lhs.getDataPointer(), ptDof, permutation, ee.get_ptr(), ptDof, &lapack_info);
	rhs += ee;
	int idx = 0;
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
		pointMass* m = it.value();
		if (m->ID() == 0)
			continue;
		VEC3D p, v;
		EPD ep, ev;
		for (int i = 0; i < 3; i++){
			p(i) = ipp(idx + i) + dt2acc * rhs(idx + i);
			v(i) = ipv(idx + i) + dt * gamma * rhs(idx + i);
		}
		m->setAcceleration(VEC3D(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2)));
		idx += 3;
		for (int i = 0; i < 4; i++){
			ep(i) = ipp(idx + i) + dt2acc * rhs(idx + i);
			ev(i) = ipv(idx + i) + dt * gamma * rhs(idx + i);
		}
		m->setEA(EPD(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2), rhs(idx + 3)));
		m->setPosition(p);
		m->setVelocity(v);
		m->setEP(ep);

		m->setEV(ev);
		idx += 4;
		m->makeTransformationMatrix();

	}

	return e_norm;
}

// bool multibodyDynamics::cpuRun()
// {
// 	unsigned int part = 0;
// 	unsigned int cStep = 0;
// 	unsigned int eachStep = 0;
// 
// 	ct = dt * cStep;
// 	qDebug() << "-------------------------------------------------------------" << endl
// 		<< "| Num. Part | Sim. Time | I. Part | I. Total | Elapsed Time |" << endl
// 		<< "-------------------------------------------------------------";
// 	QTextStream::AlignRight;
// 	//QTextStream::setRealNumberPrecision(6);
// 	QTextStream os(stdout);
// 	os.setRealNumberPrecision(6);
// 	if (saveResult(ct, part)){
// 		os << "| " << qSetFieldWidth(9) << part << qSetFieldWidth(12) << ct << qSetFieldWidth(10) << eachStep << qSetFieldWidth(11) << cStep << qSetFieldWidth(15) << "0" << qSetFieldWidth(0) << " |" << endl;
// 		//std::cout << "| " << std::setw(9) << part << std::setw(12) << ct << std::setw(10) << eachStep << std::setw(11) << cStep << std::setw(15) << 0 << std::endl;
// 	}
// 	QTime tme;
// 	tme.start();
// 	cStep++;
// 
// 	while (cStep < nstep)
// 	{
// 		if (_isWait)
// 			continue;
// 		if (_abort){
// 			emit finished();
// 			return false;
// 		}
// 		ct = dt * cStep;
// 		md->runExpression(ct, dt);
// 		prediction(cStep);
// 		correction(cStep);
// 		if (!((cStep) % step)){
// 			part++;
// 			emit sendProgress(part);
// 			if (saveResult(ct, part)){
// 				os << "| " << qSetFieldWidth(9) << part << qSetFieldWidth(12) << ct << qSetFieldWidth(10) << eachStep << qSetFieldWidth(11) << cStep << qSetFieldWidth(15) << tme.elapsed() * 0.001 << qSetFieldWidth(0) << " |" << endl;
// 			}
// 			eachStep = 0;
// 		}
// 		cStep++;
// 		eachStep++;
// 	}
// 	qf_out.close();
// 	emit finished();
// 	return true;
// }

// bool multibodyDynamics::gpuRun()
//{
//	return true;
//}