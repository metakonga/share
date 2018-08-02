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
	, outCount(0)
	, n_NR_iteration(0)
	//, permutation(NULL)
{

}

multibodyDynamics::multibodyDynamics(mbd_model *_md)
	: md(_md)
	, outCount(0)
	, n_NR_iteration(0)
	//, permutation(NULL)
{
	//ground.setID(-1);
	//ground.makeTransformationMatrix();
}

multibodyDynamics::~multibodyDynamics()
{
	//qf_out.close();
	//if (lagMul) delete[] lagMul;
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
	unsigned int nac = md->artificialCoordinates().size();
	int od = md->mode2D() ? 3 : 7;
	int minus2D = md->mode2D() ? -3 : 0;
	if (nm)
	{
		mdim = nm * od + nac;
	}
	unsigned int cnt=0;
	foreach(artificialCoordinate* ac, md->artificialCoordinates())
	{
		ac->setMatrixLocation(mdim - nac + cnt);
		cnt++;
	}
	sdim;
	foreach(kinematicConstraint* kconst, md->kinConstraint())
	{
		kconst->setStartRow(sdim + nac);
		nnz += kconst->maxNNZ();
		switch (kconst->constType())
		{
		case kinematicConstraint::FIXED:		sdim += 6 + minus2D; break;
		case kinematicConstraint::REVOLUTE:		sdim += 5 + minus2D; break;
		case kinematicConstraint::TRANSLATIONAL:	sdim += 5 + minus2D; break;
		case kinematicConstraint::SPHERICAL:		sdim += 3; break;
		case kinematicConstraint::UNIVERSAL:		sdim += 4; break;
		case kinematicConstraint::CABLE:			sdim += 1; continue;
		case kinematicConstraint::GEAR:				sdim += 3; continue;
		default:
			break;
		}		
		kconst->setFirstColumn(kconst->iMass()->ID() * od);
		kconst->setSecondColumn(kconst->jMass()->ID() * od);
	}
	dof = mdim - sdim;
	if(!(md->mode2D()))
		sdim += md->pointMasses().size();
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
	cjaco.alloc(nnz + (nm) * 4, tdim - mdim, mdim);
	FULL_LEOM();
	QFile qf(model::path + model::name + "/" + "mat.txt");
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

	unsigned int i = 0;
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++)
	{
		pointMass* m = it.value();
		int idx = i * od;
		m->setAcceleration(VEC3D(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2)));
		if(md->mode2D() == false)
			m->setEA(EPD(rhs(idx + 3), rhs(idx + 4), rhs(idx + 5), rhs(idx + 6)));
		i++;
	}

	lagMul = rhs.get_ptr() + mdim;
	double* m_lag = rhs.get_ptr() + mdim;
	foreach(kinematicConstraint* kconst, md->kinConstraint())
	{
		kconst->setLagrangeMultiplierPointer(m_lag);
		m_lag += kconst->numConst();
	}
 	qf_out.setFileName(model::path + "/" + model::name + ".mrf");
	return true;
}

void multibodyDynamics::calcMassMatrix(double mul /* = 1.0 */)
{
	int cnt = 0;
	lhs.zeros();
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
		pointMass* m = it.value();
		if (md->mode2D())
		{
			lhs(cnt, cnt) = lhs(cnt + 1, cnt + 1) = mul * m->getMass();
			lhs(cnt + 2, cnt + 2) = mul * m->getPriInertia().z;
			cnt += 3;
		}
		else
		{
			lhs(cnt, cnt) = lhs(cnt + 1, cnt + 1) = lhs(cnt + 2, cnt + 2) = mul * m->getMass();
			cnt += 3;
			MAT44D LTJL = 4.0*mul*transpose(m->getEP().G(), m->getInertia() * m->getEP().G());
			for (int j(0); j < 4; j++){
				for (int k(0); k < 4; k++){
					if (LTJL(j, k))
						lhs(cnt + j, cnt + k) = LTJL(j, k);
				}
			}
			cnt += 4;
		}		
	}
}

VEC4D multibodyDynamics::calcInertiaForce(EPD& ev, MAT33D& J, EPD& ep)
{
	double GvP0 = -ev.e1*ep.e0 + ev.e0*ep.e1 + ev.e3*ep.e2 - ev.e2*ep.e3;
	double GvP1 = -ev.e2*ep.e0 - ev.e3*ep.e1 + ev.e0*ep.e2 + ev.e1*ep.e3;
	double GvP2 = -ev.e3*ep.e0 + ev.e2*ep.e1 - ev.e1*ep.e2 + ev.e0*ep.e3;
	return VEC4D(
		8*(-ev.e1*J.a00*GvP0 - ev.e2*J.a11*GvP1 - ev.e3*J.a22*GvP2),
		8*(ev.e0*J.a00*GvP0 - ev.e3*J.a11*GvP1 + ev.e2*J.a22*GvP2),
		8*(ev.e3*J.a00*GvP0 + ev.e0*J.a11*GvP1 - ev.e1*J.a22*GvP2),
		8*(-ev.e2*J.a00*GvP0 + ev.e1*J.a11*GvP1 + ev.e0*J.a22*GvP2));
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
//	md->userDefineProcess();
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
		pointMass* m = it.value();
		if (m->MassType() == pointMass::GROUND)
			continue;
		nf = m->getMass() * g + m->getCollisionForce() + m->getExternalForce() + m->getHydroForce();
		mm = calcMoment(m->getEP(), m->getCollisionMoment() + m->getExternalMoment() + m->getHydroMoment());
		rf = calcInertiaForce(m->getEV(), m->getInertia(), m->getEP());
		//rf.x += mm.x; rf.y += mm.y; rf.z += mm.z; rf.w += mm.w;
		rf += mm;
		if (md->mode2D())
		{
			VEC3D m3 = m->getCollisionMoment() + m->getExternalMoment() + m->getHydroMoment();
			(*out)(cnt++) = nf.x;
			(*out)(cnt++) = nf.y;
			(*out)(cnt++) = m3.z;
		}
		else
		{
			out->insert(cnt, POINTER3(nf), POINTER4(rf), 3, 4);
			cnt += 7;
		}
	}
	foreach(forceElement* fe, md->forceElements())
	{
		fe->calcForce(out);
	}
// 	double f = 40.0 * sin(2.0 * M_PI * simulation::ctime);
// 	(*out)(1) += f;
}

void multibodyDynamics::sparseConstraintJacobian()
{
	int sr, ic, jc, i = 0;
	pointMass *ib, *jb;
	//VEC3D dij;
	EPD ep;
	VEC3D dij;
	cjaco.zeroCount();
	if (!(md->mode2D()))
	{
		sr = tdim - md->pointMasses().size();
		for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
			pointMass* m = it.value();
			EPD ep2 = 2 * m->getEP();
			cjaco.extraction(sr++, i * 7 + 3, ep2.Pointer(), 4);
			i++;
		}
	}
	for (kConstIterator it = md->kinConstraint().begin(); it != md->kinConstraint().end(); it++){
		kinematicConstraint* kconst = it.value();
		sr = kconst->startRow() + (md->pointMasses().size()) * 7;
		ib = kconst->iMass(); ic = kconst->iColumn();
		jb = kconst->jMass(); jc = kconst->jColumn();
		switch (kconst->constType())
		{
		case kinematicConstraint::FIXED:
			if (ib->MassType() != pointMass::GROUND)
			{
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, ic + i) = -1;
				ep = ib->getEP();// m->getParameterv(ib);
				cjaco.extraction(sr + 0, ic + 3, POINTER(B(ep, -kconst->sp_i())), MAT3X4);
				cjaco.extraction(sr + 3, ic + 3, POINTER(transpose(jb->toGlobal(kconst->h_j()), B(ep, kconst->g_i()))), VEC4);
				cjaco.extraction(sr + 4, ic + 3, POINTER(transpose(jb->toGlobal(kconst->h_j()), B(ep, kconst->f_i()))), VEC4);
				cjaco.extraction(sr + 5, ic + 3, POINTER(transpose(jb->toGlobal(kconst->f_j()), B(ep, kconst->f_i()))), VEC4);
			}
			if (jb->MassType() != pointMass::GROUND)
			{
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, jc + i) = 1;
				ep = jb->getEP();
				cjaco.extraction(sr + 0, jc + 3, POINTER(B(ep, kconst->sp_j())), MAT3X4);
				cjaco.extraction(sr + 3, jc + 3, POINTER(transpose(ib->toGlobal(kconst->g_i()), B(ep, kconst->h_j()))), VEC4);
				cjaco.extraction(sr + 4, jc + 3, POINTER(transpose(ib->toGlobal(kconst->f_i()), B(ep, kconst->h_j()))), VEC4);
				cjaco.extraction(sr + 5, jc + 3, POINTER(transpose(ib->toGlobal(kconst->f_i()), B(ep, kconst->f_j()))), VEC4);
			}
			break;
		case kinematicConstraint::REVOLUTE:
			if (md->mode2D())
			{
				kconst->constraintJacobian2D(cjaco);
				break;
			}
				
			if (ib->MassType() != pointMass::GROUND)
			{
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, ic + i) = -1;
				ep = ib->getEP();// m->getParameterv(ib);
				cjaco.extraction(sr + 0, ic + 3, POINTER(B(ep, -kconst->sp_i())), MAT3X4);
				cjaco.extraction(sr + 3, ic + 3, POINTER(transpose(jb->toGlobal(kconst->h_j()), B(ep, kconst->f_i()))), VEC4);
				cjaco.extraction(sr + 4, ic + 3, POINTER(transpose(jb->toGlobal(kconst->h_j()), B(ep, kconst->g_i()))), VEC4);			
			}
			if (jb->MassType() != pointMass::GROUND)
			{
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, jc + i) = 1;
				ep = jb->getEP();
				VEC3D spj = jb->toGlobal(kconst->sp_j());
				cjaco.extraction(sr + 0, jc + 3, POINTER(B(ep, kconst->sp_j())), MAT3X4);
				cjaco.extraction(sr + 3, jc + 3, POINTER(transpose(ib->toGlobal(kconst->f_i()), B(ep, kconst->h_j()))), VEC4);
				cjaco.extraction(sr + 4, jc + 3, POINTER(transpose(ib->toGlobal(kconst->g_i()), B(ep, kconst->h_j()))), VEC4);
			}
			//std::cout << *sjc << std::endl;
			break;
		case kinematicConstraint::SPHERICAL:
			if (ib->MassType() != pointMass::GROUND)
			{
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, ic + i) = -1;
				ep = ib->getEP();// m->getParameterv(ib);
				cjaco.extraction(sr + 0, ic + 3, POINTER(B(ep, -kconst->sp_i())), MAT3X4);
			}
			if (jb->MassType() != pointMass::GROUND)
			{
				//if (!ib->ID()) ib = &ground;
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, jc + i) = 1;
				ep = jb->getEP();
				cjaco.extraction(sr + 0, jc + 3, POINTER(B(ep, kconst->sp_j())), MAT3X4);
			}
			break;
		case kinematicConstraint::UNIVERSAL:
			if (ib->MassType() != pointMass::GROUND)
			{
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, ic + i) = -1;
				ep = ib->getEP();// m->getParameterv(ib);
				cjaco.extraction(sr + 0, ic + 3, POINTER(B(ep, -kconst->sp_i())), MAT3X4);
				cjaco.extraction(sr + 3, ic + 3, POINTER(transpose(jb->toGlobal(kconst->h_j()), B(ep, kconst->h_i()))), VEC4);
			}
			if (jb->MassType() != pointMass::GROUND)
			{
				//if (!ib->ID()) ib = &ground;
				for (unsigned i(0); i < 3; i++) cjaco(sr + i, jc + i) = 1;
				ep = jb->getEP();
				cjaco.extraction(sr + 0, jc + 3, POINTER(B(ep, kconst->sp_j())), MAT3X4);
				cjaco.extraction(sr + 3, jc + 3, POINTER(transpose(ib->toGlobal(kconst->h_i()), B(ep, kconst->h_j()))), VEC4);
			}
			break;
		case kinematicConstraint::TRANSLATIONAL:
			//if (!ib) ib = &ground;
			//else if (!jb) jb = &ground;
// 			if (jb->ID()) dij = jb->getPosition() + jb->toGlobal(kconst->sp_j());
// 			if (ib->ID()) dij -= ib->getPosition() + ib->toGlobal(kconst->sp_i());
			dij = ( jb->getPosition() + jb->toGlobal(kconst->sp_j()) ) - ( ib->getPosition() + ib->toGlobal(kconst->sp_i()) );
			if (ib->MassType() != pointMass::GROUND)
			{
				//if(!(jb->ID())) jb = &ground;
				cjaco.extraction(sr + 0, ic + 3, POINTER(transpose(jb->toGlobal(kconst->h_j()), B(ib->getEP(), kconst->f_i()))), VEC4);
				cjaco.extraction(sr + 1, ic + 3, POINTER(transpose(jb->toGlobal(kconst->h_j()), B(ib->getEP(), kconst->g_i()))), VEC4);
				cjaco.extraction(sr + 2, ic + 0, POINTER((-kconst->f_i())), POINTER(transpose(dij + ib->toGlobal(kconst->sp_i()), B(ib->getEP(), kconst->f_i()))), VEC3_4);
				cjaco.extraction(sr + 3, ic + 0, POINTER((-kconst->g_i())), POINTER(transpose(dij + ib->toGlobal(kconst->sp_i()), B(ib->getEP(), kconst->g_i()))), VEC3_4);
				cjaco.extraction(sr + 4, ic + 3, POINTER(transpose(jb->toGlobal(kconst->f_j()), B(ib->getEP(), kconst->f_i()))), VEC4);
			}
			if (jb->MassType() != pointMass::GROUND)
			{
				//if (!(ib->ID())) ib = &ground;
				cjaco.extraction(sr + 0, jc + 3, POINTER(transpose(ib->toGlobal(kconst->f_i()), B(jb->getEP(), kconst->h_j()))), VEC4);
				cjaco.extraction(sr + 1, jc + 3, POINTER(transpose(ib->toGlobal(kconst->g_i()), B(jb->getEP(), kconst->h_j()))), VEC4);
				cjaco.extraction(sr + 2, jc + 0, POINTER(ib->toGlobal(kconst->f_i())), POINTER(transpose(ib->toGlobal(kconst->f_i()), B(jb->getEP(), kconst->sp_j()))), VEC3_4);
				cjaco.extraction(sr + 3, jc + 0, POINTER(ib->toGlobal(kconst->g_i())), POINTER(transpose(ib->toGlobal(kconst->g_i()), B(jb->getEP(), kconst->sp_j()))), VEC3_4);
				cjaco.extraction(sr + 4, jc + 3, POINTER(transpose(ib->toGlobal(kconst->f_i()), B(jb->getEP(), kconst->f_j()))), VEC4);
			}
			break;
		//case kinematicConstraint::CABLE:
 			
	//		return;
// 		case kinematicConstraint::GEAR:
// 
// 			break;
		}
	}
	foreach(cableConstraint* cc, md->cableConstraints())
	{
		sr = cc->startRow() + (md->pointMasses().size()) * 7;
		cc->constraintJacobian(sr, cjaco);
	}
	foreach(gearConstraint* gc, md->gearConstraints())
	{
		sr = gc->startRow() + (md->pointMasses().size()) * 7;
		gc->constraintJacobian(sr, cjaco);
	}
// 	for (dConstIterator it = md->drivingConstraints().begin(); it != md->drivingConstraints().end(); it++){
// 		drivingConstraint* dconst = it.value();
// 		sr = dconst->startRow() + (md->pointMasses().size() - 1) * 7;
// 		ic = dconst->startColumn();
// 		//for (int i = 0; i < 7; i++) 
// 			//if (dconst->use(i)) 
// 		cjaco(sr, ic) = 1.0;
// 	}
}

void multibodyDynamics::FULL_LEOM()
{
	calcMassMatrix();
	calcForceVector();
	sparseConstraintJacobian();
	for (int i(0); i < cjaco.nnz(); i++){
		lhs(cjaco.ridx[i], cjaco.cidx[i]) = lhs(cjaco.cidx[i], cjaco.ridx[i]) = cjaco.value[i];
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
			m->getEA()
		};
		model::rs->insertPointMassResult(m->name(), pmr);
			//m->appendResultData(pmr);
	}
	foreach(kinematicConstraint* k, md->kinConstraint())
	{
		k->calculation_reaction_force(ct);
	}
 	return true;
}

unsigned int multibodyDynamics::N_NR_Iteration()
{
	unsigned int ret = n_NR_iteration;
	n_NR_iteration = 0;
	return ret;
}

bool multibodyDynamics::oneStepAnalysis(double ct, unsigned int cstep)
{
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
	unsigned int i = 0;
	calcForceVector(&pre);
	sparseConstraintJacobian();
	for (i = 0; i < cjaco.nnz(); i++){
		pre(cjaco.cidx[i]) -= cjaco.value[i] * lagMul[cjaco.ridx[i] - mdim];
	}
	pre *= alpha / (1 + alpha);
	i = 0;
	int od = md->mode2D() ? 3 : 7;
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
		pointMass* m = it.value();
// 		if (m->ID() == 0)
// 			continue;
		unsigned int idx = i * od;
		for (int j = 0; j < 3; j++){
			ipp(idx + j) = m->getPosition()(j)+m->getVelocity()(j)* simulation::dt + rhs(idx + j) * dt2accp;
			ipv(idx + j) = m->getVelocity()(j)+rhs(idx + j) * dt2accv;
		}
		if (!(md->mode2D()))
		{
			idx += 3;
			for (int j = 0; j < 4; j++){
				double vt = m->getEV()(j)* simulation::dt;
				double vvt = 1.0 + vt;
				double t = m->getEP()(j)+m->getEV()(j)* simulation::dt;
				ipp(idx + j) = m->getEP()(j)+m->getEV()(j)* simulation::dt + rhs(idx + j) * dt2accp;
				ipv(idx + j) = m->getEV()(j)+rhs(idx + j) * dt2accv;
			}
		}		
		i++;
	}

	i = 0;
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
		pointMass* m = it.value();
// 		if (m->ID() == 0)
// 			continue;
		VEC3D p, v;
		EPD ep, ev;
		unsigned int idx = i * od;
		for (int j = 0; j < 3; j++)
		{
			p(j) = ipp(idx + j) + rhs(idx + j) * dt2acc;
			v(j) = ipv(idx + j) + rhs(idx + j) * dt * gamma;
		}
		m->setPosition(p);
		m->setVelocity(v);
		if (!(md->mode2D()))
		{
			idx += 3;
			for (int j = 0; j < 4; j++){
				ep(j) = ipp(idx + j) + rhs(idx + j) * dt2acc;
				ev(j) = ipv(idx + j) + rhs(idx + j) * dt * gamma;
			}
			m->setEP(ep);
			m->setEV(ev);
			m->makeTransformationMatrix();
		}
		else
		{
			m->makeTransformationMatrix2D();
		}
		i++;
	}
}

void multibodyDynamics::calcMassSystemJacobian(double mul)
{
	if (md->mode2D())
		return;
	int sr = 0;
	EPD e;
	EPD edd;
	MAT33D inertia;
	MAT44D data;
	unsigned int i = 0;
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
		pointMass* m = it.value();
		if (m->MassType() == pointMass::GROUND)
			continue;
		sr = i * 7 + 3;
		e = m->getEP();
		edd = m->getEA();
		inertia = m->getInertia();
		data = 4.0 * mul * (-transpose(e.G(), inertia * edd.G()) + opMiner(inertia * (e.G() * edd)));
		lhs.plus(sr, sr, POINTER(data), MAT4x4);
		i++;
	}
}

void multibodyDynamics::calcForceSystemJacobian()
{
	if (md->mode2D())
		return;
	int sr = 0;
	int i = 0;
	EPD e;
	EPD ed;
	MAT33D inertia;
	MAT44D data;
	for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++)
	{
		pointMass* m = it.value();
		sr = i * 7 + 3;
		e = m->getEP();
		ed = m->getEV();
		inertia = m->getInertia();
		data = -8.0 * gamma * dt * transpose(ed.G(), inertia * ed.G());
		lhs.plus(sr, sr, POINTER(data), MAT4x4);
		data = -8.0 * beta * dt * dt * (transpose(ed.G(), inertia * e.G()) + opMiner(inertia * (ed.G() * e)));
		lhs.plus(sr, sr, POINTER(data), MAT4x4);
		i++;
	}
//  	foreach(forceElement* fe, md->forceElements())
//  	{
//  		fe->derivate(lhs, -simulation::dt * gamma);
//  		fe->derivate_velocity(lhs, -beta * simulation::dt * simulation::dt);
//  	}
}

MAT44D multibodyDynamics::D(VEC3D& a, VEC3D& lag)
{
	double one = 2 * (a.x*lag.x + a.y*lag.y + a.z*lag.z);
	VEC3D upper = transpose(a, tilde(lag));
	VEC3D left = a.cross(lag);
	MAT33D mid(
		-a.z*lag.z - a.y*lag.y - lag.z*a.z - lag.y*a.y + one,								 a.y*lag.y + lag.y*a.x,								    a.z*lag.x + lag.z*a.x,
									   a.x*lag.y + lag.x*a.y, -a.z*lag.z - a.x*lag.x - lag.z*a.z - lag.x*a.x + one,									a.z*lag.y + lag.z*a.y,
									   a.x*lag.z + lag.x*a.z,								 a.y*lag.z + lag.y*a.z,	-a.y*lag.y - a.x*lag.x - lag.y*a.y - lag.x*a.x + one);
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
	VEC3D v3, v3g, v3f, v;
	pointMass *ib = NULL;
	pointMass *jb = NULL;
	if (!(md->mode2D()))
	{
		int mr = tdim - mdim - md->pointMasses().size();// tdim - md->pointMasses().size();;
		for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
			pointMass* m = it.value();
			double d = m->getEP().dot();
			rhs[mr++] = divbeta * (d - 1.0);
		}
	}	

	for (kConstIterator it = md->kinConstraint().begin(); it != md->kinConstraint().end(); it++){
		kinematicConstraint *kconst = it.value();
		ib = kconst->iMass();// if (!ib) ib = &ground;
		jb = kconst->jMass();// if (!jb) jb = &ground;
		switch (kconst->constType()){
		case kinematicConstraint::FIXED:
			v3 = jb->getPosition() + jb->toGlobal(kconst->sp_j()) - ib->getPosition() - ib->toGlobal(kconst->sp_i());
			rhs[kconst->startRow() + 0] = divbeta * v3.x;
			rhs[kconst->startRow() + 1] = divbeta * v3.y;
			rhs[kconst->startRow() + 2] = divbeta * v3.z;
			v3 = jb->toGlobal(kconst->h_j());
			//v = ib->toGlobal(kconst->g_i());
			rhs[kconst->startRow() + 3] = divbeta * v3.dot(ib->toGlobal(kconst->g_i()));
			rhs[kconst->startRow() + 4] = divbeta * v3.dot(ib->toGlobal(kconst->f_i()));
			rhs[kconst->startRow() + 5] = divbeta * jb->toGlobal(kconst->f_j()).dot(ib->toGlobal(kconst->f_i()));
			break;
		case kinematicConstraint::REVOLUTE:
			if (md->mode2D())
			{
				kconst->constraintEquation2D(divbeta, rhs);
				break;
			}				
			v3 = jb->getPosition() + jb->toGlobal(kconst->sp_j()) - ib->getPosition() - ib->toGlobal(kconst->sp_i());
			rhs[kconst->startRow() + 0] = divbeta * v3.x;
			rhs[kconst->startRow() + 1] = divbeta * v3.y;
			rhs[kconst->startRow() + 2] = divbeta * v3.z;
			v3 = jb->toGlobal(kconst->h_j());
			//v = ib->toGlobal(kconst->g_i());
			rhs[kconst->startRow() + 3] = divbeta * v3.dot(ib->toGlobal(kconst->f_i()));
			rhs[kconst->startRow() + 4] = divbeta * v3.dot(ib->toGlobal(kconst->g_i()));
			break;
		case kinematicConstraint::UNIVERSAL:
			v3 = jb->getPosition() + jb->toGlobal(kconst->sp_j()) - ib->getPosition() - ib->toGlobal(kconst->sp_i());
			rhs[kconst->startRow() + 0] = divbeta * v3.x;
			rhs[kconst->startRow() + 1] = divbeta * v3.y;
			rhs[kconst->startRow() + 2] = divbeta * v3.z;
			v3 = jb->toGlobal(kconst->h_j());
			//v = ib->toGlobal(kconst->g_i());
			rhs[kconst->startRow() + 3] = divbeta * v3.dot(ib->toGlobal(kconst->h_i()));
			break;
		case kinematicConstraint::SPHERICAL:
			v3 = jb->getPosition() + jb->toGlobal(kconst->sp_j()) - ib->getPosition() - ib->toGlobal(kconst->sp_i());
			rhs[kconst->startRow() + 0] = divbeta * v3.x;
			rhs[kconst->startRow() + 1] = divbeta * v3.y;
			rhs[kconst->startRow() + 2] = divbeta * v3.z;
			break;
		case kinematicConstraint::TRANSLATIONAL:
			v3 = jb->toGlobal(kconst->h_j());
			v3g = ib->toGlobal(kconst->g_i());
			v3f = ib->toGlobal(kconst->f_i());
			rhs[kconst->startRow() + 0] = divbeta * v3.dot(v3f);
			rhs[kconst->startRow() + 1] = divbeta * v3.dot(v3g);
			v3 = jb->getPosition() + jb->toGlobal(kconst->sp_j()) - ib->getPosition();// -ib->toGlobal(kconst->sp_i());
			rhs[kconst->startRow() + 2] = divbeta * (v3.dot(v3f) - kconst->sp_i().dot(kconst->f_i()));
			rhs[kconst->startRow() + 3] = divbeta * (v3.dot(v3g) - kconst->sp_i().dot(kconst->g_i()));
			rhs[kconst->startRow() + 4] = divbeta * v3f.dot(jb->toGlobal(kconst->f_j()));
			break;
		}
	}
	foreach(cableConstraint *cc, md->cableConstraints())
	{
		cc->constraintEquation(cc->startRow(), rhs, divbeta);
	}
	foreach(gearConstraint *gc, md->gearConstraints())
	{
		gc->constraintEquation(gc->startRow(), rhs, divbeta);
	}
}

bool multibodyDynamics::correction(unsigned int cs)
{
	int niter= 0;
	double e_norm = 1;
	linearSolver ls(LAPACK_COL_MAJOR, tdim, 1, tdim, tdim);
	while (1){
		niter++;
		if (niter > 100)
		{
			errors::setError(errors::MBD_EXCEED_NR_ITERATION);
			return false;
		}
			
		calcForceVector(&ee);
		sparseConstraintJacobian();
		for (int i = 0; i < cjaco.nnz(); i++){
			ee(cjaco.cidx[i]) -= cjaco.value[i] * lagMul[cjaco.ridx[i] - mdim];
		}

		calcMassMatrix(divalpha);
		for (unsigned int i = 0, k = 0; i < md->pointMasses().size(); i++, k += 7){
			ee(k + 0) += -lhs(k + 0, k + 0) * rhs(k + 0);
			ee(k + 1) += -lhs(k + 1, k + 1) * rhs(k + 1);
			ee(k + 2) += -lhs(k + 2, k + 2) * rhs(k + 2);
			for (unsigned int j = 3; j < 7; j++){
				ee(k + j) += -(lhs(k + j, k + 3) * rhs(k + 3) + lhs(k + j, k + 4)*rhs(k + 4) + lhs(k + j, k + 5)*rhs(k + 5) + lhs(k + j, k + 6)*rhs(k + 6));
			}
		}
		calcMassSystemJacobian(divalpha * beta * dt * dt);
		calcForceSystemJacobian();
		//calcConstraintSystemJacobian(beta * dt * dt);
		for (int i(0); i < cjaco.nnz(); i++){
			lhs(cjaco.ridx[i], cjaco.cidx[i]) = lhs(cjaco.cidx[i], cjaco.ridx[i]) = cjaco.value[i];
		}
		for (int i(0); i < mdim; i++) ee(i) -= pre(i);
		constraintEquation();
		e_norm = ee.norm();
		int info = ls.solve(lhs.getDataPointer(), ee.get_ptr());
		rhs += ee;
		int idx = 0;
		for (massIterator it = md->pointMasses().begin(); it != md->pointMasses().end(); it++){
			pointMass* m = it.value();
			VEC3D p, v;
			EPD ep, ev;
			for (int i = 0; i < 3; i++){
				p(i) = ipp(idx + i) + dt2acc * rhs(idx + i);
				v(i) = ipv(idx + i) + dt * gamma * rhs(idx + i);
			}
			m->setAcceleration(VEC3D(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2)));
			m->setPosition(p);
			m->setVelocity(v);
			idx += 3;
			if (!(md->mode2D()))
			{
				for (int i = 0; i < 4; i++){
					ep(i) = ipp(idx + i) + dt2acc * rhs(idx + i);
					ev(i) = ipv(idx + i) + dt * gamma * rhs(idx + i);
				}
				m->setEA(EPD(rhs(idx + 0), rhs(idx + 1), rhs(idx + 2), rhs(idx + 3)));
				m->setEP(ep);
				m->setEV(ev);
				idx += 4;
				m->makeTransformationMatrix();
			}
			else
			{
				m->makeTransformationMatrix2D();
			}
		}
		if (e_norm <= 1e-5)
		{
			break;
		}
	}
	n_NR_iteration += niter;
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
	for (unsigned int i = 0, k = 0; i < md->pointMasses().size()-1; i++, k += 7){
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