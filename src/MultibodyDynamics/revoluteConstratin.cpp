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

void revoluteConstraint::constraintEquation(double m, double* rhs)
{
	if (ib->NumDOF() == DIM2)
	{
		VEC3D v3 = jb->Position() + jb->toGlobal2D(spj) - ib->Position() - ib->toGlobal2D(spi);
		rhs[srow + 0] = m * v3.x;
		rhs[srow + 1] = m * v3.y;
	}
	else
	{
		VEC3D v3;
		v3 = jb->Position() + jb->toGlobal(spj) - ib->Position() - ib->toGlobal(spi);
		rhs[srow + 0] = m * v3.x;
		rhs[srow + 1] = m * v3.y;
		rhs[srow + 2] = m * v3.z;
		v3 = jb->toGlobal(hj);
		//v = ib->toGlobal(kconst->g_i());
		rhs[srow + 3] = m * v3.dot(ib->toGlobal(fi));
		rhs[srow + 4] = m * v3.dot(ib->toGlobal(gi));
	}
}

void revoluteConstraint::constraintEquation(double *rhs, double* q, double* dq)
{
	int i = ib->ID() * 7;
	int j = jb->ID() * 7;
	VEC3D ri = i < 0 ? ib->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
	VEC3D rj = j < 0 ? jb->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
	EPD epi = i < 0 ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD epj = j < 0 ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	VEC3D v3;
	v3 = (rj + epj.toGlobal(spj)) - (ri + epi.toGlobal(spi));
	rhs[srow + 0] = v3.x;
	rhs[srow + 1] = v3.y;
	rhs[srow + 2] = v3.z;
	v3 = jb->toGlobal(hj);
	//v = ib->toGlobal(kconst->g_i());
	rhs[srow + 3] = v3.dot(epi.toGlobal(fi));
	rhs[srow + 4] = v3.dot(epi.toGlobal(gi));
}

void revoluteConstraint::constraintJacobian(SMATD& cjaco)
{

	if (ib->NumDOF() == DIM2)
	{
		if (ib->MassType() != pointMass::GROUND)
		{
			int ic = ib->ID() * 3;
			for (unsigned t(0); t < 2; t++)
				cjaco(srow + t, ic + t) = -1;
			double th = ib->Position().z;
			cjaco(srow + 0, ic + 2) = -spi.x * sin(th) - spi.y * cos(th);
			cjaco(srow + 1, ic + 2) = spi.x * cos(th) - spi.y * sin(th);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			int jc = jb->ID() * 3;
			for (unsigned t(0); t < 2; t++)
				cjaco(srow + t, jc + t) = 1;
			double th = jb->Position().z;
			cjaco(srow + 0, jc + 2) = -spj.x * sin(th) - spj.y * cos(th);
			cjaco(srow + 1, jc + 2) = spj.x * cos(th) - spj.y * sin(th);
		}
	}
	else
	{
		if (ib->MassType() != pointMass::GROUND)
		{
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, icol + i) = -1;
			EPD ep = ib->getEP();// m->getParameterv(ib);
			cjaco.extraction(srow + 0, icol + 3, POINTER(B(ep, -spi)), MAT3X4);
			cjaco.extraction(srow + 3, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ep, fi))), VEC4);
			cjaco.extraction(srow + 4, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ep, gi))), VEC4);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			for (unsigned i(0); i < 3; i++) cjaco(srow + i, jcol + i) = 1;
			EPD ep = jb->getEP();
		//	VEC3D spj = jb->toGlobal(spj);
			cjaco.extraction(srow + 0, jcol + 3, POINTER(B(ep, spj)), MAT3X4);
			cjaco.extraction(srow + 3, jcol + 3, POINTER(transpose(ib->toGlobal(fi), B(ep, hj))), VEC4);
			cjaco.extraction(srow + 4, jcol + 3, POINTER(transpose(ib->toGlobal(gi), B(ep, hj))), VEC4);
		}
	}
}

void revoluteConstraint::constraintJacobian(SMATD& cjaco, double* q, double* dq)
{
	if (ib->MassType() != pointMass::GROUND)
	{
		unsigned int id = ib->ID() * 7;
		for (unsigned i(0); i < 3; i++) cjaco(srow + i, icol + i) = -1;
		EPD ep(q[id + 3], q[id + 4], q[id + 5], q[id + 6]);// EPD ep = ib->getEP();// m->getParameterv(ib);
		cjaco.extraction(srow + 0, icol + 3, POINTER(B(ep, -spi)), MAT3X4);
		cjaco.extraction(srow + 3, icol + 3, POINTER(transpose(ep.toGlobal(hj), B(ep, fi))), VEC4);
		cjaco.extraction(srow + 4, icol + 3, POINTER(transpose(ep.toGlobal(hj), B(ep, gi))), VEC4);
	}
	if (jb->MassType() != pointMass::GROUND)
	{
		unsigned int id = jb->ID() * 7;
		for (unsigned i(0); i < 3; i++) cjaco(srow + i, jcol + i) = 1;
		//EPD ep = jb->getEP();
		EPD ep(q[id + 3], q[id + 4], q[id + 5], q[id + 6]);//	VEC3D spj = jb->toGlobal(spj);
		cjaco.extraction(srow + 0, jcol + 3, POINTER(B(ep, spj)), MAT3X4);
		cjaco.extraction(srow + 3, jcol + 3, POINTER(transpose(ep.toGlobal(fi), B(ep, hj))), VEC4);
		cjaco.extraction(srow + 4, jcol + 3, POINTER(transpose(ep.toGlobal(gi), B(ep, hj))), VEC4);
	}
}

void revoluteConstraint::derivate(MATD& lhs, double mul)
{
	VEC3D L(lm[0], lm[1], lm[2]);
	MAT44D Di = -D(spi, L) + lm[3] * D(fi, jb->toGlobal(hj)) + lm[4] * D(gi, jb->toGlobal(hj));
	MAT44D Dj = D(spj, L) + lm[3] * D(hj, ib->toGlobal(fi)) + lm[4] * D(hj, ib->toGlobal(gi));

	lhs.plus(icol + 3, icol + 3, POINTER(Di), MAT4x4, mul);
	lhs.plus(jcol + 3, jcol + 3, POINTER(Dj), MAT4x4, mul);
	
	Di = lm[3] * transpose(B(ib->getEP(), fi), B(jb->getEP(), hj));
	Di += lm[4] * transpose(B(ib->getEP(), gi), B(jb->getEP(), hj));
	Dj = lm[3] * transpose(B(jb->getEP(), hj), B(ib->getEP(), fi));
	Dj += lm[4] * transpose(B(jb->getEP(), hj), B(ib->getEP(), gi));

	lhs.plus(icol + 3, jcol + 3, POINTER(Di), MAT4x4, mul);
	lhs.plus(jcol + 3, icol + 3, POINTER(Dj), MAT4x4, mul);

}

void revoluteConstraint::differentialEquation(double *rhs, double* q, double* dq, double t)
{
	bool ig = ib->MassType() == pointMass::GROUND;
	bool jg = jb->MassType() == pointMass::GROUND;
	int i = ig ? -1 : ib->ID() * 7;
	int j = jg ? -1 : jb->ID() * 7;
	EPD pi = ig ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD pj = jg ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	EPD dpi = ig ? ib->getEV() : EPD(dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	EPD dpj = jg ? jb->getEV() : EPD(dq[j + 3], dq[j + 4], dq[j + 5], dq[j + 6]);
	VEC3D gamma = spherical_differential(dpi, dpj);
	rhs[srow + 0] = gamma.x;
	rhs[srow + 1] = gamma.y;
	rhs[srow + 2] = gamma.z;
	double d1_0 = dot_1_differential(hi, fj, pi, pj, dpi, dpj);
	double d1_1 = dot_1_differential(hi, gj, pi, pj, dpi, dpj);
	rhs[srow + 3] = d1_0;
	rhs[srow + 4] = d1_1;
}

void revoluteConstraint::calculation_reaction_force(double ct)
{
// 	SMATD jaco;
// 	VECD out;
// 	out.alloc(14);
// 	jaco.alloc(54, 5, 14);
// 	constraintJacobian(jaco);
// 	for (unsigned int i = 0; i < jaco.nnz(); i++)
// 	{
// 		out(jaco.cidx[i]) += jaco.value[i] * lm[jaco.ridx[i] - srow];
// 	}
// 	resultStorage::reactionForceData rfd =
// 	{
// 		ct,
// 		VEC3D(out(0), out(1), out(2)),
// 		VEC4D(out(3), out(4), out(5), out(6)),
// 		VEC3D(out(7), out(8), out(9)),
// 		VEC4D(out(10), out(11), out(12), out(13))
// 	};
// 	model::rs->insertReactionForceResult(nm, rfd);// .push_back(rfd);
}