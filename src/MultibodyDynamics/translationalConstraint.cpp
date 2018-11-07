#include "translationalConstraint.h"
#include <QDebug>

translationalConstraint::translationalConstraint()
	: kinematicConstraint()
{

}

translationalConstraint::translationalConstraint(
	mbd_model* _md, QString& _nm, kinematicConstraint::Type kt,
	pointMass* ip, VEC3D& _spi, VEC3D& _fi, VEC3D& _gi, 
	pointMass* jp, VEC3D& _spj, VEC3D& _fj, VEC3D& _gj)
	: kinematicConstraint(_md, _nm, kt, ip, _spi, _fi, _gi, jp, _spj, _fj, _gj)
{

}

translationalConstraint::~translationalConstraint()
{

}

void translationalConstraint::calculation_reaction_force(double ct)
{
// 	SMATD jaco;
// 	VECD out;
// 	out.alloc(14);
// 	jaco.alloc(52, 5, 14);
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

void translationalConstraint::constraintEquation(double m, double* rhs)
{
	if (ib->NumDOF() == DIM2)
	{

	}
	else
	{
		VEC3D v3 = jb->toGlobal(hj);
		VEC3D v3g = ib->toGlobal(gi);
		VEC3D v3f = ib->toGlobal(fi);
		rhs[srow + 0] = m * v3.dot(v3f);
		rhs[srow + 1] = m * v3.dot(v3g);
		v3 = jb->Position() + jb->toGlobal(spj) - ib->Position();// -ib->toGlobal(kconst->sp_i());
		rhs[srow + 2] = m * (v3.dot(v3f) - spi.dot(fi));
		rhs[srow + 3] = m * (v3.dot(v3g) - spi.dot(gi));
		rhs[srow + 4] = m * v3f.dot(jb->toGlobal(fj));
	}
}

void translationalConstraint::constraintEquation(double *rhs, double* q, double* dq)
{
	int i = ib->ID() * 7;
	int j = jb->ID() * 7;
	VEC3D ri = i < 0 ? ib->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
	VEC3D rj = j < 0 ? jb->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
	EPD epi = i < 0 ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD epj = j < 0 ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	VEC3D v3 = epj.toGlobal(hj);
	VEC3D v3g = epi.toGlobal(gi);
	VEC3D v3f = epi.toGlobal(fi);
	rhs[srow + 0] = v3.dot(v3f);
	rhs[srow + 1] = v3.dot(v3g);
	v3 = rj + epj.toGlobal(spj) - ri;// jb->Position() + jb->toGlobal(spj) - ib->Position();// -ib->toGlobal(kconst->sp_i());
	rhs[srow + 2] = (v3.dot(v3f) - spi.dot(fi));
	rhs[srow + 3] = (v3.dot(v3g) - spi.dot(gi));
	rhs[srow + 4] = v3f.dot(epj.toGlobal(fj));
}

void translationalConstraint::constraintJacobian(SMATD& cjaco)
{
	if (ib->NumDOF() == DIM2)
	{

	}
	else
	{
		VEC3D dij = (jb->Position() + jb->toGlobal(spj)) - (ib->Position() + ib->toGlobal(spi));
		//qDebug() << "constraint : " << nm;
		if (ib->MassType() != pointMass::GROUND)
		{
			//qDebug() << "i body : " << ib->Name();
			
			cjaco.extraction(srow + 0, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ib->getEP(), fi))), VEC4);
			cjaco.extraction(srow + 1, icol + 3, POINTER(transpose(jb->toGlobal(hj), B(ib->getEP(), gi))), VEC4);
			cjaco.extraction(srow + 2, icol + 0, POINTER((-fi)), POINTER(transpose(dij + ib->toGlobal(spi), B(ib->getEP(), fi))), VEC3_4);
			cjaco.extraction(srow + 3, icol + 0, POINTER((-gi)), POINTER(transpose(dij + ib->toGlobal(spi), B(ib->getEP(), gi))), VEC3_4);
			cjaco.extraction(srow + 4, icol + 3, POINTER(transpose(jb->toGlobal(fj), B(ib->getEP(), fi))), VEC4);
		}
		if (jb->MassType() != pointMass::GROUND)
		{
			//qDebug() << "j body : " << jb->Name();
			cjaco.extraction(srow + 0, jcol + 3, POINTER(transpose(ib->toGlobal(fi), B(jb->getEP(), hj))), VEC4);
			cjaco.extraction(srow + 1, jcol + 3, POINTER(transpose(ib->toGlobal(gi), B(jb->getEP(), hj))), VEC4);
			cjaco.extraction(srow + 2, jcol + 0, POINTER(ib->toGlobal(fi)), POINTER(transpose(ib->toGlobal(fi), B(jb->getEP(), spj))), VEC3_4);
			cjaco.extraction(srow + 3, jcol + 0, POINTER(ib->toGlobal(gi)), POINTER(transpose(ib->toGlobal(gi), B(jb->getEP(), spj))), VEC3_4);
			cjaco.extraction(srow + 4, jcol + 3, POINTER(transpose(ib->toGlobal(fi), B(jb->getEP(), fj))), VEC4);
		}
	}
}

void translationalConstraint::constraintJacobian(SMATD& cjaco, double* q, double* dq)
{
	int i = ib->ID() * 7;
	int j = jb->ID() * 7;
	VEC3D ri = i < 0 ? ib->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
	VEC3D rj = j < 0 ? jb->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
	EPD epi = i < 0 ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD epj = j < 0 ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	VEC3D dij = (rj + epj.toGlobal(spj)) - (ri + epi.toGlobal(spi));
	//qDebug() << "constraint : " << nm;
	if (ib->MassType() != pointMass::GROUND)
	{
		//qDebug() << "i body : " << ib->Name();

		cjaco.extraction(srow + 0, icol + 3, POINTER(transpose(epj.toGlobal(hj), B(epi, fi))), VEC4);
		cjaco.extraction(srow + 1, icol + 3, POINTER(transpose(epj.toGlobal(hj), B(epi, gi))), VEC4);
		cjaco.extraction(srow + 2, icol + 0, POINTER((-fi)), POINTER(transpose(dij + epi.toGlobal(spi), B(epi, fi))), VEC3_4);
		cjaco.extraction(srow + 3, icol + 0, POINTER((-gi)), POINTER(transpose(dij + epi.toGlobal(spi), B(epi, gi))), VEC3_4);
		cjaco.extraction(srow + 4, icol + 3, POINTER(transpose(epj.toGlobal(fj), B(epi, fi))), VEC4);
	}
	if (jb->MassType() != pointMass::GROUND)
	{
		//qDebug() << "j body : " << jb->Name();
		cjaco.extraction(srow + 0, jcol + 3, POINTER(transpose(epi.toGlobal(fi), B(epj, hj))), VEC4);
		cjaco.extraction(srow + 1, jcol + 3, POINTER(transpose(epi.toGlobal(gi), B(epj, hj))), VEC4);
		cjaco.extraction(srow + 2, jcol + 0, POINTER(epi.toGlobal(fi)), POINTER(transpose(epi.toGlobal(fi), B(epj, spj))), VEC3_4);
		cjaco.extraction(srow + 3, jcol + 0, POINTER(epi.toGlobal(gi)), POINTER(transpose(epi.toGlobal(gi), B(epj, spj))), VEC3_4);
		cjaco.extraction(srow + 4, jcol + 3, POINTER(transpose(epi.toGlobal(fi), B(epj, fj))), VEC4);
	}
}

void translationalConstraint::derivate(MATD& lhs, double mul)
{
	MAT44D Dv;
	MAT34D Bv;
// 	bool ig = ib->MassType() != pointMass::GROUND;
// 	bool jg = jb->MassType() != pointMass::GROUND;
	VEC3D dij = (jb->Position() + jb->toGlobal(spj)) - (ib->Position() + ib->toGlobal(spi));

	Dv = lm[0] * D(fi, jb->toGlobal(hj));
	Dv += lm[1] * D(gi, jb->toGlobal(hj));
	Dv += lm[2] * D(fi, dij + ib->toGlobal(spi));
	Dv += lm[3] * D(gi, dij + ib->toGlobal(spi));
	Dv += lm[4] * D(fi, jb->toGlobal(fj));
	lhs.plus(icol + 3, icol + 3, POINTER(Dv), MAT4x4, mul);
	Bv = -lm[2] * B(ib->getEP(), fi) - lm[3] * B(ib->getEP(), gi);
	lhs.plus(icol, icol + 3, POINTER(Bv), MAT3X4, mul);
	lhs.plus(icol + 3, icol, POINTER(Bv), MAT4x3, mul);

	Bv = -Bv;
	lhs.plus(icol + 3, jcol, POINTER(Bv), MAT4x3, mul);
	lhs.plus(jcol, icol + 3, POINTER(Bv), MAT3X4, mul);

	Dv = lm[0] * D(hj, ib->toGlobal(fi));
	Dv += lm[1] * D(hj, ib->toGlobal(gi));
	Dv += lm[2] * D(spj, ib->toGlobal(fi));
	Dv += lm[3] * D(spj, ib->toGlobal(gi));
	Dv += lm[4] * D(fj, ib->toGlobal(fi));
	lhs.plus(jcol + 3, jcol + 3, POINTER(Dv), MAT4x4, mul);

	Dv = lm[0] * transpose(B(ib->getEP(), fi), B(jb->getEP(), hj));
	Dv += lm[1] * transpose(B(ib->getEP(), gi), B(jb->getEP(), hj));
	Dv += lm[2] * transpose(B(ib->getEP(), fi), B(jb->getEP(), spj));
	Dv += lm[3] * transpose(B(ib->getEP(), gi), B(jb->getEP(), spj));
	Dv += lm[4] * transpose(B(ib->getEP(), fi), B(jb->getEP(), fj));
	lhs.plus(icol + 3, jcol + 3, POINTER(Dv), MAT4x4, mul);

	Dv = lm[0] * transpose(B(jb->getEP(), hj), B(ib->getEP(), fi));
	Dv += lm[1] * transpose(B(jb->getEP(), hj), B(ib->getEP(), gi));
	Dv += lm[2] * transpose(B(jb->getEP(), spj), B(ib->getEP(), fi));
	Dv += lm[3] * transpose(B(jb->getEP(), spj), B(ib->getEP(), gi));
	Dv += lm[4] * transpose(B(jb->getEP(), fj), B(ib->getEP(), fi));
	lhs.plus(jcol + 3, icol + 3, POINTER(Dv), MAT4x4, mul);
}

void translationalConstraint::differentialEquation(double *rhs, double* q, double *dq, double t)
{
	bool ig = ib->MassType() == pointMass::GROUND;
	bool jg = jb->MassType() == pointMass::GROUND;
	int i = ig ? -1 : ib->ID() * 7;
	int j = jg ? -1 : jb->ID() * 7;
	VEC3D ri = ig ? ib->Position() : VEC3D(q[i + 0], q[i + 1], q[i + 2]);
	VEC3D rj = jg ? jb->Position() : VEC3D(q[j + 0], q[j + 1], q[j + 2]);
	VEC3D dri = ig ? ib->getVelocity() : VEC3D(dq[i + 0], q[i + 1], q[i + 2]);
	VEC3D drj = jg ? jb->getVelocity() : VEC3D(dq[j + 0], q[j + 1], q[j + 2]);
	EPD pi = ig ? ib->getEP() : EPD(q[i + 3], q[i + 4], q[i + 5], q[i + 6]);
	EPD pj = jg ? jb->getEP() : EPD(q[j + 3], q[j + 4], q[j + 5], q[j + 6]);
	EPD dpi = ig ? ib->getEV() : EPD(dq[i + 3], dq[i + 4], dq[i + 5], dq[i + 6]);
	EPD dpj = jg ? jb->getEV() : EPD(dq[j + 3], dq[j + 4], dq[j + 5], dq[j + 6]);
	VEC3D dij = rj + pj.toGlobal(spj) - ri - pi.toGlobal(spi);
	double d1_0 = dot_1_differential(fi, hj, pi, pj, dpi, dpj);
	double d1_1 = dot_1_differential(gi, hj, pi, pj, dpi, dpj);
	double d2_0 = dot_2_differential(fi, dri, drj, dij, pi, pj, dpi, dpj);
	double d2_1 = dot_2_differential(gi, dri, drj, dij, pi, pj, dpi, dpj);
	double d1_2 = dot_1_differential(fi, fj, pi, pj, dpi, dpj);
	rhs[srow + 0] = d1_0;
	rhs[srow + 1] = d1_1;
	rhs[srow + 2] = d2_0;
	rhs[srow + 3] = d2_1;
	rhs[srow + 4] = d1_2;
}
