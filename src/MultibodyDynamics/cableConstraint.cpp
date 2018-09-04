#include "cableConstraint.h"
#include "pointMass.h"
#include "model.h"

cableConstraint::cableConstraint(QString _name, pointMass* _fi, VEC3D& _fspi, pointMass* _fj, VEC3D& _fspj, pointMass* _si, VEC3D& _sspi, pointMass* _sj, VEC3D& _sspj)
	: kinematicConstraint()
	, fi(_fi)
	, fj(_fj)
	, si(_si)
	, sj(_sj)
	, fspi(_fspi)
	, fspj(_fspj)
	, sspi(_sspi)
	, sspj(_sspj)
	, f_c(0)
	, s_c(0)
	, ispulling(false)
	, pre_f_c(0)
{
	kinematicConstraint::nm = _name;
	kinematicConstraint::type = kinematicConstraint::CABLE;
	kinematicConstraint::nconst = 1;
	kinematicConstraint::maxnnz = 45;

	f_c = (fj->getPosition() + fj->toGlobal(fspj) - fi->getPosition() - fi->toGlobal(fspi)).length();
	s_c = (sj->getPosition() + sj->toGlobal(sspj) - si->getPosition() - si->toGlobal(sspi)).length();

	pre_f_c = f_c;
}

cableConstraint::~cableConstraint()
{

}

void cableConstraint::setApplyCableForceMass(pointMass* acm)
{
	a_cf_mass = acm;
}

// void cableConstraint::constraintEquation(int sr, double* rhs, double mul)
// {
// 	VEC3D fdij = fj->getPosition() + fj->toGlobal(fspj) - fi->getPosition() - fi->toGlobal(fspi);
// 	double fc = fdij.length() - f_c;
// 	VEC3D sdij = sj->getPosition() + sj->toGlobal(sspj) - si->getPosition() - si->toGlobal(sspi);
// 	double sc = sdij.length() - s_c;
// 	rhs[sr] = mul * (fc + sc);
// }

void cableConstraint::constraintEquation(double m, double* rhs)
{
	VEC3D fdij = fj->getPosition() + fj->toGlobal(fspj) - fi->getPosition() - fi->toGlobal(fspi);
	double fc = fdij.length() - f_c;
	VEC3D sdij = sj->getPosition() + sj->toGlobal(sspj) - si->getPosition() - si->toGlobal(sspi);
	double sc = sdij.length() - s_c;
	rhs[srow] = m * (fc + sc);
}

void cableConstraint::constraintJacobian(SMATD& cjaco)
{
	//updateCableInitLength();
	VEC3D fdij = fj->getPosition() + fj->toGlobal(fspj) - fi->getPosition() - fi->toGlobal(fspi);
	VEC3D sdij = sj->getPosition() + sj->toGlobal(sspj) - si->getPosition() - si->toGlobal(sspi);
	//int ic = 0, jc = 0;
	// 	if (fdij.length() > f_c)
	// 		ispulling = true;
	// 	else
	// 		ispulling = false;
	VEC3D D1;
	VEC4D D2;
	double div = 0;
	if (fi->MassType() != pointMass::GROUND)
	{
		//ic = srow ? fi->ID() * 7 : 0;
		div = 1.0 / fdij.length();
		D1 = div * (-fdij);
		D2 = div * transpose(-fdij, B(fi->getEP(), fspi));
		cjaco.extraction(srow, icol, POINTER(D1), POINTER(D2), VEC3_4);
	}
	if (fj->MassType() != pointMass::GROUND)
	{
		//jc = sr ? fj->ID() * 7 : 0;
		div = 1.0 / fdij.length();
		D1 = div * fdij;
		D2 = div * transpose(fdij, B(fj->getEP(), fspj));
		cjaco.extraction(srow, jcol, POINTER(D1), POINTER(D2), VEC3_4);
	}
	if (si->MassType() != pointMass::GROUND)
	{
		//ic = sr ? si->ID() * 7 : 0;
		div = 1.0 / sdij.length();
		D1 = div * (-sdij);
		D2 = div * transpose(-sdij, B(si->getEP(), sspi));
		cjaco.extraction(srow, icol, POINTER(D1), POINTER(D2), VEC3_4);
	}
	if (sj->MassType() != pointMass::GROUND)
	{
		//jc = sr ? sj->ID() * 7 : 7;
		div = 1.0 / sdij.length();
		D1 = div * sdij;
		D2 = div * transpose(sdij, B(sj->getEP(), sspj));
		cjaco.extraction(srow, jcol, POINTER(D1), POINTER(D2), VEC3_4);
	}
}

// void cableConstraint::constraintJacobian(int sr, SMATD& cjaco)
// {
// 	//updateCableInitLength();
// 	VEC3D fdij = fj->getPosition() + fj->toGlobal(fspj) - fi->getPosition() - fi->toGlobal(fspi);
// 	VEC3D sdij = sj->getPosition() + sj->toGlobal(sspj) - si->getPosition() - si->toGlobal(sspi);
// 	int ic = 0, jc = 0;
// // 	if (fdij.length() > f_c)
// // 		ispulling = true;
// // 	else
// // 		ispulling = false;
// 	VEC3D D1;
// 	VEC4D D2;
// 	double div = 0;
// 	if (fi->MassType() != pointMass::GROUND)
// 	{
// 		ic = sr ? fi->ID() * 7 : 0;
// 		div = 1.0 / fdij.length();
// 		D1 = div * (-fdij);
// 		D2 = div * transpose(-fdij, B(fi->getEP(), fspi));
// 		cjaco.extraction(sr, ic, POINTER(D1), POINTER(D2), VEC3_4);
// 	}
// 	if (fj->MassType() != pointMass::GROUND)
// 	{
// 		jc = sr ? fj->ID() * 7 : 0;
// 		div = 1.0 / fdij.length();
// 		D1 = div * fdij;
// 		D2 = div * transpose(fdij, B(fj->getEP(), fspj));
// 		cjaco.extraction(sr, jc, POINTER(D1), POINTER(D2), VEC3_4);
// 	}
// 	if (si->MassType() != pointMass::GROUND)
// 	{
// 		ic = sr ? si->ID() * 7 : 0;
// 		div = 1.0 / sdij.length();
// 		D1 = div * (-sdij);
// 		D2 = div * transpose(-sdij, B(si->getEP(), sspi));
// 		cjaco.extraction(sr, ic, POINTER(D1), POINTER(D2), VEC3_4);
// 	}
// 	if (sj->MassType() != pointMass::GROUND)
// 	{
// 		jc = sr ? sj->ID() * 7 : 7;
// 		div = 1.0 / sdij.length();
// 		D1 = div * sdij;
// 		D2 = div * transpose(sdij, B(sj->getEP(), sspj));
// 		cjaco.extraction(sr, jc, POINTER(D1), POINTER(D2), VEC3_4);
// 	}
// }

void cableConstraint::saveCableConstraintData(QTextStream& qts)
{
	qts << endl
		<< "CABLE_CONSTRAINT_DATA" << endl
		<< "NAME " << nm << endl
		<< "FIRST_I_BODY " << fi->Name() << endl
		<< "FIRST_J_BODY " << fj->Name() << endl
		<< "SECOND_I_BODY " << si->Name() << endl
		<< "SECOND_J_BODY " << sj->Name() << endl
		<< "FIRST_I_SP " << fspi.x << " " << fspi.y << " " << fspi.z << endl
		<< "FIRST_J_SP " << fspj.x << " " << fspj.y << " " << fspj.z << endl
		<< "SECOND_I_SP " << sspi.x << " " << sspi.y << " " << sspi.z << endl
		<< "SECOND_J_SP " << sspj.x << " " << sspj.y << " " << sspj.z << endl
		<< "FIRST_INIT_LENGTH " << f_c << endl
		<< "SECOND_INIT_LENGTH " << s_c << endl;
}

void cableConstraint::updateCableInitLength()
{
	f_c = (fj->getPosition() + fj->toGlobal(fspj) - fi->getPosition() - fi->toGlobal(fspi)).length();
	s_c = (sj->getPosition() + sj->toGlobal(sspj) - si->getPosition() - si->toGlobal(sspi)).length();
	pre_f_c = f_c;
}

void cableConstraint::calculation_reaction_force2body()
{
	SMATD jaco;
	VECD out;
	out.alloc(14);
	jaco.alloc(14, 1, 14);
	constraintJacobian(jaco);
	double cur_f_c = (fj->getPosition() + fj->toGlobal(fspj) - fi->getPosition() - fi->toGlobal(fspi)).length();
	if (cur_f_c > pre_f_c)
	{
		for (unsigned int i = 0; i < jaco.nnz(); i++)
		{
			out(jaco.cidx[i]) += jaco.value[i] * lm[jaco.ridx[i] - srow];
		}
		VEC3D cf = VEC3D(out(0), out(1), out(2));
		a_cf_mass->setExternalMoment(0.4 * VEC3D(0, 0, cf.length()));
	}
	else
	{
		a_cf_mass->setExternalMoment(VEC3D(0.0, 0.0, 0.0));
	}
	pre_f_c = cur_f_c;
	/*reactionForceResults.push_back(rfd)*/
}

void cableConstraint::calculation_reaction_force(double ct)
{
	SMATD jaco;
	VECD out;
	out.alloc(14);
	jaco.alloc(14, 1, 14);
	constraintJacobian(jaco);
	for (unsigned int i = 0; i < jaco.nnz(); i++)
	{
		out(jaco.cidx[i]) += jaco.value[i] * lm[jaco.ridx[i] - srow];
	}
	resultStorage::reactionForceData rfd =
	{
		ct,
		VEC3D(out(0), out(1), out(2)),
		VEC4D(out(3), out(4), out(5), out(6)),
		VEC3D(out(7), out(8), out(9)),
		VEC4D(out(10), out(11), out(12), out(13))
	};
	model::rs->insertReactionForceResult(nm, rfd);// .push_back(rfd);

	model::rs->insertGLineData(nm + "_f_0", fi->getPosition() + fi->toGlobal(fspi));
	model::rs->insertGLineData(nm + "_f_1", fj->getPosition() + fj->toGlobal(fspj));
	model::rs->insertGLineData(nm + "_s_0", si->getPosition() + si->toGlobal(sspi));
	model::rs->insertGLineData(nm + "_s_1", sj->getPosition() + sj->toGlobal(sspj));
}
