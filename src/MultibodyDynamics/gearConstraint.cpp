#include "gearConstraint.h"
#include "artificialCoordinate.h"
#include "simulation.h"

gearConstraint::gearConstraint()
	: kinematicConstraint()
	, ratio(0)
	, ki(0)
	, kj(0)
	, dthi(0)
	, dthj(0)
	, fmi(NULL)
	, fmj(NULL)
	, smi(NULL)
	, smj(NULL)
	, fac(NULL)
	, sac(NULL)
{

}

gearConstraint::gearConstraint(QString _nm, pointMass* i, kinematicConstraint* ik, pointMass* j, kinematicConstraint* jk, double r)
	: kinematicConstraint()
	, ratio(r)
	, ki(0)
	, kj(0)
	, dthi(0)
	, dthj(0)
	, fmi(NULL)
	, fmj(NULL)
	, smi(NULL)
	, smj(NULL)
	, fac(NULL)
	, sac(NULL)
{
	kinematicConstraint::nm = _nm;
	kinematicConstraint::type = kinematicConstraint::GEAR;
	kinematicConstraint::nconst = 1;
	kinematicConstraint::maxnnz = 8;
	kinematicConstraint::ib = i;
	kinematicConstraint::jb = j;
	kinematicConstraint::fi = VEC3D(0.0, 0.0, 1.0);// ik->f_i();
	kinematicConstraint::fj = VEC3D(0.0, 0.0, 1.0);// jk->f_j();
	kinematicConstraint::gi = VEC3D(0.0, -1.0, 0.0);
	kinematicConstraint::gj = VEC3D(0.0, 1.0, 0.0);
	kinematicConstraint::hi = ik->h_i();
	kinematicConstraint::hj = -kinematicConstraint::hi;

	init_f_i = VEC3D(0.0, 0.0, 1.0);
	init_f_j = VEC3D(0.0, 0.0, 1.0);

	init_g_i = VEC3D(0.0, 0.0, -1.0);// i->toGlobal(-kinematicConstraint::gi);
	init_g_j = VEC3D(0.0, 0.0, -1.0); //j->toGlobal(-kinematicConstraint::gj);

	ki = 0;
	kj = 0;
	nri = 0;
	nrj = 0;

	memset(&roti, 0, sizeof(rotationInfo));
	memset(&rotj, 0, sizeof(rotationInfo));
}

gearConstraint::gearConstraint(QString _nm)
	: kinematicConstraint(_nm)
	, fmi(NULL)
	, fmj(NULL)
	, smi(NULL)
	, smj(NULL)
	, fac(NULL)
	, sac(NULL)
{
	kinematicConstraint::type = kinematicConstraint::GEAR;
}

gearConstraint::~gearConstraint()
{

}

// void gearConstraint::constraintEquation(int sr, double* rhs, double mul)
// {
// 	//EPD pi = EPD()
// // 	double dth = i->getEP().dot(j->getEP());
//  //	double thi = 2.0 * acos(i->getEP().e0);
//  //	double thj = 2.0 * acos(j->getEP().e0);
// // 	VEC3D Gi = i->toGlobal(gi);
// // 	VEC3D Gj = j->toGlobal(gj);
// //  	double ti = i->toGlobal(gi).dot(init_f_i);
// // 	double tj = j->toGlobal(gj).dot(init_f_j);// init_g_j.dot(j->toGlobal(fj));
// // 	if (abs(tj) > 1.0)
// // 		tj = sign2(tj);
// 	//dthi = /*asin(ti);*/ asin2(ti, ki);
// 	//dthj = /*asin(tj);*/ asin2(tj, kj);
// 	//double th = acos(i->toGlobal(fi).dot(j->toGlobal(fj)));//2.0 * acos(kinematicConstraint::i->getEP().e0);
// 	//double thj = 2.0 * acos(kinematicConstraint::j->getEP().e0);
// // 	EPD p1 = EPD(cos(th), sin(th), 0, 0);
// // 	EPD p2 = EPD(cos(ratio * th), sin(ratio * th), 0, 0);
// 	//double v = ratio*dthi -  dthj;//p1.dot(p2) - cos(th - ratio * th);
// 	EPD fPi = fmi->getEP();
// 	EPD fPj = fmj->getEP();
// 	EPD sPi = smi->getEP();
// 	EPD sPj = smj->getEP();
// 	double fd = fPi.dot(fPj);
// 	double fth = 2.0 * acos(fd);
// 	double sd = sPi.dot(sPj);
// 	double sth = 2.0 * acos(sd);
// 	double v = ratio * (fth - 0) + (sth - 0);
// 	rhs[sr] = mul * v;
// 
// 	unsigned int loc_fac = fac->getMatrixLocation();
// 	unsigned int loc_sac = sac->getMatrixLocation();
// 	rhs[sr + 1] = mul * (fd - cos(0.5 * fth));
// 	rhs[sr + 2] = mul * (sd - cos(0.5 * sth));
// // 	dthi = i->Position().z;
// // 	dthj = j->Position().z;
// // 	double v = (dthi - 0) * ratio + (dthj - 0);
// // 	rhs[sr] = mul * v;// (thj //i->getEP().dot(j->getEP()) - cos(0.5 * v);
// }
// 
// void gearConstraint::constraintJacobian(int sr, SMATD& cjaco)
// {
// // 	double thi = 2.0 * acos(kinematicConstraint::i->getEP().e0);
// // 	double thj = 2.0 * acos(kinematicConstraint::j->getEP().e0);
// // 
// // 	dthi = 0.0;
// // 	dthj = 0.0;
// // 	double lowi = 0.0;
// // 	double lowi1 = 0.0;
// // 	double lowj = 0.0;
// // 	double lowj1 = 0.0;
// // // 	EPD e = EPD(1.0, 0.0, 0.0, 0.0);
// // // 	double v = (thj / ratio) - thi;// / (ratio;// thi + ratio * thj;
// // 	VEC3D GJ = j->toGlobal(gj);
// // 	VEC3D GI = i->toGlobal(gi);
// // 
// // 	
// // 	double costhi = i->toGlobal(fi).dot(g->toGlobal(fj));
// // 	double sinthi = i->toGlobal(gi).dot(g->toGlobal(fj));
// // 
// // 	double costhj = j->toGlobal(fj).dot(g->toGlobal(fi));
// // 	double sinthj = j->toGlobal(gj).dot(g->toGlobal(fi));
// // 
// // 	lowi = abs(sinthi) != 1.0 ? 1.0 / sqrt(1.0 - pow(sinthi, 2.0)) : 0.0;
// // 	lowj = abs(sinthj) != 1.0 ? 1.0 / sqrt(1.0 - pow(sinthj, 2.0)) : 0.0;
// // 	double si = sign2(costhi);
// // 	double sj = sign2(costhj);
// // 	//double dlowj = sinthj - pre_sinthj
// // 	VEC4D v4i = ratio * si * lowi * transpose(g->toGlobal(fi), B(i->getEP(), gi));// costhi * transpose(g->toGlobal(fi), B(i->getEP(), gi));// 2.0 * transpose(g->toGlobal(fj), costhi * B(i->getEP(), gi) - sinthi * B(i->getEP(), fi));
// // 	VEC4D v4j = -1.0 * sj * lowj * transpose(g->toGlobal(fi), B(j->getEP(), gj));// -2.0 * costhj * transpose(g->toGlobal(fi), B(i->getEP(), gi));
// // // 	if (abs(ti) > 0.99)
// // // 	{
// // // 		double s = sign2(ti);
// // // 		dthj = 0.5 * M_PI - (0.5 * M_PI - asin(tj));
// // 		//ti = s - s*(s*ti - 1.0);
// // //	}
// // 	double pp = i->getEP().dot(i->getEP());
// // //	dthi = asin(ti);// asin2(ti, ki);
// // 	//dthj = asin(tj);// asin2(tj, kj);
// // 	double sgi = 1.0;
// // 	double sgj = 1.0;
// // 	double ssgj = 1.0;
// // 	if (sinthi >= 0 && costhi >= 0)
// // 	{
// // 		if (ki == 3)
// // 			nri++;
// // 		dthi = asin(sinthi);
// // 		ki = 0;
// // 	}
// // 	else if (sinthi >= 0 && costhi < 0)
// // 	{
// // 		dthi = M_PI - asin(sinthi);
// // 		ki = 1;
// // 	}
// // 	else if (sinthi < 0 && costhi < 0)
// // 	{
// // 		dthi = M_PI - asin(sinthi);
// // 		ki = 2;
// // 	}
// // 	else if (sinthi < 0 && costhi >= 0)
// // 	{
// // 		if (ki == 0)
// // 			nri--;
// // 		dthi = 2.0 * M_PI + asin(sinthi);
// // 		ki = 3;
// // 	}
// // 		
// // 	if (nri)
// // 		dthi += nri*2.0* M_PI;
// // 
// // 	if (sinthj >= 0 && costhj >= 0)
// // 	{
// // 		if (kj == 3)
// // 			nrj++;
// // 		dthj = asin(sinthj);
// // 		kj = 0;
// // 	}	
// // 	else if (sinthj >= 0 && costhj < 0)
// // 	{
// // 		dthj = M_PI - asin(sinthj);
// // 		kj = 1;
// // 	}
// // 	else if (sinthj < 0 && costhj < 0)
// // 	{
// // 		dthj = M_PI - asin(sinthj);
// // 		kj = 2;
// // 	}
// // 	else if (sinthj < 0 && costhj >= 0)
// // 	{
// // 		if (kj == 0)
// // 			nrj--;
// // 		dthj = 2.0 * M_PI + asin(sinthj);
// // 		kj = 3;
// // 	}
// // 
// // 	if (nrj)
// // 		dthj += nrj*2.0* M_PI;
// // 		
// // 	unsigned int ic = sr ? i->ID() * 7 : 0;
// // 	//v4j.y *= ssgj;
// // 	//v4i.x *= ssgj;
// // 	unsigned int jc = sr ? j->ID() * 7 : 0;
// // 	cjaco.extraction(sr, ic + 3, POINTER(v4i), VEC4);
// // 	cjaco.extraction(sr, jc + 3, POINTER(v4j), VEC4);
// 	unsigned int loc_fac = fac->getMatrixLocation();
// 	unsigned int loc_sac = sac->getMatrixLocation();
// 
// // 	int ic = i->ID() * 3 + 2;
// // 	int jc = j->ID() * 3 + 2;
// 	int c = 0;
// 	EPD p;
// 	cjaco(sr, loc_fac) = ratio;
// 	cjaco(sr, loc_sac) = 1.0;
// 	if (fmi->MassType() != pointMass::GROUND)
// 	{
// 		c = fmi->ID() * 7 + 3;
// 		p = fmj->getEP();
// 		cjaco(sr + 1, c + 0) = p.e0;
// 		cjaco(sr + 1, c + 1) = p.e1;
// 		cjaco(sr + 1, c + 2) = p.e2;
// 		cjaco(sr + 1, c + 3) = p.e3;
// 		cjaco(sr + 2, c + 0) = p.e0;
// 		cjaco(sr + 2, c + 1) = p.e1;
// 		cjaco(sr + 2, c + 2) = p.e2;
// 		cjaco(sr + 2, c + 3) = p.e3;
// 	}
// 	if (fmj->MassType() != pointMass::GROUND)
// 	{
// 		c = fmj->ID() * 7 + 3;
// 		p = fmi->getEP();
// 		cjaco(sr + 1, c + 0) = p.e0;
// 		cjaco(sr + 1, c + 1) = p.e1;
// 		cjaco(sr + 1, c + 2) = p.e2;
// 		cjaco(sr + 1, c + 3) = p.e3;
// 		cjaco(sr + 2, c + 0) = p.e0;
// 		cjaco(sr + 2, c + 1) = p.e1;
// 		cjaco(sr + 2, c + 2) = p.e2;
// 		cjaco(sr + 2, c + 3) = p.e3;
// 	}
// 	if (smi->MassType() != pointMass::GROUND)
// 	{
// 		c = smi->ID() * 7 + 3;
// 		p = smj->getEP();
// 		cjaco(sr + 1, c + 0) = p.e0;
// 		cjaco(sr + 1, c + 1) = p.e1;
// 		cjaco(sr + 1, c + 2) = p.e2;
// 		cjaco(sr + 1, c + 3) = p.e3;
// 		cjaco(sr + 2, c + 0) = p.e0;
// 		cjaco(sr + 2, c + 1) = p.e1;
// 		cjaco(sr + 2, c + 2) = p.e2;
// 		cjaco(sr + 2, c + 3) = p.e3;
// 	}
// 	if (smj->MassType() != pointMass::GROUND)
// 	{
// 		c = smj->ID() * 7 + 3;
// 		p = smi->getEP();
// 		cjaco(sr + 1, c + 0) = p.e0;
// 		cjaco(sr + 1, c + 1) = p.e1;
// 		cjaco(sr + 1, c + 2) = p.e2;
// 		cjaco(sr + 1, c + 3) = p.e3;
// 		cjaco(sr + 2, c + 0) = p.e0;
// 		cjaco(sr + 2, c + 1) = p.e1;
// 		cjaco(sr + 2, c + 2) = p.e2;
// 		cjaco(sr + 2, c + 3) = p.e3;
// 	}
// 
// // 	cjaco(sr, ic) = ratio;
// // 	cjaco(sr, jc) = 1.0;
// }

void gearConstraint::constraintEquation(double m, double* rhs)
{
	dthi = ib->Position().z;
	dthj = jb->Position().z;
	double v = (dthi - 0) * ratio + (dthj - 0);
	rhs[srow] = m * v;
}

void gearConstraint::constraintJacobian(SMATD& cjaco)
{
	cjaco(srow, icol) = ratio;
	cjaco(srow, jcol) = 1.0;
}

void gearConstraint::constraintJacobian(SMATD& cjaco, double* q, double* dq)
{

}

void gearConstraint::derivate(MATD& lhs, double mul)
{

}

void gearConstraint::saveGearConstraintData(QTextStream& qts)
{
	qts << endl
		<< "GEAR_CONSTRAINT_DATA" << endl
		<< "NAME " << nm << endl
		<< "FIRST_I_BODY " << fmi->Name() << endl
		<< "FIRST_J_BODY " << fmj->Name() << endl
		<< "SECOND_I_BODY " << smi->Name() << endl
		<< "SECOND_J_BODY " << smj->Name() << endl
		<< "FIRST_SP " << fs.x << " " << fs.y << " " << fs.z << endl
		<< "SECOND_SP " << ss.x << " " << ss.y << " " << ss.z << endl
		<< "GEAR_RATIO " << ratio << endl
		<< "CONTACT_POINT " << cpoint.x << " " << cpoint.y << " " << cpoint.z << endl;
}

void gearConstraint::decisionPreviousData()
{
// 	pre_thi = dthi;
// 	pre_thj = dthj;

	///double costhj = j->toGlobal(fj).dot(g->toGlobal(fi));
	///double sinthj = j->toGlobal(gj).dot(g->toGlobal(fi));
}
