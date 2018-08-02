#include "mbd_model.h"
#include "GLWidget.h"
//#include "cube.h"
#include "pointMass.h"
#include "forceElement.h"
//#include "plane.h"
#include "fixedConstraint.h"
#include "revoluteConstraint.h"
#include "sphericalConstraint.h"
#include "translationalConstraint.h"
#include "axialRotationForce.h"
#include "universalConstraint.h"
//#include "drivingConstraint.h"
//#include "gearConstraint.h"
#include "springDamperModel.h"
#include <direct.h>
#include <QString>

mbd_model::mbd_model()
	: model()
	, model_path("")
	, ground(NULL)
	, is2D(false)
{
	ground = new pointMass(QString("ground"), pointMass::GROUND);
}

mbd_model::mbd_model(QString _name)
	: model()
	, model_path("")
	, ground(NULL)
	, is2D(false)
{
	ground = new pointMass(QString("ground"), pointMass::GROUND);
}

mbd_model::~mbd_model()
{

	if (ground) delete ground; ground = NULL;
	qDeleteAll(masses);
	qDeleteAll(consts);
	qDeleteAll(forces);
}

kinematicConstraint* mbd_model::createKinematicConstraint(
	QString _name, kinematicConstraint::Type kt,
	pointMass* i, VEC3D& spi, VEC3D& fi, VEC3D& gi,
	pointMass* j, VEC3D& spj, VEC3D& fj, VEC3D& gj)
{
	kinematicConstraint* kin;
	switch (kt)
	{
	case kinematicConstraint::FIXED:
		kin = new fixedConstraint(this, _name, kt, i, spi, fi, gi, j, spj, fj, gj);
		break;
	case kinematicConstraint::REVOLUTE:
		kin = new revoluteConstraint(this, _name, kt, i, spi, fi, gi, j, spj, fj, gj);
		break;
	case kinematicConstraint::SPHERICAL:
		kin = new sphericalConstraint(this, _name, kt, i, spi, fi, gi, j, spj, fj, gj);
		break;
	case kinematicConstraint::TRANSLATIONAL:
		kin = new translationalConstraint(this, _name, kt, i, spi, fi, gi, j, spj, fj, gj);
		break;
	case kinematicConstraint::UNIVERSAL:
		kin = new universalConstraint(this, _name, kt, i, spi, fi, gi, j, spj, fj, gj);
		break;
	}
	consts[_name] = kin;
	return kin;
}

kinematicConstraint* mbd_model::createKinematicConstraint(QTextStream& qts)
{
	QString ch;
	QString kin_name;
	QString base_name, action_name;
	VEC3D loc;
	VEC3D spi, fi, gi;
	VEC3D spj, fj, gj;
	int tp;
	qts >> ch >> kin_name >> ch >> tp
		>> ch >> loc.x >> loc.y >> loc.z;
	qts >> ch >> base_name
		>> ch >> spi.x >> spi.y >> spi.z
		>> ch >> fi.x >> fi.y >> fi.z
		>> ch >> gi.x >> gi.y >> gi.z
		>> ch >> action_name
		>> ch >> spj.x >> spj.y >> spj.z
		>> ch >> fj.x >> fj.y >> fj.z
		>> ch >> gj.x >> gj.y >> gj.z;
	pointMass *base = base_name == "ground" ? ground : masses[base_name];
	pointMass *action = action_name == "ground" ? ground : masses[action_name];
	kinematicConstraint* kin =  createKinematicConstraint(
		kin_name, (kinematicConstraint::Type)tp,
		base, spi, fi, gi,
		action, spj, fj, gj);
	kin->setLocation(loc);
	return kin;
}

kinematicConstraint* mbd_model::createCableConstraint(QString _name, pointMass* fi, VEC3D& fspi, pointMass* fj, VEC3D& fspj, pointMass* si, VEC3D& sspi, pointMass* sj, VEC3D& sspj)
{
	cableConstraint *cable = new cableConstraint(
		_name, fi, fspi, fj, fspj, si, sspi, sj, sspj);
	cables[_name] = cable;
	consts[_name] = cable;
	VEC3D sp = fi->getPosition() + fi->toGlobal(fspi);
	VEC3D ep = fj->getPosition() + fj->toGlobal(fspj);
	//GLWidget::GLObject()->createLine(_name + "_f", sp.x, sp.y, sp.z, ep.x, ep.y, ep.z);

	sp = si->getPosition() + si->toGlobal(sspi);
	ep = sj->getPosition() + sj->toGlobal(sspj);
	//GLWidget::GLObject()->createLine(_name + "_s", sp.x, sp.y, sp.z, ep.x, ep.y, ep.z);
	return cable;
}

kinematicConstraint* mbd_model::createCableConstraint(QTextStream& qts)
{
	QString pass;
	QString c_name;
	QString fi_name, fj_name, si_name, sj_name;
	VEC3D fspi, fspj, sspi, sspj;
	double fc = 0;
	double sc = 0;
	qts >> pass >> c_name
		>> pass >> fi_name
		>> pass >> fj_name
		>> pass >> si_name
		>> pass >> sj_name
		>> pass >> fspi.x >> fspi.y >> fspi.z
		>> pass >> fspj.x >> fspj.y >> fspj.z
		>> pass >> sspi.x >> sspi.y >> sspi.z
		>> pass >> sspj.x >> sspj.y >> sspj.z
		>> pass >> fc >> pass >> sc;
	pointMass *fi = fi_name == "ground" ? ground : masses[fi_name];
	pointMass *fj = fj_name == "ground" ? ground : masses[fj_name];
	pointMass *si = si_name == "ground" ? ground : masses[si_name];
	pointMass *sj = sj_name == "ground" ? ground : masses[sj_name];
	return createCableConstraint(c_name, fi, fspi, fj, fspj, si, sspi, sj, sspj);
}

kinematicConstraint* mbd_model::createGearConstraint(QString _name, pointMass* i, kinematicConstraint* ik, pointMass* j, kinematicConstraint* jk, double r)
{
	gearConstraint *gear = new gearConstraint(_name, i, ik, j, jk, r);
	gear->setGroundPointer(ground);
	gears[_name] = gear;
	consts[_name] = gear;
	return gear;
}

springDamperModel* mbd_model::createSpringDamperElement(
	QString _name,
	pointMass* i, VEC3D& bLoc, 
	pointMass* j, VEC3D& aLoc, 
	double k, double c)
{
	springDamperModel* fe = new springDamperModel(_name, this, i, bLoc, j, aLoc, k, c);
	forces[_name] = fe;
	return fe;
}

springDamperModel* mbd_model::createSpringDamperElement(QTextStream& qts)
{
	QString ch, nm;
	QString i_name, j_name;
	VEC3D iLoc, jLoc, iLocal, jLocal;
	double init_l, k, c;
	qts >> ch >> nm
		>> ch >> i_name
		>> ch >> j_name
		>> ch >> iLoc.x >> iLoc.y >> iLoc.z
		>> ch >> jLoc.x >> jLoc.y >> jLoc.z
		>> ch >> iLocal.x >> iLocal.y >> iLocal.z
		>> ch >> jLocal.x >> jLocal.y >> jLocal.z
		>> ch >> init_l >> k >> c;
	pointMass *fi = i_name == "ground" ? ground : masses[i_name];
	pointMass *fj = j_name == "ground" ? ground : masses[j_name];
	return createSpringDamperElement(nm, fi, iLoc, fj, jLoc, k, c);
}

axialRotationForce* mbd_model::createAxialRotationForce(
	QString _name, pointMass* i, pointMass* j, VEC3D loc, VEC3D u, double v)
{
	axialRotationForce* arf = new axialRotationForce(_name, this, loc, u, i, j);
	arf->setForceValue(v);
	forces[_name] = arf;
	return arf;
}

axialRotationForce* mbd_model::createAxialRotationForce(QTextStream& qts)
{
	QString ch, nm;
	QString i_name, j_name;
	VEC3D loc, u;
	double fv = 0.0;
	qts >> ch >> nm
		>> ch >> i_name
		>> ch >> j_name
		>> ch >> loc.x >> loc.y >> loc.z
		>> ch >> u.x >> u.y >> u.z
		>> ch >> fv;
	pointMass *fi = i_name == "ground" ? ground : masses[i_name];
	pointMass *fj = j_name == "ground" ? ground : masses[j_name];
	return createAxialRotationForce(nm, fi, fj, loc, u, fv);
	//arf->setForceValue(fv);
}

void mbd_model::set2D_Mode(bool b)
{
	is2D = b;
}

bool mbd_model::mode2D()
{
	return is2D;
}

rigidBody* mbd_model::createRigidBody(QString _name)
{
	rigidBody* rb = new rigidBody(_name);
	masses[_name] = rb;
// 	QMap<QString, object*>::iterator it = objs.find(_name);
// 	object* obj = NULL;
// 	mass *ms = new mass(this, _name);
// 	if (it == objs.end()){
// 		masses[NULL] = ms;
// 		return ms;
// 	}
// 	else{
// 		obj = it.value();
// 		masses[obj] = ms;
// 		obj->addPointMass(ms);
// 	}
// 
// 	switch (obj->objectType()){				//calculate Mass Center Position
// 	case CUBE:{
// 		cube *c = getChildObject<cube*>(_name);
// 		VEC3D minp = c->min_point();
// 		VEC3D maxp = c->max_point();
// 		VEC3D CuCenterp;
// 		CuCenterp.x = (maxp.x + minp.x) / 2;
// 		CuCenterp.y = (maxp.y + minp.y) / 2;
// 		CuCenterp.z = (maxp.z + minp.z) / 2;
// 		ms->setMassPoint(CuCenterp);
// 		break;
// 	}
// 	case PLANE:{
// 		plane *pl = getChildObject<plane*>(_name);// .find(_name).value();
// 		VEC3D xw = pl->XW();
// 		VEC3D w2 = pl->W2();
// 		VEC3D w3 = pl->W3();
// 		VEC3D w4 = pl->W4();
// 		VEC3D PlMidp1;
// 		PlMidp1.x = (xw.x + w2.x) / 2;
// 		PlMidp1.y = (xw.y + w2.y) / 2;
// 		PlMidp1.z = (xw.z + w2.z) / 2;
// 		VEC3D PlMidp2;
// 		PlMidp2.x = (w3.x + w4.x) / 2;
// 		PlMidp2.y = (w3.y + w4.y) / 2;
// 		PlMidp2.z = (w3.z + w4.z) / 2;
// 		VEC3D PlCenterp;
// 		PlCenterp.x = (PlMidp1.x + PlMidp2.x) / 2;
// 		PlCenterp.y = (PlMidp1.y + PlMidp2.y) / 2;
// 		PlCenterp.z = (PlMidp1.z + PlMidp2.z) / 2;
// 		ms->setMassPoint(PlCenterp);
// 		break;
// 	}
// 	case POLYGON:{
// 		polygonObject* pobj = getChildObject<polygonObject*>(_name);// <QString, polygonObject>::iterator po = pObjs.find(_name);
// 		ms->setMassPoint(pobj->getOrigin());
// 		break;
// 	}
// 	case CYLINDER:{
// 		cylinder *cy = getChildObject<cylinder*>(_name);// .find(_name).value();
// 		VEC3D goc;
// 		ms->setPosition(cy->origin());
// 		ms->setEP(cy->orientation());
// 		//goc = cy.origin().To<double>();
// 		break;
// 	}
// 	}
// 	db->addChild(database::MASS_ROOT, _name);
	return rb;
}

rigidBody* mbd_model::createRigidBody(QTextStream& qts)
{
	rigidBody* rb = new rigidBody();
	rb->openData(qts);
	masses[rb->name()] = rb;
	return rb;
}

pointMass* mbd_model::Ground()
{
	return ground;
}

// particle_system* mbd_model::makeParticleSystem(QString _name)
// {
// 	if (!ps)
// 		ps = new particle_system(_name, this);
// 	return ps;
// }

// drivingConstraint* mbd_model::makeDrivingConstraint(QString _name, kinematicConstraint* kconst, tDriving td, double val)
// {
// 	drivingConstraint* dc = new drivingConstraint(_name);
// 	dc->define(kconst, td, val);
// 	dconsts[_name] = dc;
// 	return dc;
// }

void mbd_model::exportPointMassResultData2TXT()
{
	foreach(pointMass* pm, masses)
	{
		pm->exportResultData2TXT();
	}
	
}

void mbd_model::exportReactionForceResultData2TXT()
{
	foreach(kinematicConstraint* kc, consts)
	{
		kc->exportResultData2TXT();
	}
}

void mbd_model::loadPointMassResultDataFromTXT()
{
// 	foreach(pointMass* pm, masses)
// 	{
// 		QString pname = model::path + model::name + "/" + pm->name() + ".txt";
// 		QFile qf(pname);
// 		QString ch;
// 		qf.open(QIODevice::ReadOnly);
// 		if (qf.isOpen())
// 		{
// 			QTextStream qts(&qf);
// 			for (int i = 0; i < 22; i++)
// 			{
// 				qts >> ch;
// 			}
// 			while (!qts.atEnd())
// 			{
// 				
// 			
// 				respointMassResultData pmrd;
// 				qts >> pmrd.time
// 					>> pmrd.pos.x >> pmrd.pos.y >> pmrd.pos.z
// 					>> pmrd.ep.e0 >> pmrd.ep.e1 >> pmrd.ep.e2 >> pmrd.ep.e3
// 					>> pmrd.vel.x >> pmrd.vel.y >> pmrd.vel.z
// 					>> pmrd.omega.x >> pmrd.omega.y >> pmrd.omega.z
// 					>> pmrd.acc.x >> pmrd.acc.y >> pmrd.acc.z
// 					>> pmrd.alpha.x >> pmrd.alpha.y >> pmrd.alpha.z;
// 				pm->appendResultData(pmrd);
// 			}			
// 		}
// 		qf.close();
// 	}
}

// collision* mbd_model::makeCollision(QString _name, double _rest, double _fric, double _rfric, double _coh, double _ratio, tCollisionPair tcp, tContactModel tcm, void* o1, void* o2)
// {
// 	collision* c = NULL;
// 	switch (tcp){
// 	case PARTICLES_PARTICLES:
// 		c = new collision_particles_particles(_name, this, ps, tcm);
// 		c->setContactParameter(
// 			ps->youngs(), ps->youngs(), ps->poisson(), ps->poisson(),
// 			ps->shear(), ps->shear(), _rest, _fric, _rfric, _coh, _ratio);
// 		break;
// 	case PARTICLES_PLANE:{
// 		plane *pl = (plane*)o2;
// 		c = new collision_particles_plane(_name, this, ps, pl, tcm);
// 		//	ps->addCollision(c);
// 		c->setContactParameter(
// 			ps->youngs(), pl->youngs(), ps->poisson(), pl->poisson(),
// 			ps->shear(), pl->shear(), _rest, _fric, _rfric, _coh, _ratio);
// 		break;
// 	}
// 	case PARTICLES_CYLINDER:{
// 		cylinder *cy = (cylinder*)o2;
// 		c = new collision_particles_cylinder(_name, this, ps, cy, tcm);
// 		//	ps->addCollision(c);
// 		c->setContactParameter(
// 			ps->youngs(), cy->youngs(), ps->poisson(), cy->poisson(),
// 			ps->shear(), cy->shear(), _rest, _fric, _rfric, _coh, _ratio);
// 		break;
// 	}
// 	case PARTICLES_POLYGONOBJECT:{
// 		polygonObject *po = (polygonObject*)o2;
// 		c = new collision_particles_polygonObject(_name, this, ps, po, tcm);
// 		//ps->addCollision(c);
// 		c->setContactParameter(
// 			ps->youngs(), po->youngs(), ps->poisson(), po->youngs(),
// 			ps->shear(), po->shear(), _rest, _fric, _rfric, _coh, _ratio);
// 		break;
// 	}
// 	}
// 
// 
// 	if (c)
// 		cs[_name] = c;
// 
// 	return c;
// }

void mbd_model::saveModel(QTextStream& qts)
{
	foreach(pointMass* pm, masses)
	{
		pm->saveData(qts);
	}

	foreach(kinematicConstraint* kc, consts)
	{
		if (kc->constType() == kinematicConstraint::CABLE)
			continue;
		kc->saveData(qts);
	}
	foreach(cableConstraint* cc, cables)
	{
		cc->saveCableConstraintData(qts);
	}
	foreach(forceElement* fe, forces)
	{
		fe->saveData(qts);
	}
// 	QString model_file = model_path + "/" + name + ".mde";
// 	QFile io_model(model_file);
// 	io_model.open(QIODevice::WriteOnly);
// 	QTextStream ts(&io_model);
// 	ts << "model_path " << model_path << endl;
// 	ts << "sim_type " << tsim << endl;
// 	ts << "name " << name << endl;
// 	ts << "unit " << (int)unit << endl;
// 	ts << "gravity_direction " << (int)dg << endl;
// 	ts << "gravity " << grav.x << " " << grav.y << " " << grav.z << endl;
// 	if (objs.size())
// 	{
// 		foreach(object* value, objs)
// 		{
// 			value->save_object_data(ts);
// 		}
// 	}
// 
// 	if (ps){
// 		if (ps->numParticle())
// 		{
// 			ps->saveParticleSystem(io_model);
// 		}
// 	}
// 
// 	if (cs.size())
// 	{
// 		foreach(collision* value, cs)
// 		{
// 			value->save_collision_data(ts);
// 		}
// 	}
// 	io_model.close();
// 	//comm("Model file was saved - " + model_file);
}

// void mbd_model::openModeler(GLWidget *gl, QString& file)
// {
// 	QFile pf(file);
// 	pf.open(QIODevice::ReadOnly);
// 	QTextStream in(&pf);
// 	QString ch;
// 	int ts;
// 	int _unit;
// 	int _dg;
// 	int isExist = false;
// 	unsigned int id;
// 	while (!in.atEnd())
// 	{
// 		in >> ch;
// 		if (ch == "model_path") in >> model_path;
// 		else if (ch == "name") in >> name;
// 		else if (ch == "unit") in >> _unit;
// 		else if (ch == "gravity_direction") in >> _dg;
// 		else if (ch == "sim_type") in >> ts;
// 		else if (ch == "gravity") in >> grav.x >> grav.y >> grav.z;
// 		//		break;
// 		else if (ch == "OBJECT")
// 		{
// 			in >> ch;
// 			if (ch == "CUBE") {
// 				int up;
// 				int tr, tmat;
// 				in >> id >> ch >> tr >> tmat >> up >> isExist;
// 				VEC3D min_p, max_p;
// 				in >> min_p.x >> min_p.y >> min_p.z;
// 				in >> max_p.x >> max_p.y >> max_p.z;
// 				cube *c = makeCube(ch, (tMaterial)tmat, (tRoll)tr);
// 				c->setUpdate((bool)up);
// 				c->setID(id);
// 				c->define(min_p, max_p);
// 				gl->makeCube(c);
// 				if (isExist){
// 					mass *m = makeMass(c->objectName());
// 					m->openData(in);
// 					m->setBaseGeometryType(CUBE);
// 					m->setGeometryObject(gl->getVObjectFromName(m->name()));
// 				}
// 			}
// 			else if (ch == "PLANE") {
// 				int up;
// 				int tr, tmat;
// 				in >> id >> ch >> tr >> tmat >> up >> isExist;
// 				VEC3D p1, p2, p3, p4;
// 				in >> p1.x >> p1.y >> p1.z;
// 				in >> p2.x >> p2.y >> p2.z;
// 				in >> p3.x >> p3.y >> p3.z;
// 				in >> p4.x >> p4.y >> p4.z;
// 				plane *p = makePlane(ch, (tMaterial)tmat, (tRoll)tr);
// 				p->setUpdate((bool)up);
// 				p->setID(id);
// 				p->define(p1, p2, p3, p4);
// 				gl->makePlane(p);
// 
// 				if (isExist){		//plane추가
// 					mass *m = makeMass(p->objectName());
// 					m->openData(in);
// 					m->setBaseGeometryType(PLANE);
// 					m->setGeometryObject(gl->getVObjectFromName(m->name()));
// 				}
// 			}
// 			else if (ch == "POLYGON"){
// 				int up;
// 				int tr, tmat;
// 				in >> id >> ch >> tr >> tmat >> up >> isExist;
// 				in >> ch;
// 				polygonObject* po = makePolygonObject(NO_FORMAT, ch);
// 				po->setUpdate(up);
// 				polygons.push_back(po);
// 				gl->makePolygonObject(po);
// 
// 				if (isExist){		//polygon 추가
// 					mass *m = makeMass(po->objectName());
// 					m->openData(in);
// 					m->setBaseGeometryType(POLYGON);
// 					m->setPolyGeometryObject(gl->getVPolyObjectFromName(m->name()));
// 				}
// 			}
// 			else if (ch == "CYLINDER"){
// 				double br, tpr, len;
// 				int tr, tmat, up;
// 				VEC3D org, bpos, tpos;
// 				in >> id >> ch >> tr >> tmat >> up >> isExist;
// 				in >> br >> tpr >> len;
// 				in >> org.x >> org.y >> org.z;
// 				in >> bpos.x >> bpos.y >> bpos.z;
// 				in >> tpos.x >> tpos.y >> tpos.z;
// 				cylinder *cy = makeCylinder(ch, (tMaterial)tmat, (tRoll)tr);
// 				cy->setUpdate((bool)up);
// 				cy->setID(id);
// 				cy->define(br, tpr, bpos, tpos);
// 				gl->makeCylinder(cy);
// 				if (isExist){
// 					mass *m = makeMass(cy->objectName());
// 					m->openData(in);
// 					m->setBaseGeometryType(CYLINDER);
// 					m->setGeometryObject(gl->getVObjectFromName(m->name()));
// 				}
// 			}
// 		}
// 		else if (ch == "PARTICLES")
// 		{
// 			unsigned int np;
// 			QString pfile, bo;
// 			double rho, E, pr, sh;// rest, fric, rfric, coh, sratio;
// 			in >> ch >> np;
// 			in >> ch >> pfile;
// 			in >> ch >> bo;
// 			in >> ch >> rho;
// 			in >> ch >> E;
// 			in >> ch >> pr;
// 			in >> ch >> sh;
// 			makeParticleSystem("particles");
// 			ps->setParticlesFromFile(pfile, bo, np, rho, E, pr, sh);
// 			//ps->setCollision(rest, fric, rfric, coh, sratio);
// 			gl->makeParticle(ps);
// 		}
// 		else if (ch == "STACK")
// 		{
// 			int tg;
// 			unsigned int nStack, npPerStack;
// 			double interval;
// 			in >> tg;
// 			in >> ch >> nStack;
// 			in >> ch >> interval;
// 			in >> ch >> npPerStack;
// 			ps->setGenerationMethod((tGenerationParticleMethod)tg, nStack, interval, npPerStack);
// 		}
// 		else if (ch == "CLUSTER")
// 		{
// 			int consist = 0;
// 			unsigned int nc = 0;
// 			in >> ch >> consist;
// 			ps->setParticleCluster(consist);
// 			//particle_cluster::setConsistNumber(consist);
// 			//ps->setParticleCluster(consist);
// 		}
// 		else if (ch == "COLLISION")
// 		{
// 			double rest, fric, rfric, coh, sratio;
// 			int tcm;
// 			QString nm, io, jo;
// 			in >> nm;
// 			in >> rest >> fric >> rfric >> coh >> sratio >> tcm;
// 			in >> ch >> io;
// 			in >> ch >> jo;
// 
// 			tCollisionPair cp = getCollisionPair(io != "particles" ? objs[io]->objectType() : ps->objectType(), jo != "particles" ? objs[jo]->objectType() : ps->objectType());
// 			if (io == "particles" && jo != "particles")
// 				makeCollision(nm, rest, fric, rfric, coh, sratio, cp, (tContactModel)tcm, ps, objs[jo]);
// 			else if (io != "particles" && jo == "particles")
// 				makeCollision(nm, rest, fric, rfric, coh, sratio, cp, (tContactModel)tcm, ps, objs[io]);
// 			else if (io == "particles" && jo == "particles")
// 				makeCollision(nm, rest, fric, rfric, coh, sratio, cp, (tContactModel)tcm, ps, ps);
// 			else
// 				makeCollision(nm, rest, fric, rfric, coh, sratio, cp, (tContactModel)tcm, objs[io], objs[jo]);
// 		}
// 	}
// 	pf.close();
// 	unit = (tUnit)_unit;
// 	dg = (tGravity)_dg;
// 	switch (dg)
// 	{
// 	case PLUS_X: grav = VEC3D(9.80665f, 0.f, 0.f); break;
// 	case PLUS_Y: grav = VEC3D(0.f, 9.80665f, 0.f); break;
// 	case PLUS_Z: grav = VEC3D(0.f, 0.f, 9.80665f); break;
// 	case MINUS_X: grav = VEC3D(-9.80665f, 0.f, 0.f); break;
// 	case MINUS_Y: grav = VEC3D(0.f, -9.80665f, 0.f); break;
// 	case MINUS_Z: grav = VEC3D(0.f, 0.f, -9.80665f); break;
// 	}
// }

void mbd_model::runExpression(double ct, double dt)
{
// 	if (objs.size()){
// 		foreach(object* value, objs)
// 		{
// 			if (value->expression())
// 				value->runExpression(ct, dt);
// 		}
// 	}
// 	if (dconsts.size()){
// 		foreach(drivingConstraint* value, dconsts)
// 		{
// 			value->driving(ct);
// 		}
// 	}
}

void mbd_model::updateCableInitLength()
{
	foreach(cableConstraint* cc, cables)
	{
		cc->updateCableInitLength();
	}
}

void mbd_model::updateAndInitializing()
{
// 	memset(rhs)
// 	foreach()
}

// unsigned int mbd_model::numPolygonSphere()
// {
// 	unsigned int _np = 0;
// 	if (npoly)
// 	{
// 		foreach(object* value, objs)
// 		{
// 			if (value->objectType() == POLYGON)
// 				_np += getChildObject<polygonObject*>(value->objectName())->numIndex();
// 		}
// 	}
// 	return _np;
// }

// void mbd_model::comm(QString com)
// {
// 	pte->appendPlainText(com);
// }

// void mbd_model::actionDelete(const QString& tg)
// {
// 	object* obj = objs.take(tg);
// 	if (obj)
// 		delete obj;
// }
