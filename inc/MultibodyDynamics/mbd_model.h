#ifndef MBD_MODEL_H
#define MBD_MODEL_H

#include <QMap>
#include <QFile>
#include <QTextStream>
#include <QString>
#include "kinematicConstraint.h"
#include "cableConstraint.h"
#include "gearConstraint.h"
#include "rigidBody.h"
#include "artificialCoordinate.h"
#include "model.h"
#include "algebraMath.h"

class forceElement;
class axialRotationForce;
class springDamperModel;

class mbd_model : public model
{
public:
	mbd_model();
	mbd_model(QString _name);
	virtual ~mbd_model();

	rigidBody* createRigidBody(QString _name);
	rigidBody* createRigidBody(QTextStream& qts);
	pointMass* Ground();
	kinematicConstraint* createKinematicConstraint(
		QString _name, kinematicConstraint::Type kt, 
		pointMass* i, VEC3D& spi, VEC3D& fi, VEC3D& gi, 
		pointMass* j, VEC3D& spj, VEC3D& fj, VEC3D& gj);
	kinematicConstraint* createKinematicConstraint(QTextStream& qts);
	kinematicConstraint* createCableConstraint(
		QString _name,
		pointMass* fi, VEC3D& fspi, pointMass* fj, VEC3D& fspj,
		pointMass* si, VEC3D& sspi, pointMass* sj, VEC3D& sspj);
	kinematicConstraint* createGearConstraint(
		QString _name, pointMass* i, kinematicConstraint* ik, 
		pointMass* j, kinematicConstraint* jk, double r);
	kinematicConstraint* createGearConstraint(
		QString _name, kinematicConstraint* ik, kinematicConstraint* jk,
		double cx, double cy, double cz, double r);
	kinematicConstraint* createGearConstraint(QTextStream& qts);
	kinematicConstraint* createCableConstraint(QTextStream& qts);
	springDamperModel* createSpringDamperElement(
		QString _name,
		pointMass* i, VEC3D& bLoc,
		pointMass* j, VEC3D& aLoc,
		double k, double c);
	springDamperModel* createSpringDamperElement(QTextStream& qts);
	axialRotationForce* createAxialRotationForce(
		QString _name, pointMass* i, pointMass* j, VEC3D loc, VEC3D u, double v);
	axialRotationForce* createAxialRotationForce(QTextStream& qts);
	artificialCoordinate* createArtificialCoordinate(QString _nm);

	void set2D_Mode(bool b);
	bool mode2D();
	QString& modelPath() { return model_path; }
	QString& modelName() { return name; }
	QMap<QString, kinematicConstraint*>& kinConstraint() { return consts; }
	QMap<QString, pointMass*>& pointMasses() { return masses; }
	QMap<QString, forceElement*>& forceElements() { return forces; }
	QMap<QString, cableConstraint*>& cableConstraints() { return cables; }
	QMap<QString, gearConstraint*>& gearConstraints() { return gears; }
	QMap<QString, artificialCoordinate*>& artificialCoordinates() { return acoordinates; }
	void exportPointMassResultData2TXT();
	void exportReactionForceResultData2TXT();
	void loadPointMassResultDataFromTXT();
	void saveModel(QTextStream& qts);
	void runExpression(double ct, double dt);
	void updateCableInitLength();
	void updateAndInitializing();
	virtual void userDefineProcess(){};

protected:
	pointMass *ground;
	VEC3D grav;

	QString model_path;
	QString name;
	bool is2D;

	QMap<QString, pointMass*> masses;
	QMap<QString, kinematicConstraint*> consts;
	QMap<QString, cableConstraint*> cables;
	QMap<QString, gearConstraint*> gears;
	QMap<QString, forceElement*> forces;
	QMap<QString, artificialCoordinate*> acoordinates;
};

#endif