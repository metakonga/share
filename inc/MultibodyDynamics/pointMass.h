#ifndef MASSPOINT_H
#define MASSPOINT_H

//#include "vectorTypes.h"
#include "algebraMath.h"
#include <QTextStream>

class pointMass
{
public:
	enum Type { GROUND = 0, RIGID, FLEXIBLE };
	pointMass();
	pointMass(QString& _name, Type _type);
	~pointMass();

	int ID() { return id; }
	Type MassType() { return type; }
	QString& name(){ return nm; }
	double getMass() { return ms; }
	MAT33D getInertia() { return inertia; }
	VEC3D getPosition() { return pos; }
	EPD getEP() { return ep; }
	EPD getEV() { return ev; }
	EPD getEA() { return ea; }
	VEC3D getExternalForce() { return ef; }
	VEC3D getCollisionForce() { return cf; }
	VEC3D getExternalMoment() { return em; }
	VEC3D getCollisionMoment() { return cm; }
	VEC3D getHydroForce() { return hf; }
	VEC3D getHydroMoment() { return hm; }
	VEC3D getVelocity() { return vel; }
	VEC3D getAcceleration() { return acc; }
	VEC3D getSymInertia() { return sym_iner; }
	VEC3D getPriInertia() { return prin_iner; }

	void setID(int _id) { id = _id; }
	void setMassType(Type tp) { type = tp; }
	//void setName(QString _n)
	void setMass(double _ms) { ms = _ms; }
	//void setMassPoint(VEC3D _mp) { pos = _mp; }
	void setSymetryInertia(double ixx, double iyy, double izz);
	void setDiagonalInertia(double ixy, double ixz, double iyz);
	void setInertia() { inertia.diagonal(POINTER3(prin_iner)); }
	void setPosition(VEC3D& _p) { pos = _p; }
	void setVelocity(VEC3D& _v) { vel = _v; }
	void setAcceleration(VEC3D& a) { acc = a; }
	void setHydroForce(VEC3D& _hf) { hf = _hf; }
	void setHydroMoment(VEC3D& _hm) { hm = _hm; }
	void setExternalForce(VEC3D& _f) { ef = _f; }
	void setExternalMoment(VEC3D& _m) { em = _m; }
	void addExternalForce(VEC3D& _f) { ef += _f; }
	void addExternalMoment(VEC3D& _m) { em += _m; }
	void setCollisionForce(VEC3D& _f) { cf = _f; }
	void setCollisionMoment(VEC3D& _m) { cm = _m; }
	void addCollisionForce(VEC3D& _f) { cf += _f; }
	void addCollisionMoment(VEC3D& _m) { cm += _m; }
	void setEP(EPD& _ep) { ep = _ep; }
	void setEV(EPD& _ev) { ev = _ev; }
	void setEA(EPD& _ea) { ea = _ea; }
	void saveData(QTextStream& ts) const;
	void openData(QTextStream& ts);
	void makeTransformationMatrix();
	void makeTransformationMatrix2D();
	MAT33D getTransformationMatrix();
	void exportResultData2TXT();

	VEC3D toLocal(VEC3D &v);
	VEC3D toLocal2D(VEC3D &v);
	VEC3D toGlobal(VEC3D &v);
	VEC3D toGlobal2D(VEC3D &v);
	VEC3D toDGlobal(VEC3D &v);

	double rotationAngle();

	static int Count() { return count; }

protected:
	static int count;
	Type type;
	QString nm;
	int id;
	double ms;				// pointMass of object
	VEC3D sym_iner;			// Ixx, Iyy, Izz
	VEC3D prin_iner;		// Ixy, Ixz, Iyz
	MAT33D inertia;

	double ang;
	VEC3D pos;
	VEC3D vel;
	VEC3D acc;
	VEC3D ef;
	VEC3D em;
	VEC3D cf;
	VEC3D cm;
	VEC3D hf;
	VEC3D hm;

	EPD ep;					// euler parameter
	EPD ev;					// angular velocity of ep
	EPD ea;					// angular acceleration of ep

	MAT33D A;				// transformation matrix	
};

#endif