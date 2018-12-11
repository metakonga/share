#ifndef HARD_MOVING_H
#define HARD_MOVING_H

#include <QString>
#include <QTextStream>
#include "algebraMath.h"

class pointMass;

class hardMoving
{
public:
	enum Type{ POSITION_MOVING = 0, CONSTANT_VELOCITY_MOVING };
	hardMoving(QString nm);
	~hardMoving();

	bool Define(Type tp, pointMass* tg, VEC3D dir, double et, double v);
	void RunMovement(double ct);
	void Save(QTextStream& qts);
	static bool IsOver() { return isOver; }
	static void SetOver(bool b) { isOver = b; }

private:
	static bool isOver;
	Type type;
	QString name;
	pointMass* target;
	VEC3D init_p;
	VEC3D direction;
	double value;
	double end_time;
};

#endif