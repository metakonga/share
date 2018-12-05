#include "hardMoving.h"
#include "pointMass.h"

bool hardMoving::isOver = false;

hardMoving::hardMoving(QString nm)
	: name(nm)
	, target(NULL)
	, value(0.0)
	, end_time(0.0)
{

}

hardMoving::~hardMoving()
{

}

void hardMoving::RunMovement(double ct)
{
	if (ct <= end_time)
	{
		VEC3D new_p;
		switch (type)
		{
		case CONSTANT_VELOCITY_MOVING:
			new_p = init_p + ct * value * direction;
			target->setPosition(new_p);
			//target->setVelocity(value * direction);
			break;
		}
		isOver = false;
	}
	else
	{
		isOver = true;
	}	
}

void hardMoving::Save(QTextStream& qts)
{
	qts << "ELEMENT hard_moving" << endl
		<< "NAME " << name << endl
		<< "TARGET " << target->Name() << endl
		<< "TYPE " << (int)type << endl
		<< "END_TIME " << end_time << endl
		<< "DIRECTION " << direction.x << " " << direction.y << " " << direction.z << endl
		<< "PARAMETERS " << value << endl;
}

bool hardMoving::Define(hardMoving::Type tp, pointMass* tg, VEC3D dir, double et, double v)
{
	init_p = tg->Position();
	type = tp;
	target = tg;
	direction = dir;
	value = v;
	end_time = et;
	return true;
}

