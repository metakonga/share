#include "rigidBody.h"

rigidBody::rigidBody(QString _name)
	: pointMass(_name, RIGID)
{

}

rigidBody::rigidBody()
	: pointMass(QString(""), RIGID)
{
	
}

rigidBody::~rigidBody()
{

}