#ifndef WAVEGENERATORMODEL_HPP
#define WAVEGENERATORMODEL_HPP

#include "sph_model.h"
#include "mbd_model.h"
#include "rectangle.h"
#include "simulation.h"

#include "gearConstraint.h"

class waveGeneratorModel : public mbd_model
{
public:
	waveGeneratorModel(sph_model* _sph) 
		: mbd_model("wave_generator") 
		, m_drum0(NULL)
		, m_drum1(NULL)
	{
		setUp();
	}
	virtual ~waveGeneratorModel(){}

	void setUp()
	{
		set2D_Mode(true);
		isFirstTime = false;
		m_drum0 = createRigidBody("drum0");
		m_drum0->setID(0);
		m_drum0->setMass(2);
		m_drum0->setPosition(VEC3D(4.5, 1.7, 0.0));
		m_drum0->setDiagonalInertia(31.57, 31.57, 31.57);
		m_drum0->setSymetryInertia(0, 0, 0);
		m_drum0->setInertia();

		m_drum1 = createRigidBody("drum1");
		m_drum1->setID(1);
		m_drum1->setMass(2);
		m_drum1->setPosition(VEC3D(4.7, 1.7, 0.0));
		m_drum1->setDiagonalInertia(31.57, 31.57, 31.57);
		m_drum1->setSymetryInertia(0, 0, 0);
		m_drum1->setInertia();

		m_generator0 = createRigidBody("generator0");
		m_generator0->setID(2);
		m_generator0->setMass(1);
		m_generator0->setPosition(VEC3D(5.0, 1.7, 0.0));
		m_generator0->setDiagonalInertia(2.0e-4, 2.0e-4, 2.0e-4);
		m_generator0->setSymetryInertia(0, 0, 0);
		m_generator0->setInertia();

		m_generator1 = createRigidBody("generator1");
		m_generator1->setID(3);
		m_generator1->setMass(1);
		m_generator1->setPosition(VEC3D(5.2, 1.7, 0.0));
		m_generator1->setDiagonalInertia(2.0e-4, 2.0e-4, 2.0e-4);
		m_generator1->setSymetryInertia(0, 0, 0);
		m_generator1->setInertia();

		kinematicConstraint *ground_drum0 = createKinematicConstraint(
			"ground_drum0", kinematicConstraint::REVOLUTE,
			ground, m_drum0->getPosition(), VEC3D(0, 1, 0), VEC3D(1, 0, 0),
			m_drum0, VEC3D(0, 0, 0), VEC3D(0, 1, 0), VEC3D(-1, 0, 0));
		ground_drum0->setLocation(m_drum0->getPosition());

		kinematicConstraint *ground_drum1 = createKinematicConstraint(
			"ground_drum1", kinematicConstraint::REVOLUTE,
			ground, m_drum1->getPosition(), VEC3D(0, 1, 0), VEC3D(1, 0, 0),
			m_drum1, VEC3D(0, 0, 0), VEC3D(0, 1, 0), VEC3D(-1, 0, 0));
		ground_drum1->setLocation(m_drum1->getPosition());

		kinematicConstraint *ground_generator0 = createKinematicConstraint(
			"ground_generator0", kinematicConstraint::REVOLUTE,
			ground, m_generator0->getPosition(), VEC3D(0, 1, 0), VEC3D(1, 0, 0),
			m_generator0, VEC3D(0, 0, 0), VEC3D(0, 1, 0), VEC3D(-1, 0, 0));
		ground_generator0->setLocation(m_generator0->getPosition());

		kinematicConstraint *ground_generator1 = createKinematicConstraint(
			"ground_generator1", kinematicConstraint::REVOLUTE,
			ground, m_generator1->getPosition(), VEC3D(0, 1, 0), VEC3D(1, 0, 0),
			m_generator1, VEC3D(0, 0, 0), VEC3D(0, 1, 0), VEC3D(-1, 0, 0));
		ground_generator1->setLocation(m_generator1->getPosition());

		gear_ratio = 20.0;

		damping = 5.0;

		//model::rs->defineUserResult("power0");
	//	model::rs->defineUserResult("power1");
		kinematicConstraint *gear0 = createGearConstraint(
			"gear0", m_drum0, ground_drum0, m_generator0, ground_generator0, 20.0);
// 		//gear0->ratio = 0.5;
		kinematicConstraint *gear1 = createGearConstraint(
			"gear1", m_drum1, ground_drum1, m_generator1, ground_generator1, 20.0);
		//gear0->setGroundPointer(ground);
		//m_drum0->setExternalMoment(VEC3D(0.0, 0.0, 1.0));
		//m_generator0->setExternalMoment(VEC3D(-0.002, 0.0, 0.0));
		//gear1->ratio = 
	}

	void applyGeneratorToque()
	{
// 		VEC3D T_drum0 = m_drum0->getExternalMoment();
// 		VEC3D T_drum1 = m_drum1->getExternalMoment();
	
		double td0 = m_drum0->getExternalMoment().z;
		double td1 = m_drum1->getExternalMoment().z;
// 		T_drum0 -= dFc0;
// 		T_drum1 -= dFc1;
	
		//double d = 0.5;
		double wg0 = m_generator0->getVelocity().z;
		double wg1 = m_generator1->getVelocity().z;
// 		VEC3D w_generator0 = 2.0 * m_generator0->getEP().L() * m_generator0->getEV();
// 		VEC3D w_generator1 = 2.0 * m_generator1->getEP().L() * m_generator1->getEV();
// 		if (!isFirstTime && simulation::ctime > 0.01)
// 		{
// 			damping = gear_ratio * T_drum0.length() / w_generator0.length();
// 			//dFc1 = T_drum1;
// 			isFirstTime = true;
// 			//return;
// 		}
// 		VEC3D T_generator0;
// 		VEC3D T_generator1;
		double tg0 = /*-gear_ratio * td0*/ - damping * wg0;
		double tg1 = /*-gear_ratio * td1*/ - damping * wg1;
//  		if (T_generator0.x < 0)
//  			T_generator0.x = 0;
//  		if (T_generator1.x < 0)
//  			T_generator1.x = 0;
		m_generator0->setExternalMoment(VEC3D(0, 0, tg0));
		m_generator1->setExternalMoment(VEC3D(0, 0, tg1));
		//isFirstTime = true;
	}

	void calculatePower(double ct)
	{
		double wg0 = m_generator0->getVelocity().z;
		double wg1 = m_generator1->getVelocity().z;
		//VEC3D w_generator0 = 2.0 * m_generator0->getEP().L() * m_generator0->getEV();
		double power0 = damping * pow(wg0, 2.0);
		//VEC3D w_generator1 = 2.0 * m_generator1->getEP().L() * m_generator1->getEV();
		double power1 = damping * pow(wg1, 2.0);
		time_double td0 = { ct, power0 };
		time_double td1 = { ct, power1 };
		model::rs->insertTimeDoubleResult(QString("power0"), td0);
		model::rs->insertTimeDoubleResult(QString("power1"), td1);
		/*powers1[ct] = power1;*/
	}

	void exportResult()
	{
// 		QFile qf(model::path + model::name + "/" + "power0.txt");
// 		qf.open(QIODevice::WriteOnly);
// 		QTextStream qts(&qf);
// 		
// 		foreach(time_double td, model::rs->timeDoubleResults()["power0"])
// 		{
// 			qts << td.time << " " << td.value << endl;
// 		}
// // 		QMapIterator<double, double> it0(powers0);
// // 		while (it0.hasNext())
// // 		{
// // 			it0.next();
// // 			qts << it0.key() << " " << it0.value() << endl;
// // 		}
// 		qf.close();
// 		qf.setFileName(model::path + model::name + "/" + "power1.txt");
// 		qf.open(QIODevice::WriteOnly);
// 		//qts(&qf);
// 		qts << "time " << "power " << endl;
// 		foreach(time_double td, model::rs->timeDoubleResults()["power1"])
// 		{
// 			qts << td.time << " " << td.value << endl;
// 		}
// // 		QMapIterator<double, double> it1(powers1);
// // 		while (it1.hasNext())
// // 		{
// // 			it1.next();
// // 			
// // 		}
// 		qf.close();
	}

	virtual void userDefineProcess()
	{
		applyGeneratorToque();
	}

	bool isFirstTime;
	double gear_ratio;
	double damping;

	VEC3D dFc0;
	VEC3D dFc1;

	rigidBody *m_drum0;
	rigidBody *m_drum1;

	rigidBody *m_generator0;
	rigidBody *m_generator1;

//	QMap<double, double> powers0;
//	QMap<double, double> powers1;
};

#endif