#include "resultStorage.h"
#include <QDebug>
#include <QFile>

resultStorage *RS;

resultStorage::resultStorage()
	: isAllocMemory(false)
	, type(NULL)
	, pos(NULL)
	, vel(NULL)
	, press(NULL)
	, isf(NULL)
	, color(NULL)
	, time(NULL)
	, cramp(NULL)
	, cpart(0)
	, np(0)
	, nfluid(0)
	, nfloat(0)
	, nbound(0)
	, ndummy(0)
	, maxPress(0)
	, pdata_type(3)
{
	RS = this;
}

resultStorage::~resultStorage()
{
	clearMemory();
}

resultStorage* resultStorage::RStorage()
{
	return RS;
}

unsigned int resultStorage::Np()
{
	return np;
}

unsigned int resultStorage::NFluid()
{
	return nfluid;
}

unsigned int resultStorage::NFloat()
{
	return nfloat;
}

unsigned int resultStorage::NBound()
{
	return nbound;
}

unsigned int resultStorage::NDummy()
{
	return ndummy;
}

ucolors::colorMap * resultStorage::ColorMap()
{
	return cramp;
}

// void resultStorage::defineUserResult(QString nm)
// {
// 	
// }

void resultStorage::insertTimeData(double ct)
{
	time[cpart] = ct; 
}

void resultStorage::insertGLineData(QString tg, VEC3D& b)
{
	glLine[tg].push_back(b);
}

void resultStorage::insertPartName(QString nm)
{
	pList.push_back(nm);
}

void resultStorage::insertPointMassResult(QString& nm, pointMassResultData& pmrd)
{
	pmrs[nm].push_back(pmrd);
}

void resultStorage::insertReactionForceResult(QString& nm, reactionForceData& rfd)
{
	rfrs[nm].push_back(rfd);
}

double* resultStorage::getPartPosition(unsigned int pt)
{
	return pos + (np * pdata_type * pt);
}

double* resultStorage::getPartVelocity(unsigned int pt)
{
	return vel + (np * 3 * pt);
}

double* resultStorage::getPartPressure(unsigned int pt)
{
	return press + (np * pt);
}

double* resultStorage::getPartColor(unsigned int pt)
{
	return color + (np * 4 * pt);
}

particle_type* resultStorage::getParticleType()
{
	return type;
}

bool* resultStorage::getPartFreeSurface(unsigned int pt)
{
	return isf + (np * pt);
}

void resultStorage::clearMemory()
{
	if (time) delete[] time; time = NULL;
	if (type) delete[] type; type = NULL;
	if (pos) delete[] pos; pos = NULL;
	if (vel) delete[] vel; vel = NULL;
	if (press) delete[] press; press = NULL;
	if (isf) delete[] isf; isf = NULL;
	if (color) delete[] color; color = NULL;
	if (cramp) delete cramp; cramp = NULL;
//	qDeleteAll(pmrd);
	isAllocMemory = false;
}

void resultStorage::allocMemory(unsigned int n, unsigned int npart, solver_type ts)
{
	int pn = ts == DEM ? 4 : 3;
	time = new double[npart]; memset(time, 0, sizeof(double)*npart);
	pos = new double[(n * pn) * npart]; memset(pos, 0, sizeof(double) * n * pn * npart);
	vel = new double[(n * 3) * npart]; memset(vel, 0, sizeof(double) * n * 3 * npart);	
	color = new double[(n * 4) * npart]; memset(color, 0, sizeof(double) * n * 4 * npart);
	cramp = new ucolors::colorMap(npart);
	pdata_type = pn;
	if (ts == SPH)
	{
		type = new particle_type[n]; memset(type, 0, sizeof(particle_type)*n);
		press = new double[n * npart]; memset(press, 0, sizeof(double) * n * npart);
		isf = new bool[n * npart]; memset(isf, 0, sizeof(bool) * n * npart);
	}
	else if (ts == DEM)
	{
		cramp->setTarget(ucolors::COLORMAP_VELOCITY_MAG);
	}
	
// 	QStringList pm = pmrd.keys();
// 	for (unsigned int i = 0; i < pmrd.size(); i++)
// 	{
// 		pmrd[pm.at(i)] = new pointMassResultData[npart];
// 	}
	isAllocMemory = true;
}

void resultStorage::setResultMemorySPH(unsigned int npart)
{
// 	nfluid = nf;
// 	nfloat = nft;
// 	nbound = nb;
// 	ndummy = nd;
// 	np = nf + nft + nb + nd;
// // 	if (np != _np)
// 		return;
// 	double *t_time = new double[1];				memset(t_time, 0, sizeof(double));
// 	tParticle *t_type = new tParticle[np];		memset(t_type, 0, sizeof(tParticle) * np);
// 	double *t_pos = new double[(np * 3)];		memset(t_pos, 0, sizeof(double) * np * 3);
// 	double *t_vel = new double[(np * 3)];		memset(t_vel, 0, sizeof(double) * np * 3);
// 	double *t_press = new double[np];			memset(t_press, 0, sizeof(double) * np);
// 	bool *t_isf = new bool[np];					memset(t_isf, 0, sizeof(bool) * np);
// 	double* t_color = new double[(np * 4)];		memset(t_color, 0, sizeof(double) * np * 4);
// 
// 	memcpy(t_time, time, sizeof(double));
// 	memcpy(t_type, type, sizeof(tParticle) * np);
// 	memcpy(t_pos, pos, sizeof(double) * np * 3);
// 	memcpy(t_vel, vel, sizeof(double) * np * 3);
// 	memcpy(t_press, press, sizeof(double) * np);
// 	memcpy(t_isf, isf, sizeof(bool) * np);
// 	memcpy(t_color, color, sizeof(double) * np * 4);
	total_parts = npart;
	clearMemory();
	if (!isAllocMemory)
		allocMemory(np, total_parts);
// 	if (!cramp)
// 	{
// 		delete cramp;
// 	}
	//cramp = new ucolors::colorMap(total_parts);
// 	memcpy(time, t_time, sizeof(double));
// 	memcpy(type, t_type, sizeof(tParticle) * np);
// 	memcpy(pos, t_pos, sizeof(double) * np * 3);
// 	memcpy(vel, t_vel, sizeof(double) * np * 3);
// 	memcpy(press, t_press, sizeof(double) * np);
// 	memcpy(isf, t_isf, sizeof(bool) * np);
// 	memcpy(color, t_color, sizeof(double) * np * 4);
// 
// 	delete[] t_time; t_time = NULL;
// 	delete[] t_type; t_type = NULL;
// 	delete[] t_pos; t_pos = NULL;
// 	delete[] t_vel; t_vel = NULL;
// 	delete[] t_press; t_press = NULL;
// 	delete[] t_isf; t_isf = NULL;
// 	delete[] t_color; t_color = NULL;
}

void resultStorage::setResultMemoryDEM(unsigned int npart, unsigned int _np)
{
	np = _np;
	total_parts = npart;
	clearMemory();
	if (!isAllocMemory)
		allocMemory(np, total_parts, DEM);
}

void resultStorage::insertTimeDoubleResult(QString& nm, time_double& td)
{
	tds[nm].push_back(td);// td;
}

void resultStorage::definePartDatasSPH(bool isOnlyCountNumPart, int index)
{
	unsigned int _np = np;
// 	if (!cramp)
// 	{
// 		if (!cpart)
// 			cramp = new ucolors::colorMap(1);
// 		else
// 			cramp = new ucolors::colorMap(cpart);
// 	}

	if (!isOnlyCountNumPart)
	{
		unsigned int id = index > -1 ? index : cpart;
		VEC3D min_vel = DBL_MAX;
		VEC3D max_vel = -DBL_MAX;
		double min_pre = DBL_MAX;
		double max_pre = -DBL_MAX;
		for (unsigned int j = 0; j < _np; j++)
		{
			unsigned int idx = np * 3 * id + j * 3;
			max_vel.x = abs(vel[idx + 0]) > max_vel.x ? abs(vel[idx + 0]) : max_vel.x;
			max_vel.y = vel[idx + 1] > max_vel.y ? vel[idx + 1] : max_vel.y;
			max_vel.z = vel[idx + 2] > max_vel.z ? vel[idx + 2] : max_vel.z;

			min_vel.x = 0.0;// abs(vel[idx + 0]) < min_vel.x ? abs(vel[idx + 0]) : min_vel.x;
			min_vel.y = vel[idx + 1] < min_vel.y ? vel[idx + 1] : min_vel.y;
			min_vel.z = vel[idx + 2] < min_vel.z ? vel[idx + 2] : min_vel.z;

			max_pre = press[np * id + j] > max_pre ? press[np * id + j] : max_pre;
			min_pre = press[np * id + j] < min_pre ? press[np * id + j] : min_pre;

			max_pre = press[np * id + j] > max_pre ? press[np * id + j] : max_pre;
			min_pre = press[np * id + j] < min_pre ? press[np * id + j] : min_pre;
		}
		cramp->setMinMax(id, min_vel.x, min_vel.y, min_vel.z, max_vel.x, max_vel.y, max_vel.z, min_pre, max_pre);

		unsigned int idx = np * id;

		double grad = maxPress * 0.1f;
		double t = 0.f;
		for (unsigned int j = 0; j < _np; j++){

			unsigned int cidx = (idx * 4) + (j * 4);
			unsigned int pidx = (idx * 3) + (j * 3);
			if (isf[idx + j]){
				color[cidx + 0] = 1.0f;
				color[cidx + 1] = 1.0f;
				color[cidx + 2] = 1.0f;
				color[cidx + 3] = 1.0f;
				continue;
			}
			if (type[j] == FLOATING)
			{
				color[cidx + 0] = 0.0f;
				color[cidx + 1] = 0.0f;
				color[cidx + 2] = 0.0f;
				color[cidx + 3] = 1.0f;
				continue;
			}
			if (type[j] == FLOATING_DUMMY)
			{
				color[cidx + 0] = 0.0f;
				color[cidx + 1] = 1.0f;
				color[cidx + 2] = 1.0f;
				color[cidx + 3] = 1.0f;
				continue;
			}
		

			if (cramp->target() == ucolors::COLORMAP_PRESSURE)
				cramp->getColorRamp(id, press[idx + j], &(color[cidx]));
			else if (cramp->target() == ucolors::COLORMAP_VELOCITY_X)
				cramp->getColorRamp(id, abs(vel[pidx + 0]), &(color[cidx]));
			color[cidx + 3] = 1.0f;
		}
		if (index == -1){
			cpart++;
		}
	}
}

void resultStorage::definePartDatasDEM(bool isOnlyCountNumPart, int index /*= -1*/)
{
	if (!isOnlyCountNumPart)
	{
		unsigned int id = index > -1 ? index : cpart;
		VEC3D min_vel = DBL_MAX;
		VEC3D max_vel = -DBL_MAX;
		double min_pre = DBL_MAX;
		double max_pre = -DBL_MAX;
		for (unsigned int j = 0; j < np; j++)
		{
			unsigned int idx = np * 3 * id + j * 3;
			max_vel.x = abs(vel[idx + 0]) > max_vel.x ? abs(vel[idx + 0]) : max_vel.x;
			max_vel.y = vel[idx + 1] > max_vel.y ? vel[idx + 1] : max_vel.y;
			max_vel.z = vel[idx + 2] > max_vel.z ? vel[idx + 2] : max_vel.z;

			min_vel.x = 0.0;// abs(vel[idx + 0]) < min_vel.x ? abs(vel[idx + 0]) : min_vel.x;
			min_vel.y = vel[idx + 1] < min_vel.y ? vel[idx + 1] : min_vel.y;
			min_vel.z = vel[idx + 2] < min_vel.z ? vel[idx + 2] : min_vel.z;

// 			max_pre = press[np * id + j] > max_pre ? press[np * id + j] : max_pre;
// 			min_pre = press[np * id + j] < min_pre ? press[np * id + j] : min_pre;
// 
// 			max_pre = press[np * id + j] > max_pre ? press[np * id + j] : max_pre;
// 			min_pre = press[np * id + j] < min_pre ? press[np * id + j] : min_pre;
		}
		cramp->setMinMax(id, min_vel.x, min_vel.y, min_vel.z, max_vel.x, max_vel.y, max_vel.z, 0, 0);

		unsigned int idx = np * id;

		double grad = maxPress * 0.1f;
		double t = 0.f;
		for (unsigned int j = 0; j < np; j++){

			unsigned int cidx = (idx * 4) + (j * 4);
			unsigned int pidx = (idx * 3) + (j * 3);
			VEC3D v(vel[pidx + 0], vel[pidx + 1], vel[pidx + 2]);
			if (cramp->target() == ucolors::COLORMAP_PRESSURE)
				cramp->getColorRamp(id, press[idx + j], &(color[cidx]));
			else if (cramp->target() == ucolors::COLORMAP_VELOCITY_X)
				cramp->getColorRamp(id, abs(v.x), &(color[cidx]));
			else if (cramp->target() == ucolors::COLORMAP_VELOCITY_Y)
				cramp->getColorRamp(id, abs(v.y), &(color[cidx]));
			else if (cramp->target() == ucolors::COLORMAP_VELOCITY_Z)
				cramp->getColorRamp(id, abs(v.z), &(color[cidx]));
			else if (cramp->target() == ucolors::COLORMAP_VELOCITY_MAG)
				cramp->getColorRamp(id, abs(v.length()), &(color[cidx]));
			color[cidx + 3] = 1.0f;
		}
		if (index == -1){
			cpart++;
		}
	}
}

void resultStorage::setInformationSPH(unsigned int _np, unsigned int _nf, unsigned int _nft, unsigned int _nb, unsigned int _nd)
{
	np = _np;
	nfluid = _nf;
	nfloat = _nft;
	nbound = _nb;
	ndummy = _nd;
}

double resultStorage::getPartTime(unsigned int pt)
{
	return time[pt];
}

void resultStorage::openResultFromFile(unsigned int idx)
{
	if (!isAllocMemory && np)
		allocMemory(np);
	QFile pf(rList.at(idx));
	pf.open(QIODevice::ReadOnly);
	double time = 0.f;
	pf.read((char*)&time, sizeof(double));
	pf.read((char*)type, sizeof(particle_type) * np);
	pf.read((char*)pos, sizeof(double) * np * 3);
	pf.read((char*)vel, sizeof(double) * np * 3);
	pf.read((char*)press, sizeof(double) * np);
	pf.read((char*)isf, sizeof(bool) * np);

	double grad = maxPress * 0.1f;
	double t = 0.f;
	for (unsigned int i = 0; i < nfluid; i++){
		// 		if (type[i] != FLUID)
		// 			continue;
		if (isf[i]){
			color[i * 4 + 0] = 1.0f;
			color[i * 4 + 1] = 1.0f;
			color[i * 4 + 2] = 1.0f;
			color[i * 4 + 3] = 1.0f;
			continue;
		}

		switch (type[i]){
		case FLUID:
			//t = (press[i] - 0) / grad;
			if (cramp->target() == ucolors::COLORMAP_PRESSURE)
				cramp->getColorRamp(0, press[i], &(color[i * 4]));
			else if (cramp->target() == ucolors::COLORMAP_VELOCITY_X)
				cramp->getColorRamp(0, vel[i * 3 + 0], &(color[i * 4]));
			color[i * 4 + 3] = 1.0f;
			break;
		}
	}
	pf.close();
// 	if (!isDefined)
// 	{
// 		define();
// 	}
}

void resultStorage::openSphResultFiles(QStringList& slist)
{
	rList = slist;
	clearMemory();
	if (!isAllocMemory)
		allocMemory(np, slist.size());
	unsigned int i = 0;
	foreach(QString str, slist)
	{
		int begin = str.lastIndexOf("/");
		int end = str.lastIndexOf(".");
		QString partName = str.mid(begin+1, end - begin - 1);
		insertPartName(partName);
		QString file = str;
		QFile pf(file);
		pf.open(QIODevice::ReadOnly);
		pf.read((char*)&(time[i]), sizeof(double));
		pf.read((char*)&np, sizeof(unsigned int));
		pf.read((char*)&nfluid, sizeof(unsigned int));
		pf.read((char*)&nfloat, sizeof(unsigned int));
		pf.read((char*)&nbound, sizeof(unsigned int));
		pf.read((char*)&ndummy, sizeof(unsigned int));
		pf.read((char*)type, sizeof(particle_type)*np);
		pf.read((char*)&(pos[(np * 3) * i]), sizeof(double) * np * 3);
		pf.read((char*)&(vel[(np * 3) * i]), sizeof(double) * np * 3);
		pf.read((char*)&(press[np * i]), sizeof(double) * np);
		pf.read((char*)&(isf[np * i]), sizeof(bool) * np);
		pf.close();
		definePartDatasSPH(false, i);
		i++;
		qDebug() << i;
	}
	//np = nfluid + nfloat + nbound + ndummy;

// 	if (!isDefined)
// 	{
// 		define();
// 	}
// 	setMinMaxResult();
}

void resultStorage::insertDataSPH(particle_type* tp, double* _p, double* _v, double* _prs, bool isCalcContour /*= false*/)
{

}

void resultStorage::exportEachResult2TXT(QString path)
{
	QMapIterator<QString, QList<pointMassResultData>> m_pmrd(pmrs);
	while (m_pmrd.hasNext())
	{
		m_pmrd.next();
		QString file_name = path + "/" + m_pmrd.key() + ".txt";
		QFile qf(file_name);
		qf.open(QIODevice::WriteOnly);
		QTextStream qts(&qf);
		qts << "time "
			<< "px " << "py " << "pz " << "ep0 " << "ep1 " << "ep2 " << "ep3 "
			<< "vx " << "vy " << "vz " << "wx " << "wy " << "wz "
			<< "ax " << "ay " << "az " << "apx " << "apy " << "apz " << endl;
		foreach(pointMassResultData p, m_pmrd.value())
		{
			qts << p.time
				<< " " << p.pos.x << " " << p.pos.y << " " << p.pos.z
				<< " " << p.ep.e0 << " " << p.ep.e1 << " " << p.ep.e2 << " " << p.ep.e3
				<< " " << p.vel.x << " " << p.vel.y << " " << p.vel.z
				<< " " << p.omega.x << " " << p.omega.y << " " << p.omega.z
				<< " " << p.acc.x << " " << p.acc.y << " " << p.acc.z
				<< " " << p.alpha.x << " " << p.alpha.y << " " << p.alpha.z 
				<< " " << p.ea.e0 << " " << p.ea.e1 << " " << p.ea.e2 << " " << p.ea.e3 << endl;
		}
		qf.close();
	}
	QMapIterator<QString, QList<reactionForceData>> m_rfd(rfrs);
	while (m_rfd.hasNext())
	{
		m_rfd.next();
		QString file_name = path + "/" + m_rfd.key() + ".txt";
		QFile qf(file_name);
		qf.open(QIODevice::WriteOnly);
		QTextStream qts(&qf);
		qts << "time "
			<< "fix " << "fiy " << "fiz " << "ri0 " << "ri1 " << "ri2 " << "ri3 "
			<< "fjx " << "fjy " << "fjz " << "rj0 " << "rj1 " << "rj2 " << "rj3 " << endl;
		foreach(reactionForceData p, m_rfd.value())
		{
			qts << p.time
				<< " " << p.iAForce.x << " " << p.iAForce.y << " " << p.iAForce.z
				<< " " << p.iRForce.x << " " << p.iRForce.y << " " << p.iRForce.z << " " << p.iRForce.w
				<< " " << p.jAForce.x << " " << p.jAForce.y << " " << p.jAForce.z
				<< " " << p.jRForce.x << " " << p.jRForce.y << " " << p.jRForce.z << " " << p.jRForce.w << endl;
		}
		qf.close();
	}
	QMapIterator<QString, QList<time_double>> m_tds(tds);
	while (m_tds.hasNext())
	{
		m_tds.next();
		QString file_name = path + "/" + m_tds.key() + ".txt";
		QFile qf(file_name);
		qf.open(QIODevice::WriteOnly);
		QTextStream qts(&qf);
		qts << "time " << "power " << endl;
		foreach(time_double td, m_tds.value())
		{
			qts << td.time << " " << td.value << endl;
		}
		qf.close();
	}
}

QMap<QString, QList<VEC3D>>& resultStorage::linePointResults()
{
	return glLine;
}

// QMap<QString, pointMassResultData*>& resultStorage::PointMassResultData()
// {
// 	return pmrd;
// }
// 
// void resultStorage::preparePointMassResultDataMemory(QString nm)
// {
// 	pmrd[nm] = NULL;
// }
// 
// void resultStorage::exportPointMassResultData2TXT()
// {
// // 	QString file_name = model::path + model::name + "/" + name + ".txt";
// // 	QFile qf(file_name);
// // 	qf.open(QIODevice::WriteOnly);
// // 	QTextStream qts(&qf);
// }

QStringList& resultStorage::partList()
{
	return pList;
}
