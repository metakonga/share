#ifndef MPHYSICS_NUMERIC_H
#define MPHYSICS_NUMERIC_H

#include "algebra/matrix3x3.hpp"
#include "Algebra/matrix4x3.hpp"
#include "algebra/matrix4x4.hpp"
#include "algebra/vector3.hpp"
#include "algebra/vector4.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/sparse_matrix.hpp"
#include "algebra/euler_parameter.hpp"

using namespace algebra;

// #define min(a,b) a < b ? a : b
// #define max(a,b) a > b ? a : b

#define sign(a) a <= 0 ? (a == 0 ? 0 : -1) : 1

#define POINTER(a) &a(0)
#define POINTER3(a) &(a.x)
#define POINTER4(a) &(a.x)

inline float frand() { return rand() / (float)RAND_MAX; }
inline double drand() { return rand() / (double)RAND_MAX; }

typedef vector3<float>			VEC3F;
typedef vector3<double>			VEC3D;
typedef vector3<int>			VEC3I;
typedef vector3<unsigned int>	VEC3UI;
typedef vector4<float>			VEC4F;
typedef vector4<double>			VEC4D;
typedef vector4<int>			VEC4I;
typedef vector4<unsigned int>	VEC4UI;
typedef algebra::vector<int>	VECI;
typedef algebra::vector<unsigned int> VECUI;
typedef algebra::vector<float>	VECF;
typedef algebra::vector<double>	VECD;
typedef matrix<float>			MATF;
typedef matrix<double>			MATD;
typedef matrix3x3<float>		MAT33F;
typedef matrix3x4<float>		MAT34F;
typedef matrix4x4<float>		MAT44F;
typedef matrix3x3<double>		MAT33D;
typedef matrix3x4<double>		MAT34D;
typedef matrix4x3<double>		MAT43D;
typedef matrix4x4<double>		MAT44D;
typedef euler_parameter<float>	EPF;
typedef euler_parameter<double>	EPD;
typedef xdyn::sparse_matrix<double>   SMATD;
typedef VECF*					VECF_PTR;
typedef VECD*					VECD_PTR;
typedef VECI*					VECI_PTR;
typedef VECUI*					VECUI_PTR;
typedef VEC3F*					VEC3F_PTR;
typedef VEC3D*					VEC3D_PTR;
typedef VEC4F*					VEC4F_PTR;
typedef VEC4D*					VEC4D_PTR;
typedef MAT33F*					MAT33F_PTR;
typedef MAT34F*					MAT34F_PTR;
typedef MAT44F*					MAT44F_PTR;
typedef MAT33D*					MAT33D_PTR;
typedef MAT34D*					MAT34D_PTR;
typedef MAT44D*					MAT44D_PTR;
typedef MATF*					MATF_PTR;
typedef MATD*					MATD_PTR;
typedef EPF*					EPF_PTR;
typedef EPD*					EPD_PTR;
typedef SMATD*					SMATD_PTR;

typedef vector4<float>			VEC4F;
typedef VEC4F*					VEC4F_PTR;

typedef struct
{
	VEC3D P;
	VEC3D Q;
	VEC3D R;
	VEC3D V;
	VEC3D W;
	VEC3D N;
}host_polygon_info;

typedef struct 
{
	VEC3D origin;
	VEC3D vel;
	VEC3D omega;
	EPD ep;
}host_polygon_mass_info;

inline
MAT44D transpose(const MAT34D& m4x3, MAT34D& m3x4)
{
	return MAT44D(
		m4x3.a00 * m3x4.a00 + m4x3.a10 * m3x4.a10 + m4x3.a20 * m3x4.a20, m4x3.a00 * m3x4.a01 + m4x3.a10 * m3x4.a11 + m4x3.a20 * m3x4.a21, m4x3.a00 * m3x4.a02 + m4x3.a10 * m3x4.a12 + m4x3.a20 * m3x4.a22, m4x3.a00 * m3x4.a03 + m4x3.a10 * m3x4.a13 + m4x3.a20 * m3x4.a23,
		m4x3.a01 * m3x4.a00 + m4x3.a11 * m3x4.a10 + m4x3.a21 * m3x4.a20, m4x3.a01 * m3x4.a01 + m4x3.a11 * m3x4.a11 + m4x3.a21 * m3x4.a21, m4x3.a01 * m3x4.a02 + m4x3.a11 * m3x4.a12 + m4x3.a21 * m3x4.a22, m4x3.a01 * m3x4.a03 + m4x3.a11 * m3x4.a13 + m4x3.a21 * m3x4.a23,
		m4x3.a02 * m3x4.a00 + m4x3.a12 * m3x4.a10 + m4x3.a22 * m3x4.a20, m4x3.a02 * m3x4.a01 + m4x3.a12 * m3x4.a11 + m4x3.a22 * m3x4.a21, m4x3.a02 * m3x4.a02 + m4x3.a12 * m3x4.a12 + m4x3.a22 * m3x4.a22, m4x3.a02 * m3x4.a03 + m4x3.a12 * m3x4.a13 + m4x3.a22 * m3x4.a23,
		m4x3.a03 * m3x4.a00 + m4x3.a13 * m3x4.a10 + m4x3.a23 * m3x4.a20, m4x3.a03 * m3x4.a01 + m4x3.a13 * m3x4.a11 + m4x3.a23 * m3x4.a21, m4x3.a03 * m3x4.a02 + m4x3.a13 * m3x4.a12 + m4x3.a23 * m3x4.a22, m4x3.a03 * m3x4.a03 + m4x3.a13 * m3x4.a13 + m4x3.a23 * m3x4.a23
		);
}

inline 
MAT43D transpose(MAT34D& m4x3, MAT33D& m3x3)
{
	return MAT43D(
		m4x3.a00 * m3x3.a00 + m4x3.a10 * m3x3.a10 + m4x3.a20 * m3x3.a20, m4x3.a00 * m3x3.a01 + m4x3.a10 * m3x3.a11 + m4x3.a20 * m3x3.a21, m4x3.a00 * m3x3.a02 + m4x3.a10 * m3x3.a12 + m4x3.a20 * m3x3.a22,
		m4x3.a01 * m3x3.a00 + m4x3.a11 * m3x3.a10 + m4x3.a21 * m3x3.a20, m4x3.a01 * m3x3.a01 + m4x3.a11 * m3x3.a11 + m4x3.a21 * m3x3.a21, m4x3.a01 * m3x3.a02 + m4x3.a11 * m3x3.a12 + m4x3.a21 * m3x3.a22,
		m4x3.a02 * m3x3.a00 + m4x3.a12 * m3x3.a10 + m4x3.a22 * m3x3.a20, m4x3.a02 * m3x3.a01 + m4x3.a12 * m3x3.a11 + m4x3.a22 * m3x3.a21, m4x3.a02 * m3x3.a02 + m4x3.a12 * m3x3.a12 + m4x3.a22 * m3x3.a22,
		m4x3.a03 * m3x3.a00 + m4x3.a13 * m3x3.a10 + m4x3.a23 * m3x3.a20, m4x3.a03 * m3x3.a01 + m4x3.a13 * m3x3.a11 + m4x3.a23 * m3x3.a21, m4x3.a03 * m3x3.a02 + m4x3.a13 * m3x3.a12 + m4x3.a23 * m3x3.a22
		);
}

// inline 
// VEC4D transpose(const MAT34D& m4x3, VEC3D& v3)
// {
// 
// }

inline
VEC4D transpose(VEC3D& v3, MAT34D& m3x4)
{
	return VEC4D(
		v3.x * m3x4.a00 + v3.y * m3x4.a10 + v3.z * m3x4.a20,
		v3.x * m3x4.a01 + v3.y * m3x4.a11 + v3.z * m3x4.a21,
		v3.x * m3x4.a02 + v3.y * m3x4.a12 + v3.z * m3x4.a22,
		v3.x * m3x4.a03 + v3.y * m3x4.a13 + v3.z * m3x4.a23
		);
}


inline
VEC3F transpose(MAT33F& m, VEC3F& a)
{
	return VEC3F(
		a.x*m.a00 + a.y*m.a10 + a.z*m.a20,
		a.x*m.a01 + a.y*m.a11 + a.z*m.a21,
		a.x*m.a02 + a.y*m.a12 + a.z*m.a22);
}

inline
VEC3D transpose(VEC3D& a, MAT33D& m)
{
	return VEC3D(
		a.x*m.a00 + a.y*m.a10 + a.z*m.a20,
		a.x*m.a01 + a.y*m.a11 + a.z*m.a21,
		a.x*m.a02 + a.y*m.a12 + a.z*m.a22);
}

inline 
MAT33D transpose(VEC3D& a, VEC3D& b)
{
	return MAT33D(
		a.x * b.x, a.x * b.y, a.x * b.z,
		a.y * b.x, a.y * b.y, a.y * b.z,
		a.z * b.x, a.z * b.y, a.z * b.z);
}

inline MAT34D operator*(const MAT33D& a, const MAT34D& b)
{
	return MAT34D(
		a.a00*b.a00 + a.a01*b.a10 + a.a02*b.a20, a.a00*b.a01 + a.a01*b.a11 + a.a02*b.a21, a.a00*b.a02 + a.a01*b.a12 + a.a02*b.a22, a.a00*b.a03 + a.a01*b.a13 + a.a02*b.a23,
		a.a10*b.a00 + a.a11*b.a10 + a.a12*b.a20, a.a10*b.a01 + a.a11*b.a11 + a.a12*b.a21, a.a10*b.a02 + a.a11*b.a12 + a.a12*b.a22, a.a10*b.a03 + a.a11*b.a13 + a.a12*b.a23,
		a.a20*b.a00 + a.a21*b.a10 + a.a22*b.a20, a.a20*b.a01 + a.a21*b.a11 + a.a22*b.a21, a.a20*b.a02 + a.a21*b.a12 + a.a22*b.a22, a.a20*b.a03 + a.a21*b.a13 + a.a22*b.a23
		);
}

inline
MAT34D B(EPD& e, VEC3D& s)
{
	return MAT34D(
		2 * (2 * s.x*e.e0 + e.e2*s.z - e.e3*s.y), 2 * (2 * s.x*e.e1 + e.e3*s.z + e.e2*s.y), 2 * (e.e1*s.y + e.e0*s.z), 2 * (e.e1*s.z - e.e0*s.y),
		2 * (2 * s.y*e.e0 - e.e1*s.z + e.e3*s.x), 2 * (s.y*e.e1 - e.e0*s.z), 2 * (2 * s.y*e.e2 + e.e3*s.z + e.e1*s.x), 2 * (e.e2*s.z + e.e0*s.x),
		2 * (2 * s.z*e.e0 - e.e2*s.x + e.e1*s.y), 2 * (s.z*e.e1 + e.e0*s.y), 2 * (e.e3*s.y - e.e0*s.x), 2 * (2 * s.z*e.e3 + e.e2*s.y + e.e1*s.x)
		);
}

inline
MAT44D D(VEC3D& s, VEC3D& d)
{
	return matrix4x4<double>(
		4 * (s.x*d.x + s.y*d.y + s.z*d.z), 2 * (s.y*d.z - s.z*d.y), 2 * (s.z*d.x - s.x*d.z), 2 * (s.x*d.y - s.y*d.x),
		2 * (s.y*d.z - s.z*d.y), 4 * s.x*d.x, 2 * (s.x*d.y + s.y*d.x), 2 * (s.x*d.z + s.z*d.x),
		2 * (s.z*d.x - s.x*s.z), 2 * (s.x*d.y + s.y*d.x), 4 * s.y*d.y, 2 * (s.y*d.z + s.z*d.y),
		2 * (s.x*d.y - s.y*d.x), 2 * (s.x*d.z + s.z*d.x), 2 * (s.y*d.z + s.z*d.y), 4 * s.z*d.z);
}

inline VEC3D operator*(const MAT34D& m3x4, EPD& v4)
{
	return VEC3D(
		m3x4.a00*v4.e0 + m3x4.a01*v4.e1 + m3x4.a02*v4.e2 + m3x4.a03*v4.e3,
		m3x4.a10*v4.e0 + m3x4.a11*v4.e1 + m3x4.a12*v4.e2 + m3x4.a13*v4.e3,
		m3x4.a20*v4.e0 + m3x4.a21*v4.e1 + m3x4.a22*v4.e2 + m3x4.a23*v4.e3
		);
}

inline
MAT44D operator*(MAT43D& m4x3, MAT34D& m3x4)
{
	return MAT44D(
		m4x3.a00 * m3x4.a00 + m4x3.a01 * m3x4.a10 + m4x3.a02 * m3x4.a20, m4x3.a00 * m3x4.a01 + m4x3.a01 * m3x4.a11 + m4x3.a02 * m3x4.a21, m4x3.a00 * m3x4.a02 + m4x3.a01 * m3x4.a12 + m4x3.a02 * m3x4.a22, m4x3.a00 * m3x4.a03 + m4x3.a01 * m3x4.a13 + m4x3.a02 * m3x4.a23,
		m4x3.a10 * m3x4.a00 + m4x3.a11 * m3x4.a10 + m4x3.a12 * m3x4.a20, m4x3.a10 * m3x4.a01 + m4x3.a11 * m3x4.a11 + m4x3.a12 * m3x4.a21, m4x3.a10 * m3x4.a02 + m4x3.a11 * m3x4.a12 + m4x3.a12 * m3x4.a22, m4x3.a10 * m3x4.a03 + m4x3.a11 * m3x4.a13 + m4x3.a12 * m3x4.a23,
		m4x3.a20 * m3x4.a00 + m4x3.a21 * m3x4.a10 + m4x3.a22 * m3x4.a20, m4x3.a20 * m3x4.a01 + m4x3.a21 * m3x4.a11 + m4x3.a22 * m3x4.a21, m4x3.a20 * m3x4.a02 + m4x3.a21 * m3x4.a12 + m4x3.a22 * m3x4.a22, m4x3.a20 * m3x4.a03 + m4x3.a21 * m3x4.a13 + m4x3.a22 * m3x4.a23,
		m4x3.a30 * m3x4.a00 + m4x3.a31 * m3x4.a10 + m4x3.a32 * m3x4.a20, m4x3.a30 * m3x4.a01 + m4x3.a31 * m3x4.a11 + m4x3.a32 * m3x4.a21, m4x3.a30 * m3x4.a02 + m4x3.a31 * m3x4.a12 + m4x3.a32 * m3x4.a22, m4x3.a30 * m3x4.a03 + m4x3.a31 * m3x4.a13 + m4x3.a32 * m3x4.a23
		);
}

inline VEC3D operator*(MAT33D& m, VEC3D& v)
{
	return VEC3D(m.a00 * v.x + m.a01 * v.y + m.a02 * v.z,
		m.a10 * v.x + m.a11 * v.y + m.a12 * v.z,
		m.a20 * v.x + m.a21 * v.y + m.a22 * v.z);
}

inline MAT33D operator* (MAT33D& n, MAT33D& m)
{
	return MAT33D(n.a00*m.a00 + n.a01*m.a10 + n.a02*m.a20, n.a00*m.a01 + n.a01*m.a11 + n.a02*m.a21, n.a00*m.a02 + n.a01*m.a12 + n.a02*m.a22
		, n.a10*m.a00 + n.a11*m.a10 + n.a12*m.a20, n.a10*m.a01 + n.a11*m.a11 + n.a12*m.a21, n.a10*m.a02 + n.a11*m.a12 + n.a12*m.a22
		, n.a20*m.a00 + n.a21*m.a10 + n.a22*m.a20, n.a20*m.a01 + n.a21*m.a11 + n.a22*m.a21, n.a20*m.a02 + n.a21*m.a12 + n.a22*m.a22);
}

inline MAT33D operator* (double c, MAT33D& m)
{
	return MAT33D(
		c * m.a00, c * m.a01, c * m.a02,
		c * m.a10, c * m.a11, c * m.a12,
		c * m.a20, c * m.a21, c * m.a22);
}

inline VEC3F operator*(MAT33F& m, VEC3F& v)
{
	return VEC3F(m.a00 * v.x + m.a01 * v.y + m.a02 * v.z,
		m.a10 * v.x + m.a11 * v.y + m.a12 * v.z,
		m.a20 * v.x + m.a21 * v.y + m.a22 * v.z);
}

inline
MAT44D opMiner(VEC3D& v)
{
	return MAT44D(
		0, -v.x, -v.y, -v.z,
		v.x, 0, v.z, -v.y,
		v.y, -v.z, 0, v.x,
		v.z, v.y, -v.x, 0);
}

inline MAT33D tilde(VEC3D& v)
{
	return MAT33D(
		0, -v.z, v.y,
		v.z, 0, -v.x,
		-v.y, v.x, 0);
}

inline VEC3D ep2e(EPD& ep)
{
	//double d1 = ((ep.e1 * ep.e3) + (ep.e0 * ep.e2)) / -((ep.e2 * ep.e3) - (ep.e0 * ep.e1));
	double m33 = 2.0 * (ep.e0 * ep.e0 + ep.e3 * ep.e3 - 0.5);
	if (m33 < 1.0e-15)
	{
		m33 = 0.0;
	}
	if (m33 > 1.0)
	{

		m33 = 1.0;
	}
	return VEC3D(
		atan2(((ep.e1 * ep.e3) + (ep.e0 * ep.e2)), -((ep.e2 * ep.e3) - (ep.e0 * ep.e1))),
		atan2(sqrt(1 - m33 * m33), m33),
		atan2(((ep.e1 * ep.e3) - (ep.e0 * ep.e2)), ((ep.e2 * ep.e3) + ep.e0 * ep.e1)));
}

inline MAT33D dA(EPD& ep, EPD& ev)
{

}

inline VEC3D global2local_eulerAngle(MAT33D& t, VEC3D& v)
{
	return VEC3D(
		t.a00 * v.x + t.a10 * v.y + t.a20 * v.z,
		t.a01 * v.x + t.a11 * v.y + t.a21 * v.z,
		t.a02 * v.x + t.a12 * v.y + t.a22 * v.z);
}

inline VEC3D local2global_eulerAngle(MAT33D& t, VEC3D& v)
{
	return VEC3D(
		t.a00 * v.x + t.a01 * v.y + t.a02 * v.z,
		t.a10 * v.x + t.a11 * v.y + t.a12 * v.z,
		t.a20 * v.x + t.a21 * v.y + t.a22 * v.z);
}

#endif