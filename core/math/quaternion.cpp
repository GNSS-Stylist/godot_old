/*************************************************************************/
/*  quaternion.cpp                                                       */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2022 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2022 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "quaternion.h"

#include "core/math/basis.h"
#include "core/string/print_string.h"

real_t Quaternion::angle_to(const Quaternion &p_to) const {
	real_t d = dot(p_to);
	return Math::acos(CLAMP(d * d * 2 - 1, -1, 1));
}

// get_euler_xyz returns a vector containing the Euler angles in the format
// (ax,ay,az), where ax is the angle of rotation around x axis,
// and similar for other axes.
// This implementation uses XYZ convention (Z is the first rotation).
Vector3 Quaternion::get_euler_xyz() const {
	Basis m(*this);
	return m.get_euler(Basis::EULER_ORDER_XYZ);
}

// get_euler_yxz returns a vector containing the Euler angles in the format
// (ax,ay,az), where ax is the angle of rotation around x axis,
// and similar for other axes.
// This implementation uses YXZ convention (Z is the first rotation).
Vector3 Quaternion::get_euler_yxz() const {
#ifdef MATH_CHECKS
	ERR_FAIL_COND_V_MSG(!is_normalized(), Vector3(0, 0, 0), "The quaternion must be normalized.");
#endif
	Basis m(*this);
	return m.get_euler(Basis::EULER_ORDER_YXZ);
}

void Quaternion::operator*=(const Quaternion &p_q) {
	real_t xx = w * p_q.x + x * p_q.w + y * p_q.z - z * p_q.y;
	real_t yy = w * p_q.y + y * p_q.w + z * p_q.x - x * p_q.z;
	real_t zz = w * p_q.z + z * p_q.w + x * p_q.y - y * p_q.x;
	w = w * p_q.w - x * p_q.x - y * p_q.y - z * p_q.z;
	x = xx;
	y = yy;
	z = zz;
}

Quaternion Quaternion::operator*(const Quaternion &p_q) const {
	Quaternion r = *this;
	r *= p_q;
	return r;
}

bool Quaternion::is_equal_approx(const Quaternion &p_quaternion) const {
	return Math::is_equal_approx(x, p_quaternion.x) && Math::is_equal_approx(y, p_quaternion.y) && Math::is_equal_approx(z, p_quaternion.z) && Math::is_equal_approx(w, p_quaternion.w);
}

real_t Quaternion::length() const {
	return Math::sqrt(length_squared());
}

void Quaternion::normalize() {
	*this /= length();
}

Quaternion Quaternion::normalized() const {
#ifdef MATH_CHECKS
	ERR_FAIL_COND_V_MSG(length() == 0, Quaternion(), "Normalized(): Length of the quaternion can't be zero.");
#endif
	return *this / length();
}

bool Quaternion::is_normalized() const {
	return Math::is_equal_approx(length_squared(), 1, (real_t)UNIT_EPSILON); //use less epsilon
}

Quaternion Quaternion::static_conjugate(const Quaternion& quat) {
	return Quaternion(-quat.x, -quat.y, -quat.z, quat.w);
}

Quaternion Quaternion::conjugate() const {
	return static_conjugate(*this);
}

Quaternion Quaternion::static_inverse(const Quaternion& quat) {
	return static_conjugate(quat) / quat.length_squared();
}

Quaternion Quaternion::inverse() const {
	return static_inverse(*this);
}

Quaternion Quaternion::static_slerp(const Quaternion& p_from, const Quaternion& p_to, const real_t& p_weight)
{
#ifdef MATH_CHECKS
	// Due to rounding errors(?) these sometimes (although rarely) trigger without a reason so they are now disabled.
//	ERR_FAIL_COND_V_MSG(!p_from.is_normalized(), Quaternion(), "The start quaternion must be normalized.");
//	ERR_FAIL_COND_V_MSG(!p_to.is_normalized(), Quaternion(), "The end quaternion must be normalized.");
#endif
	Quaternion to1;
	real_t omega, cosom, sinom, scale0, scale1;

	// calc cosine
	cosom = p_from.dot(p_to);

	// adjust signs (if necessary)
	if (cosom < 0.0f) {
		cosom = -cosom;
		to1.x = -p_to.x;
		to1.y = -p_to.y;
		to1.z = -p_to.z;
		to1.w = -p_to.w;
	}
	else {
		to1.x = p_to.x;
		to1.y = p_to.y;
		to1.z = p_to.z;
		to1.w = p_to.w;
	}

	// calculate coefficients

	if ((1.0f - cosom) > (real_t)CMP_EPSILON) {
		// standard case (slerp)
		omega = Math::acos(cosom);
		sinom = Math::sin(omega);
		scale0 = Math::sin((1.0 - p_weight) * omega) / sinom;
		scale1 = Math::sin(p_weight * omega) / sinom;
	}
	else {
		// "from" and "to" quaternions are very close
		//  ... so we can do a linear interpolation
		scale0 = 1.0f - p_weight;
		scale1 = p_weight;
	}
	// calculate final values
	return Quaternion(
		scale0 * p_from.x + scale1 * to1.x,
		scale0 * p_from.y + scale1 * to1.y,
		scale0 * p_from.z + scale1 * to1.z,
		scale0 * p_from.w + scale1 * to1.w);
}

Quaternion Quaternion::slerp(const Quaternion& p_to, const real_t& p_weight) const {
	return static_slerp(*this, p_to, p_weight);
}

Quaternion Quaternion::static_slerpni(const Quaternion& p_from, const Quaternion& p_to, const real_t& p_weight) {
#ifdef MATH_CHECKS
	// Due to rounding errors(?) these sometimes (although rarely) trigger without a reason so they are now disabled.
//	ERR_FAIL_COND_V_MSG(!p_from.is_normalized(), Quaternion(), "The start quaternion must be normalized.");
//	ERR_FAIL_COND_V_MSG(!p_to.is_normalized(), Quaternion(), "The end quaternion must be normalized.");
#endif
	real_t dot = p_from.dot(p_to);

	if (Math::absf(dot) > 0.9999f) {
		// "from" and "to" quaternions are very close
		//  ... so we can do a linear interpolation
		real_t startWeight = 1.0f - p_weight;
		return Quaternion(p_from.x * startWeight + p_to.x * p_weight, p_from.y * startWeight + p_to.y * p_weight, p_from.z * startWeight + p_to.z * p_weight, p_from.w * startWeight + p_to.w * p_weight).normalized();
//		_err_print_error(FUNCTION_STR, __FILE__, __LINE__, "Notification: ", "slerpni->linear");
	}

	real_t theta = Math::acos(dot),
		   sinT = 1.0f / Math::sin(theta),
		   newFactor = Math::sin(p_weight * theta) * sinT,
		   invFactor = Math::sin((1.0f - p_weight) * theta) * sinT;

	return Quaternion(invFactor * p_from.x + newFactor * p_to.x,
			invFactor * p_from.y + newFactor * p_to.y,
			invFactor * p_from.z + newFactor * p_to.z,
			invFactor * p_from.w + newFactor * p_to.w);
}

Quaternion Quaternion::slerpni(const Quaternion& p_to, const real_t& p_weight) const {
	return static_slerpni(*this, p_to, p_weight);
}

static real_t cubic_interpolate_real(const real_t p_pre_a, const real_t p_a, const real_t p_b, const real_t p_post_b, real_t p_c) {
	// This is cloned straight from animation.cpp

	real_t p0 = p_pre_a;
	real_t p1 = p_a;
	real_t p2 = p_b;
	real_t p3 = p_post_b;

	real_t t = p_c;
	real_t t2 = t * t;
	real_t t3 = t2 * t;

	return 0.5f *
		((p1 * 2.0f) +
			(-p0 + p2) * t +
			(2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
			(-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}

static Quaternion flip_to_shortest(const Quaternion& compared, const Quaternion& flip) {
	if (compared.dot(flip) < 0.0f) {
		// Return flipped
		return Quaternion(-flip.x, -flip.y, -flip.z, -flip.w);
	}
	else {
		// No flipping necessary
		return flip;
	}
}

Quaternion Quaternion::static_cubic_slerp(const Quaternion& p_q0, const Quaternion& p_q1, const Quaternion& p_q2, const Quaternion& p_q3, const real_t& p_t) {
	if (true) {
		// Flip quaternions to shortest path if necessary (p_q1 used as the reference)
		Quaternion q0 = flip_to_shortest(p_q1, p_q0);
		Quaternion q1 = p_q1;
		Quaternion q2 = flip_to_shortest(p_q1, p_q2);
		Quaternion q3 = flip_to_shortest(q2, p_q3);	// Need to compare with possibly already flipped q2

		return Quaternion(
			cubic_interpolate_real(q0.x, q1.x, q2.x, q3.x, p_t),
			cubic_interpolate_real(q0.y, q1.y, q2.y, q3.y, p_t),
			cubic_interpolate_real(q0.z, q1.z, q2.z, q3.z, p_t),
			cubic_interpolate_real(q0.w, q1.w, q2.w, q3.w, p_t)
		);
	}

	if (false) {
		// Test code to compare slerp and lerp
		// (just return lerp)
		// (This was active in the previous example)
		return static_lerp(p_q1, p_q2, p_t);
	}



	Quaternion q1 = p_q1;

	Quaternion q0 = (p_q0 + p_q1).length_squared() < (p_q0 - p_q1).length_squared() ? -p_q0 : p_q0;
	Quaternion q2 = (p_q1 + p_q2).length_squared() < (p_q1 - p_q2).length_squared() ? -p_q2 : p_q2;
	Quaternion q3 = (q2 + p_q3).length_squared() < (q2 - p_q3).length_squared() ? -p_q3 : p_q3;

// Original calculations for outA and outB (causing log to trigger error):
//	Quaternion outA = static_exp(-0.25 * (static_log(q2 * static_inverse(q1)) + static_log(q0 * static_inverse(q1)))) * q1;
//	Quaternion outB = static_exp(-0.25 * (static_log(q3 * static_inverse(q2)) + static_log(q1 * static_inverse(q2)))) * q2;

	const real_t lengthLimit = 1e-6;	// 1e-9 caused static_log to trigger "Length of the vector part can not be zero."-error

	// Following lines for calculating outA and outB added to prevent
	// static_log-function to trigger "Length of the vector part can not be zero."-error.
	// (length of the vector part is used as divider there so it can't be zero).
	// lengthLimit is quite arbitrary, current value just "seemed to work".
	// TODO: This could be optimized, inputs for log-function calls are calculated twice now.

	Quaternion outA;
	if ((q0 * static_inverse(q1)).vector_part().length_squared() < lengthLimit) {
		outA = q1;
	}
	else if ((q2 * static_inverse(q1)).vector_part().length_squared() < lengthLimit) {
		outA = static_exp(-0.25 * static_log(q0 * static_inverse(q1))) * q1;
	}
	else {
		outA = static_exp(-0.25 * (static_log(q2 * static_inverse(q1)) + static_log(q0 * static_inverse(q1)))) * q1;
	}

	Quaternion outB;
	if ((q3 * static_inverse(q2)).vector_part().length_squared() < lengthLimit) {
		outB = q2;
	}
	else if ((q1 * static_inverse(q2)).vector_part().length_squared() < lengthLimit) {
		outB = static_exp(-0.25 * (static_log(q3 * static_inverse(q2)))) * q2;
	}
	else {
		outB = static_exp(-0.25 * (static_log(q3 * static_inverse(q2)) + static_log(q1 * static_inverse(q2)))) * q2;
	}

	Quaternion outC = q2;

	return static_squad(q1, outA, outB, outC, p_t);
}

Quaternion Quaternion::static_squad(const Quaternion& q1, const Quaternion& a, const Quaternion& b, const Quaternion& c, const real_t t) {
	return static_slerpni(static_slerpni(q1, c, t), static_slerpni(a, b, t), 2.0 * t * (1.0 - t)).normalized();
//	return static_slerpni(static_slerpni(q1, c, t), static_slerpni(a, b, t), 2.0 * t * (1.0 - t)).normalized();
}

Quaternion Quaternion::cubic_slerp(const Quaternion& p_q, const Quaternion& p_prep, const Quaternion& p_postq, const real_t& p_t) const {
	return static_cubic_slerp(p_prep, *this, p_q, p_postq, p_t);
}

Quaternion Quaternion::static_log(const Quaternion& quat) {
	Vector3 v = quat.vector_part();
	real_t vLength = v.length();
#ifdef MATH_CHECKS
	ERR_FAIL_COND_V_MSG(vLength == 0, Quaternion(), "Length of the vector part can not be zero.");
#endif

	real_t s = quat.w;

	real_t scalarPart = Math::log(quat.length());
	Vector3 vecPart = v / vLength * acos(s / quat.length());

	return Quaternion(vecPart.x, vecPart.y, vecPart.z, scalarPart);
}

Quaternion Quaternion::log() const {
	return static_log(*this);
}

Quaternion Quaternion::static_exp(const Quaternion& quat)
{
	Vector3 v = quat.vector_part();
	real_t vLength = v.length();
#ifdef MATH_CHECKS
	ERR_FAIL_COND_V_MSG(vLength == 0, Quaternion(), "Length of the vector part can not be zero.");
#endif

	real_t s = quat.w;

	real_t scalarPart = cos(vLength);
	Vector3 vecPart = v / vLength * sin(vLength);

	return Math::exp(s) * Quaternion(vecPart.x, vecPart.y, vecPart.z, scalarPart);
}

Quaternion Quaternion::exp() const {
	return static_exp(*this);
}

Quaternion::operator String() const {
	return "(" + String::num_real(x, false) + ", " + String::num_real(y, false) + ", " + String::num_real(z, false) + ", " + String::num_real(w, false) + ")";
}

Vector3 Quaternion::get_axis() const {
	if (Math::abs(w) > 1 - CMP_EPSILON) {
		return Vector3(x, y, z);
	}
	real_t r = ((real_t)1) / Math::sqrt(1 - w * w);
	return Vector3(x * r, y * r, z * r);
}

real_t Quaternion::get_angle() const {
	return 2 * Math::acos(w);
}

Quaternion::Quaternion(const Vector3 &p_axis, real_t p_angle) {
#ifdef MATH_CHECKS
	ERR_FAIL_COND_MSG(!p_axis.is_normalized(), "The axis Vector3 must be normalized.");
#endif
	real_t d = p_axis.length();
	if (d == 0) {
		x = 0;
		y = 0;
		z = 0;
		w = 0;
	} else {
		real_t sin_angle = Math::sin(p_angle * 0.5f);
		real_t cos_angle = Math::cos(p_angle * 0.5f);
		real_t s = sin_angle / d;
		x = p_axis.x * s;
		y = p_axis.y * s;
		z = p_axis.z * s;
		w = cos_angle;
	}
}

// Euler constructor expects a vector containing the Euler angles in the format
// (ax, ay, az), where ax is the angle of rotation around x axis,
// and similar for other axes.
// This implementation uses YXZ convention (Z is the first rotation).
Quaternion::Quaternion(const Vector3 &p_euler) {
	real_t half_a1 = p_euler.y * 0.5f;
	real_t half_a2 = p_euler.x * 0.5f;
	real_t half_a3 = p_euler.z * 0.5f;

	// R = Y(a1).X(a2).Z(a3) convention for Euler angles.
	// Conversion to quaternion as listed in https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf (page A-6)
	// a3 is the angle of the first rotation, following the notation in this reference.

	real_t cos_a1 = Math::cos(half_a1);
	real_t sin_a1 = Math::sin(half_a1);
	real_t cos_a2 = Math::cos(half_a2);
	real_t sin_a2 = Math::sin(half_a2);
	real_t cos_a3 = Math::cos(half_a3);
	real_t sin_a3 = Math::sin(half_a3);

	x = sin_a1 * cos_a2 * sin_a3 + cos_a1 * sin_a2 * cos_a3;
	y = sin_a1 * cos_a2 * cos_a3 - cos_a1 * sin_a2 * sin_a3;
	z = -sin_a1 * sin_a2 * cos_a3 + cos_a1 * cos_a2 * sin_a3;
	w = sin_a1 * sin_a2 * sin_a3 + cos_a1 * cos_a2 * cos_a3;
}

Vector3 Quaternion::vector_part() const {
	return Vector3(x, y, z);
}

Quaternion Quaternion::static_lerp(const Quaternion& p_from, const Quaternion& p_to, const real_t& p_weight)
{
#ifdef MATH_CHECKS
	// Due to rounding errors(?) these sometimes (although rarely) trigger without a reason so they are now disabled.
//	ERR_FAIL_COND_V_MSG(!p_from.is_normalized(), Quaternion(), "The start quaternion must be normalized.");
//	ERR_FAIL_COND_V_MSG(!p_to.is_normalized(), Quaternion(), "The end quaternion must be normalized.");
#endif
	Quaternion to1;
	real_t cosom, scale0, scale1;

	// calc cosine
	cosom = p_from.dot(p_to);

	// adjust signs (if necessary)
	if (cosom < 0.0) {
		cosom = -cosom;
		to1.x = -p_to.x;
		to1.y = -p_to.y;
		to1.z = -p_to.z;
		to1.w = -p_to.w;
	}
	else {
		to1.x = p_to.x;
		to1.y = p_to.y;
		to1.z = p_to.z;
		to1.w = p_to.w;
	}

	// calculate coefficients

	//  When lerping we can just do a linear interpolation
	scale0 = 1.0 - p_weight;
	scale1 = p_weight;

	// calculate final values
	return Quaternion(
		scale0 * p_from.x + scale1 * to1.x,
		scale0 * p_from.y + scale1 * to1.y,
		scale0 * p_from.z + scale1 * to1.z,
		scale0 * p_from.w + scale1 * to1.w);
}

Quaternion Quaternion::lerp(const Quaternion& p_to, const real_t& p_weight) const {
	return static_lerp(*this, p_to, p_weight);
}

Quaternion Quaternion::static_lerpni(const Quaternion& p_from, const Quaternion& p_to, const real_t& p_weight) {
#ifdef MATH_CHECKS
	// Due to rounding errors(?) these sometimes (although rarely) trigger without a reason so they are now disabled.
//	ERR_FAIL_COND_V_MSG(!p_from.is_normalized(), Quaternion(), "The start quaternion must be normalized.");
//	ERR_FAIL_COND_V_MSG(!p_to.is_normalized(), Quaternion(), "The end quaternion must be normalized.");
#endif
		//  Lerp is just component-wise linear interpolation
	real_t startWeight = 1.0 - p_weight;
	return Quaternion(p_from.x * startWeight + p_to.x * p_weight, p_from.y * startWeight + p_to.y * p_weight, p_from.z * startWeight + p_to.z * p_weight, p_from.w * startWeight + p_to.w * p_weight).normalized();
}

Quaternion Quaternion::lerpni(const Quaternion& p_to, const real_t& p_weight) const {
	return static_lerpni(*this, p_to, p_weight);
}

#if false
static real_t hermite_interpolate_real(const real_t y0, const real_t y1, const real_t y2, const real_t y3,
	const real_t x0, const real_t x1, const real_t x2, const real_t x3,
	real_t t) {

	// Everything here is based on this wikipedia article:
	// https://en.wikipedia.org/wiki/Cubic_Hermite_spline

	real_t m0;
	real_t m1;

	// Calculate tangents using Catmull-Rom method
	if ((x2 - x0) != real_t(0.0)) {
		m0 = 0.5 * ((y2 - y0) / (x2 - x0));
	}
	else {
		m0 = 0.0;
	}

	if ((x3 - x1) != real_t(0.0)) {
		m1 = 0.5 * ((y3 - y1) / (x3 - x1));
	}
	else {
		m1 = 0.0;
	}

/*
	// "Finite difference":
	m0 = 0.5 * ((y2 -y1) / (x2 - x1) + (y1 - y0) / (x1 - x0));
	m1 = 0.5 * ((y3 - y2) / (x3 - x2) + (y2 - y1) / (x2 - x1));
*/

	return (2.0 * pow(t, 3.0) - 3.0 * pow(t, 2.0) + 1.0) * y1 +
		(pow(t, 3.0) - 2.0 * pow(t, 2.0) + t) * m0 +
		(-2.0 * pow(t, 3.0) + 3.0 * pow(t, 2.0)) * y2 +
		(pow(t, 3.0) - pow(t, 2.0)) * m1;
}
#else


static real_t hermite_interpolate_real(const real_t y0, const real_t y1, const real_t y2, const real_t y3,
	const real_t x0, const real_t x1, const real_t x2, const real_t x3,
	real_t t) {

	// Everything here is based on this wikipedia article:
	// https://en.wikipedia.org/wiki/Cubic_Hermite_spline

	real_t h00 = 2.0 * pow(t, 3) - 3.0 * pow(t, 2.0) + 1.0;
	real_t h10 = pow(t, 3.0) - 2.0 * pow(t, 2.0) + t;
	real_t h01 = -2.0 * pow(t, 3.0) + 3.0 * pow(t, 2.0);
	real_t h11 = pow(t, 3) - pow(t, 2.0);

	real_t m0;
	real_t m1;

	m0 = 0.5 * (y2 - y0);
	m1 = 0.5 * (y3 - y1);

	// "Finite difference":
	if (((x2 - x1) == 0.0) || ((x1 - x0) == 0.0)) {
		m0 = 0.0;
	}
	else {
		m0 = 0.5 * ((y2 - y1) / (x2 - x1) + (y1 - y0) / (x1 - x0)) * (x2 - x1);
	}

	if (((x3 - x2) == 0.0) || ((x2 - x1) == 0.0)) {
		m1 = 0.0;
	}
	else {
		m1 = 0.5 * ((y3 - y2) / (x3 - x2) + (y2 - y1) / (x2 - x1)) * (x2 - x1);
	}

	return h00 * y1 +
		h10 * m0 +
		h01 * y2 +
		h11 * m1;
}
#endif

Quaternion Quaternion::static_cubic_hermite_spline_interpolate(const Quaternion& p_q0, const Quaternion& p_q1, const Quaternion& p_q2, const Quaternion& p_q3,
	const real_t p_t0, const real_t p_t1, const real_t p_t2, const real_t p_t3,
	const real_t& p_t) const {

	// Note: p_t must be 0...1 (spanning to range p_t1...p_t2)

	// Flip quaternions to shortest path if necessary (p_q1 used as the reference)
	Quaternion q0 = flip_to_shortest(p_q1, p_q0);
	Quaternion q1 = p_q1;
	Quaternion q2 = flip_to_shortest(p_q1, p_q2);
	Quaternion q3 = flip_to_shortest(q2, p_q3);	// Need to compare with possibly already flipped q2

	return Quaternion(
		hermite_interpolate_real(q0.x, q1.x, q2.x, q3.x, p_t0, p_t1, p_t2, p_t3, p_t),
		hermite_interpolate_real(q0.y, q1.y, q2.y, q3.y, p_t0, p_t1, p_t2, p_t3, p_t),
		hermite_interpolate_real(q0.z, q1.z, q2.z, q3.z, p_t0, p_t1, p_t2, p_t3, p_t),
		hermite_interpolate_real(q0.w, q1.w, q2.w, q3.w, p_t0, p_t1, p_t2, p_t3, p_t)
	);
}

Quaternion Quaternion::cubic_hermite_spline_interpolate(const Quaternion& p_b, const Quaternion& p_pre_a, const Quaternion& p_post_b,
	const real_t t_a, const real_t t_b, const real_t t_pre_a, const real_t t_post_b,
	const real_t& p_t) const {

//	return static_slerp(*this, p_b, p_t);

	return static_cubic_hermite_spline_interpolate(p_pre_a, *this, p_b, p_post_b, t_pre_a, t_a, t_b, t_post_b, p_t);
}

