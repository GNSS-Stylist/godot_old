/*************************************************************************/
/*  vector3.cpp                                                          */
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

#include "vector3.h"

#include "core/math/basis.h"
#include "core/math/vector2.h"
#include "core/math/vector3i.h"
#include "core/string/ustring.h"

void Vector3::rotate(const Vector3 &p_axis, const real_t p_phi) {
	*this = Basis(p_axis, p_phi).xform(*this);
}

Vector3 Vector3::rotated(const Vector3 &p_axis, const real_t p_phi) const {
	Vector3 r = *this;
	r.rotate(p_axis, p_phi);
	return r;
}

void Vector3::set_axis(const int p_axis, const real_t p_value) {
	ERR_FAIL_INDEX(p_axis, 3);
	coord[p_axis] = p_value;
}

real_t Vector3::get_axis(const int p_axis) const {
	ERR_FAIL_INDEX_V(p_axis, 3, 0);
	return operator[](p_axis);
}

Vector3 Vector3::clamp(const Vector3 &p_min, const Vector3 &p_max) const {
	return Vector3(
			CLAMP(x, p_min.x, p_max.x),
			CLAMP(y, p_min.y, p_max.y),
			CLAMP(z, p_min.z, p_max.z));
}

void Vector3::snap(const Vector3 p_step) {
	x = Math::snapped(x, p_step.x);
	y = Math::snapped(y, p_step.y);
	z = Math::snapped(z, p_step.z);
}

Vector3 Vector3::snapped(const Vector3 p_step) const {
	Vector3 v = *this;
	v.snap(p_step);
	return v;
}

Vector3 Vector3::limit_length(const real_t p_len) const {
	const real_t l = length();
	Vector3 v = *this;
	if (l > 0 && p_len < l) {
		v /= l;
		v *= p_len;
	}

	return v;
}

Vector3 Vector3::cubic_interpolate(const Vector3 &p_b, const Vector3 &p_pre_a, const Vector3 &p_post_b, const real_t p_weight) const {
	Vector3 res = *this;
	res.x = Math::cubic_interpolate(res.x, p_b.x, p_pre_a.x, p_post_b.x, p_weight);
	res.y = Math::cubic_interpolate(res.y, p_b.y, p_pre_a.y, p_post_b.y, p_weight);
	res.z = Math::cubic_interpolate(res.z, p_b.z, p_pre_a.z, p_post_b.z, p_weight);
	return res;
}

static Vector3 hermite_interpolate_vector3(const Vector3& y0, const Vector3& y1, const Vector3& y2, const Vector3& y3,
	const real_t x0, const real_t x1, const real_t x2, const real_t x3,
	real_t t) {

	// Everything here is based on this wikipedia article:
	// https://en.wikipedia.org/wiki/Cubic_Hermite_spline

	real_t h00 = 2.0 * pow(t, 3) - 3.0 * pow(t, 2.0) + 1.0;
	real_t h10 = pow(t, 3.0) - 2.0 * pow(t, 2.0) + t;
	real_t h01 = -2.0 * pow(t, 3.0) + 3.0 * pow(t, 2.0);
	real_t h11 = pow(t, 3) - pow(t, 2.0);

	Vector3 m0;
	Vector3 m1;

	m0 = 0.5 * (y2 - y0);
	m1 = 0.5 * (y3 - y1);

	// "Finite difference":
	if (((x2 - x1) == 0.0) || ((x1 - x0) == 0.0)) {
		m0 = Vector3();
	}
	else {
		m0 = 0.5 * ((y2 - y1) / (x2 - x1) + (y1 - y0) / (x1 - x0)) * (x2 - x1);
	}

	if (((x3 - x2) == 0.0) || ((x2 - x1) == 0.0)) {
		m1 = Vector3();
	}
	else {
		m1 = 0.5 * ((y3 - y2) / (x3 - x2) + (y2 - y1) / (x2 - x1)) * (x2 - x1);
	}

	return h00 * y1 +
		h10 * m0 +
		h01 * y2 +
		h11 * m1;
}


Vector3 Vector3::cubic_hermite_spline_interpolate(const Vector3& p_b, const Vector3& p_pre_a, const Vector3& p_post_b, const real_t t_a, const real_t t_b, const real_t t_pre_a, const real_t t_post_b, const real_t p_weight) const {
	Vector3 p0 = p_pre_a;
	Vector3 p1 = *this;
	Vector3 p2 = p_b;
	Vector3 p3 = p_post_b;

	real_t t = p_weight;
	real_t t2 = t * t;
	real_t t3 = t2 * t;

	Vector3 out;
	out = hermite_interpolate_vector3(p_pre_a, *this, p_b, p_post_b, t_pre_a, t_a, t_b, t_post_b, t);
	return out;
}


Vector3 Vector3::move_toward(const Vector3 &p_to, const real_t p_delta) const {
	Vector3 v = *this;
	Vector3 vd = p_to - v;
	real_t len = vd.length();
	return len <= p_delta || len < (real_t)CMP_EPSILON ? p_to : v + vd / len * p_delta;
}

Vector2 Vector3::octahedron_encode() const {
	Vector3 n = *this;
	n /= Math::abs(n.x) + Math::abs(n.y) + Math::abs(n.z);
	Vector2 o;
	if (n.z >= 0.0f) {
		o.x = n.x;
		o.y = n.y;
	} else {
		o.x = (1.0f - Math::abs(n.y)) * (n.x >= 0.0f ? 1.0f : -1.0f);
		o.y = (1.0f - Math::abs(n.x)) * (n.y >= 0.0f ? 1.0f : -1.0f);
	}
	o.x = o.x * 0.5f + 0.5f;
	o.y = o.y * 0.5f + 0.5f;
	return o;
}

Vector3 Vector3::octahedron_decode(const Vector2 &p_oct) {
	Vector2 f(p_oct.x * 2.0f - 1.0f, p_oct.y * 2.0f - 1.0f);
	Vector3 n(f.x, f.y, 1.0f - Math::abs(f.x) - Math::abs(f.y));
	float t = CLAMP(-n.z, 0.0f, 1.0f);
	n.x += n.x >= 0 ? -t : t;
	n.y += n.y >= 0 ? -t : t;
	return n.normalized();
}

Basis Vector3::outer(const Vector3 &p_with) const {
	Vector3 row0(x * p_with.x, x * p_with.y, x * p_with.z);
	Vector3 row1(y * p_with.x, y * p_with.y, y * p_with.z);
	Vector3 row2(z * p_with.x, z * p_with.y, z * p_with.z);

	return Basis(row0, row1, row2);
}

bool Vector3::is_equal_approx(const Vector3 &p_v) const {
	return Math::is_equal_approx(x, p_v.x) && Math::is_equal_approx(y, p_v.y) && Math::is_equal_approx(z, p_v.z);
}

Vector3::operator String() const {
	return "(" + String::num_real(x, false) + ", " + String::num_real(y, false) + ", " + String::num_real(z, false) + ")";
}

Vector3::operator Vector3i() const {
	return Vector3i(x, y, z);
}
