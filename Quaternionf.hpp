/*
 * Quaternionf.hpp
 *
 *  Created on: Aug 25, 2026
 *      Author: KAVINDU
 */

#ifndef COMMON_QUATERNIONF_HPP_
#define COMMON_QUATERNIONF_HPP_

#include "Vector3f.hpp"

struct Quaternionf {
	float w, x, y, z;

	Quaternionf() : w(1.0f), x(0.0f), y(0.0f), z(0.0f) {}
	Quaternionf(float _w, float _x, float _y, float _z) : w(_w), x(_x), y(_y), z(_z) {}

	Quaternionf operator*(const Quaternionf &o) const {
		return Quaternionf(
			w * o.w - x * o.x - y * o.y - z * o.z,
			w * o.x + x * o.w + y * o.z - z * o.y,
			w * o.y - x * o.z + y * o.w + z * o.x,
			w * o.z + x * o.y - y * o.x + z * o.w
		);
	}

	Vector3f rotate(const Vector3f &v) const {
		Quaternionf qv(0, v.x, v.y, v.z);
		Quaternionf qr = (*this) * qv * conjugate();
		return Vector3f(qr.x, qr.y, qr.z);
	}

	Quaternionf conjugate() const {
		return Quaternionf(w, -x, -y, -z);
	}
};

#endif /* COMMON_QUATERNIONF_HPP_ */
