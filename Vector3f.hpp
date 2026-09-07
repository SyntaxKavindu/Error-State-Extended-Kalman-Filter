/*
 * Vector3f.hpp
 *
 *  Created on: Aug 25, 2026
 *      Author: KAVINDU
 */

#ifndef COMMON_VECTOR3F_HPP_
#define COMMON_VECTOR3F_HPP_

#include <cmath>   // sqrtf

struct Vector3f {
    float x, y, z;
    Vector3f() : x(0.0f), y(0.0f), z(0.0f) {}
    Vector3f(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

    Vector3f operator+(const Vector3f &o) const { return Vector3f(x + o.x, y + o.y, z + o.z); }
    Vector3f operator-(const Vector3f &o) const { return Vector3f(x - o.x, y - o.y, z - o.z); }
    Vector3f operator*(float s) const { return Vector3f(x * s, y * s, z * s); }
    Vector3f operator/(float s) const { return Vector3f(x / s, y / s, z / s); }

    float length() const { return sqrtf(x * x + y * y + z * z); }
    float dot(const Vector3f &o) const { return x * o.x + y * o.y + z * o.z; }
};

#endif /* COMMON_VECTOR3F_HPP_ */
