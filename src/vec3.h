#pragma once

#include <cmath>

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}

    explicit Vec3(double val) : x(val), y(val), z(val) {}

    Vec3(double _x, double _y, double _z) {
        x = _x;
        y = _y;
        z = _z;
    }
    double max_component_val() const { return x > y && x > z ? x : y > z ? y : z; }

    Vec3 abs() const { return Vec3(::abs(x), ::abs(y), ::abs(z)); }

    Vec3 norm() const { return *this * (1 / sqrt(x * x + y * y + z * z)); }

    double dot(const Vec3 &b) const { return x * b.x + y * b.y + z * b.z; }

    Vec3 cross(const Vec3 &b) const {
        return Vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }

    double luminance() const { return x * 0.2126 + y * 0.7152 + z * 0.0722; }

    Vec3 operator+(const Vec3 &b) const { return Vec3(x + b.x, y + b.y, z + b.z); }

    void operator+=(const Vec3 &b) {
        x += b.x;
        y += b.y;
        z += b.z;
    }

    Vec3 operator-(const Vec3 &b) const { return Vec3(x - b.x, y - b.y, z - b.z); }

    Vec3 operator*(double b) const { return Vec3(x * b, y * b, z * b); }

    void operator*=(double b) {
        x *= b;
        y *= b;
        z *= b;
    }

    Vec3 operator*(const Vec3 &b) const { return Vec3(x * b.x, y * b.y, z * b.z); }

    void operator*=(const Vec3 &b) {
        x *= b.x;
        y *= b.y;
        z *= b.z;
    }
};

inline Vec3 operator*(double f, const Vec3 &v) { return v * f; }
