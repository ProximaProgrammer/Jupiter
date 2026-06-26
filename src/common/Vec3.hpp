#pragma once

#include <cmath>

namespace Jupiter
{

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& b) const { return {x + b.x, y + b.y, z + b.z}; }
    Vec3 operator-(const Vec3& b) const { return {x - b.x, y - b.y, z - b.z}; }
    Vec3 operator*(double a) const { return {a * x, a * y, a * z}; }
    Vec3 operator/(double a) const { return {x / a, y / a, z / a}; }
    Vec3& operator+=(const Vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
};

inline Vec3 operator*(double a, const Vec3& v) { return v * a; }
inline double Dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Vec3 Cross(const Vec3& a, const Vec3& b)
{
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
inline double Norm(const Vec3& v) { return std::sqrt(Dot(v, v)); }
inline Vec3 Unit(const Vec3& v)
{
    const double n = Norm(v);
    return n > 0.0 ? v / n : Vec3{1.0, 0.0, 0.0};
}

inline Vec3 SphericalBasisR(double theta, double phi)
{
    const double st = std::sin(theta), ct = std::cos(theta);
    const double sp = std::sin(phi), cp = std::cos(phi);
    return {st*cp, st*sp, ct};
}

inline Vec3 SphericalBasisTheta(double theta, double phi)
{
    const double st = std::sin(theta), ct = std::cos(theta);
    const double sp = std::sin(phi), cp = std::cos(phi);
    return {ct*cp, ct*sp, -st};
}

inline Vec3 SphericalBasisPhi(double phi)
{
    const double sp = std::sin(phi), cp = std::cos(phi);
    return {-sp, cp, 0.0};
}

inline Vec3 SphericalVelocityToCartesian(double theta, double phi, double vr, double vt, double vp)
{
    return SphericalBasisR(theta, phi)*vr + SphericalBasisTheta(theta, phi)*vt + SphericalBasisPhi(phi)*vp;
}

inline Vec3 PositionCartesian(double r, double theta, double phi)
{
    return SphericalBasisR(theta, phi) * r;
}

} // namespace Jupiter
