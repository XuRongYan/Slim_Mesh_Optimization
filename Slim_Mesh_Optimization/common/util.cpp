#include "util.h"

#include <fstream>
#include <cfloat>
#include <iostream>

namespace common
{
inline double clamp(float a, float low, float high){
    return (a > low && a < high) ? a : ((a <= low) ? low : high);
}

inline int dblcmp(float a, float b)
{
    return a > b ? 1 : (fabs(a - b) < eps ? 0 : -1);
}

int barycentric(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                const Eigen::Vector3f &c, const Eigen::Vector3f &p,
                Eigen::Vector3f &bary)
{
    Eigen::Vector3f u = b - a;
    Eigen::Vector3f v = c - a;
    Eigen::Vector3f w = p - a;

    Eigen::Vector3f vw = v.cross(w);
    Eigen::Vector3f vu = v.cross(u);

    // test sign of r
    if (vw.dot(vu) < 0)
        return 0;

    Eigen::Vector3f uw = u.cross(w);
    Eigen::Vector3f uv = -vu;

    // test sign of t
    if (uw.dot(uv) < 0)
        return 0;

    double denom = uv.norm();
    bary[1] = vw.norm() / denom;
    bary[2] = uw.norm() / denom;
    bary[0] = 1 - bary[1] - bary[2];

    return (bary[0] >= 0);
}

Eigen::Vector3f closestPointOnTriangle(const Eigen::Vector3f &p0,
                                       const Eigen::Vector3f &p1,
                                       const Eigen::Vector3f &p2,
                                       const Eigen::Vector3f &srcPos,
                                       Eigen::Vector3f &bary)
{
    const Eigen::Vector3f edge0 = p1 - p0;
    const Eigen::Vector3f edge1 = p2 - p0;
    const Eigen::Vector3f v0 = p0 - srcPos;

    const double a = edge0.dot(edge0);
    const double b = edge0.dot(edge1);
    const double c = edge1.dot(edge1);
    const double d = edge0.dot(v0);
    const double e = edge1.dot(v0);

    const double det = a*c - b*b;
    double s = b*e - c*d;
    double t = b*d - a*e;
    if ( s + t < det ){
        if ( s < 0.f ){
            if ( t < 0.f ){
                if ( d < 0.f ){
                    s = clamp( -d/a, 0.f, 1.f );
                    t = 0.f;
                }
                else{
                    s = 0.f;
                    t = clamp( -e/c, 0.f, 1.f );
                }
            }
            else{
                s = 0.f;
                t = clamp( -e/c, 0.f, 1.f );
            }
        }
        else if ( t < 0.f ){
            s = clamp( -d/a, 0.f, 1.f );
            t = 0.f;
        }
        else{
            float invDet = 1.f / det;
            s *= invDet;
            t *= invDet;
        }
    }
    else{
        if ( s < 0.f ){
            float tmp0 = b+d;
            float tmp1 = c+e;
            if ( tmp1 > tmp0 ){
                float numer = tmp1 - tmp0;
                float denom = a-2*b+c;
                s = clamp( numer/denom, 0.f, 1.f );
                t = 1-s;
            }
            else{
                t = clamp( -e/c, 0.f, 1.f );
                s = 0.f;
            }
        }
        else if ( t < 0.f ){
            if ( a+d > b+e ){
                float numer = c+e-b-d;
                float denom = a-2*b+c;
                s = clamp( numer/denom, 0.f, 1.f );
                t = 1-s;
            }
            else{
                s = clamp( -e/c, 0.f, 1.f );
                t = 0.f;
            }
        }
        else{
            float numer = c+e-b-d;
            float denom = a-2*b+c;
            s = clamp( numer/denom, 0.f, 1.f );
            t = 1.f - s;
        }
    }
    bary[0] = 1 - s - t;
    bary[1] = s;
    bary[2] = t;
    return p0 + s * edge0 + t * edge1;
}

int line_triangle_intersect(const Eigen::Vector3f &o, const Eigen::Vector3f &d,
                            const Eigen::Vector3f &p0, const Eigen::Vector3f &p1,
                            const Eigen::Vector3f &p2, Eigen::Vector3f &bary, Eigen::Vector3f &pt)
{
    const double eps = 5e-6;
    const Eigen::Vector3f e1 = p1 - p0;
    const Eigen::Vector3f e2 = p2 - p0;
    const Eigen::Vector3f p = d.cross(e2);
    double delta = p.dot(e1);

    if(delta >-eps && delta < eps)
        return 0;

    delta = 1.0f / delta;
    const Eigen::Vector3f s = o - p0;
    const double u = s.dot(p) * delta;
    if(u < 0.0f || u > 1.0f)
        return 0;

    const Eigen::Vector3f q = s.cross(e1);
    const double v = d.dot(q) * delta;
    if(v < 0.0f || v > 1.0f)
        return 0;

    const double w = 1 - u - v;
    if(w < 0.0f || w > 1.0f)
        return 0;

    bary << w,u,v;

    const double t = e2.dot(q) * delta;
    pt = o + d * t;
    return 1;
}

int line_plane_intersect(const Eigen::Vector3f &p, const Eigen::Vector3f &n,
                         const Eigen::Vector3f &p0, const Eigen::Vector3f &p1,
                         Eigen::Vector3f &q)
{
    float a = n.dot(p - p0);
    float b = n.dot(p1 - p0);
    if(abs(b) > 1e-6){
        float s = a / b;
        if(s > 1.0f || s < 0.0f){
            return 0;
        }
        else if(s < 1e-6){
            q = p0;
            return 1;
        }
        else if(s < 1.0){
            q = p0 + (p1 - p0) * s;
            return 3;
        }else{
            q = p1;
        }
    }
    return 0;
}

float distPointLineSquared(const Eigen::Vector3f& _p,
                           const Eigen::Vector3f& _v0,
                           const Eigen::Vector3f& _v1,
                           Eigen::Vector3f*       _min_v )
{
    Eigen::Vector3f d1(_p-_v0), d2(_v1-_v0),min_v(_v0);;
    float t = d1.dot(d2)/ d2.squaredNorm();

    if (t >  1.0)
        d1 = _p - (min_v = _v1);
    else if (t > 0.0)
        d1 = _p - (min_v = _v0 + d2*t);
    if(_min_v)
        *_min_v = min_v;
    return d1.squaredNorm();
}

float distPointTriangleSquared(const Eigen::Vector3f& _p,
                               const Eigen::Vector3f& _v0,
                               const Eigen::Vector3f& _v1,
                               const Eigen::Vector3f& _v2,
                               Eigen::Vector3f& _nearestPoint,
                               bool _stable)
{
    Eigen::Vector3f v0v1 = _v1 - _v0;
    Eigen::Vector3f v0v2 = _v2 - _v0;
    Eigen::Vector3f n = v0v1.cross(v0v2); // not normalized !
    float d = n.squaredNorm();

    // Check if the triangle is degenerated
    if (d < FLT_MIN && d > -FLT_MIN) {
        if (_stable) {
            const float l0 = v0v1.squaredNorm();
            const float l1 = v0v2.squaredNorm();
            const float l2 = (_v2 - _v1).squaredNorm();
            if (l0 > l1 && l0 > l2) {
                return distPointLineSquared(_p, _v0, _v1, &_nearestPoint);
            }
            else if (l1 > l0 && l1 > l2) {
                return distPointLineSquared(_p, _v0, _v2, &_nearestPoint);
            }
            else {
                return distPointLineSquared(_p, _v1, _v2, &_nearestPoint);
            }
        } else {
            std::cerr << "distPointTriangleSquared: Degenerated triangle !\n";
            return -1.0;
        }
    }
    float invD = 1.0f / d;

    // these are not needed for every point, should still perform
    // better with many points against one triangle
    Eigen::Vector3f v1v2 = _v2 - _v1;
    float inv_v0v2_2 = 1.0f / v0v2.squaredNorm();
    float inv_v0v1_2 = 1.0f / v0v1.squaredNorm();
    float inv_v1v2_2 = 1.0f / v1v2.squaredNorm();

    Eigen::Vector3f v0p = _p - _v0;
    Eigen::Vector3f t = v0p.cross(n);
    float s01, s02, s12;
    float a = t.dot(v0v2) * -invD;
    float b = t.dot(v0v1) * invD;

    if (a < 0)
    {
        // Calculate the distance to an edge or a corner vertex
        s02 = v0v2.dot(v0p) * inv_v0v2_2;
        if (s02 < 0.0)
        {
            s01 = v0v1.dot(v0p) * inv_v0v1_2;
            if (s01 <= 0.0) {
                v0p = _v0;
            } else if (s01 >= 1.0) {
                v0p = _v1;
            } else {
                v0p = _v0 + v0v1 * s01;
            }
        } else if (s02 > 1.0) {
            s12 = v1v2.dot(_p - _v1) * inv_v1v2_2;
            if (s12 >= 1.0) {
                v0p = _v2;
            } else if (s12 <= 0.0) {
                v0p = _v1;
            } else {
                v0p = _v1 + v1v2 * s12;
            }
        } else {
            v0p = _v0 + v0v2 * s02;
        }
    } else if (b < 0.0) {
        // Calculate the distance to an edge or a corner vertex
        s01 = v0v1.dot(v0p) * inv_v0v1_2;
        if (s01 < 0.0)
        {
            s02 = v0v2.dot(v0p ) * inv_v0v2_2;
            if (s02 <= 0.0) {
                v0p = _v0;
            } else if (s02 >= 1.0) {
                v0p = _v2;
            } else {
                v0p = _v0 + v0v2 * s02;
            }
        } else if (s01 > 1.0) {
            s12 = v1v2.dot( _p - _v1 ) * inv_v1v2_2;
            if (s12 >= 1.0) {
                v0p = _v2;
            } else if (s12 <= 0.0) {
                v0p = _v1;
            } else {
                v0p = _v1 + v1v2 * s12;
            }
        } else {
            v0p = _v0 + v0v1 * s01;
        }
    } else if (a+b > 1.0) {
        // Calculate the distance to an edge or a corner vertex
        s12 = v1v2.dot( _p - _v1 ) * inv_v1v2_2;
        if (s12 >= 1.0) {
            s02 = v0v2.dot(v0p) * inv_v0v2_2;
            if (s02 <= 0.0) {
                v0p = _v0;
            } else if (s02 >= 1.0) {
                v0p = _v2;
            } else {
                v0p = _v0 + v0v2*s02;
            }
        } else if (s12 <= 0.0) {
            s01 = v0v1.dot(v0p) * inv_v0v1_2;
            if (s01 <= 0.0) {
                v0p = _v0;
            } else if (s01 >= 1.0) {
                v0p = _v1;
            } else {
                v0p = _v0 + v0v1 * s01;
            }
        } else {
            v0p = _v1 + v1v2 * s12;
        }
    } else {
        // Calculate the distance to an interior point of the triangle
        _nearestPoint = _p - n*(n.dot(v0p) * invD);
        return (_nearestPoint - _p).squaredNorm();
    }

    _nearestPoint = v0p;

    return (_nearestPoint - _p).squaredNorm();
}


int big_angle(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
              const Eigen::Vector3f &n,float &angle)
{
    const double dot = a.dot(b);
    const double det = n.dot(a.cross(b));
    angle = std::atan2(det, dot);
    return angle <= 0;
}

int big_angle(const Eigen::Vector2f &a, const Eigen::Vector2f &b, float &angle)
{
    const double dot = a.dot(b);
    const double det = a[0]*b[1] - a[1]*b[0];
    angle = std::atan2(det, dot);
    return angle >= M_PI;
}

int solve_tri_diagnal(const size_t n,const float *l,
                      const float *d,const float *u, float *x)
{
    const float eps = 1e-6;
    Eigen::VectorXf dd = Eigen::Map<const Eigen::VectorXf>(d, n);
    for(size_t i = 1;i < n;++i){
        if(fabs(dd[i-1]) < eps)
            return 0;
        float fac = l[i] / dd[i-1];
        dd[i] -= fac * u[i-1];
        x[i] -= fac * x[i-1];
    }
    if(fabs(dd[n-1]) < eps){
        return 0;
    }
    x[n-1] /= dd[n-1];
    for(int i = n - 2; i >= 0;--i){
        if(fabs(dd[i]) < eps)
            return 0;
        x[i] = (x[i] - u[i] * x[i+1]) / dd[i];
    }
    return 1;
}

int solve_tri_diagnal_symmetric(const size_t n,const float *l,
                                float *d, float *x)
{
    const float eps = 1e-6;
    for(size_t i = 1;i < n;++i){
        if(fabs(d[i-1]) < eps)
            return 0;
        float fac = l[i] / d[i-1];
        d[i] -= fac * l[i];
        x[i] -= fac * x[i-1];
    }
    if(fabs(d[n-1]) < eps){
        return 0;
    }
    x[n-1] /= d[n-1];
    for(int i = n - 2; i >= 0;--i){
        if(fabs(d[i]) < eps)
            return 0;
        x[i] = (x[i] - l[i+1] * x[i+1]) / d[i];
    }
    return 1;
}

}
