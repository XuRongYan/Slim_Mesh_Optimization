#ifndef COMMON_UTIL_H
#define COMMON_UTIL_H

#include <Eigen/Dense>

#include "conf.h"

namespace common{
static const float eps = 1.1e-6;

int COMMON_API barycentric(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                           const Eigen::Vector3f &c, const Eigen::Vector3f &p,
                           Eigen::Vector3f &bary);

int COMMON_API line_triangle_intersect(const Eigen::Vector3f &o, const Eigen::Vector3f &d,
                                       const Eigen::Vector3f &p0, const Eigen::Vector3f &p1,
                                       const Eigen::Vector3f &p2, Eigen::Vector3f &bary, Eigen::Vector3f &pt);

int COMMON_API line_plane_intersect(const Eigen::Vector3f &p, const Eigen::Vector3f &n,
                                    const Eigen::Vector3f &p0, const Eigen::Vector3f &p1,
                                    Eigen::Vector3f &q);

Eigen::Vector3f COMMON_API closestPointOnTriangle(const Eigen::Vector3f &p0,
                                                  const Eigen::Vector3f &p1,
                                                  const Eigen::Vector3f &p2,
                                                  const Eigen::Vector3f &srcPos,
                                                  Eigen::Vector3f &bary);

float COMMON_API distPointTriangleSquared(const Eigen::Vector3f& _p,
                                          const Eigen::Vector3f& _v0,
                                          const Eigen::Vector3f& _v1,
                                          const Eigen::Vector3f& _v2,
                                          Eigen::Vector3f& _nearestPoint,
                                          bool _stable = false);

int COMMON_API  big_angle(const Eigen::Vector2f &a, const Eigen::Vector2f &b, float &angle);
int COMMON_API  big_angle(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                          const Eigen::Vector3f &n, float &angle);

//size(l) = size(d) = size(u), l[0] = u[n-1] = 0
int COMMON_API solve_tri_diagnal(const size_t n,const float *l,
                                 const float *d,const float *u,float *x);

int COMMON_API solve_tri_diagnal_symmetric(const size_t n,const float *l,
                                           float *d, float *x);
}

#endif
