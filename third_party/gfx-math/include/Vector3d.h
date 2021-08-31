//
//  vector3d.h
//  gfx-math
//
//  Created by Daosheng Mu on 9/26/20.
//  Copyright © 2020 Daosheng Mu. All rights reserved.
//

#ifndef GFX_MATH_VECTOR3D_H
#define GFX_MATH_VECTOR3D_H

#include <math.h>
#include <cfloat>
#include <vector>

namespace gfx_math {

template <class Type>
class Vector3D {

public:
  Vector3D() = default;
  
  Vector3D(Type aX, Type aY, Type aZ) :
    x(aX), y(aY), z(aZ) {
    
  }
  
  void Set(Type aX, Type aY, Type aZ) {
    x = aX; y = aY; z = aZ;
  }
  
  Type DotProduct(const Vector3D<Type>& aInput) const {
    return x * aInput.x + y * aInput.y + z * aInput.z;
  }

  Vector3D<Type> CrossProduct(const Vector3D<Type>& aInput) const {
    return Vector3D(y * aInput.z - aInput.y * z,
                    z * aInput.x - aInput.z * x,
                    x * aInput.y - aInput.x * y);
  }
  
  bool operator == (const Vector3D<Type> & aIn) const {
    return (x == aIn.x && y == aIn.y && z == aIn.z);
  }
  
  bool operator != (const Vector3D<Type> & aRhs) const {
    return !(*this == aRhs);
  }
  
  Vector3D<Type>& operator /= (const Type& aScale) {
    x /= aScale;
    y /= aScale;
    z /= aScale;
    return Vector3D<Type>(this);
  }
  
  Vector3D<Type>& operator *= (const Type& aScale) {
    x *= aScale;
    y *= aScale;
    z *= aScale;
    return *this;
  }
  
  friend Vector3D<Type> operator * (const Type& aScale, const Vector3D<Type> &aV) {
    return Vector3D<Type>(aV.x * aScale, aV.y * aScale, aV.z * aScale);
  }
  
  friend Vector3D<Type> operator / (const Type& aDividend, const Vector3D<Type> &aV) {
    return Vector3D<Type>(aDividend / aV.x, aDividend / aV.y, aDividend / aV.z);
  }
  
  Vector3D<Type> operator * (const Type& aScale) const {
    return Vector3D<Type>(x * aScale, y * aScale, z * aScale);
  }
  
  Vector3D<Type> operator + (const Vector3D<Type>& aInput) const {
    return Vector3D<Type>(x + aInput.x, y + aInput.y, z + aInput.z);
  }
  
  Vector3D<Type> operator - (const Vector3D<Type>& aInput) const {
    return Vector3D<Type>(x - aInput.x, y - aInput.y, z - aInput.z);
  }
  
  void Normalize() {
    Type len = Length();
    
    if (len != 0) {
      x /= len;
      y /= len;
      z /= len;
    } else {
      x = y = z = 0;
    }
  }
  
  Type Length() const {
    Type sq = x * x + y * y + z * z;
    if (sq == 0) {
      return 0;
    }
    return sqrt(sq);
  }

  Type x, y, z;
};

// Good references, https://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld018.htm

/**
  We according to the plane formula aX + bY + cZ + d = 0 and p = p0 + dir * t
  For looking for a point from direction r on the plane, we can suppose Dot(N, p) + d = 0.
  and Dot(N, p0 + dir * t) + d = 0. If t < 0, it means the p is behind p0.
 */
// Return: false: ray no interset
//         true: ray intersect with plane.
template <class Type>
bool RayIntersectWithPlane(const Vector3D<Type>& aDir,
                           const Vector3D<Type>& aP0,
                           const Vector3D<Type>& aPlaneNrm,
                           Type aPlaneDist,
                           Vector3D<Type>& aInterPoint) {
  Type denom = aPlaneNrm.DotProduct(aDir);
  if (fabs(denom) <= FLT_EPSILON) { // prevent 0, perpendicular, no intersect.
    return false;
  }
  
  Type t = -(aPlaneDist + aPlaneNrm.DotProduct(aP0)) / denom;
  if (t < 0) {
    return false;
  }

  aInterPoint = aP0 + aDir * t;
  return true;
}

/**
  We use the formula of b^2 - 4ac to determine if there is any intersect point in a sphere.
  p = p0 + dir * t.  Dist(p - origin) = r => (p - origin) ^ 2 = r ^ 2.
  Based on quadratic function, at^2 + bt + c  = 0.
  a = 1, b = 2 * dir * (p0 - origin), c = (p0 - origin)^2 - radius ^ 2.
  x = (-b +- sqrt(b^2 - 4ac)) / 2a, discr = b ^ 2 - 4ac
  if discr == 0, it means we only have one intersection, the p = -b / 2a.
  Two intersects, (−b +- sqrt(discr)) / 2a
  https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
 */

template <class Type>
bool RayIntersectWithSphere(const Vector3D<Type>& aDir,
                            const Vector3D<Type>& aP0,
                            const Vector3D<Type>& aSphereOrigin,
                            Type aRadius,
                            std::vector<Vector3D<Type>>& aInterPoint) {
  Vector3D<Type> L = aP0 - aSphereOrigin;
  Type a = aDir.DotProduct(aDir);
  Type b = 2 * aDir.DotProduct(L);
  Type c = L.DotProduct(L) - aRadius * aRadius;
  
  assert(aInterPoint.size() == 0);
  
  Type discr = b * b - 4 * a * c;
  if (discr < 0) {
    return false;
  } else if (discr == 0) {  // −b/2a, the only intersect point.
    Type t0;
    t0 = -0.5 * b / a; // ray on the sphere and dist is radius.
    Vector3D<Type> p1 = aP0 + aDir * t0;
    aInterPoint.push_back(p1);
    return true;
  }
   
  // Two intersects, (−b +- sqrt(discr))/2a
  Type t0, t1;
  Type q = (b > 0) ?   // To fix catastrophic cancellation or round off error.
      -0.5 * (b + sqrt(discr)) :
      -0.5 * (b - sqrt(discr));
  t0 = q / a;
  t1 = c / q;

  if (t0 > t1) {
    std::swap(t0, t1);
  }
  
  if (t0 >= 0) {
    Vector3D<Type> p1 = aP0 + aDir * t0;
    aInterPoint.push_back(p1);
  }

  Vector3D<Type> p2 = aP0 + aDir * t1;
  aInterPoint.push_back(p2);
  return true;
}

/**
  First of all, getting the normal from the three points of a triangle.
  Finding the intersect point p by using ax + by + cz + d = 0 and p = p0 + dir * t.
  Getting d by using -d = dot(N, v0) because v0 is on the plane
  Then, finding a point on the plane (x, y, z) = p. t = -(dot(N, P0) + d) / dot(N, dir)
  So, p = p0 + (dir * t).
  Finally, checking if p in the boundary of triangle v0, v1, v2 by edge0 = v1- v0, edge1 = p- v0
  Cross product edge0 and edge1 => dir , if dot(dir, N) > 0, it means they point to the same direction,
  so p is at the left side of edge0, Continue to try other two edges of this triangle.
    
 https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
 */

template <class Type>
bool RayIntersectWithTriangle(const Vector3D<Type>& aDir, const Vector3D<Type>& aP0,
                              const Vector3D<Type>& v0, const Vector3D<Type>& v1, const Vector3D<Type>& v2,
                              Vector3D<Type>& aInterPoint) {
  // compute plane's normal
  Vector3D<Type> v01 = v1 - v0;
  Vector3D<Type> v02 = v2 - v0;
  // no need to normalize
  Vector3D<Type> N = v01.CrossProduct(v02);

  // Step 1: finding P

  // check if ray and plane are parallel?
  Type NdotRayDirection = N.DotProduct(aDir);
  // TODO: if (det < kEpsilon) return false; // backface culling.
  if (fabs(NdotRayDirection) < FLT_EPSILON) // almost 0
     return false; // they are parallel so they don't intersect!

  // compute d parameter using equation 2
  float d = -N.DotProduct(v0);
  // compute t (equation 3)
  Type t = -(N.DotProduct(aP0) + d) / NdotRayDirection;
  // check if the triangle is in behind the ray
  if (t < 0) {
    return false; // the triangle is behind
  }

  // compute the intersection point using equation 1
  aInterPoint = aP0 + (aDir * t);

  // Step 2: inside-outside test
  Vector3D<Type> C; // vector perpendicular to triangle's plane

  // edge 0
  Vector3D<Type> edge0 = v1 - v0;
  Vector3D<Type> vp0 = aInterPoint - v0;
  C = edge0.CrossProduct(vp0);
  if (N.DotProduct(C) < 0) {
    return false; // P is on the right side
  }

  // edge 1
  Vector3D<Type> edge1 = v2 - v1;
  Vector3D<Type> vp1 = aInterPoint - v1;
  C = edge1.CrossProduct(vp1);
  if (N.DotProduct(C) < 0) {
    return false; // P is on the right side
  }

  // edge 2
  Vector3D<Type> edge2 = v0 - v2;
  Vector3D<Type> vp2 = aInterPoint - v2;
  C = edge2.CrossProduct(vp2);
  if (N.DotProduct(C) < 0) {
    return false; // P is on the right side;
  }

  return true;
}

/**
 Try to find the intersect point on AABB from ray. p = p0 + dir * t
 According to the direction of the ray to get t0.x = (bb0.x - orig.x) / dir.x
 and see which would be tmin and tmax, Be sure the dir is normalized.
 If t0x > t1y || t0y > t1x, it means no intersection. Then, continue to check t0z, t1z
 Finally, we would try to get the first intersect point by checking tmin or tmax depend on if it is > 0.
 */

// Ray-Box
//
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
template <class Type>
bool RayIntersectWithAABBox(const Vector3D<Type>& aDir, const Vector3D<Type>& aP0,
                            const Vector3D<Type>& aBoundMin, const Vector3D<Type>& aBoundMax,
                            Vector3D<Type>& aInterPoint) {
  Type tmin, tmax, tymin, tymax, tzmin, tzmax;
  
//  if (r.dir.x >= 0) {
//      tmin = (min.x - r.orig.x) / r.dir.x;
//      tmax = (max.x - r.orig.x) / r.dir.x;
//  }
//  else {
//      tmin = (max.x - r.orig.x) / r.dir.x;
//      tmax = (min.x - r.orig.x) / r.dir.x;
//  }
  
  // Using invdir in order to fix the direction of ray is 0 it causes a division by zero.
  Vector3D<Type> invdir = 1 / aDir;
  Vector3D<Type> bounds[2] = {aBoundMin, aBoundMax};
  bool signX = invdir.x < 0, signY = invdir.y < 0, signZ = invdir.z < 0;
  
  tmin = (bounds[signX].x - aP0.x) * invdir.x;
  tmax = (bounds[1 - signX].x - aP0.x) * invdir.x;
  tymin = (bounds[signY].y - aP0.y) * invdir.y;
  tymax = (bounds[1 - signY].y - aP0.y) * invdir.y;
  
  if ((tmin > tymax) || (tymin > tmax)) {
    return false;
  }
  
  if (tymin > tmin) {
    tmin = tymin;
  }
  if (tymax < tmax) {
    tmax = tymax;
  }
  
  tzmin = (bounds[signZ].z - aP0.z) * invdir.z;
  tzmax = (bounds[1 - signZ].z - aP0.z) * invdir.z;

  if ((tmin > tzmax) || (tzmin > tmax)) {
    return false;
  }

  if (tzmin > tmin) {
    tmin = tzmin;
  }
  if (tzmax < tmax) {
    tmax = tzmax;
  }

  Type t = tmin;
  if (t < 0) {
     t = tmax;

    if (t < 0) {
      return false;
    }
  }
  
  aInterPoint = aP0 + (aDir * t);
  return true;
}

/**
  Using ax+by+cy+d = 0 to see if the point is on the plane
  If result = 0, it means on the plane, > 0 means in front of,
  < 0 means behind the plane.
 */

// Return: -1: back
//         0: on plane
//         1: front
template <class Type>
int PointIntersectWithPlane(const Vector3D<Type>& aP0,
                            const Vector3D<Type>& aPlaneNrm,
                            Type aPlaneDist) {
  Type denom = aPlaneNrm.DotProduct(aP0);
  Type distance = denom + aPlaneDist;
  
  if (fabs(distance) <= FLT_EPSILON) {
    return 0; // on the plane
  } else if (distance > 0) {
    return 1;
  } else {  // result < 0
    return -1;
  }
}

template <class Type>
Type PointDistanceWithPlane(const Vector3D<Type>& aP0,
                           const Vector3D<Type>& aPlaneNrm,
                           Type aPlaneDist) {
  Type denom = aPlaneNrm.DotProduct(aP0);
  Type len = aPlaneNrm.Length();
  if (len == 0) {
    return 0;
  }
  
  return fabs(denom + aPlaneDist) / len;
}

/**
  Get the closest point (x, y, z) from boundMax, boundMin and the sphere origin.
  if the distance of origin and (x, y, z) <= radius, it means there are intersections.
 */
template <class Type>
bool SphereIntersectWithAABBox(const Vector3D<Type>& aOrigin, Type aRadius,
                               const Vector3D<Type>& aBoundMin,
                               const Vector3D<Type>& aBoundMax) {
  Type radius2 = aRadius * aRadius;
  
  // Get box closest point from the sphere center.
  Type x = std::max(aBoundMin.x, std::min(aOrigin.x, aBoundMax.x));
  Type y = std::max(aBoundMin.y, std::min(aOrigin.y, aBoundMax.y));
  Type z = std::max(aBoundMin.z, std::min(aOrigin.z, aBoundMax.z));
  
  Type distance = (x - aOrigin.x) * (x - aOrigin.x) +
                  (y - aOrigin.y) * (y - aOrigin.y) +
                  (z - aOrigin.z) * (z - aOrigin.z);
  return distance <= radius2;
}

/**
 Compare if the bounding box A max in the middle between bounding box B min and max.
*/
template <class Type>
bool AABBoxIntersectWithAABBox(const Vector3D<Type>& aBoundAMin,
                               const Vector3D<Type>& aBoundAMax,
                               const Vector3D<Type>& aBoundBMin,
                               const Vector3D<Type>& aBoundBMax) {
  return (aBoundAMin.x <= aBoundBMax.x && aBoundAMax.x >= aBoundBMin.x) &&
          (aBoundAMin.y <= aBoundBMax.y && aBoundAMax.y >= aBoundBMin.y) &&
          (aBoundAMin.z <= aBoundBMax.z && aBoundAMax.z >= aBoundBMin.z);
}

/**
 Because the reflection direction is mirror to the incoming direction.  We can make it L - N1 = -(N - N2).
 Besides, N1 and N2 are the projected vectors from L and R to N, and N1 is the same with N2.
 Θr = Θl, so dot(R^, N^) == dot(L^, N^)  and R^ - N1^  = -(L^ - N2^),
 N1^ = (dot(R^, N^) / (|R^| |N^|) ) * N^
 N2^ = (dot(L^, N^) / (|L^| |N^|) ) * N^,  Because L^, R^, and N^ are unit vectors, their length are one.
 Hence,  R^ - dot(L^, N^) * N^ = -L^ + dot(L^, N^) * N^,  R^ = 2 * dot(L^, N^) * N^ - L.
 https://www.fabrizioduroni.it/2017/08/25/how-to-calculate-reflection-vector.html
*/
template <class Type>
Vector3D<Type> ReflectionVector(const Vector3D<Type>& aDir, const Vector3D<Type>& aNormal) {
  // 2 * dot(N^, I) * N^ - I;
  Vector3D<Type> normal(aNormal);
  normal.Normalize();
  Vector3D<Type> reflection = 2 * normal.DotProduct(aDir) * normal - aDir;
  return reflection;
}

typedef Vector3D<float> Vector3Df;

} // End of namespace gfx_math


#endif /* GFX_MATH_VECTOR3D_H */
