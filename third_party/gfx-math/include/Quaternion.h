//
//  Quaternion.h
//  gfx-math
//
//  Created by Daosheng Mu on 9/26/20.
//  Copyright © 2020 Daosheng Mu. All rights reserved.
//

#ifndef GFX_MATH_QUATERNION_H
#define GFX_MATH_QUATERNION_H

#include "Vector3d.h"

namespace gfx_math {

template <class Type>
class Quaternion {
public:
  Quaternion() {
    Identity();
  }
  
  Quaternion(Type aX, Type aY, Type aZ, Type aW) :
    x(aX), y(aY), z(aZ), w(aW) {
    
  }
  
  void Identity() {
    x = y = z = 0.0f;
    w = 1.0f;
  }
  
  Type Length() const {
    return sqrt(x * x + y * y + z * z + w * w);
  }
  
  void Normalize() {
    Type l = Length();
    if (l) {
      l = 1.0f / l;
      x *= l;
      y *= l;
      z *= l;
      w *= l;
    } else {
      x = y = z = 0.f;
      w = 1.f;
    }
  }
  
  void Conjugate() {
    x *= -1.f;
    y *= -1.f;
    z *= -1.f;
  }
  
  void Inverse() {
    Type dot = x * x + y * y + z * z + w * w;
    Type invDot = dot ? (1.0 / dot) : 0;
    
    // Inverse = Conjugate() / (Length() ^ 2)
    x = -x / invDot;
    y = -y / invDot;
    z = -z / invDot;
    w = w / invDot;
  }
  
  Quaternion<Type>& operator *= (const Quaternion<Type>& aInput) {
    x = x * aInput.x + w * aInput.x + y * aInput.z - z * aInput.y;
    y = y * aInput.w + w * aInput.y + z * aInput.x - x * aInput.z;
    z = z * aInput.w + w * aInput.z + x * aInput.y - y * aInput.x;
    w = w * aInput.w - x * aInput.x - y * aInput.y - z * aInput.z;
    
    return *this;
  }
  
  Type DotProduct(const Quaternion<Type>& aInput) {
    return x * aInput.x + y * aInput.y + z * aInput.z;
  }
  
  Type GetAngle(const Quaternion<Type>& aInput) {
    Type dot = DotProduct(aInput);

    return acos(2 * dot * dot - 1);
  }
  
  Quaternion<Type> operator * (const Quaternion<Type>& aInput) {
    Quaternion<Type> output;
    
    output[0] = x * aInput.x + w * aInput.x + y * aInput.z - z * aInput.y;
    output[1] = y * aInput.w + w * aInput.y + z * aInput.x - x * aInput.z;
    output[2] = z * aInput.w + w * aInput.z + x * aInput.y - y * aInput.x;
    output[3] = w * aInput.w - x * aInput.x - y * aInput.y - z * aInput.z;
    
    return output;
  }
  
  void RotateX(const Quaternion<Type>& aInput, Type aRadian) {
    aRadian *= 0.5;
    Type bx = sin(aRadian),
    bw = cos(aRadian);

    x = x * bw + w * bx;
    y = y * bw + z * bx;
    z = z * bw - y * bx;
    w = w * bw - x * bx;
  }
  
  void RotateY(const Quaternion<Type>& aInput, Type aRadian) {
    Type rad = aRadian * 0.5;
    Type by = sin(rad), bw = cos(rad);

    x = x * bw - z * by;
    y = y * bw + w * by;
    z = z * bw + x * by;
    w = w * bw - y * by;
  }
  
  void RotateZ(const Quaternion<Type>& aInput, Type aRadian) {
    Type rad = aRadian * 0.5;
    Type bz = sin(rad), bw = cos(rad);

    x = x * bw + y * bz;
    y = y * bw - x * bz;
    z = z * bw + w * bz;
    z = w * bw - z * bz;
  }
      
  Quaternion<Type> Slerp(const Quaternion<Type>& a, const Quaternion<Type>& b, Type t) {
    // benchmarks:
    //    http://jsperf.com/quaternion-slerp-implementations
    Type ax = a.x,
      ay = a.y,
      az = a.z,
      aw = a.w;
    Type bx = b.x,
      by = b.y,
      bz = b.z,
      bw = b.w;

    Type omega, cosom, sinom, scale0, scale1;

    // calc cosine
    cosom = ax * bx + ay * by + az * bz + aw * bw;
    // adjust signs (if necessary)
    if (cosom < 0.0) {
      cosom = -cosom;
      bx = -bx;
      by = -by;
      bz = -bz;
      bw = -bw;
    }
    // calculate coefficients
    if (1.0 - cosom > FLT_EPSILON) {
      // standard case (slerp)
      omega = acos(cosom);
      sinom = sin(omega);
      scale0 = sin((1.0 - t) * omega) / sinom;
      scale1 = sin(t * omega) / sinom;
    } else {
      // "from" and "to" quaternions are very close
      //  ... so we can do a linear interpolation
      scale0 = 1.0 - t;
      scale1 = t;
    }
    
    Quaternion<Type> output;
    
    // calculate final values
    output.x = scale0 * ax + scale1 * bx;
    output.y = scale0 * ay + scale1 * by;
    output.z = scale0 * az + scale1 * bz;
    output.w = scale0 * aw + scale1 * bw;

    return output;
  }
  
  static Quaternion<Type> SetByAxisAngle(const Vector3D<Type>& aAxis, Type aRad) {
    Type rad = aRad * 0.5;
    Type s = sin(rad);
    
    Quaternion<Type> output;
    
    output.x = s * aAxis.x;
    output.y = s * aAxis.y;
    output.z = s * aAxis.z;
    output.w = cos(rad);
    return output;
  }
  
  Type x, y, z, w;
};

typedef Quaternion<float> Quaternionf;

} // End of namespace gfx_math


#endif /* GFX_MATH_QUATERNION_H */
