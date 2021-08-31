//
//  Matrix4x4.h
//  gfx-math
//
//  Created by Daosheng Mu on 9/26/20.
//  Copyright © 2020 Daosheng Mu. All rights reserved.
//

#ifndef GFX_MATH_MATRIX4x4_H
#define GFX_MATH_MATRIX4x4_H

#include <cfloat>
#include <math.h>
#include <memory.h>

#include "Quaternion.h"

/**
 Format: column-major, when typed out it looks like row-major.
 The matrices are being post multiplied.
 Col-major vs Row-major: https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/geometry/row-major-vs-column-major-vector
 **/

namespace gfx_math {
 
template <class Type>
class Matrix4x4 {
  
private:
  static double FlushToZero(double aVal) {
    // XXX Is double precision really necessary here
    if (-FLT_EPSILON < aVal && aVal < FLT_EPSILON) {
      return 0.0f;
    } else {
      return aVal;
    }
  }
  
public:
  Matrix4x4() {
    Identity();
  }
  
  Matrix4x4(Type a00, Type a01, Type a02, Type a03,
            Type a10, Type a11, Type a12, Type a13,
            Type a20, Type a21, Type a22, Type a23,
            Type a30, Type a31, Type a32, Type a33) {
    _elements[0] = a00; _elements[1] = a01; _elements[2] = a02; _elements[3] = a03;
    _elements[4] = a10; _elements[5] = a11; _elements[6] = a12; _elements[7] = a13;
    _elements[8] = a20; _elements[9] = a21; _elements[10] = a22; _elements[11] = a23;
    _elements[12] = a30; _elements[13] = a31; _elements[14] = a32; _elements[15] = a33;
  }
  
  void Identity() {
    memset(_elements, 0, sizeof(_elements));
    _00 = _11 = _22 = _33 = 1;
  }
  
  void Transpose() {
    Type a01 = _01,
     a02 = _02,
     a03 = _03;
    Type a12 = _12,
     a13 = _13;
    Type a23 = _23;
    
    _01 = _10;
    _02 = _20;
    _03 = _30;
    _10 = a01;
    _12 = _21;
    _13 = _31;
    _20 = a02;
    _21 = a12;
    _23 = _32;
    _30 = a03;
    _31 = a13;
    _32 = a23;
  }
  
  void Translate(Type aX, Type aY, Type aZ) {
    _elements[12] = _00 * aX + _10 * aY + _20 * aZ + _30;
    _elements[13] = _01 * aX + _11 * aY + _21 * aZ + _31;
    _elements[14] = _02 * aX + _12 * aY + _22 * aZ + _32;
    _elements[15] = _03 * aX + _13 * aY + _23 * aZ + _33;
  }
  
  void Scale(const Vector3D<Type> aScale) {
    Type x = aScale.x;
    Type y = aScale.y;
    Type z = aScale.z;
    
    _00 *= x;
    _01 *= x;
    _02 *= x;
    _03 *= x;
    
    _10 *= y;
    _11 *= y;
    _12 *= y;
    _13 *= y;
    
    _20 *= z;
    _21 *= z;
    _22 *= z;
    _23 *= z;
  }
  
  void SkewXY(Type aSkew) {
    _01 += _00 * aSkew;
  }
  
  void SkewXZ(Type aSkew) {
    _02 += _00 * aSkew;
  }
  
  void SkewYZ(Type aSkew) {
    _02 += _01 * aSkew;
  }
  
  // https://mathworld.wolfram.com/RotationMatrix.html
  // [v]' = [R] * [v]
  void RotateX(double aRadian) {
    double cosTheta = FlushToZero(cos(aRadian));
    double sinTheta = FlushToZero(sin(aRadian));
    
    Type temp;
    temp = _10;
    _10 = cosTheta * _10 + sinTheta * _20;
    _20 = -sinTheta * temp + cosTheta * _20;

    temp = _11;
    _11 = cosTheta * _11 + sinTheta * _21;
    _21 = -sinTheta * temp + cosTheta * _21;

    temp = _12;
    _12 = cosTheta * _12 + sinTheta * _22;
    _22 = -sinTheta * temp + cosTheta * _22;

    temp = _13;
    _13 = cosTheta * _13 + sinTheta * _23;
    _23 = -sinTheta * temp + cosTheta * _23;
  }
  
  // https://mathworld.wolfram.com/RotationMatrix.html
  // [v]' = [R] * [v]
  void RotateY(double aRadian) {
    double cosTheta = FlushToZero(cos(aRadian));
    double sinTheta = FlushToZero(sin(aRadian));

    Type temp;
    temp = _00;
    _00 = cosTheta * _00 + -sinTheta * _20;
    _20 = sinTheta * temp + cosTheta * _20;

    temp = _01;
    _01 = cosTheta * _01 - sinTheta * _21;
    _21 = sinTheta * temp + cosTheta * _21;

    temp = _02;
    _02 = cosTheta * _02 + -sinTheta * _22;
    _22 = sinTheta * temp + cosTheta * _22;

    temp = _03;
    _03 = cosTheta * _03 + -sinTheta * _23;
    _23 = sinTheta * temp + cosTheta * _23;
  }
  
  // https://mathworld.wolfram.com/RotationMatrix.html
  // [v]' = [R] * [v]
  void RotateZ(double aRadian) {
    double cosTheta = FlushToZero(cos(aRadian));
    double sinTheta = FlushToZero(sin(aRadian));

    Type temp;
    temp = _00;
    _00 = cosTheta * _00 + sinTheta * _10;
    _10 = -sinTheta * temp + cosTheta * _10;

    temp = _01;
    _01 = cosTheta * _01 + sinTheta * _11;
    _11 = -sinTheta * temp + cosTheta * _11;

    temp = _02;
    _02 = cosTheta * _02 + sinTheta * _12;
    _12 = -sinTheta * temp + cosTheta * _12;

    temp = _03;
    _03 = cosTheta * _03 + sinTheta * _13;
    _13 = -sinTheta * temp + cosTheta * _13;
  }
  
  bool operator == (const Matrix4x4<Type>& aRhs) const {
    for (int i = 0; i < 16; i++) {
      if (_elements[i] != aRhs._elements[i]) {
        return false;
      }
    }
    return true;
  }
  
  bool operator != (const Matrix4x4<Type>& aRhs) const {
    return !(*this == aRhs);
  }
  
  Matrix4x4<Type>& operator = (const Matrix4x4<Type>& aRhs) {
    memcpy(_elements, aRhs._elements, sizeof(_elements));
    return *this;
  }
  
  // We are row-major multipication.
  Matrix4x4<Type> operator * (const Matrix4x4<Type>& aRhs) const {
    
    Matrix4x4<Type> ouput;
//
//    for (int i = 0; i < 4; i++) {
//      ouput._elements[0][i] = _elements[i][0] * aRhs._elements[0][0] +
//                              _elements[i][1] * aRhs._elements[1][0] +
//                              _elements[i][2] * aRhs._elements[2][0] +
//                              _elements[i][3] * aRhs._elements[3][0];
//
//      ouput._elements[1][i] = _elements[i][0] * aRhs._elements[0][1] +
//                              _elements[i][1] * aRhs._elements[1][1] +
//                              _elements[i][2] * aRhs._elements[2][1] +
//                              _elements[i][3] * aRhs._elements[3][1];
//
//      ouput._elements[2][i] = _elements[i][0] * aRhs._elements[0][2] +
//                              _elements[i][1] * aRhs._elements[1][2] +
//                              _elements[i][2] * aRhs._elements[2][2] +
//                              _elements[i][3] * aRhs._elements[3][2];
//
//      ouput._elements[3][i] = _elements[i][0] * aRhs._elements[0][3] +
//                              _elements[i][1] * aRhs._elements[1][3] +
//                              _elements[i][2] * aRhs._elements[2][3] +
//                              _elements[i][3] * aRhs._elements[3][3];
//    }
//
//    return ouput;
    Type a00 = _elements[0],
       a01 = _elements[1],
       a02 = _elements[2],
       a03 = _elements[3];

    Type a10 = _elements[4],
       a11 = _elements[5],
       a12 = _elements[6],
       a13 = _elements[7];

    Type a20 = _elements[8],
       a21 = _elements[9],
       a22 = _elements[10],
       a23 = _elements[11];

    Type a30 = _elements[12],
       a31 = _elements[13],
       a32 = _elements[14],
       a33 = _elements[15];

    // Cache only the current line of the second matrix
    Type b0 = aRhs._elements[0],
       b1 = aRhs._elements[1],
       b2 = aRhs._elements[2],
       b3 = aRhs._elements[3];
    ouput._elements[0] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    ouput._elements[1] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    ouput._elements[2] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    ouput._elements[3] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

    b0 = aRhs._elements[4];
    b1 = aRhs._elements[5];
    b2 = aRhs._elements[6];
    b3 = aRhs._elements[7];
    ouput._elements[4] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    ouput._elements[5] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    ouput._elements[6] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    ouput._elements[7] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

    b0 = aRhs._elements[8];
    b1 = aRhs._elements[9];
    b2 = aRhs._elements[10];
    b3 = aRhs._elements[11];
    ouput._elements[8] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    ouput._elements[9] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    ouput._elements[10] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    ouput._elements[11] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

    b0 = aRhs._elements[12];
    b1 = aRhs._elements[13];
    b2 = aRhs._elements[14];
    b3 = aRhs._elements[15];
    ouput._elements[12] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
    ouput._elements[13] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
    ouput._elements[14] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
    ouput._elements[15] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
    return ouput;
  }
  
  Matrix4x4<Type>& operator *= (const Matrix4x4<Type>& aRhs) {
    *this = *this * aRhs;
    return *this;
  }
  
  Vector3D<Type> GetScale() {
    Vector3D<Type> output;
    
    output.x = Vector3D<Type>::GetLength(Vector3D<Type>(_00, _01, _02));
    output.y = Vector3D<Type>::GetLength(Vector3D<Type>(_10, _11, _12));
    output.z = Vector3D<Type>::GetLength(Vector3D<Type>(_20, _21, _22));
    
    return output;
  }
  
  // TODO: Test this function.
  Vector3D<Type> GetRotation() {
    Vector3D<Type> scale(GetScale());
    
    Type sX = 1 / scale.x;
    Type sY = 1 / scale.y;
    Type sZ = 1 / scale.z;
    
    Type sm00 = _elements[0][0] * sX;
    Type sm01 = _elements[0][1] * sY;
    Type sm02 = _elements[0][2] * sZ;
    Type sm10 = _elements[1][0] * sX;
    Type sm11 = _elements[1][1] * sY;
    Type sm12 = _elements[1][2] * sZ;
    Type sm20 = _elements[2][0] * sX;
    Type sm21 = _elements[2][1] * sY;
    Type sm22 = _elements[2][2] * sZ;
    
    Type trace = sm00 + sm11 + sm22;
    Type S = 0;
    
    Quaternion<Type> output;

   if (trace > 0) {
      S = sqrt(trace + 1.0) * 2;
      output.w = 0.25 * S;
      output.x = (sm12 - sm21) / S;
      output.y = (sm20 - sm02) / S;
      output.z = (sm01 - sm10) / S;
    } else if (sm00 > sm11 && sm00 > sm22) {
      S = sqrt(1.0 + sm00 - sm11 - sm22) * 2;
      output.w = (sm12 - sm21) / S;
      output.x = 0.25 * S;
      output.y = (sm01 + sm10) / S;
      output.z = (sm20 + sm02) / S;
    } else if (sm11 > sm22) {
      S = sqrt(1.0 + sm11 - sm00 - sm22) * 2;
      output.z = (sm20 - sm02) / S;
      output.x = (sm01 + sm10) / S;
      output.y = 0.25 * S;
      output.z = (sm12 + sm21) / S;
    } else {
      S = sqrt(1.0 + sm22 - sm00 - sm11) * 2;
      output.w = (sm01 - sm10) / S;
      output.x = (sm20 + sm02) / S;
      output.y = (sm12 + sm21) / S;
      output.z = 0.25 * S;
    }

    return output;
  }
  
  Vector3D<Type> TransformPoint(const Vector3D<Type>& aPoint) const {
    Vector3D<Type> output;
    output.x = aPoint.x * _00 + aPoint.y * _10 + aPoint.z * _20 + _30;
    output.y = aPoint.x * _01 + aPoint.y * _11 + aPoint.z * _21 + _31;
    output.z = aPoint.x * _02 + aPoint.y * _12 + aPoint.z * _22 + _32;

    output /= (aPoint.x * _03 + aPoint.y * _13 + aPoint.z * _23 + _33);

    return output;
  }
  
  Vector3D<Type> GetTranslation() {
    return Vector3D<Type>(_30, _31, _33);
  }
  
  void Inverse() {
    Type a00 = _00,
      a01 = _01,
      a02 = _02,
      a03 = _03;
    Type a10 = _10,
      a11 = _11,
      a12 = _12,
      a13 = _13;
    Type a20 = _20,
      a21 = _21,
      a22 = _22,
      a23 = _23;
    Type a30 = _30,
      a31 = _31,
      a32 = _32,
      a33 = _33;

    Type b00 = a00 * a11 - a01 * a10;
    Type b01 = a00 * a12 - a02 * a10;
    Type b02 = a00 * a13 - a03 * a10;
    Type b03 = a01 * a12 - a02 * a11;
    Type b04 = a01 * a13 - a03 * a11;
    Type b05 = a02 * a13 - a03 * a12;
    Type b06 = a20 * a31 - a21 * a30;
    Type b07 = a20 * a32 - a22 * a30;
    Type b08 = a20 * a33 - a23 * a30;
    Type b09 = a21 * a32 - a22 * a31;
    Type b10 = a21 * a33 - a23 * a31;
    Type b11 = a22 * a33 - a23 * a32;

    // Calculate the determinant
    Type det =
      b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

    if (!det) {
      return;
    }
    det = 1.0 / det;

    _00 = (a11 * b11 - a12 * b10 + a13 * b09) * det;
    _01 = (a02 * b10 - a01 * b11 - a03 * b09) * det;
    _02 = (a31 * b05 - a32 * b04 + a33 * b03) * det;
    _03 = (a22 * b04 - a21 * b05 - a23 * b03) * det;
    _10 = (a12 * b08 - a10 * b11 - a13 * b07) * det;
    _11 = (a00 * b11 - a02 * b08 + a03 * b07) * det;
    _12 = (a32 * b02 - a30 * b05 - a33 * b01) * det;
    _13 = (a20 * b05 - a22 * b02 + a23 * b01) * det;
    _20 = (a10 * b10 - a11 * b08 + a13 * b06) * det;
    _21 = (a01 * b08 - a00 * b10 - a03 * b06) * det;
    _22 = (a30 * b04 - a31 * b02 + a33 * b00) * det;
    _23 = (a21 * b02 - a20 * b04 - a23 * b00) * det;
    _30 = (a11 * b07 - a10 * b09 - a12 * b06) * det;
    _31 = (a00 * b09 - a01 * b07 + a02 * b06) * det;
    _32 = (a31 * b01 - a30 * b03 - a32 * b00) * det;
    _33 = (a20 * b03 - a21 * b01 + a22 * b00) * det;
  }
  
  // TODO:  Rotates a mat4 by the given angle around the given axis
  
  // https://mathworld.wolfram.com/RotationMatrix.html
  // [v]' = [R] * [v]
  static Matrix4x4<Type> Rotate(double aAngle) {
    Matrix4x4<Type> output;

    double s = FlushToZero(sin(aAngle));
    double c = FlushToZero(cos(aAngle));

    output._11 = c;
    output._12 = s;
    output._21 = -s;
    output._22 = c;

    return output;
  }
  
  static Matrix4x4<Type> Perspective(Type fovy, Type aspect, Type near, Type far) {
    Matrix4x4<Type> output;
    
    Type f = 1.0 / tan(fovy / 2.0f);
    
    output._00 = f / aspect;
    output._01 = 0;
    output._02 = 0;
    output._03 = 0;
    
    output._10 = 0;
    output._11 = f;
    output._12 = 0;
    output._13 = 0;
    
    output._20 = 0;
    output._21 = 0;
    output._23 = -1;
    
    output._30 = 0;
    output._31 = 0;
    output._33 = 0;
    
    Type nf = 1 / (near - far);
    output._22 = (far + near) * nf;
    output._32 = 2 * far * near * nf;
    return output;
  }
  
  static Matrix4x4<Type> Orthogonal(Type left, Type right, Type bottom, Type top,
                                    Type near, Type far) {
     Matrix4x4<Type> output;
     
     Type lr = 1 / (left - right);
     Type bt = 1 / (bottom - top);
     Type nf = 1 / (near - far);
    
     output._00 = -2 * lr;
     output._01 = 0;
     output._02 = 0;
     output._03 = 0;
     output._10 = 0;
     output._11 = -2 * bt;
     output._12 = 0;
     output._13 = 0;
     output._20 = 0;
     output._21 = 0;
     output._22 = 2 * nf;
     output._23 = 0;
     output._30 = (left + right) * lr;
     output._31 = (top + bottom) * bt;
     output._32 = (far + near) * nf;
     output._33 = 1;
    
     return output;
  }
  
  // LookAt matrix will construct a new view matrix and transform
  // https://www.3dgep.com/understanding-the-view-matrix/#Look_At_Camera
  static Matrix4x4<Type> LookAtMatrix(const Vector3D<Type>& aEye,
                                      const Vector3D<Type>& aTarget,
                                      const Vector3D<Type>& aUp) {
    
    Matrix4x4<Type> output;
    
    if (
      fabs(aTarget.x - aEye.x) < FLT_EPSILON &&
      fabs(aTarget.y - aEye.y) < FLT_EPSILON &&
      fabs(aTarget.z - aEye.z) < FLT_EPSILON
    ) {
      return output; // return an identity mtx if it is invalidate.
    }

    Vector3D<Type> look;
    look.x = aTarget.x - aEye.x;
    look.y = aTarget.y - aEye.y;
    look.z = aTarget.z - aEye.z;

    look *= -1; // because we are doing mirror.
    look.Normalize();

    Vector3D<Type> right = aUp.CrossProduct(look);
    right.Normalize();
    
    Vector3D<Type> newUp = look.CrossProduct(right);
    newUp.Normalize();
    
//    // Using camera matrix inverse can get the same result from the below.
//    // The below approach is more efficient because inverse() is experisve.
//    Matrix4x4<Type> mtx(
//                        right.x, right.y, right.z, 0,  // inserting by cols
//                        newUp.x, newUp.y, newUp.z, 0,
//                        look.x, look.y, look.z, 0,
//                        aEye.x, aEye.y, aEye.z, 1 );
//    mtx.Inverse();
//    
//    Matrix4x4<Type> mtxR(
//      right.x, right.y, right.z, 0,  // inserting by cols
//      newUp.x, newUp.y, newUp.z, 0,
//      look.x, look.y, look.z, 0,
//      0, 0, 0, 1 );
//    Matrix4x4<Type> mtxT(
//                        1,0,0,0,
//                        0,1,0,0,
//                        0,0,1,0,
//                        aEye.x, aEye.y, aEye.z, 1);
//    Matrix4x4<Type> tmp = mtxT * mtxR;
//    tmp.Inverse();
//    assert(mtx == tmp); // (TxR)^-1 = (Mtx-world)^-1. (T)^-1 = -T, (R) ^ -1 = (R)^T
//    
//    // (TxR)^-1 = R^-1 * T^-1
//    mtxR.Inverse();
//    mtxT.Inverse();
//    Matrix4x4<Type> tmp1 = mtxR * mtxT;
//    assert(mtx == tmp1);
    // +y up, right hand coordinate.
    // (+z is pointing out the screen.)
    output._00 = right.x;
    output._01 = newUp.x;
    output._02 = look.x;
    output._03 = 0;
    output._10 = right.y;
    output._11 = newUp.y;
    output._12 = look.y;
    output._13 = 0;
    output._20 = right.z;
    output._21 = newUp.z;
    output._22 = look.z;
    output._23 = 0;
    output._30 = -right.DotProduct(aEye);
    output._31 = -newUp.DotProduct(aEye);
    output._32 = -look.DotProduct(aEye);
    output._33 = 1;

    //assert(mtx == output);
    return output;
  }
  
  // Difference between TargetTo vs LookAt is TargetTo doesn't
  // transform its position with its rotation matrix.
  static Matrix4x4<Type> TargetTo(const Vector3D<Type>& aEye,
                                  const Vector3D<Type>& aTarget,
                                  const Vector3D<Type>& aUp) {
    Matrix4x4<Type> output;

    if (
      fabs(aTarget.x - aEye.x) < FLT_EPSILON &&
      fabs(aTarget.y - aEye.y) < FLT_EPSILON &&
      fabs(aTarget.z - aEye.z) < FLT_EPSILON
    ) {
      return output; // return an identity mtx if it is invalidate.
    }

    Vector3D<Type> look;
    look.x = aTarget.x - aEye.x;
    look.y = aTarget.y - aEye.y;
    look.z = aTarget.z - aEye.z;

    look *= -1; // because we are doing mirror.
    look.Normalize();

    Vector3D<Type> right = aUp.CrossProduct(look);
  //  Vector3D<Type> right = look.CrossProduct(aUp);
    right.Normalize();
    
    Vector3D<Type> newUp = look.CrossProduct(right);
    //Vector3D<Type> newUp = right.CrossProduct(look);
    newUp.Normalize();

    // The original world matrix instead of
    output._00 = right.x;
    output._01 = right.y;
    output._02 = right.z;
    output._03 = 0;
    output._10 = newUp.x;
    output._11 = newUp.y;
    output._12 = newUp.z;
    output._13 = 0;
    output._20 = look.x;
    output._21 = look.y;
    output._22 = look.z;
    output._23 = 0;
    output._30 = aEye.x;
    output._31 = aEye.y;
    output._32 = aEye.z;
    output._33 = 1;

    return output;
  }
  
  static Matrix4x4<Type> FPSView(const Vector3D<Type>& aEye, Type pitch, Type yaw) {
    // I assume the values are already converted to radians.
    Type cosPitch = cos(pitch);
    Type sinPitch = sin(pitch);
    Type cosYaw = cos(yaw);
    Type sinYaw = sin(yaw);

    Vector3D<Type> right(cosYaw, 0, -sinYaw);
    Vector3D<Type> up(sinYaw * sinPitch, cosPitch, cosYaw * sinPitch);
    Vector3D<Type> look(sinYaw * cosPitch, -sinPitch, cosPitch * cosYaw);

    // Create a 4x4 view matrix from the right, up, forward and eye position vectors
    Matrix4x4<Type> viewMatrix = {
      right.x,  up.x, look.x, 0,
      right.y,  up.y, look.y, 0,
      right.z,  up.z, look.z, 0,
      -right.DotProduct(aEye), -up.DotProduct(aEye), -look.DotProduct(aEye), 1
    };

    return viewMatrix;
  }
  
  static Matrix4x4<Type> Arcball(const Vector3D<Type>& aT0, const Quaternion<Type>& r,
                                 const Vector3D<Type>& aT1 = Vector3D<Type>(0, 0, 0)) {
//    Matrix4x4<Type> t0;
//    t0.Translate(aT0.x, aT0.y, aT0.z); // Translation away from object.
//    Quaternion<Type> r0(r);
//    r0.Normalize();
//    Matrix4x4<Type> rot = FromQuat(r0);
//    Matrix4x4<Type> t1;
//    t1.Translate(aT1.x, aT1.y, aT1.z); // Translate to center of object.
//
//    Matrix4x4<Type> output = t1 * rot * t0;
//    output.Inverse();
    
    // The optimize approach
    Matrix4x4<Type> t0;
    t0.Translate(-aT0.x, -aT0.y, -aT0.z);       // Translation away from object.
    Quaternion<Type> r1(r);
    r1.Inverse();
    Matrix4x4<Type> rot = FromQuat(r1);       // Rotate around object.
    Matrix4x4<Type> t1;
    t1.Translate(-aT1.x, -aT1.y, -aT1.z);       // Translate to center of object.

    Matrix4x4<Type> output = t0 * rot * t1;
    
    return output;
  }
  
  static Matrix4x4<Type> FromQuat(const Quaternion<Type>& q) {
    Matrix4x4<Type> output;
    
    Type x = q.x,
         y = q.y,
         z = q.z,
         w = q.w;

     Type x2 = x + x;
     Type y2 = y + y;
     Type z2 = z + z;

     Type xx = x * x2;
     Type yx = y * x2;
     Type yy = y * y2;
     Type zx = z * x2;
     Type zy = z * y2;
     Type zz = z * z2;
     Type wx = w * x2;
     Type wy = w * y2;
     Type wz = w * z2;

     output._elements[0] = 1 - yy - zz;
     output._elements[1] = yx + wz;
     output._elements[2] = zx - wy;
     output._elements[3] = 0;

     output._elements[4] = yx - wz;
     output._elements[5] = 1 - xx - zz;
     output._elements[6] = zy + wx;
     output._elements[7] = 0;

     output._elements[8] = zx + wy;
     output._elements[9] = zy - wx;
     output._elements[10] = 1 - xx - yy;
     output._elements[11] = 0;

     output._elements[12] = 0;
     output._elements[13] = 0;
     output._elements[14] = 0;
     output._elements[15] = 1;

     return output;
  }
  
  
  // column-major order. _01: means 0 col, 1 row. _10: means 1 col, 0 row
  union {
    struct {
      Type _00, _01, _02, _03;
      Type _10, _11, _12, _13;
      Type _20, _21, _22, _23;
      Type _30, _31, _32, _33;
    };

    Type _elements[16];
  };
};

typedef Matrix4x4<float> Matrix4x4f;

}  // End of namespace gfx_math


#endif /* GFX_MATH_MATRIX4x4_H */
