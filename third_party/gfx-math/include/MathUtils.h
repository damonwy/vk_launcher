//
//  vector3d.h
//  gfx-math
//
//  Created by Daosheng Mu on 9/26/20.
//  Copyright © 2020 Daosheng Mu. All rights reserved.
//

#ifndef GFX_MATH_UTILS_H
#define GFX_MATH_UTILS_H

#include <math.h>
#include <cfloat>
#include <vector>

namespace gfx_math {

#define PI 3.14159265

template <class Type>
Type DegreesToRadians(Type aDegree) {
  return aDegree * PI / 180.0f;
}

} // End of namespace gfx_math


#endif /* GFX_MATH_UTILS_H */
