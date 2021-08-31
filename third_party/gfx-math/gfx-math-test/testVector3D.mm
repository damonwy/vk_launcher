//
//  testVector3D.m
//  gfx-math-test
//
//  Created by Daosheng Mu on 1/9/21.
//  Copyright © 2021 Daosheng Mu. All rights reserved.
//

#import <XCTest/XCTest.h>
#include "Vector3d.h"
using namespace gfx_math;

@interface testVector3D : XCTestCase

@end

@implementation testVector3D

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testVectorBasic {
  Vector3Df a(1, 0, 0);
  
  float res = a.DotProduct(Vector3Df(1, 0, 0));
  XCTAssertEqual(res, 1.0f);
}

- (void)testRayIntersection {
  Vector3Df ray(1, 0, 0);
  Vector3Df p0(-2, 0, 0);
  Vector3Df planeNrm(-1, 0, 0);
  float planeDist = 10.0f;
  Vector3Df interPoint;
  
  bool result = RayIntersectWithPlane(ray, p0, planeNrm, planeDist, interPoint);
  XCTAssertTrue(result);
  
  result = PointIntersectWithPlane(interPoint, planeNrm, planeDist);
  XCTAssertFalse(result);
  
  float dist = PointDistanceWithPlane(interPoint + Vector3Df(1, 0, 0), planeNrm, planeDist);
  XCTAssertEqual(dist, 1.0f);
  
  std::vector<Vector3Df> interPointList;
  result = RayIntersectWithSphere(ray, Vector3Df(-10, 0, 0), Vector3Df(0,0,0), 5.0f, interPointList);
  XCTAssertTrue(result);
  XCTAssertEqual(interPointList.size(), 2);
  
  interPointList.clear();
  result = RayIntersectWithSphere(ray, p0, Vector3Df(0,0,0), 5.0f, interPointList);
  XCTAssertTrue(result);
  XCTAssertEqual(interPointList.size(), 1);
  
  Vector3Df p;
  result = RayIntersectWithTriangle(Vector3Df(0,0,-1), p0, Vector3Df(-10,-1,-5),
                                    Vector3Df(10,-1,-5), Vector3Df(5,1,-5), p);
  XCTAssertTrue(result);
  XCTAssertTrue(p == Vector3Df(-2, 0, -5));
  
  result = RayIntersectWithAABBox(ray, p0, Vector3Df(2, -2, -2), Vector3Df(6, 2, 2), p);
  XCTAssertTrue(result);
  XCTAssertTrue(p == Vector3Df(2, 0, 0));
  
  result = SphereIntersectWithAABBox(Vector3Df(0,0,0), 5.0f, Vector3Df(2, -2, -2), Vector3Df(6, 2, 2));
  XCTAssertTrue(result);
  
  result = AABBoxIntersectWithAABBox(Vector3Df(0, 0, -2), Vector3Df(4, 4, 2),
                                     Vector3Df(2, -2, -2), Vector3Df(6, 2, 2));
  XCTAssertTrue(result);
}

- (void)testReflection {
  Vector3Df p = ReflectionVector(Vector3Df(0.707, 0.707, 0), Vector3Df(0, 1, 0));
  XCTAssertTrue(p == Vector3Df(-0.707, 0.707, 0));
}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
    }];
}

@end
