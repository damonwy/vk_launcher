//
//  testMatrix4x4.m
//  gfx-math-test
//
//  Created by Daosheng Mu on 1/9/21.
//  Copyright © 2021 Daosheng Mu. All rights reserved.
//

#import <XCTest/XCTest.h>

#include <iostream>
#include "Vector3d.h"
#include "Matrix4x4.h"
#include "Quaternion.h"

using namespace gfx_math;

@interface testMatrix4x4 : XCTestCase

@end

@implementation testMatrix4x4

- (void)setUp {
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

- (void)testBasic {
  Matrix4x4f matA;
  Matrix4x4f matB;
  matA.Identity();
  XCTAssertEqual(matA, matB);
  
  Matrix4x4f matC = matA * matB;
  XCTAssertEqual(matC, matA);
}

- (void)testCameraMatrix {
  Matrix4x4f matC;
  
  matC = Matrix4x4f::LookAtMatrix(Vector3Df(0, 0.5, 0.5),
                                        Vector3Df(0,0.2,-1),
                                        Vector3Df(0,1,0));
  
  Matrix4x4f matD = Matrix4x4f::TargetTo(Vector3Df(0, 0.5, 0.5),                                                                                               Vector3Df(0,0.2,-1),
                                                     Vector3Df(0,1,0));
  XCTAssertNotEqual(matC, matD); // TODO: We haven't known how to test it.
  
  Matrix4x4f matE = Matrix4x4f::FPSView(Vector3Df(0, 0.5, 0.5), 0.1, 0.2);
  
  Matrix4x4f matF = Matrix4x4f::Arcball(Vector3Df(0, 0.5, 0.5),
                                        Quaternion<float>(1,1,0,1),
                                        Vector3Df(1.0, 0.5, 0.5));
  XCTAssertNotEqual(matE, matF); // TODO: We haven't known how to test it.
}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
    }];
}

@end
