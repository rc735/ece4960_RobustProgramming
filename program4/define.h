//
//  define.h
//  Program4
//
//  Created by Matthew Magaldi on 4/20/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#ifndef define_h
#define define_h

#define FORWARD_EULER   0
#define BACKWARD_EULER  1
#define TRAPEZOIDAL     2
#define RUNGE_KUTTA0    3
#define RUNGE_KUTTA1    4

#define GAMMA           0.5
#define ERROR_TOL       1e-2
#define H               0.2e-9//0.1
#define P               4
#define VT				0.026

#define TASK_NUM		2	//1-validation, 2-RC, 3-EKV Common source amp
#define NUM_METHODS     5   //number of methods to iterate through
#define NUM_VARIABLES   2   //size of x
#define NUM_FUNCTIONS   1   //size of f

#endif /* define_h */
