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
#define ERROR_TOL       1e-7
#define P               4
#define VT				      0.026

#define NUM_METHODS     5   //number of methods to iterate through

// global variables - defined in main.cpp
extern int task_num;
extern double H;
extern int num_variables;

#endif /* define_h */
