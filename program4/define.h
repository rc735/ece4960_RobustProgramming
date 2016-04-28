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

#define TASK_NUM		    2	  //1-validation, 2-RC, 3-EKV Common source amp
#define NUM_METHODS     5   //number of methods to iterate through


#if TASK_NUM==1
#define NUM_VARIABLES   1   //size of x
#define H               0.1

#elif TASK_NUM==2
#define NUM_VARIABLES   2
#define H               0.2e-9

#elif TASK_NUM==3
#define NUM_VARIABLES   3
#define H               0.2e-9
#endif

#endif /* define_h */
