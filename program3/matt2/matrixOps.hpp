//
//  matrixOps.hpp
//  Program3
//
//  Created by Matthew Magaldi on 3/23/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#ifndef matrixOps_hpp
#define matrixOps_hpp

#include <stdio.h>


void MatrixInversion(double **, int , double **);
int GetMinor(double **, double **, int , int , int );
double CalcDeterminant( double **, int );
void MatrixVectorMultiply(double **, double *, double *x, int);

#endif /* matrixOps_hpp */
