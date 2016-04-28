//
//  methods.hpp
//  Program4
//
//  Created by Matthew Magaldi on 4/20/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#ifndef methods_hpp
#define methods_hpp

#include <stdio.h>
#include <cmath>

#include "define.h"

void f(double*, const double*, double, int);

void Feuler(double*, const double*, double, int);

void Beuler(double*, const double*, double, int);

void trapezoidal(double*, const double*, double, int);

void RK34woAdapt(double*, const double*, double, int);

void RK34wAdapt(double*, const double*, double, double &, double h, int);

double normNum(double* x, int num_x);

#endif /* methods_hpp */
