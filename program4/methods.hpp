//
//  methods.hpp
//  Program4
//
//  Created by Matthew Magaldi on 4/20/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#ifndef METHODS_H
#define METHODS_H

#include <stdio.h>
#include <cmath>
#include <iostream>

#include "define.h"
#include "param_extraction.h"

using namespace std;


void f(double*, const double*, double, int);

void Feuler(double*, double, int);

void Beuler(double*, double, int);

void trapezoidal(double*, double, int);

void RK34woAdapt(double*, double, int);

void RK34wAdapt(double*, double, double &, double h, int);

double normNum(double* x, int num_x);

#endif /* METHODS_H */
