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

double f(const double*, double, int);

double Feuler(const double*, double, int);

double Beuler(const double*, double, int);

double trapezoidal(const double*, double, int);

double RK34woAdapt(const double*, double, int);

double RK34wAdapt(const double*, double, double &, double h, int);

#endif /* methods_hpp */
