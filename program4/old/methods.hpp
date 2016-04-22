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

double f(double, double);

double Feuler(double, double);

double Beuler(double, double);

double trapezoidal(double, double);

double RK34woAdapt(double, double);

double RK34wAdapt(double, double, double &, double h);

#endif /* methods_hpp */
