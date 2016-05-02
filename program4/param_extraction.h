#ifndef PARAM_EXTRACTION_H
#define PARAM_EXTRACTION_H

#include <iostream>
#include <cmath>

#include "matrix_ops.h"
#include "methods.hpp"

using namespace std;


/**
 * newton_f - function used in newton method
 */
double newton_f(double x_i, double x_i1, double h,
                double phi_i, double phi_i1, int task);

/**
 * newtonMethod - for estimating the power law example using newtown method
 */
void newtonMethod(double * x, double t, double h, int num_x, int task,
                    int max_iter = 10);


#endif /* PARAM_EXTRACTION_H */
