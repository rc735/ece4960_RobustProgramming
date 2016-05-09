#ifndef PARAM_EXTRACTION_H
#define PARAM_EXTRACTION_H

#include <iostream>
#include <cmath>

#include "defines.h"
#include "matrix_ops.h"
#include "model_transistor.h"

using namespace std;

/**
 * generateS - for the power law; calculating one value of S if
 *              S is a polynomial expression
 *
 *  @param a      parameters to calcuate S
 *  @param x      variables to calculate S
 *  @param size   number of parameter and variable elements
 *  @return       the resulting S based on S = a0 + a1*x1 + ...
 */
double generateS(double* a, double* x, int size);

/**
 * generateS - for power law; generating an array of S values
 *  
 *  @param a      parameters to calculate S
 *  @param x      variables to calculate S
 *  @param S      results stored in this array; S
 *  @param asize  number of parameters
 *  @param xsize  number of variables
 */
void generateS(double* a, double* x, double* S, int asize, int xsize);

/**
 * randError - calculates a random error based on the percent
 *
 *  @param a          the value to be calculated with the error
 *  @param percent    the percent error to be multiplied to the value
 *  @result           outputs the new value with random error
 */
double randError(double a, int percent);

/**
 *generateV - calculates the sum of the squared error for all variables
 *
 *  @param a            the parameters to calculate V
 *  @param x            the variables to calculate V
 *  @param S_measured   the measured S values
 *  @param paramSize    number of parameter elements
 *  @param measSize     number of measurements of S
 *  @return             the squared error
 */
double generateV(double *a, double *x, double *S_measured,
                 int paramSize, int measSize);


/**
 * newtonMethod - for estimating the power law example using newtown method
 */
void newtonMethod(double *guess, double *a_final, double *x,
                  double *Smeas, int paramSize, int measSize);

/**
 * secantMethod - for estimating the power law example using secant method
 *                slope is estimated using the equation in:
 *                https://en.wikipedia.org/wiki/Secant_method
 */
void secantMethod(double *guess, double *a_final, double *x,
                  double *Smeas, int paramSize, int measSize);


#endif /* PARAM_EXTRACTION_H */
