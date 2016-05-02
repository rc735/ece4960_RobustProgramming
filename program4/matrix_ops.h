#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <iostream>
#include <cmath>

using namespace std;

/*
 *http://www.sanfoundry.com/cpp-program-find-determinant-given-matrix/
 * getDeterminant - calculates the determinant of a double array of size n
 */
//double getDeterminant(double mat[SIZE_PARAM][SIZE_PARAM], int n);

/*
 *http://www.sanfoundry.com/cpp-program-find-determinant-given-matrix/
 * getDeterminant - calculates the determinant of a double array of size n
 */
double getDeterminant(double ** mat, int n);


/**
 * inverse - calculates the inverse of a double array and outputs
 *            the result to final
 */
//void inverse(double mat[SIZE_PARAM][SIZE_PARAM], int n,
//              double final[SIZE_PARAM][SIZE_PARAM]);

/**
 * inverse - calculates the inverse of a double array and outputs
 *            the result to final
 */
void inverse(double **  mat, int n, double **final);

#endif /* MATRIX_OPS_H */
