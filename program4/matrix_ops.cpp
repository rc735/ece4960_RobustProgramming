#include "matrix_ops.h"


/*
 *http://www.sanfoundry.com/cpp-program-find-determinant-given-matrix/
 * getDeterminant - calculates the determinant of a double array of size n
 */
double getDeterminant(double ** mat, int n)
{
  double d = 0;
  int subi, subj;

  // 1 by 1 matrix
  if(n == 1)
  {
    return mat[0][0];
  }
  // 2 by 2 matrix
  else if(n == 2)
  {
    return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
  }
  
  double ** submat = new double*[n-1];
  for(int i = 0; i < n-1; i++)
  {
    submat[i] = new double[n-1];
  }

  // n by n matrix
  for(int c = 0; c < n; c++)
  {
    subi = 0;
    for(int i = 1; i < n; i++)
    {
      subj = 0;
      for(int j = 0; j < n; j++)
      {
        if(j == c)
        {
          continue;
        }
        submat[subi][subj] = mat[i][j];
        subj++;
      }
      subi++;
    }
    d = d + (pow(0-1.0, c) * mat[0][c] * getDeterminant(submat, n-1));
  }

  // deallocating sub-matrix
  for(int i = 0; i < n-1; i++)
  {
    delete[] submat[i];
    submat[i] = NULL;
  }
  delete[] submat;
  submat = NULL;

  return d;
}


/**
 * inverse - calculates the inverse of a double array and outputs
 *            the result to final
 */
void inverse(double ** mat, int n, double **final)
{
  if(n == 1)
  {
    final[0][0] = 1/mat[0][0];
  }
  else if(n == 2)
  {
    double det = getDeterminant(mat, n);
    if(det == 0)
    {
      cout << "matrix_ops.cpp: inverse() - determinant is 0" << endl;
      return;
    }
    final[0][0] = mat[1][1] / det;
    final[0][1] = (0 - mat[0][1]) / det;
    final[1][0] = (0 - mat[1][0]) / det;
    final[1][1] = mat[0][0] / det;
  }

  //https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3.C3.973_matrices
  else if(n == 3)
  {
    double det = getDeterminant(mat, n);
    if(det == 0)
    {
      cout << "matrix_ops.cpp: inverse() - determinant is 0" << endl;
      return;
    }
    final[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) / det;
    final[1][0] = (0 - (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])) / det;
    final[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) / det;
    final[0][1] = (0 - (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1])) / det;
    final[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) / det;
    final[2][1] = (0 - (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0])) / det;
    final[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) / det;
    final[1][2] = (0 - (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0])) / det;
    final[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) / det;
  }
  else
  {
    cout << "matrix_ops.cpp: inverse() - greater than the allowed size" << endl;
    return;
  }

}
