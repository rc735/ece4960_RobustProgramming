#include "matrix_ops.h"
#include "defines.h"

/*
 *http://www.sanfoundry.com/cpp-program-find-determinant-given-matrix/
 * getDeterminant - calculates the determinant of a double array of size n
 */
double getDeterminant(double mat[SIZE_PARAM][SIZE_PARAM], int n)
{
  double d = 0;
  int c, subi, i, j, subj;
  double submat[SIZE_PARAM][SIZE_PARAM];

  if(n == 2)
  {
    return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
  }

  for(c = 0; c < n; c++)
  {
    subi = 0;
    for(i = 1; i < n; i++)
    {
      subj = 0;
      for(j = 0; j < n; j++)
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
  return d;
}

/*
 *http://www.sanfoundry.com/cpp-program-find-determinant-given-matrix/
 * getDeterminant - calculates the determinant of a double array of size n
 */
double getDeterminant(double **mat, int n)
{
  double d = 0;
  int c, subi, i, j, subj;
  double submat[SIZE_PARAM][SIZE_PARAM];

  if(n == 2)
  {
    return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
  }

  for(c = 0; c < n; c++)
  {
    subi = 0;
    for(i = 1; i < n; i++)
    {
      subj = 0;
      for(j = 0; j < n; j++)
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
  return d;
}

/**
 * inverse - calculates the inverse of a double array and outputs
 *            the result to final
 */
void inverse(double mat[SIZE_PARAM][SIZE_PARAM], int n,
              double final[SIZE_PARAM][SIZE_PARAM])
{
  if(n == 2)
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


/**
 * inverse - calculates the inverse of a double array and outputs
 *            the result to final
 */
void inverse(double **mat, int n, double **final)
{
  if(n == 2)
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
