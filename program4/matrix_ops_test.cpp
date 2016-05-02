#include <iostream>

#include "matrix_ops.h"

#define rows    3
#define cols    3

int main(void)
{
  double ** mat = new double*[rows];
  double ** invMat = new double*[rows];
  for(int i = 0; i < cols; i++)
  {
    mat[i] = new double[cols];
    invMat[i] = new double[cols];
  }

  /**
  double input_mat[rows][rows] = {{1,4,5,6,9},
                                  {1,4,7,2,3},
                                  {6,4,7,3,3},
                                  {9,5,3,2,3},
                                  {1,2,3,4,5}};
  */
  double input_mat[rows][rows] = {{5,6,9},
                                  {1,2,3},
                                  {6,4,7}};
  
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < cols; j++)
    {
      mat[i][j] = input_mat[i][j];
      invMat[i][j] = 0;
    }
  }

  cout << "det = " << getDeterminant(mat, rows) << endl;
  inverse(mat, rows, invMat);
  cout << "invMat = ";
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < cols; j++)
    {
      cout << invMat[i][j] << "\t";
    }
  }
  cout << endl;

  for(int i = 0; i < cols; i++)
  {
    delete[] mat[i];
    delete[] invMat[i];
    mat[i] = NULL;
    invMat[i] = NULL;
  }
  delete[] mat;
  delete[] invMat;
  invMat = NULL;
  mat = NULL;

  return 1;
}
