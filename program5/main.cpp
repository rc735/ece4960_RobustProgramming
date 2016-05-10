#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "define.h"
#include "program3/param_extraction.h"
#include "program3/matrix_ops.h"
#include "program4/solver.hpp"


using namespace std;

//----------- GLOBAL VARIABLES ---------
int task_num;
double H;

//  MAIN()
int main(int argc, const char * argv[])
{
  //-------------------------------------------------------------
  // calculating the V1(t) and V2(t) - program 4
  // utilizing ODE solvers - Runge Kutta 34 w/ time adaptivity
  task_num = 2;     //RC circuit
  H = 1e-10;
  double timeSpan = 39e-6;
  double x0[NUM_VAR] = {0.0, 0.0};

  // stores all time steps and values of x at each time step
  matrix matrix1;

  matrix1 = solver(x0, timeSpan, RUNGE_KUTTA1, NUM_VAR);
  cout << "NUM ITERATIONS:\t" << matrix1[0].size() << endl;
  cout << "FINAL VALUES:\t";
  for(int i = 0; i < NUM_VAR + 1; i++)
  {
    cout << matrix1[i][matrix1[0].size()-1] << "\t";
  }
  cout << endl;

  char filename[20];
  ofstream myfile;
  sprintf(filename, "output_%d.txt", 1);
  myfile.open(filename);
  if(myfile.is_open())
  {
    cout << "writing to: " << filename << endl << endl;
    for(int i = 0; i < (int) matrix1[0].size(); i++)
    {
	    for (int k = 0; k < NUM_VAR+1; k++)
      {
	      myfile << matrix1[k][i] << "\t";
			}
      myfile << "\n";
    }
    myfile.close();
  }
  else
  {
    cout << "ERROR - file failed to open" << endl;
  }

  cout << "V1(t) and V2(t) have finished calculating and is located in "
          << filename << endl << endl;


  //-------------------------------------------------------------
  // parameter extraction - program 3
  cout << "Begin parameter extraction for a2/a1, tau1, tau2" << endl;

  //determine which index in matrix to start the parameter extraction
  int idx;
  for(int i = 0; i < (int) matrix1[0].size(); i++)
  {
    if(matrix1[0][i] >= 20.0005e-6)
    {
      idx = i;
      break;
    }
  }

  // V2(t) = c1(exp(-t/tau1) + c2/c1 * exp(-t/tau2))
  // a[0] = c2/c1
  // a[1] = tau1
  // a[2] = tau2
  double guess[SIZE_PARAM] = {1.5, 2e4*1e-12, 2e4*1e-12};
  double a_final[SIZE_PARAM];
  double Smeas[SIZE_MEAS], x[SIZE_MEAS];
  for(int i = 0; i < SIZE_MEAS; i ++)
  {
    x[i] = matrix1[0][i+idx] - 20.0005e-6;
    Smeas[i] = matrix1[2][i];
  }
  newtonMethod_pg3(guess, a_final, x, Smeas, SIZE_PARAM, SIZE_MEAS);
  cout << "FINAL NEWTON: ";
  for(int i = 0; i< SIZE_PARAM; i++)
  {
    cout << a_final[i] << " ";
  }
  cout << endl << endl;
  

  return 0;
}
