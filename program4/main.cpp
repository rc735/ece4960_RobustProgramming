//
//  main.cpp
//  Program4
//
//  Created by Matthew Magaldi on 4/19/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "define.h"
#include "solver.hpp"

using namespace std;

//global variables
int task_num;
double H;
int num_variables;

/**
 * main()
 *    @param argv[1]      the task number to execute
 *                            1 --> validation
 *                            2 --> RC
 *                            3 --> EKV common source amp
 */
int main(int argc, const char * argv[]) {
    // check for number of arguments
    if(argc != 2)
    {
      cout << "ERROR - choose a task to perform" << endl;
      cout << "\t(1 -> validation; 2 -> RC; 3 -> CS Amp)" << endl;
      exit(1);
    }

    // initialize global variables and initial guess "x0"
    task_num = atoi(argv[1]);
    double * xtemp;
    double timeSpan;
    if(task_num == 1)             // validation
    {
      xtemp = new double[1];
      xtemp[0] = 2.0;
      timeSpan = 7;
      H = 0.1;
      num_variables = 1;
    }
    else if(task_num == 2)        // RC
    {
      xtemp = new double[2];
      xtemp[0] = 0.0;
      xtemp[1] = 0.0;
      timeSpan = 100e-9;
      H = 0.2e-9;
      num_variables = 2;
    }
    else if(task_num == 3)        // CS amp
    {
      xtemp = new double[2];
      xtemp[0] = 2.5;
      xtemp[1] = 2.5;
      timeSpan = 100e-9;
      H = 0.2e-9;
      num_variables = 2;
    }
    else
    {
      cout << "ERROR - Invalid Task Selection" << endl;
      cout << "\t(1 -> validation; 2 -> RC; 3 -> CS Amp)" << endl;
      exit(1);
    }
    double x0[num_variables];
    for(int i = 0; i < num_variables; i++)
    {
      x0[i] = xtemp[i];
    }
    delete[] xtemp;
    xtemp = NULL;

    // stores all time steps and values of x at each time step
    matrix matrix1;

    //Create Output Files for Plotting
    cout << "Creating output text files..." << endl << endl;
    char filename[20];
    for(int j = 0; j < NUM_METHODS; j++)
    {
      matrix1 = solver(x0, timeSpan, j, num_variables);
      
      // for verification purposes
      cout << "NUM ITERATIONS:\t" << matrix1[0].size() << endl;
      cout << "FINAL VALUES:\t";
      for(int i = 0; i < num_variables+1; i++)
      {
        cout << matrix1[i][matrix1[0].size()-1] << "\t";
      }
      cout << endl;

      ofstream myfile;
      sprintf(filename, "output_%d.txt", j);
      myfile.open(filename);
      if(myfile.is_open())
      {
        cout << "writing to: " << filename << endl << endl;
        for(int i = 0; i < (int) matrix1[0].size(); i++)
        {
			    for (int k = 0; k < num_variables+1; k++) {
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
    }
    
    cout << "Finished execution of program 4" << endl;
    
    return 0;
}
