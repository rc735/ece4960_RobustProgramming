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

#include "define.h"
#include "solver.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    
    matrix matrix1;
    
#if TASK_NUM==1
    double x0[NUM_VARIABLES] = {2.0};
    double timeSpan = 7;

#elif TASK_NUM==2
    double x0[NUM_VARIABLES] = {0.0, 0.0};
    double timeSpan = 100e-9;

#elif TASK_NUM==3
    double x0[NUM_VARIABLES] = {0.0, 0.0};
    double timeSpan = 100e-9;

#else
    cout << "ERROR - choose TASK_NUM in define.h" << endl;
    exit(1);
#endif

    //Create Output Files for Plotting
    cout << "Creating output text files..." << endl << endl;
    char filename[20];
    for(int j = 0; j < NUM_METHODS; j++)
    {
      matrix1 = solver(x0, timeSpan, j, NUM_VARIABLES);
      
      // for verification purposes
      cout << "NUM ITERATIONS:\t" << matrix1[0].size() << endl;
      cout << "FINAL VALUES:\t";
      for(int i = 0; i < NUM_VARIABLES+1; i++)
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
			    for (int k = 0; k < NUM_VARIABLES + 1; k++) {
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
