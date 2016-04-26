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
    
    double x0[NUM_VARIABLES] = {2.0};
    double timeSpan = 7;

    /** ORIGINAL FILE MAKER
    matrix1 = solver(x0, timeSpan, FORWARD_EULER);
    
    //output matrix to file to graph
    cout << "creating output text file..." << endl;
    ofstream myfile;
    myfile.open ("output.txt");
    if (myfile.is_open())
    { // ok, proceed with output
        cout  << "file is open." << endl;
        for (int i=0; i < (int)matrix1[0].size(); i++) {
            myfile << matrix1[0][i] << "\t" << matrix1[1][i] << "\n";
        }
        myfile.close();
    } else {
        cout << "could not open file." << endl;
    }
    */

    cout << "Creating output text file..." << endl;
    char filename[20];
    for(int j = 0; j < NUM_METHODS; j++)
    {
      matrix1 = solver(x0, timeSpan, j, NUM_VARIABLES);

      ofstream myfile;
      sprintf(filename, "output_%d.txt", j);
      myfile.open(filename);
      if(myfile.is_open())
      {
        cout << "writing to file" << endl;
        for(int i = 0; i < (int) matrix1[0].size(); i++)
        {
          myfile << matrix1[0][i] << "\t" << matrix1[1][i] << endl;
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
