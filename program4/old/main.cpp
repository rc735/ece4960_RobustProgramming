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
#include "solver.hpp"
#include "define.h"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    
    matrix matrix1;
    
    double x0 = 2.0;
    double timeSpan = 7.0;
    matrix1 = solver(x0, timeSpan, BACKWARD_EULER);
    
    //output matrix to file to graph
    cout << "creating output text file..." << endl;
    ofstream myfile;
    myfile.open ("output.txt");
    if (myfile.is_open())
    { // ok, proceed with output
        cout  << "file is open." << endl;
        for (int i=0; i < matrix1[0].size(); i++) {
            myfile << matrix1[0][i] << "\t" << matrix1[1][i] << "\n";
        }
        myfile.close();
    } else {
        cout << "could not open file." << endl;
    }
    
    cout << "Finished execution of program 4" << endl;
    
    return 0;
}
