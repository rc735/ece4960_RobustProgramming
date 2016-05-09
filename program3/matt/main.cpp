//
//  main.cpp
//  Program3
//
//  Created by Matthew Magaldi on 3/19/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#include <iostream>
#include "validation.hpp"
#include "matrixOps.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    
    validityTest();
    
    /*
    //test the matrix operations here
    double ** hessian = 0;
    double ** inverse = 0;
    hessian = new double*[2];
    inverse = new double*[2];
    for (int i = 0; i < 2 ; i++) {
        hessian[i] = new double[2];
        inverse[i] = new double[2];
    }
    
    
    hessian[0][0] = 0.23523;
    hessian[0][1] = 2.23;
    hessian[1][0] = 1.32;
    hessian[1][1] = 0.235235;
    
    MatrixInversion(hessian, 2, inverse);
    
    for (int i=0; i < 2; i++) {
        for (int j=0; j < 2; j++) {
            cout << inverse[i][j] << endl;
        }
    }
    */
    
    return 0;
}
