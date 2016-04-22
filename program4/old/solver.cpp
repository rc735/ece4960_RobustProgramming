//
//  solver.cpp
//  Program4
//
//  Created by Matthew Magaldi on 4/20/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#include <iostream>

#include "solver.hpp"

using namespace std;


//iterative solver function
//returns a matrix with 2 rows that represent the time and data while columns represent progression
//want to return a matrix of data and time from this equation

//in order to go to more dependent variables, could expand matrix size along the rows ****
matrix solver(double x, double t_final, int methodChoice)
{
    int no_of_rows = 2;
    
    //initialize matrix of vectors
    matrix matrix1;
    matrix1.resize(no_of_rows, std::vector<double> ());
    
    //set initial condition
    matrix1[1].push_back(x);
    matrix1[0].push_back(0);
    
    //iterate
    double timei0, phi, h, h_old, error, R, cur_time = H;
    h_old = H;
    h = H;
    for (int i = 1; cur_time <= t_final; i++) {
        matrix1[0].push_back(cur_time);
        timei0 = matrix1[0][i];
        phi = 0;                 //increment value(s), now choose function
        //h = H;
        
        switch (methodChoice) {
            case FORWARD_EULER:
                cur_time += H;
                phi = Feuler(x, timei0);
                break;
            case BACKWARD_EULER:
                cur_time += H;
                phi = Beuler(x, timei0);
                break;
            case TRAPEZOIDAL:
                cur_time += H;
                phi = trapezoidal(x,timei0);
                break;
            case RUNGE_KUTTA0:
                cur_time += H;
                phi = RK34woAdapt(x, timei0);
                break;
            case RUNGE_KUTTA1:
                error = 0;
                phi = RK34wAdapt(x, timei0, error, h);
                R = ERROR_TOL/error;
                if (R > 2.0 || R < 0.5) {
                    h = GAMMA*h_old*pow((ERROR_TOL/error),1.0/3.0);
                }
                if (h <= 0) {
                    cout << "ERROR: time step is 0" << endl;
                    exit(1);
                }
                cur_time += h;
                cout << "h: " << h << endl;
                h_old = h;
                break;
            default:
                break;
        }
        //x should be a vector of depedent variables
        x = x + phi*h;
        matrix1[1].push_back(x);
    }
    
    return matrix1;
}

