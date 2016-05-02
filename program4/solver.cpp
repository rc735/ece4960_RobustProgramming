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
matrix solver(const double* x_guess, double t_final, int methodChoice, int num_x)
{
    int no_of_rows = num_x+1;
    
    //initialize matrix of vectors
    matrix matrix1;
    matrix1.resize(no_of_rows, vector<double>());
    
    //set initial condition
    matrix1[0].push_back(0);
    double x[num_x];
    double phi[num_x];
    for(int i = 0; i < num_x; i++)
    {
      x[i] = x_guess[i];
      matrix1[i+1].push_back(x_guess[i]);

      //initialize phi
      phi[i] = 0;//increment value(s), now choose function
    }
    
    //iterate
    double timei0, h, h_old, error, R, cur_time = H;
    h_old = H;
    h = H;
    for (int i = 1; cur_time <= t_final; i++) {
        matrix1[0].push_back(cur_time);
        timei0 = matrix1[0][i];
        
        switch (methodChoice) {
            case FORWARD_EULER:
                cur_time += H;
                //phi = Feuler(x, timei0, num_x);
                Feuler(x, timei0, num_x);
                break;

            case BACKWARD_EULER:
                cur_time += H;
                //phi = Beuler(x, timei0, num_x);
                Beuler(x, timei0, num_x);
                break;

            case TRAPEZOIDAL:
                cur_time += H;
                //phi = trapezoidal(x,timei0, num_x);
                trapezoidal(x, timei0, num_x);
                break;

            case RUNGE_KUTTA0:
                cur_time += H;
                //phi = RK34woAdapt(x, timei0, num_x);
                RK34woAdapt(x, timei0, num_x);
                break;

            case RUNGE_KUTTA1:
                error = 0;
                RK34wAdapt(x, timei0, error, h, num_x);

                //1e-15 is the absolute tolerance
                //cout << error << " ==> " << (ERROR_TOL * normNum(x, num_x) + 1e-15) << endl;
                //R = error / (ERROR_TOL * normNum(x, num_x) + 1e-15);
                R = error / normNum(x, num_x);
                //cout << "R =  " << R << endl;
                /*
                if (R > 2.0 || R < 0.5) {
                    //h = GAMMA*h_old*pow(R,1.0/3.0);
                    h = h/R;
                    cout << "h = " << h_old << " --> " << h << endl;
                }*/
                if(R > 1e-2)
                {
                  h = h/2.0;
                }
                else if(R < ERROR_TOL)
                {
                  h = 2.0*h;
                }
                if (h <= 0) {
                    cout << "ERROR: time step is 0" << endl;
                    exit(1);
                }
                cur_time += h;
                //cout << "h: " << h << endl;
                h_old = h;
                break;

            default:
                break;
        } 
        //x should be a vector of dependent variables
        for(int j = 0; j < num_x; j++)
        {
          //x[j] += phi[j]*h;
          matrix1[j+1].push_back(x[j]);
        }
    }
    
    return matrix1;
}

