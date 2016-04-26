//
//  methods.cpp
//  Program4
//
//  Created by Matthew Magaldi on 4/20/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#include "methods.hpp"
#include <iostream>

using namespace std;


//function f for validation -- can be anything
double f(const double* x, double t, int num_x)
{
    return 4*exp(0.8*t) - 0.5*x[0];
    //return pow(t, 4)*sin(2*t) - pow(t, 2) + 4*pow(t, 3) + (2/t)*x[0];
}


//forward euler
double Feuler(const double* x, double t, int num_x)
{
    return f(x,t, num_x);
}

//backward euler
double Beuler(const double* x, double t, int num_x)
{
    double temp_x[num_x];
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + f(x, t, num_x)*H;
    }
    return f(temp_x,t+H, num_x);
}

//trapezoidal
double trapezoidal(const double* x, double t, int num_x)
{
    double dxdt0 = f(x,t, num_x);

    double temp_x[num_x];
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + dxdt0*H;
    }
    double dxdt1 = f(temp_x,t+H, num_x);
    return (dxdt0+dxdt1)/2.0;
}

double RK34woAdapt(const double* x, double t, int num_x)
{
    //implement RK34
    double k1, k2, k3, k4, RK3, RK4;
    
    //k1
    k1 = f(x,t, num_x);

    //k2
    double temp_x[num_x];
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.5*k1*H;
    }
    k2 = f(temp_x,t+0.5*H, num_x);

    //k3
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.75*k2*H;
    }
    k3 = f(temp_x,t+0.75*H, num_x);

    //k4
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + k3*H;
    }
    k4 = f(temp_x, t+H, num_x);            //what about saving this for next time step

    RK3 = (k1+ 4*k2 + k3)/6;
    RK4 = (7*k1 + 6*k2 + 8*k3 + 3*k4)/24;
    
    return RK4;
}

double RK34wAdapt(const double* x, double t, double & error, double h, int num_x)
{
    //implement RK34
    double k1, k2, k3, k4, RK3, RK4;

    //k1
    k1 = f(x,t, num_x);

    //k2
    double temp_x[num_x];
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.5*k1*h;
    }
    k2 = f(temp_x,t+0.5*h, num_x);

    //k3
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.75*k2*h;
    }
    k3 = f(temp_x,t+0.75*h, num_x);

    //k4
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + k3*h;
    }
    k4 = f(temp_x, t+h, num_x);            //what about saving this for next time step

    RK3 = (k1+ 4*k2 + k3)/6;
    RK4 = (7*k1 + 6*k2 + 8*k3 + 3*k4)/24;
    error = RK4-RK3;    //how do i implement tolerance with time stepping
    //cout << "error: " << error << endl;
    return RK4;
}

