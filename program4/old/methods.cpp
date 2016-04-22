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
double f(double x, double t)
{
    return 4*exp(0.8*t) - 0.5*x;
    //return pow(t, 4)*sin(2*t) - pow(t, 2) + 4*pow(t, 3) + (2/t)*x;
}


//forward euler
double Feuler(double x, double t)
{
    return f(x,t);
}

//backward euler
double Beuler(double x, double t)
{
    double i = f(x,t);
    return f(x+i*H,t+H);
}

//trapezoidal
double trapezoidal(double x, double t)
{
    double dxdt0 = f(x,t);
    double dxdt1 = f(x+dxdt0*H,t+H);
    return (dxdt0+dxdt1)/2;
}

double RK34woAdapt(double x, double t)
{
    //implement RK34
    double k1, k2, k3, k4, RK3, RK4;
    k1 = f(x,t);
    k2 = f(x+0.5*k1*H,t+0.5*H);
    k3 = f(x+0.75*k2*H,t+0.75*H);
    k4 = f(x+k3*H, t+H);            //what about saving this for next time step
    RK3 = (k1+ 4*k2 + k3)/6;
    RK4 = (7*k1 + 6*k2 + 8*k3 + 3*k4)/24;
    
    return RK4;
}

double RK34wAdapt(double x, double t, double & error, double h)
{
    //implement RK34
    double k1, k2, k3, k4, RK3, RK4;
    k1 = f(x,t);
    k2 = f(x+0.5*k1*h,t+0.5*h);
    k3 = f(x+0.75*k2*h,t+0.75*h);
    k4 = f(x+k3*h, t+h);            //what about saving this for next time step
    RK3 = (k1+ 4*k2 + k3)/6;
    RK4 = (7*k1 + 6*k2 + 8*k3 + 3*k4)/24;
    error = RK4-RK3;    //how do i implement tolerance with time stepping
    //cout << "error: " << error << endl;
    return RK4;
}


