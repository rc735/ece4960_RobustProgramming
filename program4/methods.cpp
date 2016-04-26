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
void f(double* phi, const double* x, double t, int num_x)
{
  	//phi[0] = 4*exp(0.8*t) - 0.5*x[0];

  	double R = 10000; //10k Ohms or 10000 Ohms
	double C = 1e-12; //1pF or 1e-12 Farads
	double i; //calculate the current based on time in 20ns periods
	//0-1ns - linear ramp-up from 0 to 0.1mA
	//1-10ns - constant current at 0.1mA
	//10-11ns - linear ramp-down from 0.1mA to 0mA
	//11-20ns - 0mA current
	double time_in_period = fmod(t,20);
	if (time_in_period <= 1.0) {
		i = time_in_period * 0.0001;
	} else if (time_in_period > 1.0 && time_in_period <= 10.0) {
		i = 0.0001;
	} else if (time_in_period > 10.0 && time_in_period <= 11.0) {
		i = (11.0-time_in_period) * 0.0001;
	} else {
		i = 0;
	}
	
	phi[0] = -(2/(C*R))*x[0] + (1/(C*R))*x[1] + i/C;
	phi[1] = (1/(C*R))*x[0] - (2/(C*R))*x[1];
    //return 4*exp(0.8*t) - 0.5*x[0];
    //return pow(t, 4)*sin(2*t) - pow(t, 2) + 4*pow(t, 3) + (2/t)*x[0];*/
}

//forward euler
void Feuler(double* phi, const double* x, double t, int num_x)
{
    f(phi, x,t, num_x);
}

//backward euler
void Beuler(double* phi, const double* x, double t, int num_x)
{
    double temp_phi[num_x];     //maybe not needed; can be replaced by temp_x
    double temp_x[num_x];
    f(temp_phi, x, t, num_x);
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + temp_phi[i]*H;
    }
    f(phi, temp_x, t+H, num_x);
}

//trapezoidal
void trapezoidal(double* phi, const double* x, double t, int num_x)
{
    double dxdt0[num_x];
    f(dxdt0, x, t, num_x);

    double temp_x[num_x];
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + dxdt0[i]*H;
    }
    double dxdt1[num_x];
    f(dxdt1, temp_x,t+H, num_x);

    for(int i = 0; i < num_x; i++)
    {
      phi[i] = (dxdt0[i] + dxdt1[i])/2.0;
    }
}

void RK34woAdapt(double* phi, const double* x, double t, int num_x)
{
    //implement RK34
    double k1[num_x], k2[num_x], k3[num_x], k4[num_x];
    double RK3[num_x], RK4[num_x];
    
    //k1
    f(k1, x, t, num_x);

    //k2
    double temp_x[num_x];
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.5*k1[i]*H;
    }
    f(k2, temp_x, t+0.5*H, num_x);

    //k3
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.75*k2[i]*H;
    }
    f(k3, temp_x, t+0.75*H, num_x);

    //k4
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + k3[i]*H;
    }
    f(k4, temp_x, t+H, num_x);            //what about saving this for next time step

    //RK3 = (k1+ 4*k2 + k3)/6;
    //RK4 = (7*k1 + 6*k2 + 8*k3 + 3*k4)/24;
    for(int i = 0; i < num_x; i++)
    {
      RK3[i] = (k1[i] + 4*k2[i] + k3[i])/6; // TODO: should be changed to /6.0 instead
      RK4[i] = (7*k1[i] + 6*k2[i] + 8*k3[i] + 3*k4[i])/24; // TODO: should be changed to /24.0
      phi[i] = RK4[i];
    }
}

void RK34wAdapt(double* phi, const double* x, double t, double & error, double h, int num_x)
{
    //implement RK34
    double k1[num_x], k2[num_x], k3[num_x], k4[num_x];
    double RK3[num_x], RK4[num_x];

    //k1
    f(k1, x, t, num_x);

    //k2
    double temp_x[num_x];
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.5*k1[i]*h;
    }
    f(k2, temp_x, t+0.5*h, num_x);

    //k3
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + 0.75*k2[i]*h;
    }
    f(k3, temp_x, t+0.75*h, num_x);

    //k4
    for(int i = 0; i < num_x; i++)
    {
      temp_x[i] = x[i] + k3[i]*h;
    }
    f(k4, temp_x, t+h, num_x);            //what about saving this for next time step

    //RK3 = (k1+ 4*k2 + k3)/6;
    //RK4 = (7*k1 + 6*k2 + 8*k3 + 3*k4)/24;
    //error = RK4-RK3;    //how do i implement tolerance with time stepping
    error = 0;
    for(int i = 0; i < num_x; i++)
    {
      RK3[i] = (k1[i] + 4*k2[i] + k3[i])/6;     //TODO change to 6.0
      RK4[i] = (7*k1[i] + 6*k2[i] + 8*k3[i] + 3*k4[i])/24;    //TODO: change to 24.0
      phi[i] = RK4[i];
      error += RK4[i] - RK3[i];
    }
    error /= (double)num_x;
}

