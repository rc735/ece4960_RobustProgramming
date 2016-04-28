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

/**
 * generateId - utilizes the Enz-Krummenacher-Vittoz model to calculate
 *              the current, Id
 *  @param Vgs        gate-source voltage
 *  @param Vds        drain-source voltage
 *  @param Is         initial current
 *  @param k          kappa
 *  @param Vth        threshold voltage
 *  @return           the drain current
 *
 * (Source is connected to Bulk)
 */
double Id(double Vgs, double Vds,
                  double Is,  double k,   double Vth)
{
    return (Is*pow(log(1.0+exp(k*(Vgs-Vth)/(2.0*VT))),2.0)
            - Is*pow(log(1.0+exp((k*(Vgs-Vth)-Vds)/(2.0*VT))),2.0));
}

//function f for validation -- can be anything
void f(double* phi, const double* x, double t, int num_x)
{
	//task 3 and 4
	double Is = 5e-6; //5 micro amps
	double k = 0.7;	  //kappa
	double Vth = 1.5; //threshold voltage
	
	double Vdd = 5;
  	double R = 10000; //10k Ohms or 10000 Ohms
	double C = 1e-12; //1pF or 1e-12 Farads
	double i; //calculate the current based on time in 20ns periods
	//0-1ns - linear ramp-up from 0 to 0.1mA
	//1-10ns - constant current at 0.1mA
	//10-11ns - linear ramp-down from 0.1mA to 0mA
	//11-20ns - 0mA current
	double t_new = t * 1e9;	//convert time to ns
	double time_in_period = fmod(t_new,20);
	if (time_in_period <= 1.0) {
		i = time_in_period * 0.0001;
	} else if (time_in_period > 1.0 && time_in_period <= 10.0) {
		i = 0.0001;
	} else if (time_in_period > 10.0 && time_in_period <= 11.0) {
		i = (11.0-time_in_period) * 0.0001;
	} else {
		i = 0;
	}

	if (TASK_NUM == 1) {
		//task 2
		phi[0] = 4.0*exp(0.8*t) - 0.5*x[0];
	} else if (TASK_NUM == 2) {
		//task 3
		phi[0] = -x[0]/(R*C) - (x[0]-x[1])/(R*C) + i/C;
		phi[1] = -x[1]/(R*C) - (x[1]-x[0])/(R*C);
	} else if (TASK_NUM == 3) {
		//task 4
		//TODO: compute vin from i using thevenin relation
		phi[0] = -x[0]/(R*C) + i/(C);
		phi[1] = -Id(x[0],x[1],Is, k, Vth)/C - x[1]/(R*C) + Vdd/(R*C);
	} else {
		cout << "TASK_NUM invalid. Check \"define.h\" for correction." << endl;
		exit(1);
	}
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
    double error = 0;
    for(int i = 0; i < num_x; i++)
    {
      RK3[i] = (k1[i] + 4.0*k2[i] + k3[i])/6.0;
      RK4[i] = (7.0*k1[i] + 6.0*k2[i] + 8.0*k3[i] + 3.0*k4[i])/24.0;
      phi[i] = RK4[i];
      
      // error = x_3rd - x_4th
      error += pow((-5*k1[i] + 6*k2[i] + 8*k3[i] - 9*k4[i]) * H / 72.0, 2.0);
    }
    error = sqrt(error);
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
      RK3[i] = (k1[i] + 4.0*k2[i] + k3[i])/6.0;
      RK4[i] = (7.0*k1[i] + 6.0*k2[i] + 8.0*k3[i] + 3.0*k4[i])/24.0;
      phi[i] = RK4[i];

      // error = x_3rd - x_4th
      error += pow((-5*k1[i] + 6*k2[i] + 8*k3[i] - 9*k4[i]) * h / 72.0, 2.0);
    }
    error = sqrt(error);
    //cout << "error = " << error << endl;
}


double normNum(double* x, int num_x)
{
  double ans = 0;
  for(int i = 0; i < num_x; i++)
  {
    ans += pow(x[i], 2.0);
  }
  return sqrt(ans);
}
