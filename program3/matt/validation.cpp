//
//  validation.cpp
//  Program3
//
//  Created by Matthew Magaldi on 3/19/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#include "validation.hpp"
#include "matrixOps.hpp"
#include <cmath>
#include <iostream>

using namespace std;

/*
 Generate S value function
 INPUTS:
    *a - pointer to parameter array
    *x - pointer to x array
    size - number of parameters
    index - index of the S_model in terms of which x data point
    *S_model - pointer to the S model array
 
Constraint is that number of parameters must be 1 or greater
 */
void generateS(double *a, double x, int size, int index, double *S_model)
{
    if (size < 1) {
        return;
    }
    
    double S = a[0];
    for (int i = 1; i < size; i++) {
        S += a[i]*pow(x, i);
    }
    S_model[index] = S;
}


/*
 Generic V least square fitting
 INPUTS:   
    *a - pointer to parameter array
    *x - pointer to x array
    *S_measured - pointer to the measurement data array
 */
double V(double *a, double *x, double *S_measured)
{
    double V = 0;
    double S_model[10];
    int size = 2;
    for (int i = 0; i < 10; i++) {
        generateS(a, x[i], size, i, S_model);
        V += pow(S_model[i]-S_measured[i], 2);
    }
    return V;
}

/*
 Computes the Hessian or second partial derivative matrix numerically
 INPUTS:
    **hessian - pointer to the hessian matrix stored in memory
    *a - pointer to parameter array
    *x - pointer to x array
    *S_measured - pointer to data array
    size - number of parameters
 */
void computeHessianMatrix(double **hessian, double * a, double * x, double * S_measured, double * derivative, int size)
{
    //know rank of hessian matrix to be size^2
    double t = 0.001;
    double da[2] = {a[0]*t,a[1]*t};     //this is delta of all parameters - need to generalize this
    double a1[2] = {a[0], a[1]};        //this is parameters with delta for first parameter
    double a2[2] = {a[0], a[1]};        //this is parameters with delta for second
    double a3[2] = {a[0], a[1]};
    double a4[2] = {a[0], a[1]};
    double a5[2] = {a[0], a[1]};
    double a6[2] = {a[0], a[1]};
    double a7[2] = {a[0], a[1]};
    double dVda1, dVda2, d2Vda1, d2Vda2, d2Vda12;
    
    //hardcode a 2-dimensional finite difference hessian matrix to test
    a2[0] += da[0]; //(x+h,y)
    a3[0] -= da[0]; //(x-h,y)
    a4[1] += da[1]; //(x,y+k)
    a5[1] -= da[1]; //(x,y-k)
    a6[0] += da[0];
    a6[1] += da[1]; //(x+h,y+k)
    a7[0] -= da[0];
    a7[1] -= da[1]; //(x-h,y-k)
    
    dVda1 = (V(a2,x,S_measured) - V(a3, x, S_measured))/(2*da[0]);
    dVda2 = (V(a4,x,S_measured) - V(a5, x, S_measured))/(2*da[1]);
    
    d2Vda1 = (V(a2,x,S_measured) - 2*V(a1, x, S_measured) + V(a3, x, S_measured))/(pow(da[0], 2));
    d2Vda2 = (V(a4,x,S_measured) - 2*V(a1, x, S_measured) + V(a5, x, S_measured))/(pow(da[1], 2));
    
    d2Vda12 = (V(a6,x,S_measured) - V(a2, x, S_measured) - V(a4, x, S_measured) + 2*V(a1, x, S_measured) -
               V(a3, x, S_measured) - V(a5, x, S_measured) + V(a7, x, S_measured));
    
    derivative[0] = dVda1;
    derivative[1] = dVda2;
    hessian[0][0] = d2Vda1;
    hessian[0][1] = d2Vda12;
    hessian[1][0] = d2Vda12;
    hessian[1][1] = d2Vda2;
    
    /*
    for (int i=0; i < size; i++) {
        a1[i] += da[i];
        dVda1_1 = (V(a,x,S_measured)-V(a1, x, S_measured))/da[i];
        derivative[i] = dVda1_1;
        for (int j = 0; j < size; j++) {
            if (i ==j) {
                a2[j] += 2*da[j];
                dVda1_2 = (V(a1,x,S_measured)-V(a2, x, S_measured))/da[j];
                d2V = (dVda1_1-dVda1_2)/da[j];
            } else {
                a2[j] += da[j];
                aT[i] += da[i];
                aT[j] += da[j];
                dVda1_2 = (V(a2,x,S_measured)-V(aT, x, S_measured))/da[i];
                d2V = (dVda1_1-dVda1_2)/da[j];
            }
            a2[j] = a[j];
            aT[i] = a[i];
            aT[j] = a[j];
            
            //update hessian matrix with scalar estimate parameter
            hessian[i][j] = d2V;
        }
        a1[i] = a[i];
    }
     */
}


/*
 Generic Newton's method
 INPUTS:
    *a - pointer to initial parameter array guess
    *b - pointer to final parameter array (result output)
    *x - pointer to measurement x values
    *S_measured - pointer to measurement values
    *size - number of parameters
 
 Constraint are:
    -   that x and S_measured array must be of the same size
    -   that a and b (input and output parameter arrays must be of the same size indicated by 'size'
 */
double * newtonsMethod(double *a, double *b, double *x, double *S_measured, int size)
{
    // initialize loop variables
    int stop = 0;
    double da[size];
    double a_old[size];
    double a_new[size];
    double derivative[2];
    for (int i = 0; i < size; i++) {
        da[i] = a[i];
        b[i] = a[i];
        derivative[i] = 0;
    }
    
    //dynamically allocate 2-d hessian and inverse matrix in memory for iteration use
    double ** hessian = 0;
    double ** inverse = 0;
    hessian = new double*[size];
    inverse = new double*[size];
    for (int i = 0; i < size ; i++) {
        hessian[i] = new double[size];
        inverse[i] = new double[size];
    }
    
    // iterative newtons method
    while (stop == 0) {
        //go on with newton's method - use b an most up to date parameter array
        computeHessianMatrix(hessian, b, x, S_measured, derivative, size);
        MatrixInversion(hessian, size, inverse);
        MatrixVectorMultiply(inverse, derivative, da, size);
        
        //now we perform a line search using da and a_old to find minimization linear search variable t
        double t = 2;
        double min_t = 1;
        double temp[2];
        double V_new, V_old;
        V_old = V(a, x, S_measured);
        for (int i = 0; i < 100; i++) {
            t = t/2;
            temp[0] = b[0] - t*da[0];
            temp[1] = b[1] - t*da[1];
            V_new = V(temp,x,S_measured);
            if (abs(V_new) < abs(V_old)) {
                //find the t that minimizes V
                min_t = t;
                V_old = V_new;
            }
        }
        
        // if step is small enough, we can quit
        double norm;
        for (int i=0; i < size; i++) {
            norm += pow((min_t*da[i]),2);
        }
        norm = sqrt(norm);
        if (norm < 1e-5) {
            stop = 1;
        }

        int count = 0;
        for(int i = 0; i < size; i++)
        {
          if(da[i] < 1e-8)
          {
            count++;
          }
        }
        if(count >= size)
          stop = 1;
        
        // otherwise continue and set a_new = a_old - t*da
        for (int i = 0; i < size; i++) {
            b[i] = b[i] - min_t*da[i];
        }
        cout << "t is: " << min_t << " and new guess is: " << b[0] << ", " << b[1] << endl;
    }
    
    //clean up memory
    for (int i = 0; i < size; i++) {
        delete [] inverse[i];
        delete [] hessian[i];
    }
    delete [] inverse;
    delete [] hessian;
    hessian = 0, inverse = 0;
    
    return b;
}



/*
 Validation of parameter extraction program for y=c0*(x^m)
 */
void validityTest()
{
    double S_model[10];
    double S_measured[10];
    double a[2] = {log(10), -0.5};  //parameters
    double b[2];
    double x[10];                   //x values
    double x_measured[10];
    int size = 2;
    double noise = 0;
    
    //Generate 10 samples of S_measured with random noises between 10-20%
    for (int i = 0; i < 10; i++) {
        x[i] = log(i+1);
    }
    
    for (int i = 0; i < 10; i++) {
        generateS(a, x[i], size, i, S_model);
        S_measured[i] = randError(S_model[i], (int)noise);
        //x_measured[i] = randError(x[i], noise);
    }
    
    double guess[2] = {log(12), -0.53};
    newtonsMethod(guess, b, x, S_measured, size);
    
    cout << "b[0] is : " << b[0] << " and b[1] is: " << b[1] << endl;
    
}


/*
 Adds random noise of noise% to x
 */
double randError(double a, int percent)
{
    if(percent == 0)
        return a;
    return a * (double) (1.0 - percent/100.0 + (random() % (percent*200))/10000.0);
}


