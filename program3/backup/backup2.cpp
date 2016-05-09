#include <iostream>
#include <cmath>
#include <vector>

#include "matrix_ops.h"

#define SIZE_PARAM      2
#define SIZE_MEAS       10

#define PERCENT_ERR     0

#define PERTURB_NUM     1e-5
#define THRESHOLD       1e-8

using namespace std;

double generateS(double* a, double* x, int size)
{
  if(size <= 0)
    return 0;

  double S = a[0];
  for(int i = 1; i < size; i++)
  {
    S += a[i] * pow(x[i-1], i);
  }
  return S;
}

void generateS(double* a, double* x, double* S, int asize, int xsize)
{
  if(asize <= 0 || xsize <= 0)
    return;

  for(int j = 0; j < xsize; j++)
  {
    S[j] = a[0];
    for(int i = 1; i < asize; i++)
    {
      S[j] += a[i] * pow(x[j], i);
    }
  }

}

double randError(double a, int percent)
{
  if(percent == 0)
  {
    return a;
  }
  return a * (double)(1.0 - percent/100.0 + (rand() % (percent*200))/10000.0);
}

double generateV(double *a, double *x, double *S_measured,
                 int paramSize, int measSize)
{
  double V = 0;
  for(int i = 0; i < measSize; i++)
  {
    V += pow(generateS(a, &x[i], paramSize) - S_measured[i], 2);
  }
  return V;
}

void newtonMethod(double *guess, double *a_final, double *x,
                  double *Smeas, int paramSize, int measSize)
{
  if(guess == NULL || x == NULL || Smeas == NULL)
    return;

  double da[paramSize];
  double a_temp[paramSize], a_temp2[paramSize];
  for(int i = 0; i < paramSize; i++)
  {
    a_final[i] = guess[i];
    a_temp[i]  = guess[i];
    a_temp2[i] = guess[i];
    da[i] = guess[i] * PERTURB_NUM;
  }

  int stop = 0;
  double dV1, dV2;
  double J[paramSize][paramSize] = {0};    //Jacobian
  double a_temp[paramSize];
  while(stop == 0)
  {
    // calculating the Jacobian
    for(int i = 0; i < paramSize; i++)
    {
      // for the diagonal elements
      a_temp[i] += da[i] * 0.5;
      dV1 = (generateV(a_temp, x, Smeas, paramSize, measSize)
                    - generateV(a_final, x, Smeas, paramSize, measSize))
                    / (0.5 * da[i]);
      a_temp[i]  -= da[i];
      a_temp2[i] -= da[i] * 0.5;
      dV2 = (generateV(a_final, x, Smeas, paramSize, measSize)
                    - generateV(a_temp, x, Smeas, paramSize, measSize))
                    / (0.5 * da[i]);
      J[i][i] = (dV1-dV2)/(da);

      // for the non-diagonal elements
      for(int j = 0; j < paramSize; i++)
      {
        if(i != j)
        {
          a_temp2[j] += da[j] * 0.5;
          dV1 = (generateV(a_temp2, x, Smeas, paramSize, measSize)
                        - generateV(a_final, x, Smeas, paramSize, measSize))
                        / (0.5 * da[i]);
          a_temp2[j] -= da[j];
          dV2 = (generateV(a_final, x, Smeas, paramSize, measSize)
                        - generateV(a_temp2, x, Smeas, paramSize, measSize))
                        / (0.5 * da[i]);
          J[j][i] = (dV1-dV2)/(da);
          a_temp2[j] += da[j] * 0.5;
        }
      }
      a_temp[i]  += da[i] * 0.5;
        a_temp2[i] += da[i] * 0.5;
    }

    // calculating the delta for parameter


  }
/*
  //copy guess into own variable
  double a_right[paramSize];
  for(int i = 0; i < paramSize; i++)
  {
    a_final[i] = guess[i];
    a_right[i] = guess[i];
  }
  
  // iterate through every parameter and find optimized value
  double V_curr, V_prev, dV, da;
  for(int i = 0; i < paramSize; i++)
  {
    da = a_final[i] * 1e-5;
    a_right[i] += da;
    V_prev = 0;
   
    while(abs(da) > 1e-8)
    {
      cout << da << " ---------- ";
      //first derivative
      V_curr = generateV(a_final, x, Smeas, paramSize, measSize);


      if(V_prev == 0 )
      {
        cout << "CALC V_PREV ";
        V_prev = generateV(a_right, x, Smeas, paramSize, measSize);
      }
      dV = (V_curr - V_prev)/(da);
      cout << V_curr << " " << V_prev << " " << dV << " ";

      //next value
      da = 0 - V_curr/dV;
      if(da > 0.5)
      {
        da = 0.5;
      }
      else if(da < -0.5)
      {
        da = -0.5;
      }
      a_final[i] = a_final[i] + da;
      V_prev = V_curr;

      cout << endl;
    }
    a_right[i] = a_final[i];
    cout << "============" << a_final[i] << endl;
  }
  */
}


/*
void secantMethod(double *a, double *x, double *Smeas, int size)
{
  if(a == NULL || x == NULL || Smeas == NULL)
    return;

  for(int i = 0; i < SIZE_MEAS; i++)
  {
    int stop = 0;
    double da = a[i]* 0.000001;
    double a1[2] = {a[0], a[1]};
    double a2[2] = {a[0], a[1]};
    double cur = a[i], next = 0;
    while (stop == 0)
    {
      double dV1, dV2, d2V;
      a1[i] = cur + da;
      a2[i] = cur - da;
      
      //calculate dV and d2V at initial a value
      dV1 = (V(a, x, S_measured) - V(a1, x, S_measured))/da;
      dV2 = (V(a, x, S_measured) - V(a2, x, S_measured))/da;
      d2V = (dV1 - dV2)/(2*da);
      next = cur - dV1/d2V;
      cout << "cur: " << cur << " and next: " << next << endl;

      // check if under threshold
      if (abs(next-cur) < 1e-8)
      {
        b[i] = next;
        stop = 1;
      }
      cur = next;
    }
  }
}*/



int main(void)
{
  double a[SIZE_PARAM] = {log(10), -0.5};
  double x[SIZE_MEAS];
  //cout << log(10) << endl;

  for(int i = 0; i < SIZE_MEAS; i++)
  {
    x[i] = log(i+1);
  }

  double Smeas[SIZE_MEAS];
  double xmeas[SIZE_MEAS];
  generateS(a, x, Smeas, SIZE_PARAM, SIZE_MEAS);

  // introduce random error
  for(int i = 0; i < SIZE_MEAS; i++)
  {
    xmeas[i] = randError(x[i], PERCENT_ERR);
    Smeas[i] = randError(Smeas[i], PERCENT_ERR);
    //Smeas[i] = randError(generateS(a, &x[i], SIZE_PARAM), PERCENT_ERR);
    cout << "Smeas: " << Smeas[i] << endl;
  }
  
  double guess[SIZE_PARAM] = {log(10), -0.5};
  double a_final[SIZE_PARAM];
  newtonMethod(guess, a_final, x, Smeas, SIZE_PARAM, SIZE_MEAS);
  for(int i = 0; i < SIZE_PARAM; i++)
  {
    cout << a_final[i] << " ";
  }
  cout << endl;

  return 0;
}
