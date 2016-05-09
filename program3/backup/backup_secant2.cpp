#include <iostream>
#include <cmath>
#include <vector>

#include "matrix_ops.h"
#include "defines.h"

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
    S[j] = generateS(a, &(x[j]),asize);
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
  double a_temp[paramSize], a_temp2[paramSize], a_temp3[paramSize];
  for(int i = 0; i < paramSize; i++)
  {
    a_final[i] = guess[i];
    a_temp[i]  = guess[i];
    a_temp2[i] = guess[i];
    a_temp3[i] = guess[i];
  }

  int stop = 0, count;
  double dV1, dV2;
  double J[SIZE_PARAM][SIZE_PARAM], invJ[SIZE_PARAM][SIZE_PARAM];    //Jacobian
  double delt_a[paramSize];
  int topCnt = 0;
  while(stop == 0)
  {
    // calculating the Jacobian
    //cout << a_final[0] << endl;
    for(int i = 0; i < paramSize; i++)
    {
      // for the diagonal elements
      //cout << a_final[i] << endl;
      da[i] = a_final[i] * PERTURB_NUM;
      a_temp[i] += da[i] * 0.5;
      
      dV1 = (generateV(a_temp, x, Smeas, paramSize, measSize)
                    - generateV(a_final, x, Smeas, paramSize, measSize))
                    / (0.5 * da[i]);
      a_temp[i]  -= da[i];
      a_temp2[i] += da[i] * 0.5;
      //dV2 = (generateV(a_final, x, Smeas, paramSize, measSize)
      //              - generateV(a_temp, x, Smeas, paramSize, measSize))
      dV2 = (generateV(a_final, x, Smeas, paramSize, measSize)
                    - generateV(a_temp, x, Smeas, paramSize, measSize))
                    / (0.5 * da[i]);
      J[i][i] = (dV1-dV2)/(da[i]);
      cout << (dV1) << " " << (dV2) << " " << (dV1-dV2) << "        " << (da[i]) << endl;

      // for the non-diagonal elements
      for(int j = 0; j < paramSize; j++)
      {
        if(i != j)
        {
          a_temp2[j] += da[j] * 0.5 ;
          a_temp3[j] += da[j] * 0.5;
          dV2 = (generateV(a_temp2, x, Smeas, paramSize, measSize)
                        - generateV(a_temp3, x, Smeas, paramSize, measSize))
                        / (0.5 * da[i]);
          a_temp2[j] -= da[j] ;
          a_temp3[j] -= da[j];
          dV1 = (generateV(a_temp3, x, Smeas, paramSize, measSize)
                        - generateV(a_temp2, x, Smeas, paramSize, measSize))
                        / (0.5 * da[i]);
          J[i][j] = (dV2-dV1)/(da[j]);
          //cout << (dV1) << " " << (dV2) << " " << (dV1-dV2) << "        " << (da[j]) << endl;

          //reset values back to a_final
          a_temp2[j] = a_final[j];
          a_temp3[j] = a_final[j];
        }
      }
      a_temp[i]  = a_final[i];
      a_temp2[i] = a_final[i];
    }
    cout << J[0][0] << " " << J[0][1] << " " << J[1][0] << " " << J[1][1] << endl;
    // calculating the new parameters
    //calculate the inverse of Jacobian
    if(getDeterminant(J, paramSize) == 0)
    {
      cout << "newton: jacobian is 0" << endl;
      stop = 1;
    }
    else{
    inverse(J, paramSize, invJ);

    //matrix multiplication
    count = 0;
    for(int i = 0; i < paramSize; i++)
    {
      delt_a[i] = 0;
      for(int j = 0; j < paramSize; j++)
      {
        delt_a[i] += invJ[i][j] * a_final[i];
        //cout << invJ[i][j] << " ";
      }
      a_final[i] -= delt_a[i];
      a_temp[i]  = a_final[i];
      a_temp2[i] = a_final[i];
      a_temp3[i] = a_final[i];

      if(abs(delt_a[i]) < 1e-7)//THRESHOLD)
      {
        count++;
      }
      cout <<  delt_a[i] << "     ";
    }
    cout << endl << endl;

    if(count >= paramSize)
    {
      stop = 1;
    }

  }
    topCnt++;
  }
  cout << "newton: number of iterations = " << topCnt << endl;
}



void secantMethod(double *guess, double *a_final, double *x,
                  double *Smeas, int paramSize, int measSize)
{
  if(a_final == NULL || x == NULL || Smeas == NULL)
    return;

  // initializing variables and assigning values
  double da[paramSize];
  double a_temp1[paramSize], a_temp2[paramSize];
  double dV_curr[paramSize], dV_prev[paramSize],
          delt_a[paramSize], a_prev[paramSize];
  for(int i = 0; i < paramSize; i++)
  {
    a_final[i] = guess[i];
    a_temp1[i] = guess[i];
    a_temp2[i] = guess[i];
    
    da[i] = a_final[i] * PERTURB_NUM;
    a_prev[i] = a_final[i] + da[i];
  }

  // generating prev guess before initial guess
  for(int i = 0; i < paramSize; i++)
  {
    a_prev[i] += da[i];
    dV_prev[i] = (generateV(a_prev, x, Smeas, paramSize, measSize)
                    - generateV(a_final, x, Smeas, paramSize, measSize))
                    / (da[i]);
    cout << generateV(a_final, x, Smeas, paramSize, measSize) << " "
          << generateV(a_prev, x, Smeas, paramSize, measSize) << endl;
    a_prev[i] = a_final[i] + da[i];
  }
  cout << dV_prev[0] << " " << dV_prev[1] << endl << endl;;

  
  // Secant Method - iterative parameter estimation
  int stop = 0, count;
  int topCnt = 0;
  double min, V_temp, delt_a_final[paramSize];
  while(stop == 0)
  {
    count = 0;
    for(int i = 0; i < paramSize; i++)
    {
      cout << "++++++++++++ " << i << endl;
      //calculate slope
      da[i] = a_final[i] * PERTURB_NUM;
      a_temp1[i] += da[i] ;
      dV_curr[i] = (generateV(a_temp1, x, Smeas, paramSize, measSize)
                    - generateV(a_final, x, Smeas, paramSize, measSize))
                    / (da[i]);
      a_temp1[i] = a_final[i];
      //cout << generateV(a_final, x, Smeas, paramSize, measSize) << " "
      //      << generateV(a_temp1, x, Smeas, paramSize, measSize) << endl;

      //calculate updated parameters (delta_a)
      delt_a[i] = dV_curr[i] * (a_final[i] - a_prev[i])
                              / (dV_curr[i] - dV_prev[i]);

      //linear search
      min = 10000;
      for(double j = 1.0; j >= 0.25; j -= 0.25)
      {
        a_temp2[i] -= j * delt_a[i];
        V_temp = generateV(a_temp2, x, Smeas, paramSize, measSize);
        a_temp2[i] = a_final[i];
        cout << V_temp << " ";
        if(abs(V_temp) < abs(min))
        {
          min = V_temp;
          delt_a_final[i] = j * delt_a[i];
        }
      }
      cout << endl;
      //cout << delt_a_final[i] << endl;

      if(abs(delt_a_final[i]) < THRESHOLD)
      {
        count++;
      }
    }
    cout << "SLOPE:        " << dV_curr[0] << " " << dV_curr[1] << endl;
    cout << "PREV:         " << dV_prev[0] << " " << dV_prev[1] << endl;
    cout << "DELT_A:       " << delt_a_final[0] << " " << delt_a_final[1] << endl;

    cout << "FINAL PARAMS: ";
    for(int i = 0; i < paramSize; i++)
    {
      a_prev[i] = a_final[i];
      a_final[i] -= delt_a_final[i];
      a_temp1[i] = a_final[i];
      dV_prev[i] = dV_curr[i];

      cout << a_final[i] << " ";
    }
    cout << endl << endl;

    if(count >= paramSize)
    {
      stop = 1;
    }
    topCnt++;
  }
  cout << "secant: number of iterations = " << topCnt << endl;
}



int main(void)
{
  double a[SIZE_PARAM] = {log(10), -0.5};
  double x[SIZE_MEAS];

  for(int i = 0; i < SIZE_MEAS; i++)
  {
    x[i] = log(i+2);
  }

  double Smeas[SIZE_MEAS], S[SIZE_MEAS];
  double xmeas[SIZE_MEAS];
  generateS(a, x, S, SIZE_PARAM, SIZE_MEAS);

  // introduce random error
  for(int i = 0; i < SIZE_MEAS; i++)
  {
    xmeas[i] = randError(x[i], PERCENT_ERR);
    Smeas[i] = randError(S[i], PERCENT_ERR);
    //Smeas[i] = randError(generateS(a, &x[i], SIZE_PARAM), PERCENT_ERR);
    cout << "x: " << x[i] << " --> xmeas: " << xmeas[i] << endl;
    cout << "S: " << S[i] << " --> Smeas: " << Smeas[i] << endl;
  }
  cout << generateV(a, xmeas, Smeas, SIZE_PARAM, SIZE_MEAS) << endl;
  
  double guess[SIZE_PARAM] = {log(14), -0.7};
  double a_final[SIZE_PARAM];
  //newtonMethod(guess, a_final, x, Smeas, SIZE_PARAM, SIZE_MEAS);
  cout << "FINAL NEWTON: ";
  for(int i = 0; i < SIZE_PARAM; i++)
  {
    cout << a_final[i] << " ";
  }
  cout << endl;

  secantMethod(guess, a_final, xmeas, Smeas, SIZE_PARAM, SIZE_MEAS);
  cout << "FINAL SECANT: ";
  for(int i = 0; i < SIZE_PARAM; i++)
  {
    cout << a_final[i] << " ";
  }
  cout << endl;
  cout << generateV(a_final, xmeas, Smeas, SIZE_PARAM, SIZE_MEAS) << endl;

/**
#if SIZE_PARAM==2
#elif SIZE_PARAM==3
  double test[3][3] = {{1, 2, 3}, {0, 4, 5}, {1, 0, 6}};
  double inv[3][3];
  inverse(test, 3, inv);
  cout << inv[0][0] << " " << inv[0][1] << " " << inv[0][2] << endl;
  cout << inv[1][0] << " " << inv[1][1] << " " << inv[1][2] << endl;
  cout << inv[2][0] << " " << inv[2][1] << " " << inv[2][2] << endl;
#endif
*/

  return 0;
}
