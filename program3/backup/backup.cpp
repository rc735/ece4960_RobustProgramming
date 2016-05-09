#include <iostream>
#include <cmath>
#include <vector>

#define SIZE_PARAM      2
#define SIZE_MEAS       10

#define PERCENT_ERR     10

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

void generateS(double* a, double* x, double S, int size)
{
  if(size <= 0)
    return 0;

  double S = a[0];
  for(int i = 1; i < size; i++)
  {
    S += a[i] * pow(x[i-1], i);
  }

}

double randError(double a, int percent)
{
  return a * (double)(1.0 - percent/100.0 + (rand() % (percent*200))/10000.0);
}

double V(double *a, double *x, double *S_measured)
{
  double V = 0;
  for(int i = 0; i < SIZE_MEAS; i++)
  {
    V += pow(generateS(a, &x[i], SIZE_PARAM) - S_measured[i], 2);
  }
  return V;
}

void newtonMethod(double *a, double *a_final, double *x,
                  double *Smeas, int paramSize)
{
  if(a == NULL || x == NULL || Smeas == NULL)
    return;

  for(int i = 0; i < paramSize; i++)
  {
    int stop = 0;
    double da = a[i]* 0.0001;
    double a1[2] = {a[0], a[1]};
    double a2[2] = {a[0], a[1]};
    double cur = a[i], next = 0;
    while (stop == 0)
    {
      cout << cur << endl;
      double dV1, dV2, d2V;
      a1[i] = cur + da;
      a2[i] = cur - da;

      //cout << a1[i] << " " << a2[i] << " " << da << endl; 
      //cout << V(a, x, Smeas) << " " << V(a1, x , Smeas) << " " << V(a2, x, Smeas) << endl;
      
      //calculate dV and d2V at initial a value
      dV1 = (V(a, x, Smeas) - V(a1, x, Smeas))/da;
      dV2 = (V(a, x, Smeas) - V(a2, x, Smeas))/da;
      d2V = (dV1 - dV2)/(2*da);
      next = cur - dV1/d2V;
      //cout << dV1 << " " << dV2 << " " << d2V << endl;
      //cout << dV1/d2V << endl;
      //cout << "cur: " << cur << " and next: " << next << endl;

      // check if under threshold
      if (abs(next-cur) < 1e-8)
      {
        a_final[i] = next;
        stop = 1;
      }
      cur = next;
      cout << endl;
    }
  }
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
  for(int i = 0; i < SIZE_MEAS; i++)
  {
    xmeas[i] = randError(x[i], PERCENT_ERR);
    Smeas[i] = randError(generateS(a,&(x[i]),SIZE_PARAM), PERCENT_ERR);
    cout << "Smeas: " << Smeas[i] << endl;
  }
  
  double guess[SIZE_PARAM] = {log(15), -0.45};
  double a_final[SIZE_PARAM];
  newtonMethod(guess, a_final, x, Smeas, SIZE_PARAM);

  return 0;
}
