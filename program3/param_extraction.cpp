#include "param_extraction.h"


/**
 * generateS - for the power law; calculating one value of S if
 *              S is a polynomial expression
 *
 *  @param a      parameters to calcuate S
 *  @param x      variables to calculate S
 *  @param size   number of parameter and variable elements
 *  @return       the resulting S based on S = a0 + a1*x1 + ...
 */
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

/**
 * generateS - for power law; generating an array of S values
 *  
 *  @param a      parameters to calculate S
 *  @param x      variables to calculate S
 *  @param S      results stored in this array; S
 *  @param asize  number of parameters
 *  @param xsize  number of variables
 */
void generateS(double* a, double* x, double* S, int asize, int xsize)
{
  if(asize <= 0 || xsize <= 0)
    return;

  for(int j = 0; j < xsize; j++)
  {
    S[j] = generateS(a, &(x[j]),asize);
  }

}

/**
 * randError - calculates a random error based on the percent
 *
 *  @param a          the value to be calculated with the error
 *  @param percent    the percent error to be multiplied to the value
 *  @result           outputs the new value with random error
 */
double randError(double a, int percent)
{
  if(percent == 0)
  {
    return a;
  }
  return a * (double)(1.0 - percent/100.0 + (rand() % (percent*200))/10000.0);
}

/**
 *generateV - calculates the sum of the squared error for all variables
 *
 *  @param a            the parameters to calculate V
 *  @param x            the variables to calculate V
 *  @param S_measured   the measured S values
 *  @param paramSize    number of parameter elements
 *  @param measSize     number of measurements of S
 *  @return             the squared error
 */
double generateV(double *a, double *x, double *S_measured,
                 int paramSize, int measSize)
{
  double V = 0;
  for(int i = 0; i < measSize; i++)
  {
    V += pow(generateS(a, &x[i], paramSize) - S_measured[i], 2.0);
  }
  return V;
}

/**
 *generateV - calculates the sum of the squared error for all variables
 *
 *  @param a            the parameters to calculate V
 *  @param Vds          the variables to calculate V
 *  @param Vgd          the varaibles to calculate V
 *  @param S_measured   the measured S values
 *  @param paramSize    number of parameter elements
 *  @param measSize     number of measurements of S
 *  @param task         for different task in program assignment
 *                          3 --> S_model = Id
 *                          4 --> S_model = Id / Id_meas
 *                          5 --> S_model = log10(Id)
 *  @return             the squared error
 */
double generateV(double *a, double *Vds, double *Vgs, double *S_measured,
                 int paramSize, int measSize, int task)
{
  double V = 0;

  if(task == 3)
  {
    for(int i = 0; i < measSize; i++)
    {
      //Vgs, Vds, Is, k, Vth
      V += pow(generateId(Vgs[i], Vds[i], a[0], a[1], a[2]) - S_measured[i], 2.0);
    }
  }
  else if(task == 4)
  {
    for(int i = 0; i < measSize; i++)
    {
      //Vgs, Vds, Is, k, Vth
      V += pow(generateId(Vgs[i], Vds[i], a[0], a[1], a[2])/S_measured[i] - 1, 2.0);
    }
  }
  else if(task == 5)
  {
    for(int i = 0; i < measSize; i++)
    {
      //Vgs, Vds, Is, k, Vth
      V += pow(log10(generateId(Vgs[i], Vds[i], a[0], a[1], a[2])) 
                - log10(S_measured[i]), 2.0);
      //if(generateId(Vgs[i], Vds[i], a[0], a[1], a[2]) <= 0)
      //  cout << "asdfasdfasdf " << endl;
    }
  }
  return V;
}

/**
 * newtonMethod - for estimating the power law example using newtown method
 */
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
  //double dV1, dV2;
  //double H[SIZE_PARAM][SIZE_PARAM], invH[SIZE_PARAM][SIZE_PARAM];    //Hessian
  double delt_a[paramSize];
  int topCnt = 0;

  double **H    = new double*[paramSize];
  double **invH = new double*[paramSize];
  double J[paramSize];
  for(int i = 0; i < paramSize; i++)
  {
    H[i]    = new double[paramSize];
    invH[i] = new double[paramSize];
  }

  while(stop == 0)
  {
    // estimating new perturb values to calculate Hessian
    for(int i = 0; i < paramSize; i++)
    {
      da[i] = a_final[i] * PERTURB_NUM;
    }

    // calculating the Hessian using central difference (finite difference)
    // http://www.math.unl.edu/~s-bbockel1/833-notes/node23.html
    for(int i = 0; i < paramSize; i++)
    {
      // for the diagonal elements
      a_temp[i]  += da[i];
      a_temp2[i] -= da[i];
      H[i][i] = (generateV(a_temp2, x, Smeas, paramSize, measSize)
                  - 2.0*generateV(a_final, x, Smeas, paramSize, measSize)
                  + generateV(a_temp, x, Smeas, paramSize, measSize))
                  / (pow(da[i],2.0));
      J[i] = (generateV(a_temp, x, Smeas, paramSize, measSize)
                - generateV(a_temp2, x, Smeas, paramSize, measSize)) / (2.0*da[i]);
      a_temp[i]  = a_final[i];
      a_temp2[i] = a_final[i];

      // for non-diagonal elements
      for(int j = 0; j < paramSize; j++)
      {
        if(i != j)
        {
          // f(x+h, y+h)
          a_temp[i] += da[i];
          a_temp[j] += da[j];
          H[i][j] = generateV(a_temp, x, Smeas, paramSize, measSize);

          // f(x+h, y-h)
          a_temp[j] -= da[j] * 2.0;
          H[i][j] -= generateV(a_temp, x, Smeas, paramSize, measSize);

          // f(x-h, y-h)
          a_temp[i] -= da[i] * 2.0;
          H[i][j] += generateV(a_temp, x, Smeas, paramSize, measSize);

          // f(x-h, y+h)
          a_temp[j] += da[j] * 2.0;
          H[i][j] -= generateV(a_temp, x, Smeas, paramSize, measSize);

          // 4h^2
          H[i][j] /= (4.0 * da[i] * da[j]);

          // resetting parameters for next iteration
          a_temp[i] = a_final[i];
          a_temp[j] = a_final[j];
          
        }
      }
    }
    //cout << "HESSIAN: ";
    //cout << H[0][0] << " " << H[0][1] << " " << H[1][0] << " " << H[1][1] << endl;

    // calculating the new parameters
    //calculate the inverse of Hessian
    if(getDeterminant(H, paramSize) == 0 || isnan(H[0][0]))
    {
      cout << "newton: det(H) is invalid" << endl;
      stop = 1;
    }
    else{
      inverse(H, paramSize, invH);
      //cout << "INVERSE: ";
      //cout << invH[0][0] << " " << invH[0][1] << " " 
      //      << invH[1][0] << " " << invH[1][1] << endl;
      //cout << "JACOBIAN: ";
      //cout << J[0] << " " << J[1] << endl;

      //matrix multiplication
      count = 0;
      for(int i = 0; i < paramSize; i++)
      {
        delt_a[i] = 0;
        for(int j = 0; j < paramSize; j++)
        {
          delt_a[i] += invH[i][j] * J[j];
        }
        a_final[i] -= delt_a[i];
        a_temp[i]  = a_final[i];
        a_temp2[i] = a_final[i];
      
        //linear search
        /**
        min = 10000;
        for(double j = 1.0; j >= 0.25; j -= 0.25)
        {
          a_temp1[i] -= j * delt_a[i];
          V_temp = generateV(a_temp1, x, Smeas, paramSize, measSize);
          a_temp1[i] = a_final[i];
          cout << V_temp << " ";
          if(abs(V_temp) < abs(min))
          {
            min = V_temp;
            delt_a_final[i] = j * delt_a[i];
          }
        }*/

        if(abs(delt_a[i]) < THRESHOLD)
        {
          count++;
        }
      }
      //cout << "PARAMS: " << a_final[0] << " " << a_final[1];
      //cout << endl << endl;

      if(count >= paramSize)
      {
        stop = 1;
      }
    }
    topCnt++;
  }
  cout << "newton: number of iterations = " << topCnt << endl;
  
  // deallocating assigned memory
  for(int i = 0; i < paramSize; i++)
  {
    delete [] H[i];
    delete [] invH[i];
  }
  delete [] H;
  delete [] invH;
  H    = NULL;
  invH = NULL;
}


/**
 * secantMethod - for estimating the power law example using secant method
 *                slope is estimated using the equation in:
 *                https://en.wikipedia.org/wiki/Secant_method
 */
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
    //cout << generateV(a_final, x, Smeas, paramSize, measSize) << " "
    //      << generateV(a_prev, x, Smeas, paramSize, measSize) << endl;
    a_prev[i] = a_final[i] + da[i];
  }
  //cout << dV_prev[0] << " " << dV_prev[1] << endl << endl;;

  
  // Secant Method - iterative parameter estimation
  int stop = 0, count;
  int topCnt = 0;
  double delt_a_final[paramSize];
//  double min, V_temp          // for linear search
  while(stop == 0)
  {
    count = 0;
    for(int i = 0; i < paramSize; i++)
    {
      //cout << "++++++++++++ " << i << endl;

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
      /*
      min = 10000;
      for(double j = 1.0; j >= 0.25; j -= 0.25)
      {
        a_temp1[i] -= j * delt_a[i];
        V_temp = generateV(a_temp1, x, Smeas, paramSize, measSize);
        a_temp1[i] = a_final[i];
        cout << V_temp << " ";
        if(abs(V_temp) < abs(min))
        {
          min = V_temp;
          delt_a_final[i] = j * delt_a[i];
        }
      }
      cout << endl;*/

      //store delta_a results into final
      delt_a_final[i] = delt_a[i];

      if(abs(delt_a_final[i]) < THRESHOLD)
      {
        count++;
      }
    }

    /*
    cout << "SLOPE:        " << dV_curr[0] << " " << dV_curr[1] << endl;
    cout << "PREV:         " << dV_prev[0] << " " << dV_prev[1] << endl;
    cout << "DELT_A:       " << delt_a_final[0] << " " << delt_a_final[1] << endl;
    cout << "FINAL PARAMS: ";*/

    //resetting values and updating parameters for next iteration
    for(int i = 0; i < paramSize; i++)
    {
      a_prev[i]  = a_final[i];
      a_final[i] -= delt_a_final[i];
      a_temp1[i] = a_final[i];
      dV_prev[i] = dV_curr[i];

      //cout << a_final[i] << " ";
    }
    //cout << endl << endl;

    //checking if all parameters reach threshold values
    if(count >= paramSize)
    {
      stop = 1;
    }
    topCnt++;
  }
  cout << "SECANT: num of iterations = " << topCnt << endl;
}


