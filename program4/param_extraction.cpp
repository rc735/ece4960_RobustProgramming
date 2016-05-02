#include "param_extraction.h"

/**
 * newton_f - function used in newton method
 *    http://www.math.pitt.edu/~sussmanm/2071Spring09/lab03/index.html
 */
double newton_f(double x_i, double x_i1, double h,
                double phi_i, double phi_i1, int task)
{
  if(task == BACKWARD_EULER)
  {
    return (x_i + h*phi_i1 - x_i1);
  }
  else if(task == TRAPEZOIDAL)
  {
    return (x_i + h*(phi_i + phi_i1)/2.0 - x_i1);
  }
  else
  {
    cout << "newton_f: invalid task" << endl;
    return 0;
  }
};


/**
 * newtonMethod - for estimating the power law example using newtown method
 */
void newtonMethod(double * x, double t, double h,
                  int num_x, int task, int max_iter)
{
  if(x == NULL || max_iter < 1 || num_x < 1)
    return;

  double perturb_num  = 1e-4;
  double thresh       = 1e-7;

  double dx[num_x];
  double x_temp[num_x], x_temp2[num_x];
  double x_next[num_x];
  for(int i = 0; i < num_x; i++)
  {
    x_temp[i]  = x[i];
    x_temp2[i] = x[i];
    
    if(x[i] == 0)
    {
      dx[i] = perturb_num;
    }
    else
    {
      dx[i] = x[i] * perturb_num;
    }
    x_next[i]  = x[i] + dx[i];
  }

  int stop = 0, count;
  double delt_x[num_x];
  int topCnt = 0;

  // original phi at (t,x)
  double phi_i[num_x];
  f(phi_i, x, t, num_x);

  double phi[num_x], phi2[num_x];
  double **J    = new double*[num_x];
  double **invJ = new double*[num_x];
  for(int i = 0; i < num_x; i++)
  {
    J[i]    = new double[num_x];
    invJ[i] = new double[num_x];
  }

  while(stop == 0)
  {
    // estimating new perturb values to calculate Jacobian
    /*
    for(int i = 0; i < num_x; i++)
    {
      dx[i] = x[i] * perturb_num;

    }*/

    // generate Jacobian
    for(int i = 0; i < num_x; i++)
    {
      // for the diagonal elements
      x_temp[i]  += dx[i];
      x_temp2[i] -= dx[i];

      f(phi,  x_temp,  t+h, num_x);
      f(phi2, x_temp2, t+h, num_x);
      //J[i][i] = ((x[i] + h*phi[i] - x_temp[i])
      //          - (x[i] + h*phi2[i] - x_temp2[i])) / (2.0*dx[i]);
      J[i][i] = (newton_f(x[i], x_temp[i], h, phi_i[i], phi[i], task)
                  - newton_f(x[i], x_temp2[i], h, phi_i[i], phi2[i], task))
                  / (2.0*dx[i]);

      x_temp[i]  = x[i];
      x_temp2[i] = x[i];

      // for non-diagonal elements
      for(int j = 0; j < num_x; j++)
      {
        if(i != j)
        {
          x_temp[j]  += dx[j];
          x_temp2[j] -= dx[j];

          f(phi,  x_temp,  t+h, num_x);
          f(phi2, x_temp2, t+h, num_x);
          //J[i][j] = ((x[i] + h*phi[i] - x_temp[i])
          //          - (x[i] + h*phi2[i] - x_temp2[i])) / (2.0*dx[j]);
          J[i][j] = (newton_f(x[i], x_temp[i], h, phi_i[i], phi[i], task)
                      - newton_f(x[i], x_temp2[i], h, phi_i[i], phi2[i], task))
                      / (2.0*dx[j]);

          // resetting parameters for next iteration
          x_temp[j] = x[j];
          x_temp[j] = x[j];
          
        }
      }
    }

    // calculating the new parameters
    //calculate the inverse of Jacobian
    if(getDeterminant(J, num_x) == 0 || isnan(J[0][0]))
    {
      cout << "newton: det(J) is invalid" << endl;
      stop = 1;
    }
    else
    {
      inverse(J, num_x, invJ);

      //matrix multiplication
      count = 0;
      f(phi, x_next, t+h, num_x);
      for(int i = 0; i < num_x; i++)
      {
        delt_x[i] = 0;
        for(int j = 0; j < num_x; j++)
        {
          //delt_x[i] += invJ[i][j] * (x[j] + h*phi[j] - x_next[j]);
          delt_x[i] += invJ[i][j]
                        * newton_f(x[j], x_next[j], h, phi_i[j], phi[j], task);
        }
        if(abs(delt_x[i]) < thresh)
        {
          count++;
        }
      }

      //update new values
      for(int i = 0; i < num_x; i++)
      {
        x_next[i] -= delt_x[i];
        x_temp[i]  = x_next[i];
        x_temp2[i] = x_next[i];
      }

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

      if((count >= num_x) || (topCnt > max_iter-2))
      {
        stop = 1;
      }
    }
    topCnt++;
  }
  
  // deallocating assigned memory and updating x with x_next
  for(int i = 0; i < num_x; i++)
  {
    //update x with x_next
    x[i] = x_next[i];

    delete [] J[i];
    delete [] invJ[i];
    J[i]    = NULL;
    invJ[i] = NULL;
  }
  delete [] J;
  delete [] invJ;
  J    = NULL;
  invJ = NULL;
}


