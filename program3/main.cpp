#include <iostream>
#include <cmath>
#include <vector>

#include <fstream>
#include <sstream>

#include "defines.h"
#include "matrix_ops.h"
#include "param_extraction.h"

using namespace std;

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
    //cout << "x: " << x[i] << " --> xmeas: " << xmeas[i] << endl;
    //cout << "S: " << S[i] << " --> Smeas: " << Smeas[i] << endl;
  }
  //cout << generateV(a, xmeas, Smeas, SIZE_PARAM, SIZE_MEAS) << endl;
  
  double guess[SIZE_PARAM] = {log(19), -1};
  double a_final[SIZE_PARAM] = {0, 0};

  // NEWTON METHOD
  newtonMethod(guess, a_final, xmeas, Smeas, SIZE_PARAM, SIZE_MEAS);
  cout << "FINAL NEWTON: ";
  for(int i = 0; i < SIZE_PARAM; i++)
  {
    cout << a_final[i] << " ";
  }
  cout << endl << endl;

  //SECANT METHOD
  cout << "=======================================" << endl;
  secantMethod(guess, a_final, xmeas, Smeas, SIZE_PARAM, SIZE_MEAS);
  cout << "FINAL SECANT: ";
  for(int i = 0; i < SIZE_PARAM; i++)
  {
    cout << a_final[i] << " ";
  }
  cout << endl;
  //cout << "||V|| = ";
  //cout << generateV(a_final, xmeas, Smeas, SIZE_PARAM, SIZE_MEAS) << endl;

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
/*
  //---------------------------------------------------------
  // transconductance.txt
  //---------------------------------------------------------
  cout << endl;
  cout << "++++++++++++++++++ TRANSCONDUCTANCE.txt +++++++++++++++++++" << endl;

  //open the txt file
  ifstream trans_file("transconductance.txt");

  //setup the variables for parsing the txt file
  string str;
  vector<double> Vds, Vgs, Ids;
  double temp1=-1, temp2=-1, temp3=-1;
  int count = 0;

  //reading the txt file
  while(getline(trans_file, str))
  {
    stringstream convertor(str);
    convertor >> temp1 >> temp2 >> temp3;
    
    if((temp1 != -1) && (temp2 != -1) && (temp3 != -1))
    {
      Vds.push_back(temp1);
      Vgs.push_back(temp2);
      Ids.push_back(temp3);
      count++;
    }
  }

  //transfer vector values to array
  double trans_Vds[count], trans_Vgs[count], trans_Ids[count];
  for(int i = 0; i < count; i++)
  {
    //variables
    trans_Vds[i] = Vds[i];
    trans_Vgs[i] = Vgs[i];

    //Smeas
    trans_Ids[i] = Ids[i];
  }

  //                       Ids,  k, Vth
  double trans_guess[3] = {1e-7, 2.0, 1.0};     // does not work for {1e-7, 1, 1}
  double trans_final[3], trans_temp[3];

  //NEWTON METHOD - TASK 3 ------------
  cout << "=============== TASK 3 - NEWTON ================" << endl;
  newtonMethod(trans_guess, trans_final, trans_Vds, trans_Vgs, trans_Ids,
                3, count, 3);
  for(int i = 0; i < 3; i++)
  {
    trans_temp[i] = trans_final[3];
    cout << trans_final[i] << " ";
  }
  cout << endl << endl;
  
  //SECANT METHOD - TASK 3 ------------
  cout << "=============== TASK 3 =================" << endl;
  secantMethod(trans_guess, trans_final, trans_Vds, trans_Vgs, trans_Ids,
                3, count, 3);
  cout << "FINAL SECANT (TASK 3): ";
  for(int i = 0; i < 3; i++)
  {
    trans_temp[i] = trans_final[3];
    cout << trans_final[i] << " ";
  }
  cout << endl;

  cout << "||V|| = ";
  cout << generateV(trans_guess, trans_Vds, trans_Vgs, trans_Ids, 3, count, 3)
            << " --> ";
  cout << generateV(trans_final, trans_Vds, trans_Vgs, trans_Ids, 3, count, 3) << endl;

  cout << "delta_S_Is = ";
  trans_temp[0] += trans_temp[0]*PERTURB_NUM;
  cout << (generateId(trans_Vgs[0], trans_Vds[0],
                      trans_temp[0], trans_temp[1], trans_temp[2])/
          generateId(trans_Vgs[0], trans_Vds[0],
                      trans_final[0], trans_final[1], trans_final[2]))/
          (trans_temp[0]*PERTURB_NUM / trans_temp[0]) << endl;
  trans_temp[0] = trans_final[0];

  cout << "delta_S_k = ";
  trans_temp[1] += trans_temp[1]*PERTURB_NUM;
  cout << (generateId(trans_Vgs[0], trans_Vds[0],
                      trans_temp[0], trans_temp[1], trans_temp[2])/
          generateId(trans_Vgs[0], trans_Vds[0],
                      trans_final[0], trans_final[1], trans_final[2]))/
          (trans_temp[1]*PERTURB_NUM / trans_temp[1]) << endl;
  trans_temp[1] = trans_final[1];

  cout << "delta_S_Vth = ";
  trans_temp[2] += trans_temp[2]*PERTURB_NUM;
  cout << (generateId(trans_Vgs[0], trans_Vds[0],
                      trans_temp[0], trans_temp[1], trans_temp[2])/
          generateId(trans_Vgs[0], trans_Vds[0],
                      trans_final[0], trans_final[1], trans_final[2]))/
          (trans_temp[2]*PERTURB_NUM / trans_temp[2]) << endl;
  trans_temp[2] = trans_final[2];
  
**  
  //SECANT METHOD - TASK 4 ------------
  cout << endl;
  cout << "================ TASK 4 ================" << endl;
  secantMethod(trans_guess, trans_final, trans_Vds, trans_Vgs, trans_Ids,
                3, count, 4);
  cout << "FINAL SECANT (TASK 4): ";
  for(int i = 0; i < 3; i++)
  {
    trans_temp[i] = trans_final[3];
    cout << trans_final[i] << " ";
  }
  cout << endl;

  cout << "||V|| = ";
  cout << generateV(trans_guess, trans_Vds, trans_Vgs, trans_Ids, 3, count, 4)
            << " --> ";
  cout << generateV(trans_final, trans_Vds, trans_Vgs, trans_Ids, 3, count, 4) << endl;

  cout << "delta_S_Is = ";
  trans_temp[0] += trans_temp[0]*PERTURB_NUM;
  cout << (generateId(trans_Vgs[0], trans_Vds[0],
                      trans_temp[0], trans_temp[1], trans_temp[2])/
          generateId(trans_Vgs[0], trans_Vds[0],
                      trans_final[0], trans_final[1], trans_final[2]))/
          (trans_temp[0]*PERTURB_NUM / trans_temp[0]) << endl;
  trans_temp[0] = trans_final[0];

  cout << "delta_S_k = ";
  trans_temp[1] += trans_temp[1]*PERTURB_NUM;
  cout << (generateId(trans_Vgs[0], trans_Vds[0],
                      trans_temp[0], trans_temp[1], trans_temp[2])/
          generateId(trans_Vgs[0], trans_Vds[0],
                      trans_final[0], trans_final[1], trans_final[2]))/
          (trans_temp[1]*PERTURB_NUM / trans_temp[1]) << endl;
  trans_temp[1] = trans_final[1];

  cout << "delta_S_Vth = ";
  trans_temp[2] += trans_temp[2]*PERTURB_NUM;
  cout << (generateId(trans_Vgs[0], trans_Vds[0],
                      trans_temp[0], trans_temp[1], trans_temp[2])/
          generateId(trans_Vgs[0], trans_Vds[0],
                      trans_final[0], trans_final[1], trans_final[2]))/
          (trans_temp[2]*PERTURB_NUM / trans_temp[2]) << endl;
  trans_temp[2] = trans_final[2];



  //SECANT METHOD - TASK 5 ------------
  **
  cout << "=======================================" << endl;
  secantMethod(trans_guess, trans_final, trans_Vds, trans_Vgs, trans_Ids,
                3, count, 5);
  cout << "FINAL SECANT (TASK 5): ";
  for(int i = 0; i < 3; i++)
  {
    cout << trans_final[i] << " ";
  }
  cout << endl;

  cout << "||V|| = ";
  cout << generateV(trans_guess, trans_Vds, trans_Vgs, trans_Ids, 3, count, 5)
            << " --> ";
  cout << generateV(trans_final, trans_Vds, trans_Vgs, trans_Ids, 3, count, 5) << endl;
  *

  //---------- - TASK 7 ------------
  cout << endl;
  cout << "============== TASK 7 =================" << endl;

  double minV1 = 10000;
  double Is[8] = {1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5};
  double k[7]  = {1, 1.5, 2, 2.5, 3, 3.5, 4};
  double Vth[10] = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
  double params7[3], result1, params7_final[3];
  for(int i = 0; i < 8; i++)
  {
    for(int j = 0; j < 7; j++)
    {
      for(int a = 0; a < 10; a++)
      {
        params7[0] = Is[i];
        params7[1] = k[j];
        params7[2] = Vth[a];
        result1 = generateV(params7, trans_Vds, trans_Vgs, trans_Ids,
                              3, count, 3);
        if(result1 < minV1)
        {
          params7_final[0] = params7[0];
          params7_final[1] = params7[1];
          params7_final[2] = params7[2];
          minV1 = result1;
          //cout << result1 << endl;
        }
      }
    }
  }
  cout << "FULL PARAMETERS SEARCH: ";
  for(int i = 0; i < 3; i++)
  {
    cout << params7_final[i] << " ";
  }
  cout << endl;


  //---------------------------------------------------------
  // outputCharacteristics.txt
  //---------------------------------------------------------
  cout << endl << endl;
  cout << "+++++++++++++++ OUTPUTCHARACTERISTICS.txt +++++++++++++++++++" << endl;

  //open the txt file
  ifstream outputC_file("outputCharacteristics.txt");

  //reset variables for parsing txt file
  temp1=-1, temp2=-1, temp3=-1;
  count = 0;
  Vgs.clear();
  Vds.clear();
  Ids.clear();
  
  //read the txt file
  while(getline(outputC_file, str))
  {
    stringstream convertor(str);
    convertor >> temp1 >> temp2 >> temp3;
    
    if((temp1 != -1) && (temp2 != -1) && (temp3 != -1))
    {
      Vgs.push_back(temp1);
      Vds.push_back(temp2);
      Ids.push_back(temp3);
      count++;
    }
  }

  //transfer the vector to array
  double outputC_Vds[count], outputC_Vgs[count], outputC_Ids[count];
  for(int i = 0; i < count; i++)
  {
    outputC_Vds[i] = Vds[i];
    outputC_Vgs[i] = Vgs[i];
    outputC_Ids[i] = Ids[i];
    //cout << outputC_Vgs[i] << " " << outputC_Vds[i] << " " << outputC_Ids[i] << endl;
  }
  
  //                       Ids,  k, Vth
  double outputC_guess[3] = {1e-7, 2.0, 1.0};     // does not work for {1e-7, 1, 1}
  double outputC_final[3], outputC_temp[3];
  
  //SECANT METHOD - TASK 3 ------------
  cout << "=============== TASK 3 =================" << endl;
  secantMethod(outputC_guess, outputC_final, outputC_Vds, outputC_Vgs, outputC_Ids,
                3, count, 3);
  cout << "FINAL SECANT (outputC - task3): ";
  for(int i = 0; i < 3; i++)
  {
    outputC_temp[i] = outputC_final[i];
    cout << outputC_final[i] << " ";
  }
  cout << endl;

  cout << "||V|| = ";
  cout << generateV(outputC_guess, outputC_Vds, outputC_Vgs, outputC_Ids, 3, count, 3)
          << " --> ";
  cout << generateV(outputC_final, outputC_Vds, outputC_Vgs, outputC_Ids, 3, count, 3) << endl;

  cout << "||delta|| = ";
  cout << sqrt(pow(outputC_final[0]*PERTURB_NUM,2.0)/pow(outputC_final[0],2.0)
                + pow(outputC_final[1]*PERTURB_NUM,2.0)/pow(outputC_final[1],2.0)
                + pow(outputC_final[2]*PERTURB_NUM,2.0)/pow(outputC_final[2],2.0)) << endl;
  
  cout << "delta_S_Is = ";
  outputC_temp[0] += outputC_temp[0]*PERTURB_NUM;
  cout << (generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_temp[0], outputC_temp[1], outputC_temp[2])/
          generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_final[0], outputC_final[1], outputC_final[2]))/
          (outputC_temp[0]*PERTURB_NUM / outputC_temp[0]) << endl;
  outputC_temp[0] = outputC_final[0];
  cout << generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_final[0], outputC_final[1], outputC_final[2]) << endl;
  cout << outputC_temp[0] << endl;

  cout << "delta_S_k = ";
  outputC_temp[1] += outputC_temp[1]*PERTURB_NUM;
  cout << (generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_temp[0], outputC_temp[1], outputC_temp[2])/
          generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_final[0], outputC_final[1], outputC_final[2]))/
          (outputC_temp[1]*PERTURB_NUM / outputC_temp[1]) << endl;
  outputC_temp[1] = outputC_final[1];

  cout << "delta_S_Vth = ";
  outputC_temp[2] += outputC_temp[2]*PERTURB_NUM;
  cout << (generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_temp[0], outputC_temp[1], outputC_temp[2])/
          generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_final[0], outputC_final[1], outputC_final[2]))/
          (outputC_temp[2]*PERTURB_NUM / outputC_temp[2]) << endl;
  outputC_temp[2] = outputC_final[2];



  //SECANT METHOD - TASK 4 ------------
  cout << "=============== TASK 4 =================" << endl;
  secantMethod(outputC_guess, outputC_final, outputC_Vds, outputC_Vgs, outputC_Ids,
                3, count, 4);
  cout << "FINAL SECANT (outputC - task4): ";
  for(int i = 0; i < 3; i++)
  {
    cout << outputC_final[i] << " ";
  }
  cout << endl;

  cout << "||V|| = ";
  cout << generateV(outputC_guess, outputC_Vds, outputC_Vgs, outputC_Ids, 3, count, 4)
          << " --> ";
  cout << generateV(outputC_final, outputC_Vds, outputC_Vgs, outputC_Ids, 3, count, 4) << endl;
  
  cout << "delta_S_Is = ";
  outputC_temp[0] += outputC_temp[0]*PERTURB_NUM;
  cout << (generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_temp[0], outputC_temp[1], outputC_temp[2])/
          generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_final[0], outputC_final[1], outputC_final[2]))/
          (outputC_temp[0]*PERTURB_NUM / outputC_temp[0]) << endl;
  outputC_temp[0] = outputC_final[0];

  cout << "delta_S_k = ";
  outputC_temp[1] += outputC_temp[1]*PERTURB_NUM;
  cout << (generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_temp[0], outputC_temp[1], outputC_temp[2])/
          generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_final[0], outputC_final[1], outputC_final[2]))/
          (outputC_temp[1]*PERTURB_NUM / outputC_temp[1]) << endl;
  outputC_temp[1] = outputC_final[1];

  cout << "delta_S_Vth = ";
  outputC_temp[2] += outputC_temp[2]*PERTURB_NUM;
  cout << (generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_temp[0], outputC_temp[1], outputC_temp[2])/
          generateId(outputC_Vgs[0], outputC_Vds[0],
                      outputC_final[0], outputC_final[1], outputC_final[2]))/
          (outputC_temp[2]*PERTURB_NUM / outputC_temp[2]) << endl;
  outputC_temp[2] = outputC_final[2];

  

  //------------ TASK 7 ------------
  cout << endl;
  cout << "============== TASK 7 =================" << endl;

  double minV2 = 10000;
  double result2;
  for(int i = 0; i < 8; i++)
  {
    for(int j = 0; j < 7; j++)
    {
      for(int a = 0; a < 10; a++)
      {
        params7[0] = Is[i];
        params7[1] = k[j];
        params7[2] = Vth[a];
        result2 = generateV(params7, outputC_Vds, outputC_Vgs, outputC_Ids,
                              3, count, 3);
        if(result2 < minV2)
        {
          params7_final[0] = params7[0];
          params7_final[1] = params7[1];
          params7_final[2] = params7[2];
          minV2 = result2;
          //cout << result2 << " " << params7[0] << " "
          //      << params7[1] << " " << params7[2] << endl;
        }
      }
    }
  }
  cout << "FULL PARAMETERS SEARCH: ";
  for(int i = 0; i < 3; i++)
  {
    cout << params7_final[i] << " ";
  }
  cout << endl;*/

  return 0;
}
