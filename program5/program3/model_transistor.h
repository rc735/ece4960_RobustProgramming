#ifndef MODEL_TRANSISTORS_H
#define MODEL_TRANSISTORS_H

#include <iostream>
#include <cmath>

#include "defines.h"

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
double generateId(double Vgs, double Vds,
                  double Is,  double k,   double Vth);

#endif /* MODEL_TRANSISTORS_H */
