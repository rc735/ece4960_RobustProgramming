#include "model_transistor.h"

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
                  double Is,  double k,   double Vth)
{
  return (Is*pow(log(1.0+exp(k*(Vgs-Vth)/(2.0*VT))),2.0)
          - Is*pow(log(1.0+exp((k*(Vgs-Vth)-Vds)/(2.0*VT))),2.0));
}
