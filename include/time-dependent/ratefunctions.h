//=======================================================
// include guard
#ifndef __RATEFUNCTIONS_TD_H_INCLUDED__
#define __RATEFUNCTIONS_TD_H_INCLUDED__

#include "globalstructs.h"
#include <gsl/gsl_math.h>


class RateFunctions
{
public:
  static double constant(double t, void *params);
  static double linear(double t, void *params);
  static double logistic(double t, void *params);
  //double polynom(double t, void *params);
  //double piecewise(double t, void *params);
  //double custom(double t, void *params);
};

#endif // __RATEFUNCTIONS_TD_H_INCLUDED__
