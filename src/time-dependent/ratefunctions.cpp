#include "ratefunctions.h"


double RateFunctions::constant(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;
  double rate = params->rate;
  return rate;
}

double RateFunctions::linear(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;
  double rate = params->rate;
  double a = params->coefs[0];
  double b = params->coefs[1];
  double min = params->coefs[2];

  double y = rate * GSL_MAX(a + t * b, min);
  return y;
}

double RateFunctions::logistic(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;
  double rate = params->rate;
  double A = params->coefs[0];
  double K = params->coefs[1];
  double B = params->coefs[2];
  double nu = params->coefs[3];
  double Q = params->coefs[4];
  double M = params->coefs[5];

  double y = rate * (A + (K - A) / pow(1 + Q * exp(-B * (t-M)), 1 / nu) );
  return y;
}
