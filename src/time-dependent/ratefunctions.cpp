/*
 * =====================================================================================
 *
 *       Filename:  ratefunctions.cpp
 *
 *    Description: Contains different time-dependent functions that may be
 *                 implemented to adjust rates over time.
 *
 *        Version:  1.0
 *        Created:  08/24/2016 16:50:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thomas McDonald (), mcdonald@jimmy.harvard.edu
 *   Organization:  DFCI
 *
 * =====================================================================================
 */

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

  double y = rate * (a + t * b);
  y = GSL_MAX(y, min);
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


double RateFunctions::Gompertz(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;
  double rate = params->rate;
  double asymptote = params->coefs[0];
  double alpha = params->coefs[1];
  double beta = params->coefs[2];

  double y = rate * (asymptote + beta * exp(- alpha * t));
  return y;
}
