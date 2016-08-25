/*
 * =====================================================================================
 *
 *       Filename:  rvfunctions.h
 *
 *    Description: Header for functions for generating random variables
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

#ifndef __RVFUNCTIONS_TD_H_INCLUDED__
#define __RVFUNCTIONS_TD_H_INCLUDED__

#include <gsl/gsl_randist.h>
#include <cmath>
#include "globalstructs.h"


extern GlobalParameters gp;

double GenerateFitness(FitnessParameters fit_params);
double GenerateMutationProb(MutationParameters mut_params);
int GeneratePunctuation(PunctuationParameters punct_params);

#endif // __RVFUNCTIONS_TD_H_INCLUDED__
