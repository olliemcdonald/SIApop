//=======================================================
// include guard
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
