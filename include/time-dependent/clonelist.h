//=======================================================
// include guard
#ifndef __CLONELIST_TD_H_INCLUDED__
#define __CLONELIST_TD_H_INCLUDED__

// dependencies
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include "globalstructs.h"
#include "rvfunctions.h"
#include "ratefunctions.h"


extern GlobalParameters gp;
extern RateFunctionsPtr func_array[];
extern gsl_integration_workspace *workspace;


class CloneList
{
private:
  struct clone *root;
  struct clone *deadroot; // to keep track of dead clones
  struct clone *currdeadnode;
  struct clone *currnode;

public:
  double tot_rate, tot_rate_homog, tot_rate_integ;
  int num_clones;
  int num_mutations;
  int tot_cell_count;
  double int_result_b, int_error_b, int_result_d, int_error_d;

  CloneList() { init(); };
  void init();

  // Abstract base class of functions
  class AdvanceStateFunction
  {
  public:
    AdvanceStateFunction() {};
    virtual ~AdvanceStateFunction() {};
    virtual void operator()(double curr_time, double next_time) = 0;
  };

  class AdvanceStateNoParams : public AdvanceStateFunction
  {
  public:
    AdvanceStateNoParams(CloneList& cl_) : cl(cl_)
    {
    }
    ~AdvanceStateNoParams(){};
    CloneList& cl;
    void operator()(double curr_time, double next_time);
  };

  class AdvanceStateFitMut : public AdvanceStateFunction
  {
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
  public:
    AdvanceStateFitMut(CloneList& cl_, FitnessParameters fit_params_, MutationParameters mut_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_)
    {
    }
    ~AdvanceStateFitMut(){};
    CloneList& cl;
    void operator()(double curr_time, double next_time);
  };

  class AdvanceStatePunct : public AdvanceStateFunction
  {
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
    PunctuationParameters punct_params;
  public:
    AdvanceStatePunct(CloneList& cl_, FitnessParameters fit_params_,
      MutationParameters mut_params_, PunctuationParameters punct_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_),punct_params(punct_params_)
    {
    }
    ~AdvanceStatePunct(){};
    CloneList& cl;
    void operator()(double curr_time, double next_time);
  };

  class AdvanceStateEpi : public AdvanceStateFunction
  {
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
    EpistaticParameters epi_params;
  public:
    AdvanceStateEpi(CloneList& cl_, FitnessParameters fit_params_,
      MutationParameters mut_params_, EpistaticParameters epi_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_),epi_params(epi_params_)
    {
    }
    ~AdvanceStateEpi(){};
    CloneList& cl;
    void operator()(double curr_time, double next_time);
  };

  // Next Step Functions
  double AdvanceTime(double curr_time);
  void InsertNode(clone* newnode, int number_mutations);
  void InsertAncestor(clone* ancestor);

  // Linked List Manipulation Functions
  void ChangeAncestorAllele(clone* thisnode, bool add_daughter);
  void CloneSort(clone* sortnode, bool is_birth);
  void CutNodeOut(clone* zeronode);
  void DeleteNode();
  void TreeTrim(double threshold, int max_pop);

  // Output Functions
  void Traverse(std::ofstream &F, int sim_number, bool count_alleles);
  void Traverse(std::ofstream &F, int sim_number, double obs_time, bool ancestry, bool count_alleles);
  void SampleAndTraverse(std::ofstream &F, int run, int sample_size, int nsamples);
  void DeleteList();
};

#endif // __CLONELIST_TD_H_INCLUDED__