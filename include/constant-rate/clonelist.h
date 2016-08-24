//=======================================================
// include guard
#ifndef __CLONELIST_H_INCLUDED__
#define __CLONELIST_H_INCLUDED__

// dependencies
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include "globalstructs.h"
#include "rvfunctions.h"


extern GlobalParameters gp;

class CloneList
{
private:
  struct clone *root;
  struct clone *deadroot; // to keep track of dead clones
  struct clone *currdeadnode;
  struct clone *currnode;

public:
  double tot_rate;
  int num_clones;
  int num_mutations;
  int tot_cell_count;

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
  public:
    AdvanceStateFitMut(CloneList& cl_, FitnessParameters fit_params_, MutationParameters mut_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_)
    {
    }
    ~AdvanceStateFitMut(){};
    CloneList& cl;
    void operator()(double curr_time, double next_time);
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
  };

  class AdvanceStatePunct : public AdvanceStateFunction
  {
  public:
    AdvanceStatePunct(CloneList& cl_, FitnessParameters fit_params_,
      MutationParameters mut_params_, PunctuationParameters punct_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_),punct_params(punct_params_)
    {
    }
    ~AdvanceStatePunct(){};
    CloneList& cl;
    void operator()(double curr_time, double next_time);
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
    PunctuationParameters punct_params;
  };

  class AdvanceStateEpi : public AdvanceStateFunction
  {
  public:
    AdvanceStateEpi(CloneList& cl_, FitnessParameters fit_params_,
      MutationParameters mut_params_, EpistaticParameters epi_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_),epi_params(epi_params_)
    {
    }
    ~AdvanceStateEpi(){};
    CloneList& cl;
    void operator()(double curr_time, double next_time);
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
    EpistaticParameters epi_params;
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

#endif // __CLONELIST_H_INCLUDED__
