//=======================================================
// include guard
#ifndef __GLOBALSTRUCTS_H_INCLUDED__
#define __GLOBALSTRUCTS_H_INCLUDED__

#include <string>

struct GlobalParameters
{
  int tot_life;
  int max_pop;
  int ancestors;
  int ancestor_clones;
  int num_sims;
  int num_samples;
  int sample_size;
  double detection_threshold;
  double observation_frequency;
  double start_time;
  bool allow_extinction;
  double birth_rate;
  double death_rate;
  double mutation_prob;
  bool trace_ancestry;
  bool count_alleles;
  gsl_rng* rng;
  double seed;
};


// Clone Structure and Linked List Class containing list of clones
struct clone{
  std::string clone_id;
  int subclone_count;
  int cell_count;
  int allele_count;
  int mut_count;
  int driver_count;
  bool is_driver;
  double birth_rate;
  double death_rate;
  double mut_prob;
  double clone_time;

  // linking other nodes
  struct clone *parent;
  struct clone *nextnode;
  struct clone *prevnode;
};

struct FitnessParameters
{
  bool is_randfitness;
  double alpha_fitness;
  double beta_fitness;
  double pass_prob;
  double upper_fitness;
  double lower_fitness;
};

struct MutationParameters
{
  bool is_mutator;
  double alpha_mutation;
  double beta_mutation;
  double pass_prob;
};

struct PunctuationParameters
{
    bool is_punctuated;
    double punctuated_prob;
    double poisson_param;
    double punctuated_multiplier;
    double punctuated_advantageous_prob;
};

struct EpistaticParameters
{
    bool is_epistasis;
    //bool epistatic_model;
    //double epistatic_driver_prob;
    int epistatic_mutation_thresh;
    double epistatic_multiplier;
};

#endif // __GLOBALSTRUCTS_H_INCLUDED__
