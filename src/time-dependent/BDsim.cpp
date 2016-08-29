/*
 * =============================================================================
 *
 *       Filename:  BDsim.cpp
 *
 *    Description: Time-Dependent Birth-Death-Mutation process simulation for
 *                 infinite-allele model with random fitness contributions using
 *                 Gillespie Algorithm. Imports data, runs SSA and outputs to
 *                 designated folder.
 *
 *        Version:  1.0
 *        Created:  08/24/2016 16:50:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thomas McDonald (), mcdonald@jimmy.harvard.edu
 *   Organization:  DFCI
 *
 * =============================================================================
 */
#include <iostream>
#include <fstream>
#include <chrono>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <gsl/gsl_randist.h>

#include "globalstructs.h"
#include "clonelist.h"
#include "parameterlist.h"
#include "ratefunctions.h"

// structure contains all global parameters used in multiple source files
GlobalParameters gp;
// Function array for different time-dependent function types
RateFunctionsPtr rate_function_array[] = {&RateFunctions::constant,
  &RateFunctions::linear, &RateFunctions::logistic,
  &RateFunctions::Gompertz, &RateFunctions::custom};
gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
// Pointer to Function class which will point to an instance of one based on parameters
CloneList::NewCloneFunction* NewClone;

int main(int argc, char *argv[])
{
  // timer for total simulation time - FOR TESTING
  auto t1 = std::chrono::high_resolution_clock::now();

  // track total error from integration - FOR TESTING
  gp.tot_error = 0;

  //  declaring random number generator and setting seed
  gp.rng = gsl_rng_alloc(gsl_rng_mt19937);
  gp.seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

  /*
    VARIABLE INPUT AND CONVERSION
  */
  // parsing arguments in command line by searching for argument options
  const char *input_params = NULL;
  const char *ancestor_file = NULL;
  // default output is current directory
  const char *output_folder = "./";

  for (int i = 1; i < argc; i++) {
    if (i + 1 != argc)
    { // Check that haven't finished parsing arguments
        if (strcmp(argv[i], "-in") == 0) { // input file
          input_params = argv[i + 1];
        } else if (strcmp(argv[i], "-out") == 0) { // output file
          output_folder = argv[i + 1];
        } else if (strcmp(argv[i], "-anc") == 0) { // ancestor file
          ancestor_file = argv[i + 1];
        }
    }
  }


  // declare and open output stream for simulation statistics
  char fn[100];
  std::ofstream sim_stats;
  sprintf(fn,"%s/sim_stats.txt", output_folder);
  sim_stats.open(fn);
  sim_stats.setf(std::ios::fixed);
  sim_stats.precision(8);


  // declare and initialize parameter list for simulation
  ParameterList params;
  params.init();

  // parsing through the input file and converting/adding to parameter list
  if( input_params != NULL )
  {
    std::string s = input_params;
    std::ifstream infile(input_params);
    if(infile.is_open())
    {
        while(getline(infile, s))
        {
            // allow comments
            if(s[0] == '#') continue;
            params.SplitAndFill(s);
        }
    }
  }

  for(std::map<std::string, std::string>::iterator it=params.begin(); it!=params.end(); ++it)
  {
    sim_stats << it->first << ", " << it->second << "\n";
  }

  // convert all parameters imported from file into respective values in gp
  params.convert("tot_life", gp.tot_life);
  params.convert("max_pop", gp.max_pop);
  params.convert("start_time", gp.start_time);
  params.convert("ancestors", gp.ancestors);
  params.convert("ancestor_clones", gp.ancestor_clones);
  params.convert("num_sims", gp.num_sims);
  params.convert("num_samples", gp.num_samples);
  params.convert("sample_size", gp.sample_size);
  params.convert("detection_threshold", gp.detection_threshold);
  params.convert("observation_frequency", gp.observation_frequency);
  if (gp.observation_frequency == 0) gp.observation_frequency = gp.tot_life;

  std::vector<double> observation_times;
  params.ParseVector("observation_times", observation_times);
  if (observation_times.size() == 0)
  {
    int num_obs =  ceil(gp.tot_life / gp.observation_frequency);
    for (int i = 1; i <= num_obs; i++)
    {
      observation_times.push_back(i * gp.observation_frequency);
    }
  }
  if(observation_times.back() != gp.tot_life) observation_times.push_back(gp.tot_life);


  params.convert("allow_extinction", gp.allow_extinction);
  params.convert("trace_ancestry", gp.trace_ancestry);
  params.convert("count_alleles", gp.count_alleles);
  params.convert("birth_rate", gp.birth_rate);
  params.convert("death_rate", gp.death_rate);
  params.convert("mutation_prob", gp.mutation_prob);
  params.convert("B_max", gp.B_max);
  params.convert("D_max", gp.D_max);

  FitnessParameters fit_params;
  params.convert("alpha_fitness", fit_params.alpha_fitness);
  params.convert("beta_fitness", fit_params.beta_fitness);
  params.convert("pass_prob", fit_params.pass_prob);
  params.convert("upper_fitness", fit_params.upper_fitness);
  params.convert("lower_fitness", fit_params.lower_fitness);
  fit_params.is_randfitness = false;
  if( ((fit_params.alpha_fitness > 0 && fit_params.beta_fitness > 0) ||
      (fit_params.upper_fitness != fit_params.lower_fitness)) &&
      (fit_params.pass_prob < 1) )
  {
    fit_params.is_randfitness = true;
  }

  MutationParameters mut_params;
  params.convert("alpha_mutation", mut_params.alpha_mutation);
  params.convert("beta_mutation", mut_params.beta_mutation);
  params.convert("pass_prob", mut_params.pass_prob);
  mut_params.is_mutator = false;
  if( (mut_params.alpha_mutation > 0 && mut_params.beta_mutation > 0) &&
      (mut_params.pass_prob < 1) )
  {
    mut_params.is_mutator = true;
  }

  PunctuationParameters punct_params;
  params.convert("punctuated_prob", punct_params.punctuated_prob);
  params.convert("poisson_param", punct_params.poisson_param);
  params.convert("punctuated_multiplier", punct_params.punctuated_multiplier);
  params.convert("punctuated_advantageous_prob", punct_params.punctuated_advantageous_prob);
  punct_params.is_punctuated = punct_params.punctuated_prob > 0;

  EpistaticParameters epi_params;
  params.convert("epistatic_mutation_thresh", epi_params.epistatic_mutation_thresh);
  params.convert("epistatic_multiplier", epi_params.epistatic_multiplier);
  epi_params.is_epistasis = epi_params.epistatic_mutation_thresh > 0;
  if (epi_params.epistatic_multiplier == 1.0)
  {
    epi_params.is_epistasis = false;
  }


  TimeDependentParameters td_birth_params;
  params.convert("birth_function", td_birth_params.type);
  params.ParseVector("td_birth_params", td_birth_params.coefs);
  td_birth_params.rate = gp.birth_rate;

  TimeDependentParameters td_death_params;
  params.convert("death_function", td_death_params.type);
  params.ParseVector("td_death_params", td_death_params.coefs);
  td_death_params.rate = gp.death_rate;

  /*
    END OF VARIABLE INPUT AND CONVERSION
  */

  // declare and open other output streams for time and end of sim clone list
  std::ofstream clonedata;
  sprintf(fn, "%s/clonedata.txt", output_folder);
  clonedata.open(fn);
  clonedata.setf(std::ios::fixed);
  clonedata.precision(12);
  if(gp.count_alleles)
  {
    clonedata << "run\tunique_id\tnumcells\tcount_alleles\tbirthrate\tdeathrate\tmutprob\t"
      "initialtime\tsubclone_count\tnum_mut\tnum_drivers\tis_driver" << "\n";
  }
  else
  {
    clonedata << "run\tunique_id\tnumcells\tbirthrate\tdeathrate\tmutprob\t"
      "initialtime\tsubclone_count\tnum_mut\tnum_drivers\tis_driver" << "\n";
  }

  std::ofstream timedata;
  sprintf(fn, "%s/timedata.txt", output_folder);
  timedata.open(fn);
  timedata.setf(std::ios::fixed);
  timedata.precision(12);
  if(gp.trace_ancestry)
  {
    if(gp.count_alleles)
    {
      timedata << "run\ttime\tunique_id\tnumcells\tallelefreq\tgrowth_rate\tinitialtime\t"
        "parent_growth_rate\tparent_initialtime" << "\n";
    }
    else
    {
      timedata << "run\ttime\tunique_id\tnumcells\tgrowth_rate\tinitialtime\t"
        "parent_growth_rate\tparent_initialtime" << "\n";
    }
  }
  else
  {
    if(gp.count_alleles)
    {
      timedata << "run\ttime\tunique_id\tnumcells\tallelefreq\tgrowth_rate\tinitialtime" << "\n";
    }
    else
    {
      timedata << "run\ttime\tunique_id\tnumcells\tgrowth_rate\tinitialtime" << "\n";
    }
  }

  // Open output stream for sampling data
  std::ofstream sample_data;
  if(gp.sample_size > 0 & gp.num_samples > 0)
  {
    sprintf(fn, "%s/sampledata.txt", output_folder);
    sample_data.open(fn);
    sample_data.setf(std::ios::fixed);
    sample_data << "run\tsample_number\tunique_id\tnumber_observed\n";
  }


  // set RNG seed
  gsl_rng_set(gp.rng, gp.seed);

  // simulation variables
  double avg_sim_endtime = 0;
  int count_detect = 0;
  double current_time;
  int curr_observation;
  double rand_next_time;
  int count_extinct = 0;

  // Beginning of simulation that has "sim" number of runs.
  for (int sim = 1; sim <= gp.num_sims; sim++)
  {
    // initialize time to zero
    current_time = gp.start_time;

    // initialize current obseration output to 0
    curr_observation = 0;

    // Define CloneList population and initialize variables;
    CloneList population;
    population.init();

    // Determine Advance function class to use based on the parameters
    if( punct_params.is_punctuated )
    {
      NewClone = new CloneList::NewClonePunct(population, fit_params, mut_params, punct_params);
    }
    else if( fit_params.is_randfitness || mut_params.is_mutator )
    {
      if ( epi_params.is_epistasis )
      {
        NewClone = new CloneList::NewCloneEpi(population, fit_params, mut_params, epi_params);
      }
      else
      {
        NewClone = new CloneList::NewCloneFitMut(population, fit_params, mut_params);
      }
    }
    else /*if (gp.is_custom_model)
    {
      NewClone = new CloneList::NewCloneCustom(popoulation);
    }
    else*/
    {
      NewClone = new CloneList::NewCloneNoParams(population);
    }


    if( ancestor_file == NULL ) // if no ancestor file exists
    {
      /*
        total rate for time-homogeneous population should be the max for the
        birth and death functions times their associated birth rates times
        the number of ancestors
      */
      population.tot_rate_homog = (gp.B_max * gp.birth_rate + gp.D_max * gp.death_rate) * gp.ancestors * gp.ancestor_clones;
      population.tot_cell_count = gp.ancestors * gp.ancestor_clones;

      // go through all ancestors and initialize clones for each
      for(int ance_clone_count = 1; ance_clone_count <= gp.ancestor_clones; ance_clone_count++)
      {
        // Create Ancestor Node
        struct clone* ancestor;
        ancestor = new struct clone;

        ancestor->cell_count = gp.ancestors;
        ancestor->allele_count = gp.ancestors;
        ancestor->mut_prob = gp.mutation_prob;
        ancestor->clone_time = current_time;
        ancestor->subclone_count = 0;
        ancestor->mut_count = 0;
        ancestor->driver_count = 0;
        ancestor->is_driver = false;

        // input default types
        ancestor->birth_params = td_birth_params;
        ancestor->death_params = td_death_params;
        (ancestor->B).function = rate_function_array[td_birth_params.type];
        (ancestor->D).function = rate_function_array[td_death_params.type];
        (ancestor->B).params = &(ancestor->birth_params);
        (ancestor->D).params = &(ancestor->death_params);

        // initializing the rate accumulation values
        ancestor->birth_rate = 0;
        ancestor->death_rate = 0;

        population.InsertAncestor(ancestor);
      }
    }
    else // ancestor file exists to read from
    {
      std::cout << "Reading ancestor file...";

      population.tot_rate_homog = 0;
      population.tot_rate = 0;
      population.tot_cell_count = 0;

      // import first line of file to a vector of keys
      std::vector<std::string> ancestor_keys;
      std::vector<std::string>::iterator it;

      std::string a = ancestor_file;
      std::ifstream ancfile(ancestor_file);

      if(ancfile.is_open())
      {
        getline(ancfile, a);
        std::istringstream ss( a );
        std::string s2;
        while(ss >> s2)
        {
          ancestor_keys.push_back(s2);
        }

        // Loop through all lines of the file
        while(getline(ancfile, a))
        {
          // create a map between column names (from vector) and values
          std::map<std::string, std::string> ancestor_map;
          std::istringstream ss( a );
          std::string s2;

          for(it = ancestor_keys.begin(); it < ancestor_keys.end(); it++)
          {
            ss >> s2;
            ancestor_map[*it];
            ancestor_map[*it] = s2;
          }
            // Move map values to structures
            struct clone* ancestor;
            ancestor = new struct clone;

            ancestor->clone_id = !ancestor_map["unique_id"].empty() ? ancestor_map["unique_id"] : "";
            ancestor->cell_count = !ancestor_map["numcells"].empty() ? stoi(ancestor_map["numcells"]) : 0;
            ancestor->allele_count = ancestor->cell_count;
            ancestor->birth_params.rate = !ancestor_map["birthrate"].empty() ? stof(ancestor_map["birthrate"]) : td_birth_params.rate;
            ancestor->death_params.rate = !ancestor_map["deathrate"].empty() ? stof(ancestor_map["deathrate"]) : td_death_params.rate;
            ancestor->mut_prob = !ancestor_map["mutprob"].empty() ? stof(ancestor_map["mutprob"]) : gp.mutation_prob;
            ancestor->clone_time = !ancestor_map["initialtime"].empty() ? stof(ancestor_map["initialtime"]) : current_time;
            ancestor->subclone_count = !ancestor_map["subclone_count"].empty() ? stoi(ancestor_map["subclone_count"]) : 0;
            ancestor->mut_count = !ancestor_map["num_mut"].empty() ? stoi(ancestor_map["num_mut"]) : 0;
            ancestor->driver_count = !ancestor_map["num_drivers"].empty() ? stoi(ancestor_map["num_drivers"]) : 0;
            ancestor->is_driver = !ancestor_map["is_driver"].empty() ? stoi(ancestor_map["is_driver"]) : false;

            ancestor->birth_params.type = !ancestor_map["birth_function"].empty() ? stoi(ancestor_map["birth_function"]) : td_birth_params.type;
            ancestor->death_params.type = !ancestor_map["death_function"].empty() ? stoi(ancestor_map["death_function"]) : td_death_params.type;

            if(!ancestor_map["bf_params"].empty()) // If the ancestor file contains a parameter list, parse and add to function
            {
              std::string bf_input = ancestor_map["bf_params"];
              std::istringstream ss(bf_input);
              std::string token;
              double value;
              while(getline(ss, token, ','))
              {
                value = stof(token);
                ancestor->birth_params.coefs.push_back(value);
              }
              // Add params to birth function
              (ancestor->B).params = &(ancestor->birth_params);

            }
            else // Use defaults if not present
            {
              ancestor->birth_params = td_birth_params;
              (ancestor->B).params = &(td_birth_params);
            }
            (ancestor->B).function = rate_function_array[ancestor->birth_params.type];


            if(!ancestor_map["df_params"].empty())
            {
              std::string df_input = ancestor_map["df_params"];
              std::istringstream ss(df_input);
              std::string token;
              double value;
              while(getline(ss, token, ','))
              {
                value = stof(token);
                ancestor->death_params.coefs.push_back(value);
              }
              (ancestor->D).params = &(ancestor->death_params);
            }
            else
            {
              ancestor->death_params = td_death_params;
              (ancestor->D).params = &(td_death_params);
            }
            (ancestor->D).function = rate_function_array[ancestor->death_params.type];


            ancestor->birth_rate = ancestor->cell_count * ancestor->birth_params.rate;
            ancestor->death_rate = ancestor->cell_count * ancestor->death_params.rate;

            population.InsertAncestor(ancestor);

            population.tot_rate_homog = population.tot_rate_homog +
                (ancestor->birth_params.rate * gp.B_max + ancestor->death_params.rate * gp.D_max) *
                ancestor->cell_count;
            population.tot_cell_count = population.tot_cell_count + ancestor->cell_count;

            ancestor_map.clear();
        }
      }
      std::cout << "Ancestor File Read...";
    }

    std::cout << "Output Ancestor Population...";
    population.Traverse(timedata, sim, current_time, gp.trace_ancestry, gp.count_alleles);
    std::cout << "Initial Traverse Done\n";

    // Begin single simulation with while loop that exists when hit max time, max pop, or extinction
    while ( (population.tot_cell_count < gp.max_pop) &&
            (population.tot_cell_count > 0) &&
            (current_time < gp.tot_life) )
    {
      // Get next event time by advancing with adaptive thinning
      rand_next_time = population.AdvanceTime(current_time);

      // Advance Simulation State (choose next event)
      population.AdvanceState(current_time, rand_next_time);

      // update current_time
      current_time = current_time + rand_next_time;

      // Method to output data at designated observation times
      while(current_time > observation_times[curr_observation])
      {
        if( (current_time < observation_times[curr_observation + 1]) ||
            (observation_times.size() == curr_observation + 1) )
        {
          population.Traverse(timedata, sim, observation_times[curr_observation], gp.trace_ancestry, gp.count_alleles);
          if((observation_times.size() == curr_observation + 1))
          {
            break;
          }
        }
        curr_observation++;
      }

      /*
      if((population.tot_cell_count % 50000) == 0)
      {
        std::cout << "Time: " << current_time << "\t size: " << population.tot_cell_count << "\n";
      }
      //*/
    }

    // nonextinction checker - if set to false and goes extinction, restart that sim
    if( (population.tot_cell_count == 0) && (gp.allow_extinction == false) )
    {
      population.DeleteList();
      count_extinct++;
      sim = sim - 1;
      continue;
    }

    // If simulation made it to the gp.max_pop in the alotted time
    if( (population.tot_cell_count >= gp.max_pop) &&
        (current_time < gp.tot_life) )
    {
      count_detect = count_detect + 1;
      avg_sim_endtime = avg_sim_endtime + (current_time) / (double)gp.num_sims;
    }

    // Final Timed Output
    population.Traverse(timedata, sim, current_time, gp.trace_ancestry, gp.count_alleles);
    // Sampling from population
    if(gp.sample_size > 0 & gp.num_samples > 0)
    {
      population.SampleAndTraverse(sample_data, sim, gp.sample_size, gp.num_samples);
    }
    // Trim tree if threshold is higher. Otherwise, Traverse
    population.TreeTrim(gp.detection_threshold, gp.max_pop);
    // Output of end state with clone info
    std::cout << "Traversing and outputting run " << sim << "\n";
    population.Traverse(clonedata, sim, gp.count_alleles);
    std::cout << "Traversal Done\n";

    population.DeleteList();

  }

  //avg_sim_endtime = avg_sim_endtime * (double)num_sims / (double)count_detect;

  gsl_rng_free(gp.rng);
  delete NewClone;

  sim_stats << "avg_sim_endtime, " << avg_sim_endtime << "\n" <<
               "count_detect, " << count_detect << "\n" <<
               "count_extinct, " << count_extinct  << "\n" <<
               "tot_integration_error, " << gp.tot_error;

  clonedata.close();
  timedata.close();
  sample_data.close();
  sim_stats.close();

  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout <<
      std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() <<
      "milliseconds\n";

  return 0;
}
