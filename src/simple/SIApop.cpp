/*
 * =============================================================================
 *
 *       Filename:  SIApop.cpp
 *
 *    Description: Birth-Death-Mutation process simulation for infinite-allele
 *                 model with random fitness contributions using Gillespie
 *                 Algorithm. Imports data, runs SSA and outputs to designated
 *                 folder.
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

// structure contains all global parameters used in multiple source files
GlobalParameters gp;

int main(int argc, char *argv[])
{

  // for timing total simulation - FOR TESTING
  auto t1 = std::chrono::high_resolution_clock::now();

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

  for (int i = 1; i < argc; i++)
  {
    if (i + 1 != argc)
    { // Check that we haven't finished parsing arguments
        if (strcmp(argv[i], "-in") == 0) // input file
        {
          input_params = argv[i + 1];
        } else if (strcmp(argv[i], "-out") == 0) // output location
        {
          output_folder = argv[i + 1];
        } else if (strcmp(argv[i], "-anc") == 0) // ancestor file
        {
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
  int count_extinct = 0;


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
  params.convert("ancestors", gp.ancestors);
  params.convert("ancestor_clones", gp.ancestor_clones);
  params.convert("num_sims", gp.num_sims);
  params.convert("num_samples", gp.num_samples);
  params.convert("sample_size", gp.sample_size);
  params.convert("detection_threshold", gp.detection_threshold);
  params.convert("allow_extinction", gp.allow_extinction);
  params.convert("birth_rate", gp.birth_rate);
  params.convert("death_rate", gp.death_rate);

  /*
    END OF VARIABLE INPUT AND CONVERSION
  */


  // declare and open other output streams for time and end of sim clone list
  std::ofstream clonedata;
  sprintf(fn, "%s/clonedata.txt", output_folder);
  clonedata.open(fn);
  clonedata.setf(std::ios::fixed);
  clonedata.precision(12);
  clonedata << "run\tunique_id\tnumcells\tbirthrate\tdeathrate\n";


  // Open output stream for sampling data
  std::ofstream sample_data;
  if(gp.sample_size > 0 & gp.num_samples > 0)
  {
    sprintf(fn, "%s/sampledata.txt", output_folder);
    sample_data.open(fn);
    sample_data.setf(std::ios::fixed);
    sample_data << "run\tsample_number\tunique_id\tnumber_obs\n";
  }

  // set RNG seed
  gsl_rng_set(gp.rng, gp.seed);


  // Beginning of simulation that has "sim" number of runs
  for (int sim = 1; sim <= gp.num_sims; sim++)
  {
    // Define CloneList population and initialize variables;
    CloneList population;
    population.init();
    population.tot_cell_count = 0;

    if( ancestor_file == NULL ) // if no ancestor file exists
    {
      // go through all ancestors and initialize clones for each
      for(int ance__clone_count = 1; ance__clone_count <= gp.ancestor_clones; ance__clone_count++)
      {
        // Create Ancestor Node
        struct clone* ancestor;
        ancestor = new struct clone;

        ancestor->cell_count = 0;
        ancestor->birth_rate = gp.birth_rate;
        ancestor->death_rate = gp.death_rate;

        double meangrowth = exp((ancestor->birth_rate - ancestor->death_rate) * gp.tot_life);
        double alpha = (ancestor->death_rate * meangrowth - ancestor->death_rate) /
                       (ancestor->birth_rate * meangrowth - ancestor->death_rate);
        double beta = (ancestor->birth_rate * meangrowth - ancestor->birth_rate) /
                       (ancestor->birth_rate * meangrowth - ancestor->death_rate);

        // Update ancestors to time
        for(int i = 0; i < gp.ancestors; ++i)
        {
          int not_zero = gsl_ran_bernoulli(gp.rng, 1 - alpha);
          if(not_zero == 1)
          {
            int new_cells = gsl_ran_geometric(gp.rng, 1 - beta);
            ancestor->cell_count = ancestor->cell_count + new_cells;
            population.tot_cell_count = population.tot_cell_count + new_cells;
          }
        }

        population.InsertAncestor(ancestor);
      }
    }
    else // ancestor file exists to read from
    {
      std::cout << "Reading ancestor file...";

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
            ancestor->cell_count = 0;
            int num_ancestors = !ancestor_map["numcells"].empty() ? stoi(ancestor_map["numcells"]) : 1;
            ancestor->birth_rate = !ancestor_map["birthrate"].empty() ? stof(ancestor_map["birthrate"]) : gp.birth_rate;
            ancestor->death_rate = !ancestor_map["deathrate"].empty() ? stof(ancestor_map["deathrate"]) : gp.death_rate;


            double meangrowth = exp((ancestor->birth_rate - ancestor->death_rate) * gp.tot_life);
            double alpha = (ancestor->death_rate * meangrowth - ancestor->death_rate) /
                           (ancestor->birth_rate * meangrowth - ancestor->death_rate);
            double beta = (ancestor->birth_rate * meangrowth - ancestor->birth_rate) /
                           (ancestor->birth_rate * meangrowth - ancestor->death_rate);

            // Update ancestors to time
            for(int i = 0; i < num_ancestors; ++i)
            {
              int not_zero = gsl_ran_bernoulli(gp.rng, 1 - alpha);
              if(not_zero == 1)
              {
                int new_cells = gsl_ran_geometric(gp.rng, 1 - beta);
                ancestor->cell_count = ancestor->cell_count + new_cells;
                population.tot_cell_count = population.tot_cell_count + new_cells;
              }
            }

            population.InsertAncestor(ancestor);

            ancestor_map.clear();
        }
      }
    }

    // nonextinction checker - if set to false and goes extinction, restart that sim
    if( (population.tot_cell_count == 0) && (gp.allow_extinction == false) )
    {
      population.DeleteList();
      count_extinct++;
      gp.num_sims++; // increase number of sims
      std::cout << "Population went extinct. Restarting.\n";
      continue;
    }

    // Sampling from population
    if(gp.sample_size > 0 & gp.num_samples > 0)
    {
      population.SampleAndTraverse(sample_data, sim, gp.sample_size, gp.num_samples);
    }
    // Output of end state with clone info
    std::cout << "Traversing and outputting run " << sim << "\n";
    population.Traverse(clonedata, sim);
    std::cout << "Traversal Done\n";

    population.DeleteList();

  }

  //avg_sim_endtime = avg_sim_endtime * (double)gp.num_sims / (double)count_detect;

  gsl_rng_free(gp.rng);

  sim_stats << "count_extinct, " << count_extinct  << "\n";

  clonedata.close();
  sample_data.close();
  sim_stats.close();

  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout <<
      std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() <<
      "milliseconds\n";

  return 0;
}
