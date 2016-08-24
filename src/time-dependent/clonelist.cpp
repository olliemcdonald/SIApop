#include "clonelist.h"

// Constructor
void CloneList::init()
{
    root = NULL;
    currnode = NULL;
    deadroot = NULL;
    currdeadnode = NULL;
    num_clones = 0; // number of types
    num_mutations = 0; // number of mutations (types and deadtypes)
    tot_rate = 0; // total rate (birth + death)
    tot_rate_homog = 0;
    tot_rate_integ = 0;
    tot_cell_count = 0; // total number of cells in all types
}

// InsertNode attaches the newly created node (in advancestate) to the end of the list
void CloneList::InsertNode(clone* newnode, int number_mutations)
{
  num_clones++;

  if( root == NULL )
  {
    newnode->nextnode = NULL;
    newnode->parent = NULL;
    newnode->prevnode = NULL;

    root = newnode;
    currnode = newnode;
  }
  else
  {
    // Declare extra pieces to append to id
    std::string add_id;
    if(number_mutations == 1)
    {
      add_id = add_id + '>' + std::to_string(num_clones);
    }
    else
    {
      for(int mutation_counter = 1; mutation_counter <= number_mutations; mutation_counter ++)
      {
        add_id = add_id + '>' + std::to_string(num_clones) + '.' + std::to_string(mutation_counter);
      }
    }
    // append new id to list of ancestors
    newnode->clone_id = newnode->clone_id + add_id;
    newnode->nextnode = NULL;
    newnode->parent = currnode;

    // move to end of list to attach newnode to the end
    while (currnode->nextnode != NULL)
    {
      currnode = currnode->nextnode;
    }

    newnode->prevnode = currnode;
    currnode->nextnode = newnode;
  }
}


/*
  Insert the ancestors into the linked list. If ancestors are defined
  in a separate file it parses the file and inserts. Otherwise it
  creates nodes based on the number of clones and individuals defined
  in input file.
*/
void CloneList::InsertAncestor(clone* ancestor)
{
  if( ancestor->clone_id.empty() )
  {
    num_clones++;
    ancestor->clone_id = std::to_string(num_clones) + ".a";
  } else
  {
    ancestor->clone_id = ancestor->clone_id + ".a";
  }

  ancestor->nextnode = NULL;
  ancestor->parent = NULL;

  if( root == NULL )
  {
    ancestor->prevnode = NULL;
    root = ancestor;
    currnode = ancestor;
  }
  else
  {
    while (currnode->nextnode != NULL)
    {
      currnode = currnode->nextnode;
    }

    ancestor->prevnode = currnode;
    currnode->nextnode = ancestor;
  }
}


/*
  ChangeAncestorAllele occurs after determining which event occurs. It travels through
  the ancestry using the parent link in the list to add or subtract a clone's
  numallele value for each ancestor. This generates a count for the number
  of individuals a particular allele is represented by (allele is the last
  term of unique_id)
*/
void CloneList::ChangeAncestorAllele(clone* thisnode, bool add_daughter)
{
  if(add_daughter) // if added a daughter add to each allele_count
  {
    while(thisnode != NULL)
    {
      // add 1 to all allele_count for each ancestor
      thisnode->allele_count++;
      thisnode = thisnode->parent;
    }
  }
  else // else subtract from each allele_count
  {
    while(thisnode != NULL)
    {
      // add 1 to all allele_count for each ancestor
      thisnode->allele_count--;
      thisnode = thisnode->parent;
    }
  }
}


/*
  Method to sort the current node based on the number of cells in it
  relative to its prevnode and nextnode. If death, move toward nextnode
  if birth, move to prevnode direction second argument is true if birth.
*/
void CloneList::CloneSort(clone* sortnode, bool is_birth)
{
  if(is_birth && sortnode->prevnode != NULL)
  {
    if( sortnode->cell_count * (sortnode->birth_rate + sortnode->death_rate) >=
      sortnode->prevnode->cell_count * (sortnode->prevnode->birth_rate + sortnode->prevnode->death_rate) )
    {
      struct clone* newprevnode;
      struct clone* newnextnode;

      newnextnode = sortnode->prevnode;
      newprevnode = newnextnode->prevnode;
      // attach prevnode and nextnode
      if(sortnode->nextnode != NULL)
      {
        sortnode->nextnode->prevnode = sortnode->prevnode;
      }
      sortnode->prevnode->nextnode = sortnode->nextnode;
      // insert clone between newprevnode and newnextnode
      sortnode->prevnode = newprevnode;
      sortnode->nextnode = newnextnode;

      if(newprevnode != NULL)
      {
        newprevnode->nextnode = sortnode;
      } else
      {
        // reroot to the current node
        root = sortnode;
      }
      newnextnode->prevnode = sortnode;
    }
  } else if(!is_birth && sortnode->nextnode != NULL)
  {
    // move newprevnode to the right
    if( sortnode->cell_count * (sortnode->birth_rate + sortnode->death_rate) <=
      sortnode->nextnode->cell_count  * (sortnode->nextnode->birth_rate + sortnode->nextnode->death_rate) )
    {
      struct clone* newprevnode;
      struct clone* newnextnode;

      newprevnode = sortnode->nextnode;
      newnextnode = newprevnode->nextnode;
      // attach prevnode and nextnode
      if(sortnode->prevnode != NULL)
      {
        sortnode->prevnode->nextnode = sortnode->nextnode;
      } else
      {
        root = sortnode->nextnode;
      }
      sortnode->nextnode->prevnode = sortnode->prevnode;

      // insert clone between newprevnode and newnextnode
      sortnode->prevnode = newprevnode;
      sortnode->nextnode = newnextnode;
      if(newnextnode != NULL)
      {
        newnextnode->prevnode = sortnode;
      }
      newprevnode->nextnode = sortnode;
    }
  }
}

double CloneList::AdvanceTime(double curr_time)
{
  struct clone *pnode;
  double rand_next_time;

  while(true)
  {
    tot_rate = 0;
    rand_next_time = gsl_ran_exponential(gp.rng, 1 / tot_rate_homog);
    pnode = root;

    while(pnode)
    {
      //    time dependent rate to test
      tot_rate = tot_rate + (GSL_FN_EVAL(&(pnode->B), curr_time + rand_next_time) +
        GSL_FN_EVAL(&(pnode->D), curr_time + rand_next_time)) * pnode->cell_count;
      pnode = pnode->nextnode;
    }

    double u_thin = gsl_ran_flat(gp.rng, 0, 1);
    double beta_ratio = tot_rate / tot_rate_homog;

    if(u_thin <= beta_ratio)
    {
      break;
    }
  }

  return rand_next_time;
}

void CloneList::AdvanceStateNoParams::operator()(double curr_time, double next_time)
{
  double summand = 0;
  struct clone *pnode;
  cl.tot_rate_integ = 0;
  bool flag = false;

  pnode = cl.root;
  // Collect total rates over time between last and next event
  while(pnode)
  {
    // integrals for determining next clone
    gsl_integration_qags(&(pnode->B), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_b), &(cl.int_error_b));
    pnode->birth_rate = cl.int_result_b;
    gp.tot_error = gp.tot_error + cl.int_error_b;

    gsl_integration_qags(&(pnode->D), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_d), &(cl.int_error_d));
    pnode->death_rate = cl.int_result_d;
    gp.tot_error = gp.tot_error + cl.int_error_d;


    cl.tot_rate_integ = cl.tot_rate_integ + (pnode->birth_rate + pnode->death_rate) * pnode->cell_count;
    pnode = pnode->nextnode;
  }

  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate_integ);
  double rand_mut_occur;
  // Put pnode back to the root
  pnode = cl.root;

  while ( (pnode) && !(flag) )
  {
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {
      rand_mut_occur = gsl_ran_flat(gp.rng, 0, 1);

      if (rand_mut_occur <= pnode->mut_prob)
      {
        // Creation of new clone
        struct clone *new_mut_node;
        new_mut_node = new struct clone;
        new_mut_node->clone_id = pnode->clone_id;
        new_mut_node->driver_count = pnode->driver_count;
        new_mut_node->subclone_count = 0;
        new_mut_node->cell_count = 1;
        new_mut_node->allele_count = 0;
        new_mut_node->clone_time = curr_time + next_time;
        new_mut_node->mut_count = pnode->mut_count + 1;
        new_mut_node->is_driver = false;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        // update rates
        new_mut_node->mut_prob = pnode->mut_prob;
        // updating time-dependent parameters
        new_mut_node->birth_params = pnode->birth_params;
        new_mut_node->death_params = pnode->death_params;
        (new_mut_node->B).function = func_array[new_mut_node->birth_params.type];
        (new_mut_node->B).params = &(new_mut_node->birth_params);
        (new_mut_node->D).function = func_array[new_mut_node->death_params.type];
        (new_mut_node->D).params = &(new_mut_node->death_params);
        // end of update rates

        // Reset birth and death rates (recalc at beginning of next AdvanceState)
        new_mut_node->birth_rate = 0;
        new_mut_node->death_rate = 0;

        cl.tot_rate_homog = cl.tot_rate_homog + new_mut_node->birth_params.rate * gp.B_max +
          new_mut_node->death_params.rate * gp.D_max;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations++;

        cl.currnode = pnode;
        pnode = new_mut_node;
        cl.InsertNode(pnode, 1);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        // Add to homogeneous process that undergoes thinning
        cl.tot_rate_homog = cl.tot_rate_homog + pnode->birth_params.rate * gp.B_max +
          pnode->death_params.rate * gp.D_max;
      }

      // Reprioritize clone based on size by moving to left
      cl.CloneSort(pnode, true);

      flag = true;
      break;

    }

    if(rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) + (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      cl.tot_cell_count = cl.tot_cell_count - 1;
      if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, false);
      cl.tot_rate_homog = cl.tot_rate_homog - pnode->birth_params.rate * gp.B_max -
        pnode->death_params.rate * gp.D_max;

      // Clean up of clones with zero to speed up later runs
      if(pnode->cell_count == 0)
      {
        // remove pnode by attaching it to deadlist
        cl.CutNodeOut(pnode);
        flag = true;
        break;
      }

      // sort by moving to the right until fits
      cl.CloneSort(pnode, false);
      flag = true;
      break;
    }

    summand = summand + (pnode->cell_count) * ((pnode->birth_rate + pnode->death_rate));
    pnode = pnode->nextnode;

  }

  if (!flag)
  {
    std::cout << "error: step not completed" << "\n";
    exit(0);
  }
}

void CloneList::AdvanceStateFitMut::operator()(double curr_time, double next_time)
{
  double summand = 0;
  struct clone *pnode;
  cl.tot_rate_integ = 0;
  bool flag = false;

  pnode = cl.root;
  // Collect total rates over time between last and next event
  while(pnode)
  {
    // integrals for determining next clone
    gsl_integration_qags(&(pnode->B), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_b), &(cl.int_error_b));
    pnode->birth_rate = cl.int_result_b;
    gp.tot_error = gp.tot_error + cl.int_error_b;

    gsl_integration_qags(&(pnode->D), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_d), &(cl.int_error_d));
    pnode->death_rate = cl.int_result_d;
    gp.tot_error = gp.tot_error + cl.int_error_d;

    cl.tot_rate_integ = cl.tot_rate_integ + (pnode->birth_rate + pnode->death_rate) * pnode->cell_count;
    pnode = pnode->nextnode;
  }

  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate_integ);
  double rand_mut_occur;
  // Put pnode back to the root
  pnode = cl.root;

  while ( (pnode) && !(flag) )
  {
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {
      rand_mut_occur = gsl_ran_flat(gp.rng, 0, 1);

      if (rand_mut_occur <= pnode->mut_prob)
      {
        // Creation of new clone
        struct clone *new_mut_node;
        new_mut_node = new struct clone;
        new_mut_node->clone_id = pnode->clone_id;
        new_mut_node->driver_count = pnode->driver_count;
        new_mut_node->subclone_count = 0;
        new_mut_node->cell_count = 1;
        new_mut_node->allele_count = 0;
        new_mut_node->clone_time = curr_time + next_time;
        new_mut_node->mut_count = pnode->mut_count + 1;
        new_mut_node->is_driver = false;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        // update rates
        // updating time-dependent parameters
        new_mut_node->birth_params = pnode->birth_params;
        new_mut_node->death_params = pnode->death_params;

        bool did_count_driver = false;
        // generation of additive rate to the fitness
        if(fit_params.is_randfitness)
        {
          double additional_rate = GenerateFitness(fit_params);

          if (additional_rate > 0)
          {
            new_mut_node->driver_count++;
            did_count_driver = true;
            new_mut_node->is_driver = true;
          }
          new_mut_node->birth_params.rate = fmax(0, additional_rate + new_mut_node->birth_params.rate);
        }

        if(mut_params.is_mutator)
        {
          double additional_mut_prob = GenerateMutationProb(mut_params);
          if(additional_mut_prob > 0)
          {
            if(!did_count_driver) new_mut_node->driver_count++;
            new_mut_node->is_driver = true;
          }
          new_mut_node->mut_prob = fmin(1, pnode->mut_prob + additional_mut_prob);
        }
        else
        {
          new_mut_node->mut_prob = pnode->mut_prob;
        }

        (new_mut_node->B).function = func_array[new_mut_node->birth_params.type];
        (new_mut_node->D).function = func_array[new_mut_node->death_params.type];
        (new_mut_node->B).params = &(new_mut_node->birth_params);
        (new_mut_node->D).params = &(new_mut_node->death_params);
        // end of update rates

        // Reset birth and death rates (recalc at beginning of next AdvanceState)
        new_mut_node->birth_rate = 0;
        new_mut_node->death_rate = 0;

        cl.tot_rate_homog = cl.tot_rate_homog + new_mut_node->birth_params.rate * gp.B_max +
          new_mut_node->death_params.rate * gp.D_max;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations++;

        cl.currnode = pnode;
        pnode = new_mut_node;
        cl.InsertNode(pnode, 1);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        // Add to homogeneous process that undergoes thinning
        cl.tot_rate_homog = cl.tot_rate_homog + pnode->birth_params.rate * gp.B_max +
          pnode->death_params.rate * gp.D_max;
      }

      // Reprioritize clone based on size by moving to left
      cl.CloneSort(pnode, true);

      flag = true;
      break;

    }

    if(rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) + (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      cl.tot_cell_count = cl.tot_cell_count - 1;
      if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, false);
      cl.tot_rate_homog = cl.tot_rate_homog - pnode->birth_params.rate * gp.B_max -
        pnode->death_params.rate * gp.D_max;

      // Clean up of clones with zero to speed up later runs
      if(pnode->cell_count == 0)
      {
        // remove pnode by attaching it to deadlist
        cl.CutNodeOut(pnode);
        flag = true;
        break;
      }

      // sort by moving to the right until fits
      cl.CloneSort(pnode, false);
      flag = true;
      break;
    }

    summand = summand + (pnode->cell_count) * ((pnode->birth_rate + pnode->death_rate));
    pnode = pnode->nextnode;

  }

  if (!flag)
  {
    std::cout << "error: step not completed" << "\n";
    exit(0);
  }
}

void CloneList::AdvanceStatePunct::operator()(double curr_time, double next_time)
{
  double summand = 0;
  struct clone *pnode;
  cl.tot_rate_integ = 0;
  bool flag = false;

  pnode = cl.root;
  // Collect total rates over time between last and next event
  while(pnode)
  {
    // integrals for determining next clone
    gsl_integration_qags(&(pnode->B), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_b), &(cl.int_error_b));
    pnode->birth_rate = cl.int_result_b;
    gp.tot_error = gp.tot_error + cl.int_error_b;

    gsl_integration_qags(&(pnode->D), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_d), &(cl.int_error_d));
    pnode->death_rate = cl.int_result_d;
    gp.tot_error = gp.tot_error + cl.int_error_d;


    cl.tot_rate_integ = cl.tot_rate_integ + (pnode->birth_rate + pnode->death_rate) * pnode->cell_count;
    pnode = pnode->nextnode;
  }

  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate_integ);
  double rand_mut_occur;
  // Put pnode back to the root
  pnode = cl.root;

  while ( (pnode) && !(flag) )
  {
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {
      rand_mut_occur = gsl_ran_flat(gp.rng, 0, 1);

      if (rand_mut_occur <= pnode->mut_prob)
      {
        int number_mutations = 1;

        // Creation of new clone
        struct clone *new_mut_node;
        new_mut_node = new struct clone;
        new_mut_node->clone_id = pnode->clone_id;
        new_mut_node->driver_count = pnode->driver_count;
        new_mut_node->subclone_count = 0;
        new_mut_node->cell_count = 1;
        new_mut_node->allele_count = 0;
        new_mut_node->clone_time = curr_time + next_time;
        new_mut_node->mut_count = pnode->mut_count + 1;
        new_mut_node->is_driver = false;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        // update rates
        // generation of punctuated number of mutations
        double rand_punct = gsl_ran_flat(gp.rng, 0, 1);
        double rand_advantage = 0;
        if(rand_punct < punct_params.punctuated_prob)
        {
          number_mutations = GeneratePunctuation(punct_params);
          rand_advantage = gsl_ran_flat(gp.rng, 0, 1);
        }

        // updating time-dependent parameters
        new_mut_node->birth_params = pnode->birth_params;
        new_mut_node->death_params = pnode->death_params;

        bool did_count_driver = false;
        // generation of additive rate to the fitness
        if(fit_params.is_randfitness)
        {
          double additional_rate = GenerateFitness(fit_params);
          if (additional_rate > 0)
          {
            new_mut_node->driver_count++;
            did_count_driver = true;
            new_mut_node->is_driver = true;
            if(rand_punct < punct_params.punctuated_prob)
            {
              additional_rate = additional_rate * punct_params.punctuated_multiplier;
              // the additional rate can go to the birth or the death rate
              if( rand_advantage < punct_params.punctuated_advantageous_prob )
              {
                new_mut_node->birth_params.rate = fmax(0, additional_rate + new_mut_node->birth_params.rate);
              }
              else
              {
                new_mut_node->death_params.rate = fmax(0, additional_rate + new_mut_node->death_params.rate);
              }
            }
          }
        }

        if(mut_params.is_mutator)
        {
          double additional_mut_prob = GenerateMutationProb(mut_params);
          if(additional_mut_prob > 0)
          {
            if(!did_count_driver) new_mut_node->driver_count++;
            new_mut_node->is_driver = true;
          }
          new_mut_node->mut_prob = fmin(1, pnode->mut_prob + additional_mut_prob);
        }
        else
        {
          new_mut_node->mut_prob = pnode->mut_prob;
        }

        (new_mut_node->B).function = func_array[new_mut_node->birth_params.type];
        (new_mut_node->B).params = &(new_mut_node->birth_params);
        (new_mut_node->D).function = func_array[new_mut_node->death_params.type];
        (new_mut_node->D).params = &(new_mut_node->death_params);
        // end of update rates

        // Reset birth and death rates (recalc at beginning of next AdvanceState)
        new_mut_node->birth_rate = 0;
        new_mut_node->death_rate = 0;

        cl.tot_rate_homog = cl.tot_rate_homog + new_mut_node->birth_params.rate * gp.B_max +
          new_mut_node->death_params.rate * gp.D_max;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations = cl.num_mutations + number_mutations;

        cl.currnode = pnode;
        pnode = new_mut_node;
        cl.InsertNode(pnode, number_mutations);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        // Add to homogeneous process that undergoes thinning
        cl.tot_rate_homog = cl.tot_rate_homog + pnode->birth_params.rate * gp.B_max +
          pnode->death_params.rate * gp.D_max;
      }

      // Reprioritize clone based on size by moving to left
      cl.CloneSort(pnode, true);

      flag = true;
      break;

    }

    if(rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) + (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      cl.tot_cell_count = cl.tot_cell_count - 1;
      if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, false);
      cl.tot_rate_homog = cl.tot_rate_homog - pnode->birth_params.rate * gp.B_max -
        pnode->death_params.rate * gp.D_max;

      // Clean up of clones with zero to speed up later runs
      if(pnode->cell_count == 0)
      {
        // remove pnode by attaching it to deadlist
        cl.CutNodeOut(pnode);
        flag = true;
        break;
      }

      // sort by moving to the right until fits
      cl.CloneSort(pnode, false);
      flag = true;
      break;
    }

    summand = summand + (pnode->cell_count) * ((pnode->birth_rate + pnode->death_rate));
    pnode = pnode->nextnode;

  }

  if (!flag)
  {
    std::cout << "error: step not completed" << "\n";
    exit(0);
  }
}

void CloneList::AdvanceStateEpi::operator()(double curr_time, double next_time)
{
  double summand = 0;
  struct clone *pnode;
  cl.tot_rate_integ = 0;
  bool flag = false;

  pnode = cl.root;
  // Collect total rates over time between last and next event
  while(pnode)
  {
    // integrals for determining next clone
    gsl_integration_qags(&(pnode->B), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_b), &(cl.int_error_b));
    pnode->birth_rate = cl.int_result_b;
    gp.tot_error = gp.tot_error + cl.int_error_b;

    gsl_integration_qags(&(pnode->D), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(cl.int_result_d), &(cl.int_error_d));
    pnode->death_rate = cl.int_result_d;
    gp.tot_error = gp.tot_error + cl.int_error_d;


    cl.tot_rate_integ = cl.tot_rate_integ + (pnode->birth_rate + pnode->death_rate) * pnode->cell_count;
    pnode = pnode->nextnode;
  }

  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate_integ);
  double rand_mut_occur;
  // Put pnode back to the root
  pnode = cl.root;

  while ( (pnode) && !(flag) )
  {
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {
      rand_mut_occur = gsl_ran_flat(gp.rng, 0, 1);

      if (rand_mut_occur <= pnode->mut_prob)
      {
        // Creation of new clone
        struct clone *new_mut_node;
        new_mut_node = new struct clone;
        new_mut_node->clone_id = pnode->clone_id;
        new_mut_node->driver_count = pnode->driver_count;
        new_mut_node->subclone_count = 0;
        new_mut_node->cell_count = 1;
        new_mut_node->allele_count = 0;
        new_mut_node->clone_time = curr_time + next_time;
        new_mut_node->mut_count = pnode->mut_count + 1;
        new_mut_node->is_driver = false;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        // update rates
        // updating time-dependent parameters
        new_mut_node->birth_params = pnode->birth_params;
        new_mut_node->death_params = pnode->death_params;

        bool did_count_driver = false;
        // generation of additive rate to the fitness
        if(fit_params.is_randfitness)
        {
          double additional_rate = GenerateFitness(fit_params);
          if (additional_rate > 0)
          {
            new_mut_node->driver_count++;
            did_count_driver = true;
            new_mut_node->is_driver = true;
            if(new_mut_node->mut_count == epi_params.epistatic_mutation_thresh)
            {
              additional_rate = additional_rate * epi_params.epistatic_multiplier;
            }
            new_mut_node->birth_params.rate = fmax(0, additional_rate + new_mut_node->birth_params.rate);
          }
        }

        if(mut_params.is_mutator)
        {
          double additional_mut_prob = GenerateMutationProb(mut_params);
          if(additional_mut_prob > 0)
          {
            if(!did_count_driver) new_mut_node->driver_count++;
            new_mut_node->is_driver = true;
          }
          new_mut_node->mut_prob = fmin(1, pnode->mut_prob + additional_mut_prob);
        }
        else
        {
          new_mut_node->mut_prob = pnode->mut_prob;
        }

        (new_mut_node->B).function = func_array[new_mut_node->birth_params.type];
        (new_mut_node->B).params = &(new_mut_node->birth_params);
        (new_mut_node->D).function = func_array[new_mut_node->death_params.type];
        (new_mut_node->D).params = &(new_mut_node->death_params);
        // end of update rates

        // Reset birth and death rates (recalc at beginning of next AdvanceState)
        new_mut_node->birth_rate = 0;
        new_mut_node->death_rate = 0;

        cl.tot_rate_homog = cl.tot_rate_homog + new_mut_node->birth_params.rate * gp.B_max +
          new_mut_node->death_params.rate * gp.D_max;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations++;

        cl.currnode = pnode;
        pnode = new_mut_node;
        cl.InsertNode(pnode, 1);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        // Add to homogeneous process that undergoes thinning
        cl.tot_rate_homog = cl.tot_rate_homog + pnode->birth_params.rate * gp.B_max +
          pnode->death_params.rate * gp.D_max;
      }

      // Reprioritize clone based on size by moving to left
      cl.CloneSort(pnode, true);

      flag = true;
      break;

    }

    if(rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) + (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      cl.tot_cell_count = cl.tot_cell_count - 1;
      if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, false);
      cl.tot_rate_homog = cl.tot_rate_homog - pnode->birth_params.rate * gp.B_max -
        pnode->death_params.rate * gp.D_max;

      // Clean up of clones with zero to speed up later runs
      if(pnode->cell_count == 0)
      {
        // remove pnode by attaching it to deadlist
        cl.CutNodeOut(pnode);
        flag = true;
        break;
      }

      // sort by moving to the right until fits
      cl.CloneSort(pnode, false);
      flag = true;
      break;
    }

    summand = summand + (pnode->cell_count) * ((pnode->birth_rate + pnode->death_rate));
    pnode = pnode->nextnode;

  }

  if (!flag)
  {
    std::cout << "error: step not completed" << "\n";
    exit(0);
  }
}

/*
  Final out to file after running through a simulation
*/
void CloneList::Traverse(std::ofstream &F, int sim_number, bool count_alleles = false)
{
  struct clone *pnode;
  if(count_alleles)
  {
    for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->allele_count << "\t" <<
           pnode->birth_params.rate << "\t" <<
           pnode->death_params.rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }

    for (pnode = deadroot; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->allele_count << "\t" <<
           pnode->birth_params.rate << "\t" <<
           pnode->death_params.rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }
  }
  else
  {
    for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->birth_params.rate << "\t" <<
           pnode->death_params.rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }

    for (pnode = deadroot; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->birth_params.rate << "\t" <<
           pnode->death_params.rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }
  }
}

/*
  For traversing and outputting at a designated observation time not the end
*/
void CloneList::Traverse(std::ofstream &F, int sim_number, double obs_time, bool ancestry = false, bool count_alleles = false)
{
  struct clone *pnode;
  if(ancestry)
  {
    if(count_alleles)
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" << pnode->allele_count << "\t" <<
             pnode->birth_params.rate - pnode->death_params.rate << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << pnode->parent->birth_params.rate - pnode->parent->death_params.rate << "\t" <<
               pnode->parent->clone_time << "\n";
        }
      }
    }
    else
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" <<
             pnode->birth_params.rate - pnode->death_params.rate << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << pnode->parent->birth_params.rate - pnode->parent->death_params.rate << "\t" <<
               pnode->parent->clone_time << "\n";
        }
      }
    }
  }
  else
  {
    if(count_alleles)
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" << pnode->allele_count << "\t" <<
             pnode->birth_params.rate - pnode->death_params.rate << "\t" <<
             pnode->clone_time << "\n";
      }
    }
    else
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" <<
             pnode->birth_params.rate - pnode->death_params.rate << "\t" <<
             pnode->clone_time << "\n";
      }
    }
  }
}

void CloneList::SampleAndTraverse(std::ofstream &F, int run, int sample_size, int nsamples)
{
  // loop through to repeat with all samples
  for(int sample_counter = 1; sample_counter <= nsamples; sample_counter++)
  {
    int samples_to_place = sample_size;
    int cells_left = tot_cell_count;

    // assign cells from samples_to_place to nodes using binomial
    // probabilities (same as multinomial r.v. but don't need to
    // know all beforehand
    struct clone *pnode = root;
    while(samples_to_place > 0)
    {
      // calculate prob and simulate number to sample from this clone
      double prob = (double)pnode->cell_count / (double)cells_left;
      int samples_placed = gsl_ran_binomial(gp.rng, prob, samples_to_place);

      // only write if sampled any
      if(samples_placed > 0)
      {
        F << run << "\t" << sample_counter << "\t" <<
          pnode->clone_id << "\t" << samples_placed << "\n";
      }

      samples_to_place = samples_to_place - samples_placed;
      cells_left = cells_left - pnode->cell_count;

      pnode = pnode->nextnode;
    }
  }
}

/*
  Run through all clones and check that they have a large enough popsize -
  remove those that don't by putting their count into their parent pop
  then deleting the node
*/
void CloneList::TreeTrim(double threshold, int max_pop)
{
  //Move to last clone
  while(currnode->nextnode != NULL)
  {
    currnode = currnode->nextnode;
  }

  double cell_cutoff = threshold * max_pop;
  // Starting with current node and working back until at root
  while( currnode->prevnode != NULL )
  {
    if( currnode->cell_count < cell_cutoff )
    {
      // add current clones numbers to that of parent
      currnode->parent->cell_count = currnode->parent->cell_count + currnode->cell_count;
      // remove the node and update pointers
      DeleteNode();
    }
    else
    {
      currnode = currnode->prevnode;
    }

  }
}

/*
  When a node hits 0 count, move to dead side of list
  if at end of list, no need to link the next node (NULL value), Otherwise
  cut out value
*/
void CloneList::CutNodeOut(clone* zeronode)
{
  if(zeronode->nextnode != NULL)
  {
    zeronode->nextnode->prevnode = zeronode->prevnode;
  }
  /* if at beginning of list, need to reroot the list to the next node
     otherwise, cut out node*/
  if(zeronode->prevnode != NULL)
  {
    zeronode->prevnode->nextnode = zeronode->nextnode;
  }
  else
  {
    root = zeronode->nextnode;
  }

  /* if dead size not rooted, root with first zero node,
     otherwise add to end of list*/
  if( deadroot == NULL)
  {
    deadroot = zeronode;
    currdeadnode = zeronode;
  }
  else
  {
    zeronode->prevnode = currdeadnode;
    currdeadnode->nextnode = zeronode;
    currdeadnode = zeronode;
    zeronode->nextnode = NULL;
  }
}

// This function used for TreeTrim to remove nodes that are too small
void CloneList::DeleteNode()
{
  struct clone *tmpcurrnode;
  tmpcurrnode = currnode;

  if(currnode->nextnode == NULL)
  {
    currnode = currnode->prevnode;
    currnode->nextnode = NULL;
  }
  else
  {
    struct clone *tmpnextnode;
    tmpnextnode = currnode->nextnode;
    currnode = currnode->prevnode; // set previous node as current
    currnode->nextnode = tmpnextnode; // set the nextnode of the current node to skip over tmpnode
    tmpnextnode->prevnode = currnode;
  }

  delete tmpcurrnode;
  num_clones = num_clones - 1; // increase numtypes (have 1 type in pop)
}

// Delete all members of the linked list
void CloneList::DeleteList()
{
  struct clone *pnode;
  while (root != NULL)
  {
    pnode = root;
    root = root->nextnode;
    delete pnode;
  }
}
