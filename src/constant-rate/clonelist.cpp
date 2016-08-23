#include "clonelist.h"


// Constructor
void CloneList::init()
{
    root = NULL;
    deadroot = NULL;
    currnode = NULL;
    currdeadnode = NULL;
    num_clones = 0; // number of types
    num_mutations = 0; // number of mutations (types and deadtypes)
    tot_rate = 0; // total rate (birth + death)
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
    ancestor->nextnode = NULL;
    ancestor->parent = NULL;
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
  AddToParent occurs after determining which event occurs. It travels through
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
      } else{
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
      } else{
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
  double rand_next_time;
  rand_next_time = gsl_ran_exponential(gp.rng, 1 / tot_rate);
  return rand_next_time;
}

void CloneList::AdvanceStateNoParams::operator()(double curr_time, double next_time)
{
  double summand = 0;
  struct clone *pnode;
  bool flag = false;

  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate);
  double rand_mut_occur;
  // Put pnode back at the root
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
        new_mut_node->mut_count = pnode->mut_count;
        new_mut_node->is_driver = false;

        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        new_mut_node->mut_count++;
        new_mut_node->mut_prob = pnode->mut_prob;

        // updating rate parameters
        new_mut_node->birth_rate = pnode->birth_rate;
        new_mut_node->death_rate = pnode->death_rate;

        cl.tot_rate = cl.tot_rate + new_mut_node->birth_rate + new_mut_node->death_rate;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations++;

        cl.currnode = pnode;
        pnode = new_mut_node;

        // add to all allele counts for ancestors
        cl.InsertNode(new_mut_node, 1);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        cl.tot_rate = cl.tot_rate + (pnode->birth_rate) + (pnode->death_rate);
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
      cl.tot_rate = cl.tot_rate - pnode->birth_rate - pnode->death_rate;

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

    summand = summand + (pnode->cell_count) * (pnode->birth_rate + pnode->death_rate);
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
  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate);
  double rand_mut_occur;

  struct clone *pnode;
  pnode = cl.root;
  bool flag = false;

  while ( (pnode) && !(flag) )
  {
    // Condition for new birth
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {

      rand_mut_occur = gsl_ran_flat(gp.rng, 0, 1);
      // Condition to determine if mutation occurs in new daughter
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
        new_mut_node->death_rate = pnode->death_rate;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;


        bool did_count_driver = false;
        // generation of additive rate to the fitness
        if(fit_params.is_randfitness){
          double additional_rate = GenerateFitness(fit_params);
          if (additional_rate > 0)
          {
            new_mut_node->driver_count++;
            did_count_driver = true;
            new_mut_node->is_driver = true;
          }
          new_mut_node->birth_rate = fmax(0, additional_rate + pnode->birth_rate);
        }
        else
        {
          new_mut_node->birth_rate = pnode->birth_rate;
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

        cl.tot_rate = cl.tot_rate + new_mut_node->birth_rate + new_mut_node->death_rate;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations++;

        cl.currnode = pnode;
        pnode = new_mut_node;
        cl.InsertNode(new_mut_node, 1);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        cl.tot_rate = cl.tot_rate + (pnode->birth_rate) + (pnode->death_rate);
      }

      // Reprioritize clone based on size
      cl.CloneSort(pnode, true);

      flag = true;
      break;
    }
    if( rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) +
      (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      cl.tot_cell_count = cl.tot_cell_count - 1;
      if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, false);
      cl.tot_rate = cl.tot_rate - pnode->birth_rate - pnode->death_rate;

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

    summand = summand + (pnode->cell_count) * (pnode->birth_rate + pnode->death_rate);
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
  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate);

  struct clone *pnode;
  pnode = cl.root;
  bool flag = false;

  while ( (pnode) && !(flag) )
  {
    // Condition for new birth
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {

      double rand_mut_occur = gsl_ran_flat(gp.rng, 0, 1);
      // Condition to determine if mutation occurs in new daughter
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
        new_mut_node->death_rate = pnode->death_rate;
        new_mut_node->birth_rate = pnode->birth_rate;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        // generation of punctuated number of mutations
        double rand_punct = gsl_ran_flat(gp.rng, 0, 1);
        double rand_advantage = 0;
        if(rand_punct < punct_params.punctuated_prob)
        {
          number_mutations = GeneratePunctuation(punct_params);
          rand_advantage = gsl_ran_flat(gp.rng, 0, 1);
        }

        bool did_count_driver = false;

        // generation of additive rate to the fitness
        if(fit_params.is_randfitness)
        {
          double additional_rate = GenerateFitness(fit_params);
          if (additional_rate > 0)
          {
            new_mut_node->driver_count++;
            new_mut_node->is_driver = true;
            did_count_driver = true;
            if(rand_punct < punct_params.punctuated_prob)
            {
              additional_rate = additional_rate * punct_params.punctuated_multiplier;
              // the additional rate can go to the birth or the death rate
              if( rand_advantage < punct_params.punctuated_advantageous_prob)
              {
                new_mut_node->birth_rate = fmax(0, additional_rate + new_mut_node->birth_rate);
              }
              else
              {
                new_mut_node->death_rate = fmax(0, additional_rate + pnode->death_rate);
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

        cl.tot_rate = cl.tot_rate + new_mut_node->birth_rate + new_mut_node->death_rate;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations = cl.num_mutations + number_mutations;

        cl.currnode = pnode;
        pnode = new_mut_node;
        cl.InsertNode(new_mut_node, number_mutations);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        cl.tot_rate = cl.tot_rate + (pnode->birth_rate) + (pnode->death_rate);
      }

      // Reprioritize clone based on size
      cl.CloneSort(pnode, true);

      flag = true;
      break;
    }
    if( rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) +
      (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      cl.tot_cell_count = cl.tot_cell_count - 1;
      if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, false);
      cl.tot_rate = cl.tot_rate - pnode->birth_rate - pnode->death_rate;

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

    summand = summand + (pnode->cell_count) * (pnode->birth_rate + pnode->death_rate);
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
  double rand_next_event = gsl_ran_flat(gp.rng, 0, cl.tot_rate);
  double rand_mut_occur;

  struct clone *pnode;
  pnode = cl.root;
  bool flag = false;

  while ( (pnode) && !(flag) )
  {
    // Condition for new birth
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {

      rand_mut_occur = gsl_ran_flat(gp.rng, 0, 1);
      // Condition to determine if mutation occurs in new daughter
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
        new_mut_node->death_rate = pnode->death_rate;
        new_mut_node->birth_rate = pnode->birth_rate;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        bool did_count_driver = false;

        // generation of additive rate to the fitness
        if(fit_params.is_randfitness){
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
          }
          new_mut_node->birth_rate = fmax(0, additional_rate + pnode->birth_rate);
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

        cl.tot_rate = cl.tot_rate + new_mut_node->birth_rate + new_mut_node->death_rate;

        cl.tot_cell_count++;
        cl.num_clones++;
        cl.num_mutations++;

        cl.currnode = pnode;
        pnode = new_mut_node;
        cl.InsertNode(new_mut_node, 1);
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
      }
      else
      {
        pnode->cell_count++;
        cl.tot_cell_count++;
        if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, true);
        cl.tot_rate = cl.tot_rate + (pnode->birth_rate) + (pnode->death_rate);
      }

      // Reprioritize clone based on size
      cl.CloneSort(pnode, true);

      flag = true;
      break;
    }
    if( rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) +
      (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      cl.tot_cell_count = cl.tot_cell_count - 1;
      if(gp.count_alleles) cl.ChangeAncestorAllele(pnode, false);
      cl.tot_rate = cl.tot_rate - pnode->birth_rate - pnode->death_rate;

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

    summand = summand + (pnode->cell_count) * (pnode->birth_rate + pnode->death_rate);
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << pnode->parent->birth_rate - pnode->parent->death_rate << "\t" <<
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << pnode->parent->birth_rate - pnode->parent->death_rate << "\t" <<
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
             pnode->clone_time << "\n";
      }
    }
  }
}

void CloneList::SampleAndTraverse(std::ofstream &F, int sim_number, int sample_size, int nsamples)
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
        F << sim_number << "\t" << sample_counter << "\t" <<
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

  if(currnode->nextnode == NULL){

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
