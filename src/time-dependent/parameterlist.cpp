#include "parameterlist.h"



void ParameterList::ParseVector(const std::string& s, std::vector<double>& v)
{
    std::istringstream s2((*this)[s]);
    double tmp;

    while(s2 >> tmp) v.push_back(tmp);
}

void ParameterList::splitAndFill(const std::string& s)
{
    std::string::size_type pos1 = s.find_first_of(" ,=\t");
    std::string::size_type pos2 = s.find_first_not_of(" ,=\t", pos1+1);
    if (pos2 == std::string::npos) return;

    std::string k = s.substr(0, pos1);
    std::string val =  s.substr(pos2);

    if (!val.compare("NA")) return;

    (*this)[k] = val;
}

void ParameterList::init()
{
    //default values
    insert(std::make_pair("tot_life", "40000"));
    insert(std::make_pair("max_pop", "10000"));
    insert(std::make_pair("start_time", "0"));
    insert(std::make_pair("ancestors", "1"));
    insert(std::make_pair("ancestor_clones", "1"));
    insert(std::make_pair("num_sims", "1"));
    insert(std::make_pair("allow_extinction", "1"));
    insert(std::make_pair("num_samples", "0"));
    insert(std::make_pair("sample_size", "0"));
    insert(std::make_pair("detection_threshold", "0"));
    insert(std::make_pair("birth_rate", "1.5"));
    insert(std::make_pair("death_rate", "1"));
    insert(std::make_pair("mutation_prob", "0.00"));
    insert(std::make_pair("alpha_fitness", "0"));
    insert(std::make_pair("beta_fitness", "0"));
    insert(std::make_pair("pass_prob", "1"));
    insert(std::make_pair("upper_fitness", "0"));
    insert(std::make_pair("lower_fitness", "0"));
    insert(std::make_pair("alpha_mutation", "0"));
    insert(std::make_pair("beta_mutation", "0"));
    insert(std::make_pair("trace_ancestry", "0"));
    insert(std::make_pair("count_alleles", "0"));
    insert(std::make_pair("punctuated_prob", "0"));
    insert(std::make_pair("poisson_param", "1"));
    insert(std::make_pair("punctuated_fitness_multiplier", "1"));
    insert(std::make_pair("punctuated_advantageous_prob", "0.1"));
    //insert(std::make_pair("epistatic_driver_prob", "0"));
    insert(std::make_pair("epistatic_mutation_thresh", "1"));
    insert(std::make_pair("epistatic_multiplier", "1"));
    insert(std::make_pair("B_max", "1.00")); // Max value for birth rate
    insert(std::make_pair("D_max", "1.00")); // Max value for death rate
    insert(std::make_pair("birth_function", "0"));
    insert(std::make_pair("death_function", "0"));
    insert(std::make_pair("td_birth_params", "1 0 1"));
    insert(std::make_pair("td_death_params", "1 0 1"));
}
