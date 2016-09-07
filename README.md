SIApop
========

SIApop simulates birth-death-mutation processes with mutations having random fitnesses to simulate clonal evolution.

Implementation is in Terminal with the following command:

    ./SIApop -in ./inputfile.txt -anc ./ancestors.txt -out ./results/


Features
--------

- Easily customizable to fit new models
- Import or use default initial conditions

Requirements
------------

- C++ compiler (clang or g++)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
  * Available as a Homebrew package for OS X users

    OR
  * Download the latest [GSL library](http://ftpmirror.gnu.org/gsl/)


Installation
------------

Install SIApop in Terminal by locating the SIApop directory and running:

    make

Two executables should be created in SIApop/bin/, reflecting time-homogeneous and inhomogeneous simulations

    /SIApop/bin/SIApop
    /SIApop/bin/SIApop-td


Using SIApop
------------
To run in Terminal, navigate to /SIApop/ and type

    ./SIApop
  or

    ./SIApop-td

A simple birth-death process should run using the default parameters and outputted to the current directory. Using a different set of parameters requires specifying an input file. Other options include specifying an ancestor file and an output directory. These options are added to the Terminal command with the following flags:

  * -in inputfile.txt
  * -out output/directory/
  * -anc ancestors.txt

A fully functioning program run including creating a results folder would look like:

    mkdir results
    ./SIApop -in ./inputfile.txt -out ./results -anc ./ancestors.txt

Input File
----------
The input file is a 2-column tab-delimited file containing the following arguments. A hashbang (#) indicates a comment and the line is ignored. The variable name is in the left column and the value is in the right. Defaults are provided to create a minimal working model, so no inputs are required.

| Variable Name       | Variable Type | Description |
| ------------------- | ------------- | ----------- |
| tot_life            | double        | total lifetime of branching process |
| max_pop             | double        | maximum population to end process at|
| start_time          | double        | starting time (typically 0)|
| ancestors           | int           | number of individuals per ancestor clone (if ancestor file not provided)|
| ancestor_clones     | int           | number of initial types |
| num_sims            | int           | total number of simulations of the same process|
| allow_extinction    | boolean (1/0) | 1 if allow a simulation to go extinct |
| detection_threshold | double        | minimum proportion of the total population such that a clone is output|
| num_samples         | int           | number of samples to take|
| sample_size         | int           | size of each sample|
| birth_rate          | double > 0    |  starting birth rate |
| death_rate          | double > 0    | starting death rate |
| mutation_prob       | double [0, 1] | default mutation probability for new clone|
| trace_ancestry      | bool (1/0)    | Track info on parent of each clone |
| count_alleles       | bool (1/0)    | adds/subtracts and individual to allele_count of individual and all ancestors |

  FITNESS DISTRIBUTION PARAMETERS

| Variable Name | Variable Type | Description |
| ------------- | ------------- | ----------- |
| alpha_fitness | double > 0             | exponential distribution parameter for positive side of fitness distribution |
| beta_fitness  | double > 0             | exponential distribution parameter for negative side of fitness distribution|
| pass_prob     | double [0,1]           | probability that additional fitness of new mutant is 0|
| upper_fitness | double                 | upper bound to fitness distribution|
| lower_fitness | double <= upper_fitness | lower bound to fitness distribution|



MUTATION DISTRIBUTION PARAMETERS

| Variable Name  | Variable Type | Description |
| -------------- | ------------- | ----------- |
| alpha_mutation | double > 0    | alpha parameter for Beta distribution for additional mutation probability in new mutant clone|
| beta_mutation  | double > 0    | beta parameter for Beta distribution for additional mutation probability in new mutant clone|

PUNCTUATED EQUILIBRIUM PARAMETERS

| Variable Name                 | Variable Type | Description |
| ----------------------------- | ------------- | ----------- |
| punctuated_prob               | double [0,1]  | probability of mutation burst |
| poisson_param                 | double > 0    | rate parameter for zero-truncated Poisson distribution number of mutations in burst |
| punctuated_fitness_multiplier | double        | amount to multiply additional fitness by |
| punctuated_advantageous_prob  | double [0,1]  | probability that burst affects birth rate instead of death rate |


EPISTATIC PARAMETERS

| Variable Name               | Variable Type | Description |
| --------------------------- | ------------- | ----------- |
|epistatic_mutation_threshold | int > 0       | number of mutation required before burst in fitness due to epistasis|
|epistatic_multiplier         | double        | amount to multiply fitness contribution in new clone by due to epistasis occurring|


TIME-DEPENDENT PARAMETERS

| Variable Name  | Variable Type                       | Description |
| -------------- | ----------------------------------- | ----------- |
|birth_function  | 0, 1, 2, 3, 4                       | see below |
|death_function  | 0, 1, 2, 3, 4                       | see below |
|td_birth_params | vector of doubles (space-delimited) | see below |
|td_death_params | vector of doubles (space-delimited) | see below |

Ancestor File
-------------

The ancestor file is a tab-delimited file with the same structure as the output.
The first line contains variable names and each line contains information for a
single clone to serve as an ancestor population. The only requirement for this
file is a column containing the number of cells. If not provided, the program
will look at the arguments ancestors and ancestor_clones described above to run
with nonunique clones. If those are not provided a default of a single ancestor
individual is used. The following table describes the possible variables for the
ancestor file. Some parameters below are included since they are present in the
output and store information when continuing a previous simulation.

PARAMETERS

| Variable Name | Variable Type | Description |
| ------------  | ------------- | ----------- |
| unique_id     | string        | id for each ancestor |
| numcells      | int           | the number of cells for the ancestor |
| mutprob       | double [0,1]  | the probability of initiating a new clone given a birth occurs |

TIME-DEPENDENT PARAMETERS

| Variable Name  | Variable Type                       | Description |
| -------------- | ----------------------------------- | ----------- |
|birth_function  | 0, 1, 2, 3, 4                       | see below |
|death_function  | 0, 1, 2, 3, 4                       | see below |
|bf_params       | vector of doubles (space-delimited) | see below |
|df_params       | vector of doubles (space-delimited) | see below |


Time-Dependent Rate Functions
------------------------------
The functions all are predefined and allow the user to provide a list of
parameters. The parameters should be listed the same way regardless of which
function is used in the form of “x1,x2,x3,x4,x5…” in the tab-delimeted file
(see example). An extra function is included called "custom" to allow the user
to define a unique function, but recompiling the program is required after.
The curves are parameterized as follows:

| Function Parameter  | Function Name     | Mathematical Description |
| ------------------- | ----------------- | ------------------------ |
| 0                   | constant          | ![constant](README/README-image-1.png) |             |
| 1                   | linear            | ![linear](README/README-image-2.png)  |
| 2                   | logistic          | ![logistic](README/README-image-3.png)  |
| 3                   | Gompertz growth   | ![Gompertz](README/README-image-4.png)  |
| 4                   | Custom            | Include own parameters  |


Example
--------

Examples of input and ancestor files for both the time-homogeneous (SIApop) and
the time-inhomogeneous (SIApop-td) are provided. The following is a complete
step-by-step workthrough and results for a time-inhomogeneous process with 5
ancestors having different rates. The files are found in

    /examples/inhomogeneous/workthrough:

1. In Terminal, navigate to folder with the exectuables, and run the following command:

        ./SIApop-td -out ./examples/inhomogeneous/workthrough -in ./examples/inhomogeneous/workthrough/inputfile.txt -anc ./examples/inhomogeneous/workthrough/ancestors.txt

2.	Three new files should appear with named clonedata.txt, sim_stats.txt, timedata.txt. These contain information about all clones at the end of the simulations (clonedata.txt), information about the simulation (sim_stats.txt), and time-course data about clone counts (timedata.txt). Example output is included from a previous run with the prefix “ex_”.

3.	Open R and change code to point to source directory. View population information by importing data and running the provided function on your own.



Contribute
----------

- Issue Tracker: github.com/SIApop/SIApop/issues
- Source Code: github.com/SIApop/SIApop

Support
-------

If you are having issues, please let us know.

License
-------

The project is licensed under the ... license.
