A simplified _Plasmodium falciparum_ epidemiological model with detailed descriptions of parasite genome recombination dynamics.

#### Installation

`python setup.py install` from this directory

Or to do active development on the project:

`python setup.py develop`

#### Details

`SNP` holds the properties of a single nucleotide polymorphism.

`Genome` is a discretized representation of SNPs on chromosomes.

`Infection` holds the properties of a human infected with one or more parasite strains, including the dynamics of parasite propagation through the mosquito vector into a new human host.

`HumanCohort` and `HumanIndividual` are the properties of humans tracked either as cohorts (for the susceptible population) or as individuals (to track individual infections).

`Population` is an object containing a list of individuals and their dynamics.

`Migrating` and `MigrationInfo` are, respectively, the properties of an individual migration event and the characteristics of migration rates for a specific population.

`Simulation` is the overall handler of setting up and parameterizing the initial conditions and dynamics of a parasite genetic-epidemiological simulation.

#### Dependencies

Some of the example output-processing scripts may require the pandas/numpy/matplotlib packages.
