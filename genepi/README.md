A simplified P.falciparum epidemiological model with detailed descriptions of   parasite genome recombination dynamics.

#### Installation

`python setup.py install` from this directory

#### Details

`SNP` holds the properties of a single nucleotide polymorphism.

`Genome` is a discretized representation of SNPs on chromosomes.

`Infection` holds the properties of an individual infected with one or more parasite strains, including the dynamics of parasite propagation through the mosquito vector into a new human host.

`Population` is an object containing a list of infections and their dynamics.

`MigratingIndividual` and `MigrationInfo` are, respectively, the properties of an individual migration event and the characteristics of migration rates for a specific population.

`Simulation` is the overall handler of setting up and parameterizing the initial conditions and dynamics of a parasite genetic-epidemiological simulation.

#### Dependencies

Some of the example output-processing scripts may require the pandas/numpy/matplotlib packages.
