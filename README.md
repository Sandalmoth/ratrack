# ratrack
Ratrack is a tool for deriving (changing) growth rates from population counts. It's intended use is to analyse how cell growth rate changes in a long term suspension cell culture undergoing some kind of adaptive process which gradually changes growth rate over time. For input, it requires only population size measurements, and can work with as little as four samples.

## Requirements and Installation
ratrack is implemented in python and C++. It relies on snakemake for automation. A conda environment file `environmennt.yml` is included which holds all requirements for building the C++ sources and running the program.

Using a conda environment to manage the requirements, installation looks like this
```
git clone https://github.com/Sandalmoth/ratrack
cd ratrack
conda env create -n ratrack -f environment.yml
mkdir intermediate results code/bin
cd code/bin
cmake .. && make
```

## Basic usage
### Running the minimal example
A minimal example is provided in the data folder named `minimal.csv` and `minimal.toml`

```csv
name,time,count,sample1,sample2
minimal,0.0,19,0.01,0.02
minimal,2.0,105,0.01,0.02
minimal,4.0,403,0.01,0.02
minimal,6.0,529,0.01,0.02
minimal,8.0,591,0.01,0.02
```

```toml
[simulation_params]
carrying_capacity = 3e6

[abc_params]
# range to be considered for growth rates
rate_limits = [0.01, 3.0]
# range of the number of control points in piecewise linear function
resolution_limits = [1, 3]
parallel_simulations = 4
simulator = 'bernoulli'
birthrate_coupling_sets = [
    'minimal',
]

[[filters]]
name = 'gauss-multiplicative'
mean = 1.0
sigma = 0.05
```

Running is a matter of:
```
conda activate ratrack
snakemake results/minimal.pdf results/minimal.fit.csv
snakemake results/minimal.pdf results/minimal.fit.csv
```
Running twice is necessary as the first run sets up the input files correctly.
Optionally, run `snakemake results/minimal.pdf results/minimal.fit.csv -n -r' first, which shows a detailed list of all commands that will be executed. Snakemake has many features which will not be elaborated on here.

## Second example
There is not too much more to know, but, the example in `demo.csv` and `demo.toml` explore some other useful features. I will not reproduce the whole of the `.csv` as it is largely similar. The `.toml` however is worth looking at.

```toml
[simulation_params]
carrying_capacity = 3e6

[abc_params]
# number of particles in ABC simulation
starting_population_size = 100
# max number of generations in SMC
max_populations = 10
rate_limits = [0.01, 3.0]
resolution_limits = [2, 2]
parallel_simulations = 4
# using the stochastic simulator instead
simulator = 'rar-engine'
# we are here grouping counts that were simulated from the same timeline
# as if we counted twice, or perhaps, the experiments are such that we
# do not expect the grouped dataseries to diverge, with the differences
# being primarily between groups
birthrate_coupling_sets = [
    ['demo0.a', 'demo0.b'],
    ['demo1.a', 'demo1.b'],
    ['demo2.a', 'demo2.b'],
]

# The plot y label can be changed (though cells is default)
[plot_params]
population_measure = 'Cells'

# This data was simulated without additional sampling noise
[[filters]]
name = 'copy'
```
Running works the same as before (though replacing the `minimal` with `demo`). The main difference here is in having several timelines at once, and grouping some of them together.
