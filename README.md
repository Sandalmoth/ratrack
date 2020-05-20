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
```

```toml
```

Running is a matter of:
```
conda activate ratrack
snakemake results/minimal.pdf results/minimal.fit.csv
snakemake results/minimal.pdf results/minimal.fit.csv
```
Running twice is necessary as the first run sets up the input files correctly.
Optionally, run `snakemake results/minimal.pdf results/minimal.fit.csv -n -r' first, which shows a detailed list of all commands that will be executed. Snakemake has many features which will not be elaborated on here.
