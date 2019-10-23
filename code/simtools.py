"""
Shared tools for abc simulation and analysis
"""

import copy
import csv
# import statistics
import subprocess
import sys
from io import StringIO

import numpy as np
import toml


FORWARD_SAMPLING = 'RV'
BACKWARD_SAMPLING = 'MLE'


# simulate a lb-process using the given parameters with external software
# n - starting number of cells
# t - time when to end (time starts at 0), so simulation length
# b - birth rate, can be a number or a list. In case of list, it is evenly spread over the timeline
#     with interpolation
# d - death rate
# p - birth rate interaction reduction. part of quadratic term in logistic equation applied to
#     reducing birth rate first. Overflow goes to increasing death rate (birth rate has to be > 0)
# q - interaction death rate. Complementary part of quadratic term that just works by increasing
#     death rate.
# For the purposes of this simulation the p/q balance likely doesn't matter much, but it is
# included for completeness.
# Returns a simulated timeline, calculated using c++ software.
# (Stochastic simulation using next reaction method)
# returns 3 vectors, time, size, rate
# time - time point for this datapoint
# size - size at that timepoint
# rate - growth rate at that timepoint (nice for visualizing interpolation)
def simulate_timeline(starting_population,
                      end_time,
                      birthrates,
                      deathrate,
                      birthrate_interaction,
                      deathrate_interaction,
                      simulator,
                      verbosity=0):
    """
    Simulate a lb-process using external software
    """
    # sanity checking
    assert starting_population > 0
    assert end_time >= 0
    if simulator != 'bernoulli':
        # bernoulli simulator handles negative rates without issues
        # NOTE remember to update this if if neccessary
        for birthrate in birthrates:
            assert birthrate >= 0
    assert deathrate >= 0
    assert birthrate_interaction >= 0
    assert deathrate_interaction >= 0
    # run external
    cmd = 'code/bin/' + simulator + \
          ' -n ' + str(starting_population) + \
          ' -t ' + str(end_time) + \
          ' -b \'' + str(birthrates) + '\'' \
          ' -d ' + str(deathrate) + \
          ' -p ' + str(birthrate_interaction) + \
          ' -q ' + str(deathrate_interaction)
    if verbosity > 0:
        print(cmd)
    output = subprocess.getoutput(cmd)
    if verbosity > 1:
        print(output)
    # output from rar-engine is actually a .tsv file (printed in stdout)
    # use some trickery to parse it with csv
    buff = StringIO(output)
    time = []
    size = []
    rate = []

    rdr = csv.DictReader(buff, dialect='excel-tab')
    try:
        for line in rdr:
            time.append(float(line['time']))
            size.append(int(float(line['size'])))
            # Casting like this can maybe lose precision.
            # But python3 doesn't want to construct integers from scientific notation,
            # whereas the float-constructor handles anything TODO fix?
            rate.append(float(line['rate']))
    except ValueError:
        print('Timeline simulation does not conform to standard', file=sys.stderr)
        print(output, file=sys.stderr)
        exit(1)

    if verbosity > 2:
        print([x for x in zip(time, size, rate)])

    return np.array(time), np.array(size), np.array(rate)


def apply_noise(size, filters):
    """
    Apply list of noise filters in order.
    Implemented types:
      copy: filter does nothing
      perfect: simulates perfect sampling
      poisson: simulates random sampling (in any number of steps)
      gauss-multiplicative: gaussian noise with constant COV
      gauss-additive: gaussian noise with constant stdev
    """
    for filt in filters:
        if filt['name'] == ['copy']:
            continue
        elif filt['name'] == 'perfect':
            size *= filt['sample']
            # size = [x * y for x, y in zip(size, filt['sample'])]
        elif filt['name'] == 'poisson':
            size = np.random.poisson(size * filt['sample'])
            # size = [np.random.poisson(x * y) for x, y in zip(size, filt['sample'])]
        elif filt['name'] == 'gauss-multiplicative':
            size *= np.random.normal(filt['mean'], filt['sigma'], size.size)
            # size = [np.random.normal(filt['mean'], filt['sigma']) * x for x in size]
        elif filt['name'] == 'gauss-additive':
            size += np.random.normal(filt['mean'], filt['sigma'], size.size)
            # size = [np.random.normal(filt['mean'], filt['sigma']) + x for x in size]

    return np.round(size)


#     """
#     holder of parameters for lb process
#     """
#     def __init__(self, starting_population, time,
#                  death_interaction, simulator):
#         self.starting_population = starting_population
#         self.end_time = time
#         self.death_interaction = death_interaction
#         self.simulator = simulator
#         self.filters = []
#         # TODO add sampling methods to parameter parsing
#         self.forward_sampling_method = 'poisson-RV'
#         self.backward_sampling_method = 'poisson-MLE'
#         self.starting_population_samplings = []
#         self.starting_population_dilutions = []

#     def get_starting_population(self):
#         """
#         Starting population can be a constant, or an R.V.
#         """
#         if self.forward_sampling_method == 'poisson-MLE' and \
#            self.backward_sampling_method == 'poisson-MLE':
#             starting_population = self.starting_population
#             for sample in self.starting_population_samplings:
#                 starting_population /= sample
#             for dilution in self.starting_population_dilutions:
#                 starting_population *= dilution
#             starting_population = int(starting_population)
#             return starting_population

#         elif self.forward_sampling_method == 'poisson-RV' and \
#              self.backward_sampling_method == 'poisson-MLE':
#             starting_population = self.starting_population
#             for sample in self.starting_population_samplings:
#                 starting_population /= sample
#             for dilution in self.starting_population_dilutions:
#                 starting_population = np.random.poisson(starting_population * dilution)
#             starting_population = int(starting_population)
#             return starting_population

#         sys.exit("Invalid combination of forward and backward sampling methods")


# # defaults, see also reconstruct() where some are changed
# LB_PARAMS = LBParameters(100000, 7.0, 0.3333e-7, 'rar-engine')
# PARAMS = {}

PARAMS = {}

OBSERVED = {}

def get_samplings_dilutions(observed, observation=0):
    """
    Get the list of all samplings and dilutions done to particular observation
    """
    samplings = []
    dilutions = []
    i = 1
    while True:
        if 'sample' + str(i) in observed.keys():
            samplings.append(observed['sample' + str(i)][observation])
        else:
            break
        i += 1
    i = 1
    while True:
        if 'dilute' + str(i) in observed.keys():
            dilutions.append(observed['dilute' + str(i)][observation])
        else:
            break
        i += 1
    return samplings, dilutions


def apply_sampling(size, samplings, dilutions):
    """
    apply sampling methods to simulated data to make it comparable to observations
    """

    if BACKWARD_SAMPLING == 'MLE':
        for dilution in dilutions:
            dilution = np.array(dilution)
            size /= dilution
            # size = [x / y for x, y in zip(size, dilution)]
    else:
        sys.exit("Unsupported backward sampling method")

    if FORWARD_SAMPLING == 'MLE':
        for sample in samplings:
            sample = np.array(sample)
            size *= sample
            # size = [x * y for x, y in zip(size, sample)]
    elif FORWARD_SAMPLING == 'RV':
        for sample in samplings:
            sample = np.array(sample)
            size = np.random.poisson(size*sample)
            # size = [np.random.poisson(x * y) for x, y in zip(size, sample)]
    else:
        sys.exit("Unsupported forward sampling method")

    return np.round(size)


def parse_observations(infile):
    """
    Parse csv of observations
    Creates a dictionary of dictionaries
    The inner ones hold single wells/colonies/whatever
    the outer ones holds that data coupled with their names
    """

    observed = {}

    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for obs in rdr:
            id_string = '.'.join([obs[x] for x in ['name', 'generation', 'well']])
            if id_string not in observed:
                observed[id_string] = {str(x): [] for x in rdr.fieldnames}
            for k in rdr.fieldnames:
                entry = obs[k]
                if entry:
                    if k in ['name', 'well']:
                        observed[id_string][k].append(str(entry))
                    elif k in ['generation', 'birthrate_group', 'deathrate_group']:
                        observed[id_string][k].append(int(entry))
                    else:
                        observed[id_string][k].append(float(entry))
                else:
                    if k in ['count', 'dead']:
                        observed[id_string][k].append(None)
                    else:
                        observed[id_string][k].append(1.0)

    #         for k, v in observed[id_string].items():
    #             if type(v[0]) in [int, float]:
    #                 observed[k] = np.array(v)

    # print(observed)

    global OBSERVED
    OBSERVED = copy.deepcopy(observed)
    return observed


def parse_params(paramfile, observed=None):
    """
    parse toml parameter file and observed data for lb-process parameters that are not the birthrate
    """
    # set model parameters
    # starting population is mean of all t=0 observations
    global PARAMS
    PARAMS = toml.load(paramfile)
    if PARAMS['simulation_params']['starting_cell_count'] == 'calculate':
        if observed is None:
            sys.exit("Cannot compute starting_cell_count: 'MLE' without observations")
        PARAMS['starting_population'] = {}
        for id_string, obs in observed.items():
            samplings, dilutions = get_samplings_dilutions(obs)
            pop = obs['count'][0]
            if FORWARD_SAMPLING == 'MLE' and BACKWARD_SAMPLING == 'MLE':
                for sample in samplings:
                    pop /= sample
                for dilution in dilutions:
                    pop *= dilution
                PARAMS['starting_population'][id_string] = lambda x=int(pop): x
            if FORWARD_SAMPLING == 'RV' and BACKWARD_SAMPLING == 'MLE':
                def f(x=pop, y=copy.deepcopy(samplings), z=copy.deepcopy(dilutions)):
                    # print(x, y, z)
                    for sample in y:
                        x /= sample
                    for dilution in z:
                        x = np.random.poisson(x * dilution)
                    return int(x)
                PARAMS['starting_population'][id_string] = f

    else:
        PARAMS['starting_population'] = lambda x=int(PARAMS['simulation_params']['starting_cell_count']): x
    # no need to simulate longer than observed segment
    PARAMS['end_time'] = {}
    if PARAMS['simulation_params']['end_time'] == 'max_observed':
        for id_string, obs in observed.items():
            if observed is None:
                sys.exit("Cannot compute end_time: 'max_observed' without observations")
            PARAMS['end_time'][id_string] = lambda x=max(obs['time']): x
    else:
        PARAMS['end_time'][id_string] = lambda x=float(PARAMS['end_time']): x

