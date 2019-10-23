"""
Plotting functions for various causes
"""


import statistics

import click
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as patches
from pyabc import History
from scipy.stats import linregress
from scipy.stats import f_oneway
import pandas as pd

import simtools


COLORS = ['k', 'r', 'b']


@click.group()
def main():
    """
    Plotting functions for examining output and the ABC fitting process
    """
    pass


# Credits goes to the nice violinplot customization example found at
# https://matplotlib.org/examples/statistics/customized_violin_demo.html on 18-03-29
# (includes also some code in plot function below)
def adjacent_values(vals, q1, q3):
    """
    used in violinplot customization
    """
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


@main.command()
@click.option('-p', '--paramfile', type=click.Path())
@click.option('-o', '--obsfile', type=click.Path())
@click.option('-d', '--dbfile', type=click.Path())
@click.option('--run-id', type=int, default=1)
@click.option('--save', type=click.Path(), default=None)
def abc_info(paramfile, obsfile, dbfile, run_id, save):
    """
    Plots for examining ABC fitting process
    """

    db_path = 'sqlite:///' + dbfile
    abc_history = History(db_path)
    abc_history.id = run_id

    observed = simtools.parse_observations(obsfile)
    simtools.parse_params(paramfile, observed)

    ### PLOTS SHOWING MODEL PROBABILITIES ###
    num_models = abc_history.nr_of_models_alive(0)
    max_points_in_models = max([abc_history.get_distribution(m=x, t=0)[0].shape[1] for x in range(num_models)])

    axs = abc_history.get_model_probabilities().plot.bar()
    axs.set_ylabel("Probability")
    axs.set_xlabel("Generation")
    resolutions = list(range(simtools.PARAMS['abc_params']['resolution_limits'][0],
                             simtools.PARAMS['abc_params']['resolution_limits'][1] + 1))
    axs.legend(resolutions,
               title="Reconstruction resolution")

    if save is not None:
        # first time, construct the multipage pdf
        pdf_out = PdfPages(save)
        pdf_out.savefig()
    else:
        plt.show()

    ### ABC SIMULATION DIAGNOSTICS ###
    fig, ax = plt.subplots(nrows=3, sharex=True)

    t_axis = list(range(abc_history.max_t + 1))

    populations = abc_history.get_all_populations()
    populations = populations[populations.t >= 0]

    ax[0].plot(t_axis, populations['particles'])
    ax[1].plot(t_axis, populations['epsilon'])
    ax[2].plot(t_axis, populations['samples'])

    ax[0].set_title('ABC parameters per generation')
    ax[0].set_ylabel('Particles')
    ax[1].set_ylabel('Epsilon')
    ax[2].set_ylabel('Samples')
    ax[-1].set_xlabel('Generation (t)')
    ax[0].xaxis.set_major_locator(MaxNLocator(integer=True))

    fig.set_size_inches(8, 5)

    if save is not None:
        pdf_out.savefig()
    else:
        plt.show()


    ### PARAMETERS OVER TIME ###
    fig, axs = plt.subplots(nrows=max_points_in_models, sharex=True, sharey=True)

    t_axis = np.arange(abc_history.max_t + 1)
    print(t_axis)
    # parameters = ['birthrate.s0.d', 'birthrate.s0.r0']
    all_parameters = [list(abc_history.get_distribution(m=m, t=0)[0].columns)
                  for m in range(num_models)]
    # abc_data, __ = abc_history.get_distribution(m=m, t=generation)
    parameters = []
    for x in all_parameters:
        for y in x:
            parameters.append(y)
    parameters = list(set(parameters))
    parameters = sorted(parameters, key=lambda x: x[-1])
    print(parameters)

    for m in range(num_models):

        qs1 = {param: [np.nan for __ in t_axis] for param in parameters}
        medians = {param: [np.nan for __ in t_axis] for param in parameters}
        qs3 = {param: [np.nan for __ in t_axis] for param in parameters}

        for i, generation in enumerate(t_axis):
            abc_data, __ = abc_history.get_distribution(m=m, t=generation)
            data = {x: np.array(abc_data[x]) for x in parameters if x in abc_data}
            for k, v in data.items():
                t_q1, t_m, t_q3 = np.percentile(
                    v, [25, 50, 75]
                )
                qs1[k][i] = t_q1
                medians[k][i] = t_m
                qs3[k][i] = t_q3


        for i, param in enumerate(parameters):
            # if len(medians[param]) == 0:
            if not medians[param]:
                continue
            print(t_axis, medians[param])
            axs[i].plot(t_axis, medians[param], color=COLORS[m])
            axs[i].fill_between(t_axis, qs1[param], qs3[param], color=COLORS[m], alpha=0.2)

            axs[i].set_ylabel(param[10:])

        axs[-1].set_xlabel('Generation (t)')

    if save is not None:
        pdf_out.savefig()
    else:
        plt.show()


@main.command()
@click.option('-p', '--paramfile', type=click.Path())
@click.option('-o', '--obsfile', type=click.Path())
@click.option('-d', '--dbfile', type=click.Path())
@click.option('--run-id', type=int, default=1)
@click.option('--save', type=click.Path(), default=None)
def result_single(paramfile, obsfile, dbfile, run_id, save):
    """
    Plot the result of a single fitting
    """

    db_path = 'sqlite:///' + dbfile
    abc_history = History(db_path)
    abc_history.id = run_id

    observed = simtools.parse_observations(obsfile)
    id_str = next(iter(observed))
    simtools.parse_params(paramfile, observed)

    # violin plot of results
    max_gen = abc_history.max_t

    num_models_total = abc_history.nr_of_models_alive(0)
    num_models_final = abc_history.nr_of_models_alive(max_gen)
    max_point_in_models = max([abc_history.get_distribution(m=x, t=max_gen)[0].shape[1]
                               for x in range(num_models_final)])

    # fig, axs = plt.subplots(ncols=num_models_final, sharey=True, sharex=True)
    # fig.set_size_inches(num_models_final*3, 3)

    if save is not None:
        # first time, construct the multipage pdf
        pdf_out = PdfPages(save)

    for j in range(num_models_total):
        model_prob = abc_history.get_model_probabilities()[j][max_gen]
        print(model_prob)
        if model_prob == 0.0:
            continue
        fig, axs = plt.subplots()
        fig.set_size_inches(4, 3)
        end_time = simtools.PARAMS['end_time'][id_str]()
        print(end_time)

        df, w = abc_history.get_distribution(m=j, t=max_gen)
        print(df)
        print(df.columns)
        # abc_data = [sorted(df['birthrate.b' + str(x)]) for x in range(df.shape[1])]
        time_axis = np.linspace(0, end_time, len(list(df.columns)))

        for x in list(df.columns):
            print(x)
            print(df[x])
        abc_data = [sorted(df[x]) for x in list(df.columns)]
        print(abc_data)

        violinparts = axs.violinplot(abc_data, positions=time_axis,
                                        widths=end_time/(max_point_in_models + 1)*0.8,
                                        showmeans=False, showmedians=False, showextrema=False)
        for part in violinparts['bodies']:
            part.set_facecolor('lightgrey')
            part.set_alpha(1)
        quartile1, medians, quartile3 = np.percentile(abc_data, [25, 50, 75], axis=1)
        whiskers = np.array([
            adjacent_values(sorted_array, q1, q3)
            for sorted_array, q1, q3 in zip(abc_data, quartile1, quartile3)])
        whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
        axs.scatter(time_axis, medians, marker='.', color='white', s=30, zorder=3)
        axs.vlines(time_axis, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
        axs.vlines(time_axis, quartile1, quartile3, color='k', linestyle='-', lw=5)

        birthrate = [statistics.median(x) for x in abc_data]
        axs.plot(time_axis, birthrate, color='k')
        axs.set_xlabel('Time [days]')
        axs.set_ylabel('Growth rate [doublings/day]')

        # axs.set_ylim(0, simtools.PARAMS['abc_params']['rate_limits'][1])

        plt.tight_layout()

        if save is not None:
            pdf_out.savefig()
        else:
            plt.show()


    # fit against timeline
    for j in range(num_models_total):
        model_prob = abc_history.get_model_probabilities()[j][max_gen]
        if model_prob == 0.0:
            continue
        fig, axs = plt.subplots()
        fig.set_size_inches(4, 3)
        end_time = simtools.PARAMS['end_time'][id_str]()

        df, w = abc_history.get_distribution(m=j, t=max_gen)
        time_axis = np.linspace(0, end_time, len(list(df.columns)))

        samplings = [simtools.get_samplings_dilutions(observed[id_str], x)[0]
                     for x, __ in enumerate(observed[id_str]['time'])]
        dilutions = [simtools.get_samplings_dilutions(observed[id_str], x)[1]
                     for x, __ in enumerate(observed[id_str]['time'])]

        samplings = list(zip(*samplings))
        dilutions = list(zip(*dilutions))

        abc_data = [sorted(df[x]) for x in list(df.columns)]
        for k, v in observed.items():
            measured = np.array(v['count'])
            for s in samplings:
                measured /= s
            for d in dilutions:
                measured *= d

            axs.scatter(v['time'], measured, marker='.', color='k')

        print(samplings, dilutions)

        simulations = None

        i = 0
        for index, row in df.iterrows():
            # print(index, row)
            time, size, rate = simtools.simulate_timeline(
                simtools.PARAMS['starting_population'][id_str](),
                simtools.PARAMS['end_time'][id_str](),
                list(row),
                0,
                0,
                simtools.PARAMS['simulation_params']['deathrate_interaction'],
                # 'bernoulli',
                'rar-engine',
                verbosity=1
            )

            if simulations is None:
                simulations = np.zeros((len(size), len(df)))

            simulations[:, i] = size
            i += 1

        qt1, qt2, qt3 = np.quantile(simulations, (0.05, 0.5, 0.95), axis=1)
        print(qt2)

        # axs.plot(time, qt1)
        axs.plot(time, qt2, color='k')
        # axs.plot(time, qt3)
        axs.fill_between(time, qt1, qt3, zorder=-1, color='lightgray')

        axs.set_xlabel('Time [days]')
        axs.set_ylabel('Population count measure')

        plt.tight_layout()

        if save is not None:
            pdf_out.savefig()
        else:
            plt.show()


    if save is not None:
        pdf_out.close()

if __name__ == '__main__':
    main()
