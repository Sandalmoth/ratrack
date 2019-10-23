"""
csv manipulation tools
"""

import csv

import click
import toml


@click.group()
def main():
    """
    tools for manipulating csv files during data preprocessing
    """
    pass


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
def zero_time_split(infile, outfile):
    """
    Set minimum time to 0 and adjust other times to match
    """
    with open(infile, 'r') as in_csv:
        with open(outfile, 'w') as out_csv:
            rdr = csv.DictReader(in_csv)
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            times = []
            for line in rdr:
                times.append(float(line['time']))
            min_time = min(times)
            # go back to start for copying with modification
            in_csv.seek(0)
            rdr.__next__() # skip header
            wtr.writeheader()
            for line in rdr:
                line['time'] = float(line['time']) - min_time
                wtr.writerow(line)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
def zero_time_longform(infile, outfile):
    """
    set minimum time to 0 (and adjust others accordingly)
    treats each set with identical name, well and generation as one dataset
    """
    min_times = {}
    with open(infile, 'r') as in_csv:
        with open(outfile, 'w') as out_csv:
            rdr = csv.DictReader(in_csv)
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            # times = []
            for line in rdr:
                # times.append(float(line['time']))
                id_string = '.'.join([line[x] for x in ['name', 'generation', 'well']])
                if id_string not in min_times:
                    min_times[id_string] = []
                min_times[id_string].append(float(line['time']))
            # min_time = min(times)
            min_times = {k: min(v) for k, v in min_times.items()}
            # go back to start for copying with modification
            in_csv.seek(0)
            rdr.__next__() # skip header
            wtr.writeheader()
            for line in rdr:
                id_string = '.'.join([line[x] for x in ['name', 'generation', 'well']])
                line['time'] = float(line['time']) - min_times[id_string]
                wtr.writerow(line)



@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfolder', type=click.Path())
def split_longform(infile, outfolder):
    """
    decompose a single longform csv into one small aptly named file for each cell colony
    """
    data = {}
    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for line in rdr:
            idx = line['name'] + line['well'] + 'g' + line['generation']
            if idx not in data:
                data[idx] = []
            data[idx].append(line)

    for dataset in data.values():
        filename = outfolder + '.'.join([dataset[0]['name'],
                                         dataset[0]['well'],
                                         'g' + dataset[0]['generation'],
                                         'csv'])
        with open(filename, 'w') as out_csv:
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            wtr.writeheader()
            for line in dataset:
                del line['name']
                del line['well']
                del line['generation']
                wtr.writerow(line)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
@click.option('-n', 'name', type=str)
@click.option('-p', 'position', type=int)
@click.option('-v', 'init_val', type=str)
def add_column(infile, outfile, name, position, init_val):
    """
    add a column named 'name' to a csv
    position if given will be the position of the new column
    """
    with open(infile, 'r') as in_csv:
        with open(outfile, 'w') as out_csv:
            rdr = csv.DictReader(in_csv)
            fieldnames = rdr.fieldnames[:]
            if name in fieldnames:
                print("column with given name already in csv")
                return
            if position is not None:
                fieldnames.insert(position, name)
            else:
                fieldnames.append(name)
            wtr = csv.DictWriter(out_csv, fieldnames=fieldnames)
            wtr.writeheader()
            for line in rdr:
                line[name] = init_val
                wtr.writerow(line)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-p', 'paramfile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
def define_groups(infile, paramfile, outfile):
    """
    define the coupling groups
    """

    data = {}

    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for l in rdr:
            id_string = '.'.join([l[x] for x in ['name', 'generation', 'well']])
            if id_string not in data:
                data[id_string] = {str(x): [] for x in rdr.fieldnames}
            for k, v in l.items():
                data[id_string][k].append(v)

    names = sorted(list({x['name'][0] for __, x in data.items()}))
    generations = sorted(list({x['generation'][0] for __, x in data.items()}))
    wells = sorted(list({x['well'][0] for __, x in data.items()}))

    print(names, generations, wells)

    params = toml.load(paramfile)

    birthrate_coupling = params['abc_params']['birthrate_coupling']

    running_birthrate_group = 0

    for id_string, obs in data.items():
        if birthrate_coupling == 'all':
            obs['birthrate_group'] = 0
        elif birthrate_coupling == 'name':
            obs['birthrate_group'] = [
                obs['name'][0] in x for x in params['abc_params']['birthrate_coupling_sets']
            ].index(True)
        elif birthrate_coupling == 'generation':
            coupling_groups = len(params['abc_params']['birthrate_coupling_sets'])
            obs['birthrate_group'] = [
                obs['generation'][0] in x for x in params['abc_params']['birthrate_coupling_sets']
            ].index(True) + \
            coupling_groups*names.index(obs['name'][0])
        elif birthrate_coupling == 'well':
            coupling_groups = len(params['abc_params']['birthrate_coupling_sets'])
            obs['birthrate_group'] = [
                obs['well'][0] in x for x in params['abc_params']['birthrate_coupling_sets']
            ].index(True) + \
            coupling_groups*generations.index(obs['generation'][0]) + \
            len(generations)*coupling_groups*names.index(obs['name'][0])
        elif birthrate_coupling == 'none':
            obs['birthrate_group'] = running_birthrate_group
            running_birthrate_group += 1
        else:
            print('Unknown birthrate_coupling:', birthrate_coupling)

    print(data)
    with open(outfile, 'w') as out_csv:
        fieldnames = list(data[list(data.keys())[0]].keys())
        wtr = csv.DictWriter(out_csv, fieldnames=fieldnames)
        wtr.writeheader()
        print(fieldnames)
        for __, obs in data.items():
            print(obs)
            for i, __ in enumerate(obs['time']):
                row = {k: obs[k][i] for k in fieldnames if k[-5:] != 'group'}
                row['birthrate_group'] = obs['birthrate_group']
                wtr.writerow(row)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-p', 'paramfile', type=click.Path())
@click.option('-o', 'outfilebase', type=click.Path())
def split_by_group(infile, paramfile, outfilebase):
    """
    decompose a single longform csv based on the groups
    also decomposes the paramfile to match
    """
    data = {}
    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for line in rdr:
            idx = line['birthrate_group']
            if idx not in data:
                data[idx] = []
            data[idx].append(line)

    params = toml.load(paramfile)
    print(params)

    for key, dataset in data.items():
        filename = outfilebase + '.b' + key + '.csv'
        pfn = outfilebase + '.b' + key + '.toml'
        params['abc_params']['birthrate_coupling'] = 'all'
        params['abc_params']['birthrate_coupling_sets'] = []
        with open(filename, 'w') as out_csv:
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            wtr.writeheader()
            for line in dataset:
                wtr.writerow(line)
        with open(pfn, 'w') as out_toml:
            toml.dump(params, out_toml)





if __name__ == '__main__':
    main()
