


rule define_groups:
    input:
        obs = "data/{filename}.csv",
        par = "data/{filename}.toml"
    output:
        temp("intermediate/{filename}.groups.csv")
    shell:
        """
        python3 code/csvtools.py define-groups \
            -i {input.obs} \
            -p {input.par} \
            -o {output}
        """


rule zero_time:
    input:
        "intermediate/{filename}.groups.csv"
    output:
        temp("intermediate/{filename}.groups.zero.csv")
    shell:
        """
        python3 code/csvtools.py zero-time-longform \
            -i {input} \
            -o {output}
        """


rule split_by_group:
    input:
        obs = "intermediate/{filename}.groups.zero.csv",
        par = "data/{filename}.toml"
    output:
        dynamic("intermediate/{filename}.g{k}.data.csv"),
        dynamic("intermediate/{filename}.g{k}.toml"),
        # groupinfo = "intermediate/{filename}.groupinfo.toml",
    shell:
        """
        python3 code/csvtools.py split-by-group \
            -i {input.obs} \
            -p {input.par}
        """


rule reconstruct:
    output:
        "intermediate/{filename}.g{k}.db"
    input:
        obs = "intermediate/{filename}.g{k}.data.csv",
        par = "intermediate/{filename}.g{k}.toml"
    run:
        shell(" \
             python3 code/abc.py reconstruct \
                 -p {input.par} \
                 -o {input.obs} \
                 -d {output} \
        ")


rule abc_plots:
    input:
        db = "intermediate/{filename}.g{k}.db",
        obs = "intermediate/{filename}.g{k}.data.csv",
        par = "intermediate/{filename}.g{k}.toml"
    output:
        "intermediate/{filename}.g{k}.abc.pdf"
    shell:
        """
        python3 code/plots.py abc-info \
            -p {input.par} \
            -o {input.obs} \
            -d {input.db} \
            --save {output}
        """


rule fit_plots:
    input:
        db = "intermediate/{filename}.g{k}.db",
        obs = "intermediate/{filename}.g{k}.data.csv",
        par = "intermediate/{filename}.g{k}.toml"
    output:
        "intermediate/{filename}.g{k}.fit.pdf"
    shell:
        """
        python3 code/plots.py result-single \
            -p {input.par} \
            -o {input.obs} \
            -d {input.db} \
            --save {output}
        """


rule fit_tables:
    input:
        db = "intermediate/{filename}.g{k}.db",
        obs = "intermediate/{filename}.g{k}.data.csv",
        par = "intermediate/{filename}.g{k}.toml"
    output:
        "intermediate/{filename}.g{k}.fit.csv"
    shell:
        """
        python3 code/plots.py tabulate-single \
            -p {input.par} \
            -o {input.obs} \
            -d {input.db} \
            -c {output}
        """


def report_inputs(wildcards):
    import toml
    params = toml.load('data/' + wildcards.filename + '.toml')
    groups = list(range(len(params['abc_params']['birthrate_coupling_sets'])))
    sources = {}
    for k in groups:
        sources['g' + str(k) + '.abc'] = 'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.abc.pdf'
        sources['g' + str(k) + '.fit'] = 'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.fit.pdf'
        sources['g' + str(k) + '.db'] = 'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.db'
    return sources



def table_inputs(wildcards):
    import toml
    params = toml.load('data/' + wildcards.filename + '.toml')
    groups = list(range(len(params['abc_params']['birthrate_coupling_sets'])))
    sources = {}
    for k in groups:
        # sources['g' + str(k) + '.db'] = 'intermediate/' \
        #     + wildcards.filename + '.g' + str(k) + '.db'
        sources['g' + str(k) + '.csv'] = 'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.fit.csv'
    return sources



def num_to_aa(n):
    """
    convert an integer to a base 25 letter number
    NOTE: this function isn't actually correct
    but it's good enough for this use-case
    """
    if n == 0:
        return 'A'
    aa = ''
    while n > 0:
        aa += chr(65 + n%25)
        n //= 25
    return aa


rule produce_table:
    input:
        unpack(table_inputs)
    output:
        "results/{filename}.fit.csv"
    run:
        print(input)
        tables = ' '.join(['-i ' + x for x in input])
        shell('python3 code/csvtools.py merge -o ' + str(output) + ' ' + tables)


rule produce_report:
    input:
        unpack(report_inputs)
    output:
        "results/{filename}.pdf"
    run:
        # identify most likely model to put into report
        import toml
        import re
        import numpy as np
        from pyabc import History
        params = toml.load('data/' + wildcards.filename + '.toml')
        groups = list(range(len(params['abc_params']['birthrate_coupling_sets'])))
        sources = {x: {} for x in groups}
        for group in groups:
            sources[group]['names'] = params['abc_params']['birthrate_coupling_sets'][group]
            db_path = 'sqlite:///' + input['g' + str(group) + '.db']
            abc_history = History(db_path)
            axs = abc_history.get_model_probabilities()
            final = np.array(axs[-1:])
            final = final[final > 0]
            sources[group]['model_fraction'] = np.max(final)
            sources[group]['best_model'] = np.where(final == sources[group]['model_fraction'])[0][0]
            sources[group]['num_models'] = len(final)

        # collect sources into output file
        files = ' '.join([num_to_aa(x) + '=' + input['g' + str(x) + '.fit']
                          for x in groups])
        pages = ' '.join([' '.join([num_to_aa(y) \
                                    + str(x*sources[y]['num_models'] + sources[y]['best_model'] + 1)
                          for x in range(2)]) for y in groups])
        command = 'pdftk ' + files + ' cat ' + pages + ' output ' + str(output)
        shell(command)


