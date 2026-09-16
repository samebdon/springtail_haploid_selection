import pandas as pd
import numpy as np
import functools
from string import ascii_lowercase
from operator import itemgetter

def load_busco_results(buscofile, genomefile):
    busco_colnames = ["busco_id", "status", "seq_code", "start", "stop"]
    chromosomes = pd.read_csv(genomefile, sep='\t', usecols=[0], header=None)[0].to_list()
    df = pd.read_csv(
        buscofile,
        sep="\t",
        usecols=[0, 1, 2, 3, 4],
        header=None,
        names=busco_colnames,
        skiprows=[0, 1, 2],
    )

    df = df[df["status"] != "Missing"]
    if "|" in df["seq_code"].iloc[0]:
        df["seq"] = (
            df["seq_code"]
            .str.split("|")
            .apply(itemgetter(1))
            .str.split(":")
            .apply(itemgetter(0))
        )
    else:
        df["seq"] = df["seq_code"]
    df.drop(labels=["seq_code"], axis=1, inplace=True)
    df = df[df["seq"].isin(chromosomes)]
    return df[["busco_id", "status", "seq", "start", "stop"]]


def get_labels(seqs):
    labels = {}
    for i, seq in enumerate(seqs):
        labels[seq] = int(i + 1)
    labels[np.nan] = np.nan
    return labels

def label_colours_by_ref(df):
    seqs = []
    for genome in df.filter(like='seq').columns:
        seqs = seqs + sorted(df[genome].dropna().unique())
    labels = get_labels(seqs)
    df['colour'] = np.nan
    for i, genome in enumerate(df.filter(like='seq').columns):
        if (i+1) == len(df.filter(like='seq').columns):
            break
        for seq in df[genome]:
            df.loc[(df[genome] == seq) & (df["colour"] != df["colour"]), 'colour'] = labels[seq]

def label_colours_by_busco(df):
    df["combs"] = df.filter(regex="status").astype(str).agg("_".join, axis=1)
    combs = np.sort(df["combs"].unique())
    labels = [n for n in range(len(combs))]
    global colour_dict
    colour_dict = {}
    for comb, label in zip(combs, labels):
        colour_dict[label] = comb
        df.loc[df["combs"] == comb, "colour"] = label
    df.drop(labels=["combs"], axis=1, inplace=True)


def create_liftover_from_busco(buscofile_list, genomefile_list):

    results_dfs = []
    for i, file in enumerate(buscofile_list):
        results_dfs.append(load_busco_results(file, genomefile_list[i]))

    it = iter(ascii_lowercase)

    for i, df in enumerate(results_dfs, start=1):
        df_label = next(it)
        df.rename(
            columns={
                col: "{}_{}".format(col, df_label)
                for col in ("seq", "start", "stop", "status")
            },
            inplace=True,
        )
    merge = functools.partial(pd.merge, how="outer", on="busco_id")
    df = functools.reduce(merge, results_dfs)

    if args["--busco_colours"]:
        label_colours_by_busco(df)
    else:
        label_colours_by_ref(df)

    #label_colours_by_ref(df)

    cols = ["colour"] + [
        label
        for label in df.columns
        if any(x in label for x in ["seq", "start", "stop"])
    ]

    df[cols].to_csv("liftover.tsv", sep="\t", index=False, header=False, na_rep='NA')

if __name__ == "__main__":

    file_1 = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/synteny/busco_results/fol_ang/run_arthropoda_odb10/full_table.tsv'
    file_2 = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/synteny/busco_results/smi_aqu/run_arthropoda_odb10/full_table.tsv'
    file_3 = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/synteny/busco_results/all_fus/run_arthropoda_odb10/full_table.tsv'

    gen1 = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/synteny/fol_ang.genomefile.tsv'
    gen2 = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/synteny/smi_aqu.genomefile.tsv'
    gen3 = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/synteny/all_fus.genomefile.tsv'


    buscofile_list = [file_1, file_2, file_3]
    genomefile_list = [gen1, gen2, gen3]
    create_liftover_from_busco(buscofile_list, genomefile_list)