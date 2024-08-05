#!/usr/bin/env python3

"""agg_sfs.py

Usage:
agg_sfs.py -i <DIR> -x <STR> -o <STR> [-h]

Options:
-i, --inputs <DIR>	       Directory of linkage group SFSs
-x, --x_chroms <STR>	   Comma separated list of X-linked linkage groups
-o, --output <STR>	       Output file prefix
"""

# Example Command
# python agg_sfs.py -i inputs -x OX359249,OX359250 -o allacma_fusca

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import sys
from docopt import docopt


def create_df(dir_name, x_chroms):
    input_dir = "./" + str(dir_name) + "/"
    files = os.listdir(input_dir)
    files = [file for file in files if "sfs" in file]
    files = [input_dir + file for file in files]

    df_list = []
    for file in files:
        chromosome = file.split(".")[2].split("_")[0]
        degen = file.split(".")[1].split("/")[2].split("_")[0]
        df = pd.read_csv(file, sep="\t", header=0, names=["count"])
        df["chromosome"] = chromosome
        df["degeneracy"] = degen
        df["iton"] = [i + 1 for i in range(df.shape[0])]
        if chromosome in str(x_chroms).split(","):
            df["linkage"] = "X"
        else:
            df["linkage"] = "A"
        df_list.append(df)
    df = pd.concat(df_list)

    return df


def norm_sfs(sfs):
    # sfs = sfs.to_numpy(dtype=float)[1:]
    sfs /= np.sum(sfs)
    return sfs


def get_exp_sfs(nsamp):
    sfs_exp_un = [1 / i for i in range(1, 2 * nsamp)]
    sfs_exp_un /= np.sum(sfs_exp_un)
    flipped = np.flip(sfs_exp_un)
    sfs_exp = sfs_exp_un + flipped  # sum them together for the 'folding'
    sfs_exp = sfs_exp[0:nsamp]  # only look at the minor allele counts
    # sfs_exp[-1] = sfs_exp[-1]/2 # because we've counted the entry for this twice
    return sfs_exp


def plot_sfs(a_sfs_df, x_sfs_df):
    plt.rcParams.update({"font.size": 15})
    sns.set_palette("muted")

    fig, axs = plt.subplots(nrows=2, figsize=[10, 10], sharex=True)

    a_ax = axs[0]
    x_ax = axs[1]

    sns.barplot(data=a_sfs_df, x="itons", y="sfs", hue="degeneracy", ax=a_ax)

    sns.barplot(data=x_sfs_df, x="itons", y="sfs", hue="degeneracy", ax=x_ax)

    a_ax.set_title("Autosome", loc="right")
    a_ax.set_xlabel("")
    a_ax.set_ylabel("")
    # a_ax.set_ylim([0,0.35])

    x_ax.set_title("X", loc="right")
    x_ax.set_ylabel("")
    x_ax.get_legend().remove()
    # x_ax.set_ylim([0,0.35])

    fig.text(0.06, 0.5, "Normalised SFS", ha="center", va="center", rotation="vertical")
    plt.savefig(f"{str(args['--output'])}.pdf", bbox_inches="tight")


if __name__ == "__main__":
    args = docopt(__doc__)

    df = create_df(args["--inputs"], args["--x_chroms"])

    sfs_list = []

    for degen, linkage in [("0D", "A"), ("4D", "A"), ("0D", "X"), ("4D", "X")]:
        sfs = (
            df[(df["degeneracy"] == degen) & (df["linkage"] == linkage)]
            .groupby("iton")["count"]
            .sum()
            .to_numpy(dtype=float)
        )
        sfs_list.append(np.array(sfs))

    zero_A_sfs_norm, four_A_sfs_norm, zero_X_sfs_norm, four_X_sfs_norm = (
        norm_sfs(sfs) for sfs in sfs_list
    )

    nsamp = len(zero_A_sfs_norm)

    sfs_exp = get_exp_sfs(nsamp)

    degeneracy = ["exp"] * nsamp + ["0D"] * nsamp + ["4D"] * nsamp
    itons = np.tile(np.arange(1, nsamp + 1), 3)

    a_sfs_df = pd.DataFrame(
        {
            "degeneracy": degeneracy,
            "itons": itons,
            "sfs": np.concatenate((sfs_exp, zero_A_sfs_norm, four_A_sfs_norm)),
        }
    )

    x_sfs_df = pd.DataFrame(
        {
            "degeneracy": degeneracy,
            "itons": itons,
            "sfs": np.concatenate((sfs_exp, zero_X_sfs_norm, four_X_sfs_norm)),
        }
    )

    plot_sfs(a_sfs_df, x_sfs_df)
    df.to_csv(f"{str(args['--output'])}.sfs.chrom.tsv", sep="\t", index=False)
    df.groupby(['linkage','iton']).sum().reset_index().to_csv(
        f"{str(args['--output'])}.sfs.total.tsv", sep="\t", index=False
    )
