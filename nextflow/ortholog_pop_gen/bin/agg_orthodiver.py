#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage: 
 agg_orthodiver.py -i <dir> [-h]

[Options]
 -i, --input_dir <dir>                     Directory of orthodiver output files.
 -h, --help                                Show this message.

"""

# Example command:
# python agg_orthodiver.py -i orthodiver_results

import os
import pandas as pd
import numpy as np
from docopt import docopt
from tqdm import tqdm

class OrthogroupObj(object):
    def __init__(
        self, name
    ):
        self.name = name
        self.Zd_pi_A = []
        self.Zd_pi_B = []
        self.Zd_dxy = []
        self.Fd_pi_A = []
        self.Fd_pi_B = []
        self.Fd_dxy = []
        self.metrics = {}

    def append_att(self, att, val):
        self.__dict__[att].append(val)

    def calculate_averages(self):
        self.metrics["Zd_pi_A"] = np.mean(self.Zd_pi_A)
        self.metrics["Zd_pi_B"] = np.mean(self.Zd_pi_B)
        self.metrics["Zd_dxy"] = np.mean(self.Zd_dxy)
        self.metrics["Fd_pi_A"] = np.mean(self.Fd_pi_A)
        self.metrics["Fd_pi_B"] = np.mean(self.Fd_pi_B)
        self.metrics["Fd_dxy"] = np.mean(self.Fd_dxy)

    def calculate_metrics(self):
        self.metrics["Zd_pi_mean"] = (
            self.metrics["Zd_pi_A"] + self.metrics["Zd_pi_B"]
        ) / 2
        self.metrics["Fd_pi_mean"] = (
            self.metrics["Fd_pi_A"] + self.metrics["Fd_pi_B"]
        ) / 2
        self.metrics["Zd_da"] = self.metrics["Zd_dxy"] - self.metrics["Zd_pi_mean"]
        self.metrics["Fd_da"] = self.metrics["Fd_dxy"] - self.metrics["Fd_pi_mean"]

        if self.metrics["Fd_dxy"] == 0:
        	self.metrics["dnds"] = np.nan
        else:
        	self.metrics["dnds"] = self.metrics["Zd_dxy"] / self.metrics["Fd_dxy"]

        if (self.metrics["Zd_pi_mean"] == 0) or (self.metrics["Fd_pi_mean"] == 0):
            self.metrics["alpha"] = 1
        else:
            self.metrics["alpha"] = 1 - (
                (self.metrics["Fd_dxy"] * self.metrics["Zd_pi_mean"])
                / (self.metrics["Zd_dxy"] * self.metrics["Fd_pi_mean"])
            )


def parse_files(file_list, att_list):
    for f in file_list:
        with open(f) as fh:
            next(fh)
            x = 0
            for line in fh:
                if "NA" in line:
                    continue
                orthogroup, pi_A, pi_B, dxy, pi_tot = line.rstrip("\n").split("\t")
                orthogroup_obj_dict[orthogroup].append_att(att_list[0], float(pi_A))
                orthogroup_obj_dict[orthogroup].append_att(att_list[1], float(pi_B))
                orthogroup_obj_dict[orthogroup].append_att(att_list[2], float(dxy))


if __name__ == "__main__":
    args = docopt(__doc__)
    input_dir = str(args["--input_dir"])

    Zd_files = []
    Fd_files = []

    for file_name in os.listdir(input_dir):
        if "0d_pi" in str(file_name):
            Zd_files.append(input_dir + "/" + file_name)
        if "4d_pi" in str(file_name):
            Fd_files.append(input_dir + "/" + file_name)

    init_df = pd.read_csv(Fd_files[0], sep="\t")
    orthogroups = set(init_df["# locus_id"].to_list())
    orthogroup_obj_dict = {name: OrthogroupObj(name=name) for name in orthogroups}
    sp_A, sp_B = (
        init_df.columns[1].split("(")[1].split(".")[0],
        init_df.columns[2].split("(")[1].split(".")[0],
    )

    metric_list = [
        "Zd_pi_A",
        "Zd_pi_B",
        "Zd_pi_mean",
        "Fd_pi_A",
        "Fd_pi_B",
        "Fd_pi_mean",
        "Zd_dxy",
        "Fd_dxy",
        "Zd_da",
        "Fd_da",
        "dnds",
        "alpha",
    ]

    att_list_Zd = ["Zd_pi_A", "Zd_pi_B", "Zd_dxy"]
    att_list_Fd = ["Fd_pi_A", "Fd_pi_B", "Fd_dxy"]

    parse_files(Zd_files, att_list_Zd)
    parse_files(Fd_files, att_list_Fd)

    print("Calculating metrics...")
    for orthogroup in tqdm(orthogroups, total=len(orthogroups), desc="[%] ", ncols=80):
        orthogroup_obj_dict[orthogroup].calculate_averages()
        orthogroup_obj_dict[orthogroup].calculate_metrics()

    fn = f"{sp_A}.{sp_B}.orthodiver_agg.tsv"

    with open(
        fn,
        "w",
    ) as fh:
        fh.write(
            f"Orthogroup\t0d_pi_{sp_A}\t0d_pi_{sp_B}\t0d_pi_mean\t4d_pi_{sp_A}\t4d_pi_{sp_B}\t4d_pi_mean\t0d_dxy\t4d_dxy\t0d_da\t4d_da\tdnds\talpha\n"
        )
        for orthogroup in orthogroups:
            line = [
                str(orthogroup_obj_dict[orthogroup].metrics[metric])
                for metric in metric_list
            ]
            line.insert(0, orthogroup)
            fh.write("\t".join(line) + "\n")