#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage: pair_fastas.py -i <DIR> -o <DIR> -a <STR> -b <STR> [-h]

  [Options]
    -i, --input_dir <DIR>                     Directory of fasta files
    -o, --output_dir <DIR>                    Output directory
    -a, --species_A <STR>                     Species A name
    -b, --species_B <STR>                     Species B name
    -h, --help                                Show this message

"""

# Example command:
# python pair_fastas.py -i hap_fastas_rn -o hap_fasta_pairs -a species_A -b species_B

import re
import os
from docopt import docopt
from itertools import product


if __name__ == "__main__":
    args = docopt(__doc__)

    input_dir = str(args["--input_dir"])
    output_dir = str(args["--output_dir"])
    species_A = str(args["--species_A"])
    species_B = str(args["--species_B"])

    A_files = []
    B_files = []

    for file_name in os.listdir(input_dir):
        if species_A in str(file_name):
            A_files.append(file_name)
        if species_B in str(file_name):
            B_files.append(file_name)

    file_pairs = list(product(A_files, B_files))

    for file_A, file_B in file_pairs:
        orthogroup = file_A.split(".")[0]

        sample_A = file_A.split(".")[1]
        sample_B = file_B.split(".")[1]

        with open(os.path.join(input_dir, file_A)) as file:
            lines_A = [line.rstrip() for line in file]

        with open(os.path.join(input_dir, file_B)) as file:
            lines_B = [line.rstrip() for line in file]

        for i in [1, 3]:
            lines_A[i] = re.sub('[a-z]', 'N', lines_A[i])
            lines_B[i] = re.sub('[a-z]', 'N', lines_A[i])

        if lines_A[1].count("A")+lines_A[1].count("T")+lines_A[1].count("C")+lines_A[1].count("G") < 3:
            continue

        if lines_B[1].count("A")+lines_B[1].count("T")+lines_B[1].count("C")+lines_B[1].count("G") < 3:
            continue

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        outfile = os.path.join(output_dir, f"{orthogroup}.{sample_A}.{species_A}.{sample_B}.{species_B}.unaln.fa")

        try:
            os.remove(outfile)
        except OSError:
            pass

        with open(outfile,"a",
        ) as f:
            f.writelines(line + '\n' for line in
                [
                    lines_A[0],
                    lines_A[1],
                    lines_B[0],
                    lines_B[1],
                    lines_A[2],
                    lines_A[3],
                    lines_B[2],
                    lines_B[3],
                ]
            )
