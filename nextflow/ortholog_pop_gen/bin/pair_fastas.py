#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage: pair_fastas.py -i <DIR> -o <DIR> -a <STR> -b <STR> [-h]

  [Options]
    -i, --input_dir <DIR>                     Directory of fasta files
    -o, --output_dir <DIR>                    Output directory
    -a, --species_a <STR>                     Species A name
    -b, --species_b <STR>					  Species B name
    -h, --help                                Show this message

"""

# Example command:
# python pair_fastas.py -i hap_fastas_rn -o hap_fasta_pairs -a species_A -b species_B

import docopt
import os
from itertools import combinations


if __name__ == "__main__":
    args = docopt(__doc__)

    input_dir = os.listdir(args['--input_dir'])
    output_dir_name = args['--output_dir']
    species_A = str(args['--species_A'])
    species_B = str(args['--species_B'])

    A_files = []
    B_files = []

	for file_name in input_dir:   
	    if species_A in str(file_name):
	    	A_files.append(file_name)
	    if species_B in str(file_name):
	    	B_files.append(file_name)

	print(A_files)
	print(B_files)