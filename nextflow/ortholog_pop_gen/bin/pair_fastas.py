#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage: pair_fastas.py -i <DIR> -o <DIR> [-h]

  [Options]
    -i, --input_dir <DIR>                     Directory of fasta files
    -o, --output_dir <DIR>                    Output directory
    -h, --help                                Show this message

"""

# Example command:
# python pair_fastas.py -i hap_fastas_rn -o hap_fasta_pairs

import docopt
import os
from itertools import combinations


path = args