import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join

input_path = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/toga/codon_fastas/'
output_path = '/Users/se13/workspace/projects/springtail_haploid_selection/data/results/toga/orthodiver_input/'
files = [f for f in listdir(input_path) if isfile(join(input_path, f))]

def get_lines(fp, line_numbers):
    return (x for i, x in enumerate(fp) if i in line_numbers)

def create_header(species, n, transcript):
    return '>' + species + '.1.' + str(n) + '.' + transcript + "\n"

def process_sequence(sequence):
    return sequence.lower().replace(' ','').replace('xxx','nnn')    

for file in files:
    transcript = file.split('.codon.')[0].replace('.','_')
    with open(input_path+file, 'r') as fasta:
        lines = list(get_lines(fasta, [0,1,2,3]))
        sequence_one, sequence_two = lines[1], lines[3]

    header_one_a_out = create_header('allacma_fusca', 1, transcript)
    header_one_b_out = create_header('allacma_fusca', 2, transcript)
    header_two_a_out = create_header('sminthurus_viridis', 1, transcript)
    header_two_b_out = create_header('sminthurus_viridis', 2, transcript)
    sequence_one_out = process_sequence(sequence_one)
    sequence_two_out = process_sequence(sequence_two)

    output_file = output_path+'v1.' + transcript + '.fasta'
    with open(output_file, 'w') as fh:
        fh.writelines([header_one_a_out,
            sequence_one_out,
            header_one_b_out,
            sequence_one_out,
            header_two_a_out,
            sequence_two_out,
            header_two_b_out,
            sequence_two_out,
            ])
