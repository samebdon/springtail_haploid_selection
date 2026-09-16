"""ASR.py

Usage: 
 ASR.py -i <file> [-h -g <file>]

Options:
 -i, --input <file>                       Annotation file (in GTF format)
 -o, --orthogroups <file>				  TSV of orthogroups and the name of all included genes
 -g, --gene_trees <file>				  Gene trees for each orthogroup (Required if N(--input)>2)
"""

# Example Command
# python ASR.py -i braker1.gtf -i braker2.gtf -i braker3.gtf -o orthogroups.tsv -g gene_trees.txt -s chromosome1 chromosome2

import pandas as pd
import numpy as np

if __name__ == "__main__":
    args = docopt(__doc__)
    anno_df = args['--input']
    print(anno_df)

"""
Inputs:
N annotation files
Orthogroups
Gene trees if N > 3
Label sex chromosomes and autosomes (maybe provide list of sex chromosomes?)

Outputs:
TSV
orthogroup 	node 	state
OG1			AF		X 	
OG1			SV 		A 	
OG1 		SA 		X 	
OG1			AF,SV 	X
OG2 		...

Steps:
Check number of annotation files
If 2 proceed without gene tree
If > 2 and no gene tree request gene tree
If 2 and gene tree proceed

Import annotations as data frames (only necessary columns)
Orthogroup object
Node object?

For each line in orthogroup file
Create orthogroup object
if gene tree store gene tree
Store gene names with empty locations as a dictionary (probably need to store species incase of duplicate gene names)
For each gene name in orthogroup save location to dictionary
if theres 2 genes can just write these locations to a file
if more than 2 genes, need to do ASR

Can do by parsimony with no model of rate from autosome to x and x to autosome
Would need to use a statistical framework if assume varying rates
Two state markov chain?

need to write a recursive algorithm for ASR, can look at advanced python for biologists recursion
Given the gene tree, create the nodes where state can be reconstructed
need to figure out the algorithm to do it, need to work back progressively
after leaf nodes need to address nodes in order of coalescence i think
Find an outgroup springtail
Do as far back as possible for each orthogroup, dont just need to do single copy orthologs can include all subsets of the dataset
so the number of genes per orthogroup is arbitrary
"""