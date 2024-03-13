"""get_orthologr_loc.py

Usage: 
 get_orthologr_loc.py -a <file> -o <file> [-h]

Options:
 -a, --annotation <file>                  Annotation file of query species (in GTF format)
 -o, --orthologr <file>				      TSV of orthologr output
"""

# Example Command
# python get_orthologr_loc.py -i braker.gtf -i orthologr.tsv

import pandas as pd
import numpy as np
from docopt import docopt

if __name__ == "__main__":
    args = docopt(__doc__)

    columns = [
        "sequence",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    bias_df = pd.read_csv('data/results/deseq2/allacma_fusca/v_2/allacma_fusca.gene_bias.chromosomes.tsv', sep = '\t')
    anno_df = pd.read_csv(args['--annotation'], sep = '\t', names=columns, index_col=None)
    ortho_df = pd.read_csv(args['--orthologr'], sep = '\t')
    ortho_df = ortho_df[(ortho_df['evalue'] < 1e-9)]
    ortho_df[['Geneid', 'query_Gene_transcript']]= ortho_df['query_id'].str.split('.', n=1, expand=True)
    tr_df = anno_df[anno_df['feature']=='transcript'][['sequence','attribute']]
    ortho_df['sequence'] = [tr_df[tr_df['attribute']==gene]['sequence'].to_list()[0] for gene in np.array(ortho_df['query_id'])]
    ortho_df['linkage'] = np.where((ortho_df['sequence']=='OX359249.1') | (ortho_df['sequence']=='OX359250.1'), 'X','autosome')
    ortho_bias_df = ortho_df.merge(bias_df[['Geneid', 'bias']], how = 'outer',on='Geneid')#.dropna()

    x = ortho_df[ortho_df['linkage']=='X']
    a = ortho_df[ortho_df['linkage']=='autosome']

    0.116933/0.100545

    #print(x['dN'].describe())
    #print(a['dN'].describe())
    #print(x['dS'].describe())
    #print(a['dS'].describe())
    print(x['dNdS'].describe())
    print(a['dNdS'].describe())

    print(ortho_bias_df[(ortho_bias_df['bias']=='male_biased') & (ortho_bias_df['linkage']=='X')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias']!='male_biased') & (ortho_bias_df['linkage']=='X')]['dNdS'].describe())
    print(ortho_bias_df[(ortho_bias_df['bias']=='male_biased') & (ortho_bias_df['linkage']=='autosome')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias']!='male_biased') & (ortho_bias_df['linkage']=='autosome')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias']=='male_biased')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias']!='male_biased')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias']=='female_biased')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias']=='unbiased')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias'] !='female_biased') & (ortho_bias_df['linkage']=='X')]['dNdS'].describe())
    #print(ortho_bias_df[(ortho_bias_df['bias'] !='female_biased') & (ortho_bias_df['linkage']=='autosome')]['dNdS'].describe())

    #Distribution of X and A ds
    #Distribution of X and A dn
    #Distribution of X and A dnds
    #Distribution of dnds for each bias class
    #Which distribution to show for faster male evolution on the X?
    #Dnds for each bias class on the X? or Dnds for male biased genes only A vs X?

    #plot dn/ds values for x linked male biased genes, other x linked genes, and autosomal genes
    #need to do a plot for faster X effect (do we expect a general faster X effect in springtails? maybe not because of no differences)
    #do a plot for faster evolution of male biased genes (we expect faster male evolution because of sexual selection, shouldnt change with weird genomics)
    # plot of faster x male only effect, do we expect a difference here too?
    # x is overrepresented by male biased genes