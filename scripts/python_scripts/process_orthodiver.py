import pandas as pd
import numpy as np

zero_d_f = 'v1.0d_pi_by_locus.txt'
four_d_f = 'v1.4d_pi_by_locus.txt'

zero_df = pd.read_csv(zero_d_f, sep = '\t', usecols=[0, 3])
four_df = pd.read_csv(four_d_f, sep = '\t', usecols=[0, 3])
zero_df.columns = ['transcript','0_dxy']
four_df.columns = ['transcript','4_dxy']
out_df = zero_df.merge(four_df, how = 'outer',on='transcript')
out_df['0/4_dxy'] = out_df['0_dxy'] / out_df['4_dxy']
out_df.to_csv('divergence_per_transcript.tsv', sep = '\t', index=None)