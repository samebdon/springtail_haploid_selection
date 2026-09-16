import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

file = 'afusca_dicmin.paf'

columns = [
'query_name',
'query_length',
'query_start',
'query_end',
'strand',
'target_name',
'target_length',
'target_start',
'target_end',
'matches',
'length',
'mapq',
'dvF',
'dfI'
]
df = pd.read_csv(file,sep='\t')
df.columns = columns
df['F'] = pd.to_numeric(df['dvF'].str.split(':', expand=True)[2])
df['I'] = pd.to_numeric(df['dfI'].str.split(':', expand=True)[2])
df['similarity'] = df['matches']/df['length']
df['identity'] = 1-(df['I']/df['matches'])
df['synteny_colour'] = df['query_name'].str.split('.',expand=True)[0].str.strip().str[-1]

filter_df = df.loc[(df['similarity'] >=0.8) & (df['identity']>=0.8)]

filter_df[[
'synteny_colour',
'query_name',
'query_start',
'query_end',
'target_name',
'target_start',
'target_end']].to_csv('afusca_dicmin.liftover.tsv',
	sep='\t',index=False, header=False)


#plt.hist(df['F'], bins=20, weights = df['length'])
#plt.savefig('afusca_dicmin.F.png')

#plt.hist(df['I'])
#plt.savefig('afusca_dicmin.I.png')