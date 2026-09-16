#!/usr/bin/env python3

import pandas as pd 
import numpy as np
import sys, os
from matplotlib import pyplot as plt
import seaborn as sns

infile = sys.argv[1]

#colnames = [
#	'scaf1_name',
#	'scaf1_len',
#	'AB',
#	'AE',
#	'polarity',
#	'scaf2_name',
#	'scaf2_len',
#	'BB',
#	'BE',
#	'iid',
#	'blocksum/2',
#	'255',
#	'divergence',
#	'didd',
#]

colnames = [
	'scaf1_name',
	'scaf1_len',
	'AB',
	'AE',
	'polarity',
	'scaf2_name',
	'scaf2_len',
	'BB',
	'BE',
	'iid',
	'blocksum/2',
	'255',
	'divergence',
	'didd',
]

df = pd.read_csv(infile, names=colnames, sep='\t')
scaf1_singles = df['scaf1_name'].value_counts()[df['scaf1_name'].value_counts()==1].reset_index()['index']
scaf2_singles = df['scaf2_name'].value_counts()[df['scaf2_name'].value_counts()==1].reset_index()['index']

df['dxy'] = df['divergence'].str.split(':',expand=True)[2].astype(float)
df['diddnum'] = df['didd'].str.split(':',expand=True)[2].astype(float)

 
df1 = df[df.scaf1_name.isin(scaf1_singles)]
df2 = df[df.scaf2_name.isin(scaf2_singles)]
df3 = df[(df.scaf1_name.isin(scaf1_singles)) & (df.scaf2_name.isin(scaf2_singles))]

print(f"S1 singles aln length: {df1['blocksum/2'].sum()}, S1 singles div bases: {np.sum(df1['blocksum/2']*df1['dxy'])}")
print(f"S2 singles aln length: {df2['blocksum/2'].sum()}, S2 singles div bases: {np.sum(df2['blocksum/2']*df2['dxy'])}")
#print(f"S3 singles aln length: {df3['blocksum/2'].sum()}, S3 singles div bases: {np.sum(df3['blocksum/2']*df3['dxy'])}")

print(df1['dxy'].mean())


fig, ax = plt.subplots(figsize=[8,8])
sns.jointplot(data = df1,
			x = 'blocksum/2',
			y = 'dxy',
			ax=ax)
plt.tight_layout()
plt.savefig('df1.png', dpi=300)
plt.close()

fig, ax = plt.subplots(figsize=[8,8])
sns.jointplot(data = df2,
			x = 'blocksum/2',
			y = 'dxy',
			ax=ax)
plt.tight_layout()
plt.savefig('df2.png', dpi=300)
plt.close()
#create a bed file from the dfs for intersecting

df[['scaf1_name','AB', 'AE', 'dxy']].to_csv('afusca_svir.bed',sep='\t', index=False, header=False)