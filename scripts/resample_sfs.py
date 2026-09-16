#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys, os
from allel import sfs_folded

ac_fs = [('A.0D.tsv', 'A.4D.tsv'), ('X.0D.tsv', 'X.4D.tsv')]
#ac_invar = [(5746147, 4214009), (1283252, 947167)] # all
#ac_invar = [(969552, 211663), (956121, 209216)] # male biased
#ac_invar = [(426701, 97523), (335707, 77777)] # female biased
ac_invar = [(3863825, 865746), (2701683, 609513)] # unbiased

for i, (non_syn_fn, syn_fn) in enumerate(ac_fs):
	outstr = syn_fn.split('.')[0]
	os.makedirs(f'results/{outstr}/', exist_ok=True)


	non_syn_ac_df = pd.read_csv(non_syn_fn, sep='\t', header=None, dtype='int')
	non_syn_biallelic_ac = non_syn_ac_df[[0,1]].to_numpy()
	syn_ac_df = pd.read_csv(syn_fn, sep='\t', header=None, dtype='int')
	syn_biallelic_ac = syn_ac_df[[0,1]].to_numpy()
	
	#full sfs
	non_syn_sfs = sfs_folded(non_syn_biallelic_ac).tolist()
	non_syn_sfs[0] = ac_invar[i][0]
	syn_sfs = sfs_folded(syn_biallelic_ac).tolist()
	syn_sfs[0] = ac_invar[i][1]

	# 11 females 1 male so -1 for the males on X sfs and -1 extra 0
	mac = len(non_syn_sfs)-1
	if i == 0:
		count = 2*(mac)
	elif i == 1:
		count = 2*(mac)-1

	non_syn_sfs = non_syn_sfs + [0]*mac
	syn_sfs = syn_sfs + [0]*mac

	with open(f"results/{outstr}/all.sfs.txt", 'w') as fh:
		fh.write("1 \n")
		fh.write(f"{str(count)} \n")
		fh.write(f"{' '.join(str(x) for x in non_syn_sfs)} \n")
		fh.write(f"{' '.join(str(x) for x in syn_sfs)} \n")
	
	#bootstrap
	# is this random choice doing really what i want? take a number from the total number the same number of times as the total number
	# I think so yeah
	for j in range(100):
		non_syn_bootstrap = sfs_folded(
			non_syn_biallelic_ac[np.random.choice(non_syn_biallelic_ac.shape[0], non_syn_biallelic_ac.shape[0])]
		).tolist()
		non_syn_bootstrap[0] = ac_invar[i][0]
		syn_bootstrap = sfs_folded(
			syn_biallelic_ac[np.random.choice(syn_biallelic_ac.shape[0], syn_biallelic_ac.shape[0])]
		).tolist()
		syn_bootstrap[0] = ac_invar[i][1]

		non_syn_bootstrap = non_syn_bootstrap + [0]*mac
		syn_bootstrap = syn_bootstrap + [0]*mac

		with open(f"results/{outstr}/{j+1}.sfs.txt", 'w') as fh:
			fh.write("1 \n")
			fh.write(f"{str(count)} \n")
			fh.write(f"{' '.join(str(x) for x in non_syn_bootstrap)} \n")
			fh.write(f"{' '.join(str(x) for x in syn_bootstrap)} \n")