#!/usr/bin/env python3

import allel
import numpy as np
import pandas as pd

# Read VCF
callset = allel.read_vcf('data/results/var_call/allacma_fusca/allacma_fusca.hard_filtered.sorted.vcf.gz')
genotypes = allel.GenotypeArray(callset['calldata/GT'])
positions = callset['variants/POS']
chroms = callset['variants/CHROM']
sample_names = callset['samples']

win_size = 10000  # 10 kb windows
all_results = []

for chrom in np.unique(chroms):
    print(f"Processing {chrom}...")
    
    # Subset to current chromosome
    mask_chr = chroms == chrom
    pos_chr = positions[mask_chr]
    geno_chr = genotypes.compress(mask_chr, axis=0)
    
    # Define windows for this chromosome
    chrom_max = pos_chr.max()
    windows = np.arange(0, chrom_max, win_size)
    
    for start in windows:
        end = start + win_size
        mask_win = (pos_chr >= start) & (pos_chr < end)
        g = geno_chr.compress(mask_win, axis=0)
        
        if g.shape[0] > 0:
            het_counts = g.count_het(axis=0)
            called_counts = g.count_called(axis=0)
            het_per_sample = np.where(
                called_counts > 0,
                het_counts / called_counts,
                np.nan
            )
        else:
            het_per_sample = np.full(len(sample_names), np.nan)
        
        all_results.append({
            "chrom": chrom,
            "window_start": start,
            "window_end": end,
            **dict(zip(sample_names, het_per_sample))
        })

# Combine into one DataFrame
df = pd.DataFrame(all_results)
df.to_csv("windowed_per_sample_diversity.tsv", sep='\t', index=False)