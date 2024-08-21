#!/usr/bin/env python3

"""calculatePiBed.py

Usage:
calculatePiBed.py -v <FILE> -b <FILE> -a <FILE> -g <FILE> [-h -n <STR> -o <FILE> -l <STR>]

Options:
-v, --vcf <FILE>	       VCF file
-b, --bed <FILE>	       Feature bed file (chrom, start, stop, feature)
-a, --accessible <FILE>    Accessible bed 
-n, --name <STR>           Chromosome name
-o, --output <FILE>	       Output file name
-l, --label <STR>          Result label
-g, --genome <FILE>        Genome file, left column chrom name right column total length

"""

# Example Command
# python calculatePiBed.py -v allacma_fusca.vcf.gz -b  allacma_fusca.0D.longest_isoforms.vcf.gz -o allacma_fusca.0D.theta.tsv

import allel
import pybedtools
import pandas as pd
import numpy as np
import sys
from docopt import docopt
from allel.util import asarray_ndim, mask_inaccessible


def parse_vcf(vcf, chromosome):
    query_fields = [
        "calldata/GT",
        "variants/CHROM",
        "variants/POS",
        "variants/is_snp",
        "calldata/DP",
    ]

    vcf_dict = allel.read_vcf(vcf, fields=query_fields, region=chromosome)

    is_SNP_array = vcf_dict["variants/is_snp"]

    mask_array = is_SNP_array == True

    snp_gts = vcf_dict["calldata/GT"][mask_array]
    snp_pos = vcf_dict["variants/POS"][mask_array]
    snp_dp = vcf_dict["calldata/DP"][mask_array]
    snp_ga = allel.GenotypeArray(snp_gts)

    gt_counts = []
    for i in range(snp_ga.shape[1]):
        homref = snp_ga[:, i].count_hom_ref()
        het = snp_ga[:, i].count_het()
        homalt = snp_ga[:, i].count_hom_alt()
        gt_counts.append([homref, het, homalt])
    gt_counts = np.array(gt_counts)

    ac = snp_ga.count_alleles()

    return snp_pos, ac, gt_counts, snp_ga, snp_dp


def get_accessible(bed, chrom):
    df = pd.read_csv(bed, sep="\t", names=["chrom", "start", "stop"])
    degen_pos = df[df["chrom"] == chrom]["start"].to_numpy()
    n = get_length(chrom, genome_df)
    acc_arr = np.full((n), False)
    acc_arr[degen_pos] = True
    return acc_arr

def get_length(chrom, genome_df):
    length = genome_df.loc[genome_df["chrom"] == chrom, "length"].iloc[0]
    return length


if __name__ == "__main__":
    args = docopt(__doc__)

    vcf_f = args["--vcf"]
    bed_f = args["--bed"]
    acc_bed_f = args["--accessible"]
    genome_file = args["--genome"]

    #if args["--name"]:
    #    name = str(args["--name"])
    #else:
    #    name = "NA"

    if args["--output"]:
        out_f = str(args["--output"])
    else:
        out_f = "bed_pi.tsv"

    if args["--label"]:
        result_label = str(args["--label"])
    else:
        result_label = "pi"

    results = []

    name = '.'.join(bed_f.split('.')[0:-1])

    genome_df = pd.read_csv(genome_file, sep="\t", names=["chrom", "length"])
    accessible_array = get_accessible(acc_bed_f, name)
    snp_pos, ac, gt_count_arr, snp_ga, snp_dp = parse_vcf(vcf_f, chromosome=name)

    dp_arr = np.concatenate((snp_pos.reshape(-1, 1), snp_dp), axis=1)

    np.savetxt(
        f"{name}.{result_label}.snp_dp.txt",
        dp_arr,
        delimiter="\t",
    )

    # ga_arr = np.concatenate((snp_pos.reshape(-1, 1), snp_ga), axis=1)

    # np.savetxt(
    #    f"{name}.{result_label}.snp_ga.txt",
    #    ga_arr,
    #    delimiter="\t",
    # )

    np.savetxt(
        f"{name}.{result_label}.gt_counts.txt",
        gt_count_arr,
        delimiter="\t",
        header="hom_ref\thet\thom_alt",
        comments="",
    )

    idx = allel.SortedIndex(snp_pos)
    bed = pybedtools.BedTool(bed_f)

    # name = chromosome name
    # result_label = degeneracy
    # accessible array is true/false arr made of 0d or 4d site positions on one chrom
    # so this should get the 0d or 4d folded SFS for the analysed chrom

    is_acc = asarray_ndim(accessible_array, 1, allow_none=True)
    pos, ac_is_acc = mask_inaccessible(is_acc, idx, ac)
    biallelic_ac = ac_is_acc.compress(ac_is_acc.is_biallelic()[:], axis=0)[:, :2]

    # print MAC arr and pos
    biallelic_pos = pos[ac_is_acc.is_biallelic()[:]]
    biallelic_arr = np.append(
        biallelic_ac, biallelic_pos.reshape(len(biallelic_pos), 1), axis=1
    )
    np.savetxt(
        f"{name}.{result_label}.biallelic_ac.txt",
        biallelic_arr,
        delimiter="\t",
        header="ref_allele_count\talt_allele_count\tpos",
        comments="",
    )

    fifton_pos = biallelic_arr[biallelic_arr[:, 0] == biallelic_arr[:, 1]][:, 2]
    fifton_arr = np.isin(snp_pos, fifton_pos)
    fifton_ga = snp_ga[fifton_arr]
    ## is this bit using the right pos?

    fifton_count_arr = np.array(
        [
            fifton_ga.count_hom(axis=1),
            fifton_ga.count_het(axis=1),
            fifton_ga.count_hom_alt(axis=1),
            fifton_pos,
        ]
    ).transpose()

    np.savetxt(
        f"{name}.{result_label}.fifton_gc.txt",
        fifton_count_arr,
        delimiter="\t",
        header="ref_hom_count\thet_count\talt_hom_count\tpos",
        comments="",
    )

    sfs = allel.sfs_folded(biallelic_ac)

    total = np.sum(accessible_array)
    extra_invar = total - np.sum(sfs)
    sfs[0] = sfs[0] + extra_invar
    np.savetxt(f"{result_label}.{name}.sfs.txt", sfs)

    for interval in bed:
        (chrom, start, stop, feature) = interval
        try:
            pi = allel.sequence_diversity(
                idx,
                ac,
                start=int(start) + 1,
                stop=int(stop), #is this right for indexing, i think it is?
                is_accessible=accessible_array,
            )
            results.append((chrom, feature, pi))
        except KeyError:
            results.append((chrom, feature, "NaN"))

    pd.DataFrame(results, columns=["chrom", "feature", result_label]).to_csv(
        out_f, sep="\t", index=False
    )
