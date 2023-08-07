"""getHet.py

Usage:
getHet.py -v <FILE> -b <FILE> [-h -o <FILE>]

Options:
-v, --vcf <FILE>	VCF file
-b, --bed <FILE>	bed file
-o, --output <STR>	Output file prefix

"""


# Example command
# python getHet.py data.vcf sample.lg_het.tsv

import pandas as pd
import allel
from docopt import docopt


def getCallableLength(callable_bed_f):
    bed_df = pd.read_csv(
        callable_bed_f,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "callable"],
    )
    bed_df["length"] = bed_df["end"] - bed_df["start"]
    return bed_df.groupby("chrom").sum().reset_index()[["chrom", "length"]]


def getHetCount(vcf_dict, chromosome):
    chromosome_array = vcf_dict["variants/CHROM"]

    if chromosome in chromosome_array:
        is_SNP_array = vcf_dict["variants/is_snp"]

        if isinstance(vcf_dict["variants/NUMALT"][0], int):
            numalt_array = vcf_dict["variants/NUMALT"]
            mask_array = (
                (numalt_array == 1)
                & (is_SNP_array == True)
                & (chromosome_array == chromosome)
            )
        else:
            #print("NUMALT not in VCF file, assuming 1...")
            mask_array = (is_SNP_array == True) & (chromosome_array == chromosome)

        snp_gts = vcf_dict["calldata/GT"][mask_array]
        snp_pos = vcf_dict["variants/POS"][mask_array]

        snp_ga = allel.GenotypeArray(snp_gts)
        return snp_ga.is_het().sum()

    else:
        return 0


args = docopt(__doc__)

vcf_f = args["--vcf"]
bed_f = args["--bed"]
prefix = args["--output"]

results_df = getCallableLength(bed_f)

query_fields = [
    "samples",
    "calldata/GT",
    "variants/CHROM",
    "variants/POS",
    "variants/NUMALT",
    "variants/is_snp",
]

vcf_dict = allel.read_vcf(vcf_f, fields=query_fields)

n_heterozygous_sites = []
for linkage_group in results_df["chrom"].to_list():
    n_heterozygous_sites.append(getHetCount(vcf_dict, linkage_group))
results_df["n_heterozygous_sites"] = n_heterozygous_sites
results_df["h_prop"] = results_df["n_heterozygous_sites"] / results_df["length"]

results_df.to_csv("{}.lg_het.tsv".format(prefix), sep="\t")

