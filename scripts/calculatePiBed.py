"""calculatePiBed.py

Usage:
calculatePiBed.py -v <FILE> -b <FILE> -a <FILE> [-h -n <STR> -o <FILE> -l <STR>]

Options:
-v, --vcf <FILE>	       VCF file
-b, --bed <FILE>	       Feature bed file (chrom, start, stop, feature)
-a, --accessible <FILE>    Accessible bed 
-n, --name <STR>           Chromosome name
-o, --output <FILE>	       Output file name
-l, --label <STR>          Result label

"""

# Example Command
# python calculatePiBed.py -v allacma_fusca.vcf.gz -b  allacma_fusca.0D.longest_isoforms.vcf.gz -o allacma_fusca.0D.theta.tsv

import allel
import pybedtools
import pandas as pd
import numpy as np
from docopt import docopt

def parse_vcf(vcf, chromosome):
    query_fields = [
        "samples",
        "calldata/GT",
        "variants/CHROM",
        "variants/POS",
        "variants/is_snp",
    ]

    vcf_dict = allel.read_vcf(
        vcf,
        fields=query_fields,
        region=chromosome
    )

    read_groups = vcf_dict["samples"]

    is_SNP_array = vcf_dict["variants/is_snp"]

    mask_array = is_SNP_array == True

    snp_gts = vcf_dict["calldata/GT"][mask_array]
    snp_pos = vcf_dict["variants/POS"][mask_array]
    snp_ga = allel.GenotypeArray(snp_gts)

    ac = snp_ga.count_alleles()

    return snp_pos, ac

def get_accessible(bed, chrom):
    df = pd.read_csv(bed, sep = '\t', columns = ["chrom","start","stop"])
    degen_pos = df[df['chrom']==chrom]["stop"].to_numpy()
    return degen_pos

if __name__ == "__main__":
    args = docopt(__doc__)

    vcf_f = args["--vcf"]
    bed_f = args["--bed"]
    acc_bed_f = args["--accessible"]

    if args["--name"]:
        name = str(args["--name"])
    else:
        name = "NA"

    if args["--output"]:
        out_f = str(args["--output"])
    else:
        out_f = "bed_pi.tsv"

    if args["--label"]
        result_label = str(args["--label"])
    else:
        result_label = "pi"

    results = []

    #could skip the vcf filtering step since is_accessible will only consider the right sites
    #the filtering itself is quite slow 
    #but, the vcf loading is slow and will be loaded for every chromosome
    accessible_array = get_accessible(bed_df, name)
    snp_pos, ac = parse_vcf(vcf_f, chromosome=name)
    bed = pybedtools.BedTool(bed_f)

    for interval in bed:
        (str(chrom), int(start), int(stop), str(feature)) = interval
        pi = allele.sequence_diversity(snp_pos, ac, start=start+1, stop=stop, is_accessible=accessible_array)
        results.append((chrom, feature, pi))

    pd.DataFrame(results, columns=["chrom","feature", result_label]).to_csv(
        out_f, sep="\t", index=False
    )