"""filterDegenotate.py

Usage:
filterDegenotate.py -b <FILE> -t <FILE> -o <STR> [-h]

Options:
-b, --bed <FILE>           Degeneracy bed file
-t, --target <FILE>        Fasta names
-o, --output <STR>         Output prefix

"""

from docopt import docopt
import pandas as pd
import numpy as np

if __name__ == "__main__":
    args = docopt(__doc__)

    bed_file = args['--bed']
    target_file = args['--target']
    prefix = args['--output']

    print("reading degen")
    degen_bed = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "stop",
            "transcript",
            "degeneracy",
        ],
        dtype={'chrom' : 'category',
        'transcript': 'category', ## str might be better
        'degeneracy': 'category'}
    )

    print("reading target")
    target_bed = pd.read_csv(
        target_file,
        sep="\t",
        header=None,
        names=["transcript"],
    )

    out_columns = ["chrom", "start", "stop"]

    target_degen_bed = degen_bed[degen_bed["transcript"].isin(target_bed["transcript"])]

    target_degen_bed[target_degen_bed["degeneracy"] == '0'][out_columns].to_csv(
        prefix + ".0D.bed", sep="\t", header=False, index=False
    )
    target_degen_bed[target_degen_bed["degeneracy"] == '4'][out_columns].to_csv(
        prefix + ".4D.bed", sep="\t", header=False, index=False
    )
