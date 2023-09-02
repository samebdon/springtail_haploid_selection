"""mergePi.py

Usage:
mergePi.py -z <FILE> -f <FILE> -o <STR> [-h]

Options:
-z, --zero_pi <FILE>	       Gene pi tsv
-f, --four_pi <FILE>	       Gene pi tsv
-o, --output <FILE>    		   Output file name
"""

import pandas as pd
import numpy as np
from docopt import docopt

if __name__ == "__main__":
    args = docopt(__doc__)

    zero_f = args["--zero_pi"]
    four_f = args["--four_pi"]
    out_f = str(args["--output"])

    zero_df = pd.read_csv(
        zero_f, sep="\t", dtype={"chrom": "category", "feature": "category"}
    )
    four_df = pd.read_csv(
        four_f, sep="\t", dtype={"chrom": "category", "feature": "category"}
    )

    zero_df['4D_pi'] = four_df['4D_pi']
    zero_df['0/4D_pi'] =  zero_df['0D_pi']/zero_df['4D_pi']

    zero_df.to_csv(out_f,sep="\t",index=False, na_rep='NaN')

