#!/usr/bin/env python3 

"""mean_depth.py

Usage:
mean_depth.py -b <FILE> [-h -m <INT>]

Options:
-b, --bed <FILE>	Mosdepth bed file
-m, --multiplier <INT> Multiplication factor (default: 2)

"""

import pandas as pd
import numpy as np
from docopt import docopt

def get_max_depth_mosdepth(bed_f):
    bed_df = pd.read_csv(
        bed_f,
        compression="gzip",
        sep="\t",
        names=["sequence_id", "start", "end", "depth"],
        dtype={"sequence_id": str, "start": int, "end": int, "depth": int},
    )
    bed_df = bed_df[bed_df["depth"] > 0]
    bed_df["length"] = bed_df["end"] - bed_df["start"]
    mean = np.average(bed_df["depth"], weights=bed_df["length"])
    variance = np.average((bed_df["depth"]-mean)**2, weights=bed_df["length"])
    SD = np.sqrt(variance)
    return round(mean+(2*SD), 2)

if __name__ == "__main__":
    args = docopt(__doc__)
    bed_f = args["--bed"]
    if args["--multiplier"]:
        factor = int(args["--multiplier"])
    else:
        factor = int(2)
        
    print(int(get_max_depth_mosdepth(bed_f)))
