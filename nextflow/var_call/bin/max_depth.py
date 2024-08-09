#!/usr/bin/env python3 

"""max_depth.py

Usage:
max_depth.py -b <FILE> [-h -m <INT>]

Options:
-b, --bed <FILE>	Mosdepth bed file
-m, --multiplier <INT> Multiplication factor (default: 2)

"""

import pandas as pd
import numpy as np
from docopt import docopt

def get_max_depth_mosdepth(bed_f, factor):
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
    #SD = np.std(bed_df["depth"])
    return int(mean+(factor*SD))
    
    

if __name__ == "__main__":
    args = docopt(__doc__)
    bed_f = args["--bed"]
    if args["--multiplier"]:
        factor = float(args["--multiplier"])
    else:
        factor = int(2)
        
    print(get_max_depth_mosdepth(bed_f, factor))
