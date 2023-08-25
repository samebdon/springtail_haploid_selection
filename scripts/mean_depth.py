"""mean_depth.py

Usage:
mean_depth.py -b <FILE> [-h -m <INT>]

Options:
-b, --bed <FILE>	Mosdepth bed file
-m, --multiplier <INT> Multiplication factor

"""

import pandas as pd
import numpy as np
from docopt import docopt


def get_mean_depth_mosdepth(bed_f):
    bed_df = pd.read_csv(
        bed_f,
        compression="gzip",
        sep="\t",
        names=["sequence_id", "start", "end", "depth"],
        dtype={"sequence_id": str, "start": int, "end": int, "depth": int},
    )
    bed_df = bed_df[bed_df["depth"] > 0]
    bed_df["length"] = bed_df["end"] - bed_df["start"]
    mean = round(np.average(bed_df["depth"], weights=bed_df["length"]), 2)
    return mean


if __name__ == "__main__":
    args = docopt(__doc__)
    bed_f = args["--bed"]
    if args["--multiplier"]:
        factor = int(args["--multiplier"])
    else:
        factor = int(1)
    print(factor * get_mean_depth_mosdepth(bed_f))
