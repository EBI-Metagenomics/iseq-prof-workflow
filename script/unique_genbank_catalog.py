#!/usr/bin/env python

import sys
from collections import OrderedDict

import pandas as pd
from tqdm import tqdm

dtype = OrderedDict()
dtype["Version"] = "str"
dtype["MolType"] = "category"
dtype["BasePairs"] = "int64"
dtype["Organism"] = "str"
dtype["TaxID"] = "int32"

input_file = sys.argv[1]
output_file = sys.argv[2]

names = list(dtype.keys())

with tqdm(total=3) as pbar:
    df = pd.read_csv(input_file, sep="\t", header=0, dtype=dtype)
    pbar.update(1)

    df = df.loc[df.reset_index().groupby(["Organism"])["BasePairs"].idxmax()]
    pbar.update(1)

    df.to_csv(output_file, encoding="ascii", sep="\t", index=False)
    pbar.update(1)
