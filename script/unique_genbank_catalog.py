import sys
from collections import OrderedDict

import dask.dataframe as dd

dtype = OrderedDict()
dtype["Accession"] = "str"
dtype["Version"] = "str"
# dtype["ID"] = "int64"
dtype["MolType"] = "category"
dtype["BasePairs"] = "int64"
dtype["Organism"] = "str"
dtype["TaxID"] = "int64"
# dtype["DB"] = "category"
# dtype["BioProject"] = "str"
# dtype["BioSample"] = "str"

input_file = sys.argv[1]
output_file = sys.argv[2]

names = list(dtype.keys())
df = dd.read_csv(input_file, sep="\t", header=0, dtype=dtype)
# df = pd.read_csv(input_file, sep="\t", header=0, dtype=dtype)
# df = pd.read_csv(input_file, sep="\t", header=None, names=names, dtype=dtype)

max_bp = {key: row.BasePairs for key, row in df.groupby("Organism").max().iterrows()}
df["longest"] = df.apply(
    lambda row: max_bp[row.Organism] == row.BasePairs, axis=1, meta=bool
)
# df["longest"] = df.apply(lambda row: max_bp[row.Organism] == row.BasePairs, axis=1)
df = df[df["longest"]]

# max_id = {key: row.ID for key, row in df.groupby("Organism").max().iterrows()}
# df["latest"] = df.apply(lambda row: max_id[row.Organism] == row.ID, axis=1)
# df = df[df["latest"]].copy()
df = df.drop_duplicates(subset=["Organism"], keep="last")
df = df.drop(columns=["longest"])

# df = df.reset_index(drop=True)
# df.to_feather(output_file)
df.to_csv(output_file, encoding="ascii", sep="\t", single_file=True, index=False)
