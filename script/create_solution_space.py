#!/usr/bin/env python

import hashlib
import random
import sys
from subprocess import check_call

from hmmer import read_domtbl
from hmmer_reader import fetch_metadata

domtbl = sys.argv[1]
hmmdb = sys.argv[2]
keyfile = sys.argv[3]
newdb = sys.argv[4]
seed = int(sys.argv[5])

with open(domtbl, "rb") as wfile:
    seed += int(hashlib.md5(wfile.read()).hexdigest()[:4], 16)

true_accs = set()
print("Reading domtblout file...")
for row in read_domtbl(domtbl):
    true_accs.add(row.target.accession)

print("Fetching metadata...")
all_accessions = fetch_metadata(hmmdb)["ACC"].tolist()

remain = list(set(all_accessions) - true_accs)

random.shuffle(remain)
nfalse = min(len(true_accs), len(remain))
false_accs = set(remain[i] for i in range(nfalse))

print(f"Saving {keyfile}...")
with open(keyfile, "w") as file:
    for acc in all_accessions:
        if acc in true_accs:
            file.write(acc + "\n")
        elif acc in false_accs:
            file.write(acc + "\n")

cmd = f"hmmfetch -f {hmmdb} {keyfile} > {newdb}"
print(f"Running {cmd}...")
check_call(cmd, shell=True)
