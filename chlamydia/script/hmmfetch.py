#!/usr/bin/env python

import sys
from pathlib import Path

from hmmer import HMMER

hmm_filepath = Path(sys.argv[1])
accs_filepath = Path(sys.argv[2])

hmmer = HMMER(hmm_filepath)

with open(accs_filepath, "r") as f:
    accs = [w.strip() for w in f.readlines()]

if len(set(accs)) < len(accs):
    print("Duplicated accessions.", file=sys.stderr)
    sys.exit(1)

hmm = hmmer.fetch(accs)

naccs = hmm.count("ACC   ")
if naccs != len(accs):
    msg = f"Number of profiles {naccs} is different from requested {len(accs)}."
    print(msg, file=sys.stderr)
    sys.exit(1)

print(hmm)
