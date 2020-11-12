#!/usr/bin/env python

import os
import sys
from pathlib import Path

if len(sys.argv) != 3:
    script_name = os.path.basename(__file__)
    print(f"Usage: {script_name} HMM_FILEPATH PROFILE")
    sys.exit(1)

hmmfilepath = Path(sys.argv[1])
profile = str(sys.argv[2]).strip().partition(".")[0]

ok = False
with open(hmmfilepath, "r") as f:
    for row in f:
        if row.startswith(f"ACC   {profile}"):
            row = row.strip()
            accession = row.partition("   ")[2].strip()
            print(accession)
            ok = True
            break

if not ok:
    sys.exit(1)
