#!/usr/bin/env python

import sys
from pathlib import Path

from iseq.profmark import genbank

gb_filepath = sys.argv[1]

print(genbank.get_accession(Path(gb_filepath)), end="")
