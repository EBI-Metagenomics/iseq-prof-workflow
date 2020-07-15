#!/usr/bin/env python

import sys
from pathlib import Path

from iseq.profmark import genbank

accession = sys.argv[1]
rettype = sys.argv[2]
output = Path(sys.argv[3])

genbank.download(accession, rettype, output)
