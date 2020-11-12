#!/usr/bin/env python

import sys
from pathlib import Path

from iseq_prof import GenBank

accession = sys.argv[1]
rettype = sys.argv[2]
output = Path(sys.argv[3])

GenBank.download(accession, rettype, output)
