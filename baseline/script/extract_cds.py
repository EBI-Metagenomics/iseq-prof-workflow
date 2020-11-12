#!/usr/bin/env python

import sys
from pathlib import Path

from iseq_prof import GenBank

gb = GenBank(Path(sys.argv[1]))
amino_filepath = Path(sys.argv[2])
nucl_filepath = Path(sys.argv[3])

gb.extract_cds(amino_filepath, nucl_filepath)
