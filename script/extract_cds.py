#!/usr/bin/env python

import sys
from pathlib import Path

from iseq.profmark import genbank

gb_filepath = Path(sys.argv[1])
amino_filepath = Path(sys.argv[2])
nucl_filepath = Path(sys.argv[3])

genbank.extract_cds(gb_filepath, amino_filepath, nucl_filepath)
