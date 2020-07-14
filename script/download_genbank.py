#!/usr/bin/env python

import sys

from iseq.profmark import download_genbank

accession = sys.argv[1]
rettype = sys.argv[2]

download_genbank(accession, rettype)
