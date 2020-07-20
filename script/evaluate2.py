#!/usr/bin/env python

import sys
from glob import glob
from pathlib import Path

from iseq.profmark import ProfMark

folder = Path(sys.argv[1])
acc = sys.argv[2]

hmmfile = glob(str(folder / "*.hmm"))[0]
tgtfile = folder / acc / "cds_nucl.fasta"
domtbloutfile = folder / acc / "domtblout.txt"
iseqoutfile = folder / acc / "output.gff"

profmark = ProfMark(hmmfile, tgtfile, domtbloutfile, iseqoutfile)
profmark.write_pickle(Path("profmark.pkl"))
