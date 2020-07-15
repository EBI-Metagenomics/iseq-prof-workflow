#!/usr/bin/env python

import sys
from pathlib import Path

from iseq.profmark import ProfMark

hmmer_file = Path(sys.argv[1])
target_file = Path(sys.argv[2])
domtblout_file = Path(sys.argv[3])
iseqout_file = Path(sys.argv[4])
profmark_file = Path(sys.argv[5])

profmark = ProfMark(hmmer_file, target_file, domtblout_file, iseqout_file)
profmark.write_pickle(profmark_file)
