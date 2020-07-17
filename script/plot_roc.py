import sys
from pathlib import Path

from iseq.profmark import ProfMark

pm = ProfMark.read_pickle(Path(sys.argv[1]))
pm.confusion_matrix.roc.savefig("roc.png")
