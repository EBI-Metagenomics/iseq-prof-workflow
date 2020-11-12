#!/usr/bin/env python

import sys

from fasta_reader import FASTAWriter, open_fasta

infasta = sys.argv[1]
outfasta = sys.argv[2]
n = int(sys.argv[3])

with FASTAWriter(outfasta) as writer:
    for i, target in enumerate(open_fasta(infasta)):
        if i == n:
            break
        writer.write_item(target.defline, target.sequence)
