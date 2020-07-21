#!/usr/bin/env python

import sys
from pandas import read_csv
import random
import time
from Bio import Entrez, SeqIO


def get_major(organism: str):
    return organism.split("_")[0]


def get_accession(df, organism: str):
    major = get_major(organism)
    df0 = df.loc[df.Organism.str.startswith(major)]

    index = df0.index.values.tolist()
    random.shuffle(index)

    accession = ""
    for i in index:
        time.sleep(2)

        whole_genome = False
        row = df0.loc[i]
        accession = row["Version"]

        with Entrez.esummary(db="nucleotide", id=accession) as handle:
            records = Entrez.parse(handle)
            try:
                record = next(records)
            except RuntimeError:
                continue

            if "whole genome shotgun sequence" in record["Title"]:
                whole_genome = True
            if "complete genome" in record["Title"]:
                whole_genome = True

        if not whole_genome:
            continue

        with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
            record = next(SeqIO.parse(handle, "genbank"))
            if len(record.features) <= 1:
                continue

        return accession
    return ""


Entrez.email = "horta@ebi.ac.uk"

catagfile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]
size = int(sys.argv[4])
seed = int(sys.argv[5])

random.seed(seed)

df = read_csv(catagfile, sep="\t", header=0)
# The DNA sequence for Porcine circovirus type 2 strain MLP-22 is 1726 base pairs long.
df = df[df["BasePairs"] >= 1726]

organisms = [r.strip() for r in open(infile, "r").readlines()]
random.shuffle(organisms)

accessions = []
majors = set()
for organism in organisms:
    major = get_major(organism)
    if major in majors:
        continue
    majors.add(major)
    acc = get_accession(df, organism)
    if len(acc) > 0:
        accessions.append(acc)
    if len(accessions) == size:
        break

with open(outfile, "w") as file:
    for acc in accessions:
        file.write(acc + "\n")
