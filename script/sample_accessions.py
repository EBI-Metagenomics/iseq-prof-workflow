#!/usr/bin/env python

import sys
from pandas import read_csv
import random
import time
from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC, DNAAlphabet, RNAAlphabet

taxname = {"virus": "Viruses", "archaea": "Archaea", "bacteria": "Bacteria"}


def get_major(organism: str):
    return organism.split("_")[0]


def is_alphabet_ambiguous(seq):

    if isinstance(seq.alphabet, DNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_dna.letters))
        if remains > 0:
            return True

    elif isinstance(seq.alphabet, RNAAlphabet):
        remains = len(set(str(seq)) - set(IUPAC.unambiguous_rna.letters))
        if remains > 0:
            return True

    else:
        raise ValueError("Unkown alphabet.")

    return False


def get_accession(df, organism: str, family: str):
    major = get_major(organism)
    df0 = df.loc[df.Organism.str.startswith(major)]

    index = df0.index.values.tolist()
    random.shuffle(index)

    accession = ""
    for i in index:
        time.sleep(1.2)

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

        time.sleep(0.2)
        with Entrez.efetch(
            db="nuccore", id=accession, rettype="gb", retmode="text"
        ) as handle:
            record = next(SeqIO.parse(handle, "genbank"))
            if taxname[family] not in record.annotations["taxonomy"]:
                continue

            for feature in record.features:
                if feature.type.lower() == "cds":
                    if not is_alphabet_ambiguous(feature.extract(record).seq):
                        break
            else:
                continue

        return accession
    return ""


Entrez.email = "horta@ebi.ac.uk"

catagfile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]
size = int(sys.argv[4])
seed = int(sys.argv[5])

family = infile.split(".")[0]

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
    acc = get_accession(df, organism, family)
    if len(acc) > 0:
        accessions.append(acc)
        print(f"New accession {len(accessions)}: {acc}")
    if len(accessions) == size:
        break

if len(accessions) != size:
    raise RuntimeError(f"Could not fetch {size} accessions.")

with open(outfile, "w") as file:
    for acc in accessions:
        file.write(acc + "\n")
