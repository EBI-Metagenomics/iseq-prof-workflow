#!/usr/bin/env python

import random
import sys
import time
from http.client import IncompleteRead
from typing import Callable
from urllib.error import HTTPError

from Bio import Entrez, SeqIO
from Bio.Data import IUPACData
from pandas import read_csv

taxname = {"virus": "Viruses", "archaea": "Archaea", "bacteria": "Bacteria"}
MIN_NUM_CDS = 10


def get_major(organism: str):
    return organism.split("_")[0]


def is_alphabet_ambiguous(seq):

    remains = len(set(str(seq)) - set(IUPACData.unambiguous_dna_letters))
    if remains == 0:
        return False

    remains = len(set(str(seq)) - set(IUPACData.unambiguous_rna_letters))
    if remains == 0:
        return False

    return True


def is_whole_genome(accession: str) -> bool:
    whole_genome = False

    with Entrez.esummary(db="nucleotide", id=accession) as handle:
        records = Entrez.parse(handle)
        try:
            record = next(records)
        except RuntimeError:
            print(f"ERROR: {accession} is not whole genome.", file=sys.stderr)
            return whole_genome

        if "whole genome shotgun sequence" in record["Title"]:
            whole_genome = True
        if "complete genome" in record["Title"]:
            whole_genome = True

    if not whole_genome:
        print(f"WARNING: {accession} is not whole genome.", file=sys.stderr)
    return whole_genome


def is_nice_data(accession: str, family: str):
    num_nice_cds = 0
    with Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text",
    ) as handle:
        record = next(SeqIO.parse(handle, "genbank"))
        if taxname[family] not in record.annotations["taxonomy"]:
            return False

        for feature in record.features:
            if feature.type.lower() == "cds":
                if not is_alphabet_ambiguous(feature.extract(record).seq):
                    num_nice_cds += 1
                # Avoid sequences that are not whole genome,
                # despite its title
                if num_nice_cds >= MIN_NUM_CDS:
                    return True
    msg = f"WARNING: Skip {accession} because it "
    msg += f"has less than {MIN_NUM_CDS} "
    msg += f"unambiguous CDSs ({num_nice_cds})."
    print(msg, file=sys.stderr)
    return False


def try_hard(func: Callable[..., bool]) -> bool:
    tries = 7
    while True:
        tries -= 1
        try:
            return func()
        except (HTTPError, IncompleteRead) as e:
            if tries == 0:
                raise e
            time.sleep(10)


def get_accession(df, organism: str, family: str):
    major = get_major(organism)
    df0 = df.loc[df.Organism.str.startswith(major)]

    index = df0.index.values.tolist()
    random.shuffle(index)

    accession = ""
    for i in index:
        row = df0.loc[i]
        accession = row["Version"]

        if not try_hard(lambda: is_whole_genome(accession)):
            continue

        if not try_hard(lambda: is_nice_data(accession, family)):
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
if len(df) == 0:
    raise RuntimeError(f"Empty dataframe. Check {catagfile}.")
# The DNA sequence for Porcine circovirus type 2 strain MLP-22 is 1726 base
# pairs long.
df = df[df["BasePairs"] >= 1726]

organisms = [r.strip() for r in open(infile, "r").readlines()]
print(f"Found {len(organisms)} organisms to check.")
random.shuffle(organisms)

accessions = set()
majors = set()
for organism in organisms:
    major = get_major(organism)
    if major in majors:
        continue
    majors.add(major)
    acc = get_accession(df, organism, family)
    if len(acc) > 0 and acc not in accessions:
        accessions.add(acc)
        print(f"New accession {len(accessions)}: {acc}")
    if len(accessions) == size:
        break

if len(accessions) != size:
    raise RuntimeError(f"Could not fetch {size} accessions.")

with open(outfile, "w") as file:
    for acc in accessions:
        file.write(acc + "\n")
