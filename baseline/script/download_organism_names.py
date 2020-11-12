#!/usr/bin/env python

import random
import sys
from ftplib import FTP, error_perm
from typing import List


def isdir(ftp, item):
    if item[1]["type"] == "dir":
        return True

    if item[1]["type"] == "file":
        return False

    try:
        list(ftp.mlsd(item[0]))
    except error_perm as e:
        if "is not a directory" in e.args[0]:
            return False
        raise e

    return True


def list_dirs(ftp) -> List[str]:
    return list(sorted([f[0] for f in list(ftp.mlsd())[2:] if isdir(ftp, f)]))


def fetch_organisms(ftp) -> List[str]:
    organisms: List[str] = []
    for di in list_dirs(ftp):
        major = di.split("_")[0]
        if len(major) == 0:
            continue

        if major[0] == major[0].lower():
            continue

        organisms.append(di)

    return organisms


random.seed(4872)
HOST = "ftp.ncbi.nlm.nih.gov"

ftp = FTP(HOST)
ftp.login()

orgtype = sys.argv[1]

ftp.cwd(f"/genomes/refseq/{orgtype}")

organisms = fetch_organisms(ftp)
with open(f"{orgtype}.txt", "w") as file:
    for organism in organisms:
        file.write(f"{organism}\n")
