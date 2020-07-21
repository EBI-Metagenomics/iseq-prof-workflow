#!/bin/bash

header="Accession\tVersion\tID\tMolType\tBasePairs\tOrganism\tTaxID\tDB\tBioProject\tBioSample"

db=$1
output=$2

curl -s ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/gb238.catalog.${db}.txt.gz \
   | gunzip -c  \
   | cut -d$'\t' -f2,4,5,6,7 \
   | grep $'\\(\tRNA\t\\|\tDNA\t\\)' \
   | grep --invert-match $'\tNoTaxID' \
   >> $output
