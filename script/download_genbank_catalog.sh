#!/bin/bash

header="Accession\tVersion\tID\tMolType\tBasePairs\tOrganism\tTaxID\tDB\tBioProject\tBioSample"

db=$1
output=$2

url=ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog
gbversion=$(curl --list-only -s $url/ | head -n 1 | cut -d. -f1)

# The DNA sequence for Porcine circovirus type 2 strain MLP-22
# is 1726 base pairs long.
curl -s $url/${gbversion}.catalog.${db}.txt.gz \
   | gunzip -c  \
   | cut -d$'\t' -f2,4,5,6,7 \
   | grep $'\\(\tRNA\t\\|\tDNA\t\\)' \
   | grep --invert-match $'\tNoTaxID' \
   | awk -F '\t' '{ if ($3 >= 1726) { print } }' \
   >> $output

if [ ! -s $ouput ]
then
   echo "ERROR: file $output is empty."
   exit 1
fi
