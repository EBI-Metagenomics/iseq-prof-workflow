#!/bin/bash

new_header="Version\tMolType\tBasePairs\tOrganism\tTaxID"
old_header="Accession\tVersion\tID\tMolType\tBasePairs\tOrganism\tTaxID\tDB\tBioProject\tBioSample"

echo -e $new_header > gb238.catalog.all.tsv

for db in est gss other;
do
   printf "Processing ${db}... "
   curl -s ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/gb238.catalog.${db}.txt.gz \
      | gunzip -c  \
      | cut -d$'\t' -f2,4,5,6,7 \
      | grep --invert-match $'\tmRNA\t' \
      | grep --invert-match $'\trRNA\t' \
      | grep --invert-match $'\ttRNA\t' \
      | grep --invert-match $'\tncRNA\t' \
      | grep --invert-match $'\tNoTaxID' \
      >> gb238.catalog.all.tsv
   echo "done."
done
