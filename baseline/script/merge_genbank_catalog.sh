#!/bin/bash

header="Version\tMolType\tBasePairs\tOrganism\tTaxID"

input1=$1
input2=$2
output=$3

echo -e $header > $output
cat $input1 >> $output
cat $input2 >> $output
