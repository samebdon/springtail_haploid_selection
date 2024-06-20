#!/bin/sh

## Usage ./splitBed.sh input.bed

input=$1

mkdir beds
for chr in `cut -f 1 $input | sort | uniq`;
do
	grep -w $chr $input > beds/$chr.bed
done
