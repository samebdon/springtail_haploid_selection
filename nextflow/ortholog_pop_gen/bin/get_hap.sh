#!/usr/bin/env bash

get_hap () {
for FILE in $1
do
        SAMPLE="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f-1)"
        SPECIES="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f2-2)"
        cat $FILE | grep $2 -A1 >> $3/$4.$SAMPLE.$SPECIES.fa
done
}

IN_DIR="$1"
QUERY_PROT="$2"
OUT_DIR="$3"
ORTHOGROUP="$4"

FILES_A="$($IN_DIR/*.1.cds.fasta)"
FILES_B="$($IN_DIR/*.2.cds.fasta)"

get_hap $FILES_A $QUERY_PROT $OUT_DIR $ORTHOGROUP
get_hap $FILES_B $QUERY_PROT $OUT_DIR $ORTHOGROUP