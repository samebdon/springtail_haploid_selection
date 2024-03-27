#!/usr/bin/env bash

IN_DIR="$1"
QUERY_PROT="$2"
OUT_DIR="$3"
ORTHOGROUP="$4"

FILES_A="$($IN_DIR/*.1.cds.fasta)"
FILES_B="$($IN_DIR/*.2.cds.fasta)"

for FILE in $FILES_A
do
        SAMPLE="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f-1)"
        SPECIES="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f2-2)"
        cat $FILE | grep $QUERY_PROT -A1 >> $OUT_DIR/$ORTHOGROUP.$SAMPLE.$SPECIES.fa
done

for FILE in $FILES_B
do
        SAMPLE="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f-1)"
        SPECIES="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f2-2)"
        cat $FILE | grep $QUERY_PROT -A1 >> $OUT_DIR/$ORTHOGROUP.$SAMPLE.$SPECIES.fa
done