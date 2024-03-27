#!/usr/bin/env bash

IN_DIR="$1"
QUERY_PROT="$2"
OUT_DIR="$3"
ORTHOGROUP="$4"

for FILE in $IN_DIR/*.1.cds.fasta
do
        SAMPLE="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f-1)"
        SPECIES="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f2-2)"
        cat $FILE | grep $QUERY_PROT -A1 >> $OUT_DIR/$ORTHOGROUP.$SAMPLE.$SPECIES.fa
done

for FILE in $IN_DIR/*.2.cds.fasta
do
        SAMPLE="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f-1)"
        SPECIES="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f2-2)"
        cat $FILE | grep $QUERY_PROT -A1 >> $OUT_DIR/$ORTHOGROUP.$SAMPLE.$SPECIES.fa
done