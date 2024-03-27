#!/usr/bin/env bash

IN_DIR="$1"
OUT_DIR="$2"

for FILE in $IN_DIR
do
	ORTHOGROUP="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f-1)"
	SAMPLE="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f2-2)"
	SPECIES="$(ls $FILE | cut -d'/' -f2- | cut -d'.' -f3-3)"
	GENE="$(cat $FILE | grep '>' | head -n1 | cut -d'>' -f2- | cut -d' ' -f-1)"

	HEAD_ONE=">${SPECIES}.${SAMPLE}.1.${GENE}"
	SEQ_ONE="$(cat $FILE | sed '2q;d')"
	HEAD_TWO=">${SPECIES}.${SAMPLE}.2.${GENE}"
	SEQ_TWO="$(cat $FILE | sed '4q;d')"

	OUTFILE="${OUT_DIR}/${ORTHOGROUP}.${SAMPLE}.${SPECIES}.fa"

	echo "${HEAD_ONE}" >> $OUTFILE
	echo "${SEQ_ONE}" >> $OUTFILE
	echo "${HEAD_TWO}" >> $OUTFILE
	echo "${SEQ_TWO}" >> $OUTFILE
done