#!/usr/bin/env bash

FASTA="$1"

ORTHOGROUP="$(ls $FILE | cut -d'.' -f-1)"

SP_ONE="$(cat $FILE | grep '>' | head -n1 | cut -d'>' -f2 | cut -d'.' -f-1)"
PROT_ONE="$(cat $FILE | grep '>' | head -n1 | cut -d'>' -f2 | cut -d'.' -f2-)"
SEQ_ONE="$(cat $FILE | sed '2q;d')"

SP_TWO="$(cat $FILE | grep '>' | tail -n1 | cut -d'>' -f2 | cut -d'.' -f-1)"
PROT_TWO="$(cat $FILE | grep '>' | tail -n1 | cut -d'>' -f2 | cut -d'.' -f2-)"
SEQ_TWO="$(cat $FILE | sed '4q;d')"

HEAD_ONE=">${SP_ONE}.sample_1.1.${PROT_ONE}"
HEAD_TWO=">${SP_TWO}.sample_2.1.${PROT_TWO}"
HEAD_THREE=">${SP_ONE}.sample_1.2.${PROT_ONE}"
HEAD_FOUR=">${SP_TWO}.sample_2.2.${PROT_TWO}"

OUTFILE="${ORTHOGROUP}.mafft.happed.fa"

echo "${HEAD_ONE}" >> $OUTFILE
echo "${SEQ_ONE}" >> $OUTFILE

echo "${HEAD_TWO}" >> $OUTFILE
echo "${SEQ_TWO}" >> $OUTFILE

echo "${HEAD_THREE}" >> $OUTFILE
echo "${SEQ_ONE}" >> $OUTFILE

echo "${HEAD_FOUR}" >> $OUTFILE
echo "${SEQ_TWO}" >> $OUTFILE
