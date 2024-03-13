#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage:
        choose_best_cds_isoforms.py -t <FILE> [-s <FLOAT> -l <INT>]

Options:
        -h, --help                      show this
        -t, --transdecoder_bed <FILE>   Transdecoder BED file
        -s, --score <FLOAT>             Transdecoder score [default: 0.00]
        -l, --length <INT>              Transdecoder protein length [default: 30]

Description:
    prints BED feature of best scoring isoform CDSs per gene based on a transdecoder BED file, if
    - status is 'complete'
    - score and length satisfy thresholds

"""

from __future__ import division
from docopt import docopt
import sys
from os.path import isfile
import re


def natsort(l):

    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split(_nsre, key)]

    return sorted(l, key=alphanum_key)


def readBed(infile):
    if not isfile(infile):
        sys.exit("[X] %s is not a file." % (infile))
    with open(infile) as fh:
        for l in fh:
            if not l.startswith("track"):
                field = l.split()
                isoform_id = field[0]
                protein_id = field[3].split("~")[2].split(";")[0]
                status = field[3].split(";")[2].split(":")[1].split("_")[0]
                protein_length = int(field[3].split(":")[2].split("_")[0].split("_")[0])
                cds_length = protein_length * 3
                score = float(field[3].split("=")[2].split(",")[0])
                strand = field[5]
                start = int(field[6])
                end = int(field[7])
                if status == 'complete':
                    if protein_length >= MIN_LENGTH:
                        if score >= MIN_SCORE:
                            bedObj = BedObj(isoform_id, protein_id, status, protein_length, cds_length, score, strand, start, end)
                            yield bedObj


class BedObj():
    def __init__(self, isoform_id, protein_id, status, protein_length, cds_length, score, strand, start, end):
        self.isoform_id = isoform_id
        self.gene_id = "_".join(isoform_id.split('_')[0:4])
        self.protein_id = protein_id
        self.status = status
        self.protein_length = protein_length
        self.cds_length = cds_length
        self.score = score
        self.strand = strand
        self.start = start
        self.end = end

    def bed(self):
        return "\t".join([str(x) for x in [self.isoform_id, self.start, self.end, self.gene_id, self.score, self.strand, self.cds_length, self.protein_id]])


class MainObj():
    def __init__(self):
        self.bedObjs_by_gene_id = {}
        self.parse_bed()
        self.write_bed()

    def parse_bed(self):
        for bedObj in readBed(BED):
            # print bedObj.__dict__
            try:
                self.bedObjs_by_gene_id[bedObj.gene_id].append(bedObj)
            except KeyError:
                self.bedObjs_by_gene_id[bedObj.gene_id] = [bedObj]

    def write_bed(self):
        for gene_id in natsort(self.bedObjs_by_gene_id):
            for idx, bedObj in enumerate(sorted(self.bedObjs_by_gene_id[gene_id], key=lambda x: (x.score), reverse=True)):
                if idx >= 1:
                    break
                print bedObj.bed()

if __name__ == "__main__":
    _nsre = re.compile('([0-9]+)')
    args = docopt(__doc__)

    BED = args['--transdecoder_bed']
    MIN_SCORE = float(args['--score'])
    MIN_LENGTH = int(args['--length'])

    MainObj()


