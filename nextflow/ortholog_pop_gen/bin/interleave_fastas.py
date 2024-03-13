#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage:
        interleave_fastas.py [-p <STR> -d <STR> -s <INT> -a <INT> -f] FASTAS ...

Options:
        -h --help                       show this
        -p, --out_prefix <STR>          Output prefix [default: ilv]
        -d, --delimiter <STR>           Delimiter by which the filename will be split [default: .]
        -s, --sample <INT>              Item containing sample_id after filename was split by delimiter [default: 1]
        -a, --allele <INT>              Item containing allele_id after filename was split by delimiter [default: 3]
        -f, --force                     Force interleave even when sequence_ids differ
"""

from __future__ import division
from docopt import docopt
import sys
from os.path import basename, isfile, join
from collections import Counter

def readFasta(infile):
    if not isfile(infile):
        sys.exit("[X] %s is not a file." % (infile))
    with open(infile) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if header:
                    yield header, ''.join(seqs).upper()
                header, seqs = ''.join(l[1:-1].split()[0]).split(":")[0], []  # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs).upper()


class SeqObj():
    def __init__(self, sample_id, allele_id, header, seq):
        self.sample_id = sample_id
        self.allele_id = allele_id
        self.header = header
        self.seq = seq
        self.length = len(seq)
        self.codon_counter = Counter(self.get_codon_count())

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.header == other.header
        else:
            return False

    def get_fasta(self):
        return ">%s.%s.%s\n%s" % (self.sample_id, self.allele_id, self.header, self.seq)

    def get_phylip(self):
        return "%s.%s\t%s" % (self.allele_id, self.header, self.seq)

    def id_length(self):
        return "%s : %s" % (self.header, self.length)

    def get_codon_count(self):
        codon_count = {codon: 0 for codon in CODONS}
        for i in range(0, len(self.seq), 3):
            try:
                codon_count[self.seq[i: i + 3]] += 1
            except KeyError:
                pass
        return codon_count

    def has_start(self):
        if self.seq[0:3] == 'ATG':
            return True
        return False

    def has_stop(self):
        if self.seq[-3:] in set(['TAA', 'TAG', 'TGA']):
            return True
        return False

class MainObj():
    def __init__(self):
        self.loci = []
        self.allele_ids = []
        self.sample_id = None
        self.parse_fastas()
        self.total_length = 0
        self.total_sequences = 0
        self.codon_counter_by_allele_id = {}
        self.write_outputs()

    def parse_fastas(self):
        print "[+] Parsing files ..."
        for fasta in FASTAS:
            idx = 0
            filename = basename(fasta).split(DELIMITER)
            self.sample_id = "%s.%s" % (PREFIX, filename[SAMPLE_IDX])
            self.allele_ids.append("%s.%s" % (self.sample_id, filename[ALLELE_IDX]))
            for header, seq in readFasta(fasta):
                if len(set(seq)) > 1:  # check that not only "-"
                    seqObj = SeqObj(self.sample_id, self.allele_ids[-1], header, seq)
                    try:
                        self.loci[idx].append(seqObj)
                    except IndexError:
                        self.loci.append([seqObj])
                    idx += 1

    def write_outputs(self):
        print "[+] Writing files ..."
        loci_output = ['#locus\tlength\tstart_frac\tstop_frac']
        loci_outfile = PREFIX + ".loci.txt"
        codon_output = ['allele_id\tn_loci\tn_codons\t%s' % "\t".join(CODONS)]
        codon_outfile = PREFIX + ".codon.txt"
        phylip_output = []
        phylip_outfile = PREFIX + ".phy"
        # padding = len(str(len(self.loci)))
        for idx, loci in enumerate(self.loci):
            # out_file = "%s.%s.fasta" % (out_path, str(idx).rjust(padding, '0'))
            if not FORCE:
                if not all(seqObj == loci[0] for seqObj in loci):
                    sys.exit("[X] Headers of loci %s are different :\n\t%s\nUse '--force' to ignore." % (idx, "\n\t".join(seqObj.header for seqObj in loci)))
            if not len(loci) == len(FASTAS):
                sys.exit("[X] Block %s only has %s sequences. Is one of the files truncated?" % (idx, len(loci)))
            # check for same length of sequences
            if not all(seqObj.length == loci[0].length for seqObj in loci):
                sys.exit("[X] Block %s has sequences of different lengths.\n%s" % (idx, len(loci), "\n".join([seqObj.id_length() for seqObj in loci])))
            if not loci[0].length % 3 == 0:
                off_by = 0
                if loci[0].length % 3 == 1 / 3:
                    off_by = 1
                else:
                    off_by = 2
                sys.exit("[X] Block %s contains sequences whose length is not divisible by three (off-by: %s)\n%s" % (idx, off_by, "\n".join([seqObj.id_length() for seqObj in loci])))

            self.total_length += loci[0].length
            self.total_sequences += 1

            has_start_fraction = 0.0
            has_stop_fraction = 0.0

            phylip_output.append('\n\t%s\t%s' % (len(self.allele_ids), loci[0].length))
            for seqObj in loci:
                try:
                    self.codon_counter_by_allele_id[seqObj.allele_id] += seqObj.codon_counter
                except KeyError:
                    self.codon_counter_by_allele_id[seqObj.allele_id] = seqObj.codon_counter

                if seqObj.has_start():
                    has_start_fraction += 1 / len(self.allele_ids)
                if seqObj.has_stop():
                    has_stop_fraction += 1 / len(self.allele_ids)

                phylip_output.append(seqObj.get_phylip())

            loci_output.append("%s\t%s\t%s\t%s" % (loci[0].header, loci[0].length, has_start_fraction, has_stop_fraction))

        with open(phylip_outfile, 'w') as fh:
            fh.write("%s" % "\n".join(phylip_output))

        with open(loci_outfile, 'w') as fh:
            fh.write("%s" % "\n".join(loci_output))

        for allele_id in self.allele_ids:
            codon_line = []
            codon_line.append(allele_id)
            codon_line.append(self.total_sequences)
            codon_line.append(self.total_length / 3)
            for codon in CODONS:
                codon_line.append(self.codon_counter_by_allele_id[allele_id][codon])
            codon_output.append("\t".join(str(x) for x in codon_line))

        with open(codon_outfile, 'w') as fh:
            fh.write("%s" % "\n".join(codon_output))


if __name__ == "__main__":

    args = docopt(__doc__)
    # print args
    SAMPLE_IDX = int(args['--sample'])
    ALLELE_IDX = int(args['--allele'])
    PREFIX = args['--out_prefix']
    DELIMITER = args['--delimiter']
    FORCE = args['--force']
    FASTAS = args['FASTAS']
    CODONS = ['ATA','ATC','ATT','ATG','ACA','ACC','ACG','ACT','AAC','AAT','AAA','AAG','AGC','AGT','AGA','AGG','CTA','CTC','CTG','CTT','CCA','CCC','CCG','CCT','CAC','CAT','CAA','CAG','CGA','CGC','CGG','CGT','GTA','GTC','GTG','GTT','GCA','GCC','GCG','GCT','GAC','GAT','GAA','GAG','GGA','GGC','GGG','GGT','TCA','TCC','TCG','TCT','TTC','TTT','TTA','TTG','TAC','TAT','TAA','TAG','TGC','TGT','TGA','TGG']
    MainObj()

