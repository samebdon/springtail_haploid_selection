#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
        orthodiver.py -d <DIR> -A <STR> -B <STRL> -o <STR> [-v --debug -h]


Options:
        -h, --help                              Show this
        -v, --version                           Print version
        -d, --fasta_dir <FILE>                  Directory containing all FASTA alignment files 
        -A, --taxon_A <STR>                     Taxon A
        -B, --taxon_B <STR>                     Taxon B
        -o, --outprefix <STRING>                Outprefix to use for output
        --debug                                 Print debugging messages to log file
        
Requirements:
    - FASTA alignment file names must end in *.fasta, *.fas, or *.fa
    - FASTA alignment file names must be in format '\\w+.OGNUMBER.\\.+'
    - FASTA headers must be in format 'SPECIES.SAMPLE.HAPLOTYPE.TRANSCRIPT'
    - '--taxon_A' / '--taxon_B' must be 'SPECIES.SAMPLE's

Example command:
    ./orthodiver.py -d . -A gonepteryx_cleopatra.GC_12 -B gonepteryx_rhamni.GR_112 -o test
    
"""

from docopt import docopt
import collections
import sys
import os
from tqdm import tqdm
from timeit import default_timer as timer
import logging

ALN_EXTENSIONS = set(["fasta", "fas", "fa"])
SITES = ['hetA', 'hetB', 'fixed', 'hetAB', 'invariant']
PIS = ['pi_A', 'pi_B', 'pi_between', 'pi_total']
DEGENERACIES = [0, 2, 3, 4]

# codon2aa
gen_code_dict = {   
                'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L', 
                'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S', 
                'tat': 'Y', 'tac': 'Y', 'taa': 'X', 'tag': 'X', 
                'tgt': 'C', 'tgc': 'C', 'tga': 'X', 'tgg': 'W', 
                'ctt': 'L', 'ctc': 'L', 'cta': 'L', 'ctg': 'L', 
                'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', 
                'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q', 
                'cgt': 'R', 'cgc': 'R', 'cga': 'R', 'cgg': 'R', 
                'att': 'I', 'atc': 'I', 'ata': 'I', 'atg': 'M', 
                'act': 'T', 'acc': 'T', 'aca': 'T', 'acg': 'T', 
                'aat': 'N', 'aac': 'N', 'aaa': 'K', 'aag': 'K', 
                'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R', 
                'gtt': 'V', 'gtc': 'V', 'gta': 'V', 'gtg': 'V', 
                'gct': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A', 
                'gat': 'D', 'gac': 'D', 'gaa': 'E', 'gag': 'E', 
                'ggt': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G'
                } 

# codon2degen
degen_dict = {
        'ttt': '002', 'ttc': '002', 'tta': '202', 'ttg': '202', 
        'tct': '004', 'tcc': '004', 'tca': '004', 'tcg': '004', 
        'tat': '002', 'tac': '002', 'taa': '022', 'tag': '002', 
        'tgt': '002', 'tgc': '002', 'tga': '020', 'tgg': '000', 
        'ctt': '004', 'ctc': '004', 'cta': '204', 'ctg': '204', 
        'cct': '004', 'ccc': '004', 'cca': '004', 'ccg': '004', 
        'cat': '002', 'cac': '002', 'caa': '002', 'cag': '002', 
        'cgt': '004', 'cgc': '004', 'cga': '204', 'cgg': '204', 
        'att': '003', 'atc': '003', 'ata': '003', 'atg': '000', 
        'act': '004', 'acc': '004', 'aca': '004', 'acg': '004', 
        'aat': '002', 'aac': '002', 'aaa': '002', 'aag': '002', 
        'agt': '002', 'agc': '002', 'aga': '202', 'agg': '202', 
        'gtt': '004', 'gtc': '004', 'gta': '004', 'gtg': '004', 
        'gct': '004', 'gcc': '004', 'gca': '004', 'gcg': '004', 
        'gat': '002', 'gac': '002', 'gaa': '002', 'gag': '002', 
        'ggt': '004', 'ggc': '004', 'gga': '004', 'ggg': '004'
        }

def compute_pi(sites):
    '''
    SITES = ['hetA', 'hetB', 'fixed', 'hetAB', 'invariant']
    PIS = ['pi_A', 'pi_B', 'pi_between', 'pi_total']
    pi_A = (hetA + hetAB) / total_sites
    pi_B = (hetB + hetAB) / total_sites
    pi_between =  (((hetA + hetB + hetAB) / 2) + fixed) / total_sites
    pi_total = ((4.0 / 6) * pi_between) + ((1.0 / 6) * pi_A) + ((1.0 / 6) * pi_B))
    '''
    pi = {}
    for degeneracy in DEGENERACIES:
        site_keys = ["%s_%s" % (degeneracy, site) for site in SITES] 
        pi_keys = ["%s_%s" % (degeneracy, pi) for pi in PIS]
        total_sites = sum([sites[site_key] for site_key in site_keys])
        if total_sites:
            pi[pi_keys[0]] = float("%.8f" % ((sites[site_keys[0]] + sites[site_keys[3]]) / total_sites))
            pi[pi_keys[1]] = float("%.8f" % ((sites[site_keys[1]] + sites[site_keys[3]]) / total_sites))
            pi[pi_keys[2]] = float("%.8f" % ((((sites[site_keys[0]] + sites[site_keys[1]] + sites[site_keys[3]]) / 2) + sites[site_keys[2]]) / total_sites))
            pi[pi_keys[3]] = float("%.8f" % (((4.0 / 6) * pi[pi_keys[2]]) + ((1.0 / 6) * pi[pi_keys[0]]) + ((1.0 / 6) * pi[pi_keys[1]])))
        else:
            for metric_key in pi_keys:
                pi[metric_key] = 'NA' 
    return pi

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")

def get_degeneracies(codons, variant_pos):
    return [int(degen_dict.get(codon)[variant_pos]) for codon in codons]

def get_genotypes(codons, variant_pos):
    return [codon[variant_pos] for codon in codons]

def get_variant_pos(codons):
    variant_pos = []
    for codon in codons:
        for idx, base in enumerate(codon):
            if not base == codons[0][idx]:
                variant_pos.append(idx)
    return list(set(variant_pos))

def get_variation_list(codons):
    variant_pos = set()
    for codon in codons:
        for idx, base in enumerate(codon):
            if not base == codons[0][idx]:
                variant_pos.add(idx)
    variation_list = ["".join([("*" if idx in variant_pos else "-") for idx in range(0,3)])] * 4
    return variation_list

def get_zygosity(geno1, geno2):
    if geno1 == geno2:
        return 'HOM'
    else:
        return 'HET'

def get_codons_and_sites(aligned_codons):
    aligned_codons_set = set(aligned_codons)
    codon_chars = "".join(aligned_codons_set)
    codon_type = None # valid_invariant/valid_variant/multidegenerate/multivariant/gap/N
    sites = collections.Counter() # degeneracy + invariant / degeneracy + mutype
    variation_list = []
    if '-' in codon_chars or 'n' in codon_chars:
        codon_type = 'unknown'
        variation_list = ['???' for codon in aligned_codons]
    elif len(aligned_codons_set) == 1:
        codon_type = 'valid_invariant'
        for idx in range(0, 3):
            degeneracies = get_degeneracies(aligned_codons, idx)
            site = "%s_invariant" % degeneracies.pop()
            sites[site] += 1
        variation_list = ['---' for codon in aligned_codons]
    else:
        variation_list = get_variation_list(aligned_codons)
        variant_positions = get_variant_pos(aligned_codons)
        if len(variant_positions) > 1:
            codon_type = 'multivariant'
        else:
            variant_position = variant_positions[0]
            degeneracies = get_degeneracies(aligned_codons, variant_position)
            if not len(set(degeneracies)) == 1:
                codon_type = 'multidegenerate'
            else:
                genotypes = get_genotypes(aligned_codons, variant_position)
                if len(set(genotypes)) == 2:
                    codon_type = 'valid_variant'
                    # what about non-biallelic SNPs ?!!!!!!!!
                    degeneracy = degeneracies[0]
                    A1, A2, B1, B2 = genotypes
                    A_zygosity = get_zygosity(A1, A2)
                    B_zygosity = get_zygosity(B1, B2)
                    if A_zygosity == 'HOM' and B_zygosity == 'HOM':
                        site = "%s_fixed" % degeneracy
                    elif A_zygosity == 'HET' and B_zygosity == 'HET':
                        site = "%s_hetAB" % degeneracy
                    elif A_zygosity == 'HET' and B_zygosity == 'HOM':
                        site = "%s_hetA" % degeneracy
                    elif A_zygosity == 'HOM' and B_zygosity == 'HET':
                        site = "%s_hetB" % degeneracy
                    else:
                        sys.exit("[X] A1=%s A2=%s B1=%s B2=%s [???]" % (A1, A2, B1, B2, degeneracy))
                    sites[site] += 1
                    # deal with invariant sites inside of valid_variant codon ...
                    for idx in range(0, 3):
                        if not idx == variant_position:
                            idx_degeneracies = get_degeneracies(aligned_codons, idx)
                            if len(set(idx_degeneracies)) == 1: # monodegenerate
                                site = "%s_invariant" % idx_degeneracies.pop()
                                sites[site] += 1
                else:
                    codon_type = 'multiallelic'
                    #sys.exit("[X] genotypes=%s" % genotypes)
    _degeneracies = [str(degen_dict.get(codon, '???')) for codon in aligned_codons]
    logging.debug("[#] %s A = [%s] vs B = [%s]" % (
        "Codons     :",
        "|".join([aligned_codons[0], aligned_codons[1]]),
        "|".join([aligned_codons[2], aligned_codons[3]])
        ))
    logging.debug("[+] %s A = [%s] vs B = [%s]" % ("Degeneracy :", "|".join([_degeneracies[0], _degeneracies[1]]), "|".join([_degeneracies[2], _degeneracies[3]]),))
    logging.debug("[+] %s A = [%s] vs B = [%s]" % ("Variation  :", "|".join([variation_list[0], variation_list[1]]), "|".join([variation_list[2], variation_list[3]]),))
    logging.debug("[+] -> codon_type = %r" % codon_type)
    logging.debug("[+] -> sites = %r" % dict(sites))
    return codon_type, sites

def get_codons_header(locus=False):
    if locus:
        col = '# locus_id'
    else:
        col = '# dataset_id'
    return "\t".join([
            col, 
            'length', 
            'codons', 
            'unknown',
            'multivariant', 
            'multidegenerate', 
            'valid_invariant', 
            'valid_variant'
            ])

def get_sites_header(degeneracy, taxon_by_label, locus=False):
    if locus:
        col = '# locus_id'
    else:
        col = '# dataset_id'
    return "\t".join([
            col, 
            '%s_total' % degeneracy, 
            '%s_invariant' % degeneracy, 
            '%s_hetA (%s)' % (degeneracy, taxon_by_label['A']), 
            '%s_hetB (%s)' % (degeneracy, taxon_by_label['B']), 
            '%s_fixed' % degeneracy, 
            '%s_hetAB' % degeneracy
            ])

def get_pi_header(degeneracy, taxon_by_label, locus=False):
    if locus:
        col = '# locus_id'
    else:
        col = '# dataset_id'
    return "\t".join([
        col, 
        '%s_pi_A (%s)' % (degeneracy, taxon_by_label['A']),
        '%s_pi_B (%s)' % (degeneracy, taxon_by_label['B']),
        '%s_pi_between' % degeneracy,
        '%s_pi_total' % degeneracy,
        ])

class SeqObj(object):
    def __init__(self, label, taxon, hap_id, transcript_id, sequence):
        self.label = label
        self.taxon = taxon
        self.hap_id = hap_id
        self.transcript_id = transcript_id
        self.codons = [sequence[i:i+3].lower() for i in range(0, len(sequence), 3)]
        self.length = len(sequence)
        self.header = "%s.%s.%s" % (self.taxon, self.hap_id, self.transcript_id)

class AlnObj(object):
    '''
    Contains information from one aligned locus
    - 4 FASTA sequences
    '''
    def __init__(self, locus_id):
        self.locus_id = locus_id # OG
        self.lengths = []
        self.codons_by_label = collections.defaultdict(list) # labels A/B
        self.codon_type_counter = collections.Counter()
        self.sites_counter = collections.Counter()
        self.pi_metrics = {}

    def get_pi_line(self, degeneracy):
        return "\t".join([
            self.locus_id, 
            str(self.pi_metrics['%s_pi_A' % degeneracy]),
            str(self.pi_metrics['%s_pi_B' % degeneracy]),
            str(self.pi_metrics['%s_pi_between' % degeneracy]),
            str(self.pi_metrics['%s_pi_total' % degeneracy])
            ])

    def get_sites_line(self, degeneracy):
        degeneracy_sites = ['%s_%s' % (degeneracy, site) for site in SITES]
        return "\t".join([
            self.locus_id, 
            str(sum([self.sites_counter[degeneracy_site] for degeneracy_site in degeneracy_sites])), 
            str(self.sites_counter['%s_invariant' % degeneracy]), 
            str(self.sites_counter['%s_hetA' % degeneracy]),
            str(self.sites_counter['%s_hetB' % degeneracy]),
            str(self.sites_counter['%s_fixed' % degeneracy]), 
            str(self.sites_counter['%s_hetAB' % degeneracy])
            ])

    def get_codon_line(self):
        estimated_codons = int(self.length())
        observed_codons = sum(self.codon_type_counter.values())
        if not estimated_codons == observed_codons * 3:
            sys.exit("estimated_codons [%s] != observed_codons [%s]" % (estimated_codons, observed_codons))
        return "\t".join([
            self.locus_id,
            str(estimated_codons), 
            str(observed_codons), 
            str(self.codon_type_counter['unknown']), 
            str(self.codon_type_counter['multivariant']), 
            str(self.codon_type_counter['multidegenerate']), 
            str(self.codon_type_counter['valid_invariant']), 
            str(self.codon_type_counter['valid_variant'])
            ])

    def add_seqObj(self, seqObj):
        if not seqObj.length % 3 == 0:
            sys.exit("[X] Sequence %s is not divisible by 3 (length=%s)..." % (seqObj.header, seqObj.length))
        self.codons_by_label[seqObj.label].append(seqObj.codons)
        self.lengths.append(seqObj.length)

    def yield_label_codon_list(self):
        for label, codon_list in self.codons_by_label.items():
            yield label, codon_list

    def yield_label_codon_counter(self):
        if len(set(self.lengths)) == 1: # sanity check for equal length of chars
            for label, codon_list in self.yield_label_codon_list():
                yield label, sum([collections.Counter(codons) for codons in codon_list], collections.Counter())

    def length(self):
        return self.lengths[0] # simple, since sanity check performed earlier
    
    def compute_locus_pi(self):
        codon_list_by_label = {}
        logging.debug("[#] Locus : %s" % self.locus_id)
        for label, codon_list in self.yield_label_codon_list():
            codon_list_by_label[label] = codon_list
        logging.debug("[#] A_1 : %s" % codon_list_by_label['A'][0])
        logging.debug("[#] A_2 : %s" % codon_list_by_label['A'][1])
        logging.debug("[#] B_1 : %s" % codon_list_by_label['B'][0])
        logging.debug("[#] B_2 : %s" % codon_list_by_label['B'][1])
        aligned_codons_list = [list(codons_A) + list(codons_B) for codons_A, codons_B in list(zip(zip(*codon_list_by_label['A']), zip(*codon_list_by_label['B'])))]
        for aligned_codons in aligned_codons_list:
            codon_type, sites = get_codons_and_sites(aligned_codons)
            self.codon_type_counter[codon_type] += 1
            self.sites_counter += sites
        self.pi_metrics = compute_pi(self.sites_counter)

class DataObj(object):
    def __init__(self, args):
        # input parameters
        self.fasta_dir = args['--fasta_dir']
        self.outprefix = args['--outprefix']
        
        self.label_by_taxon = {
            args['--taxon_A']: 'A',
            args['--taxon_B']: 'B',
        }
        self.taxon_by_label = {
            'A': args['--taxon_A'],
            'B': args['--taxon_B']
        }
        self.dataset_label = "%s_vs_%s" % (args['--taxon_A'], args['--taxon_B'])
        logging.info("[#] Starting analysis ...")
        for taxon, label in self.taxon_by_label.items():
            logging.info("[+]\tTaxon %r = %r" % (taxon, label))
        # data
        self.alnObjs = []
        self.sites_counter = collections.Counter()
        self.codon_type_counter = collections.Counter()
        self.pi_metrics = {}
        self.codon_counter_by_label = collections.defaultdict(collections.Counter) # codon usage?

    def get_codon_summary(self):
        length = sum([alnObj.length() for alnObj in self.alnObjs])
        observed_codons = sum(self.codon_type_counter.values())
        if not length == observed_codons * 3:
            sys.exit("data set length [%s] != observed_codons * 3 [%s]" % (length, observed_codons))
        return "\t".join([
            self.dataset_label,
            str(length), 
            str(observed_codons), 
            str(self.codon_type_counter['unknown']), 
            str(self.codon_type_counter['multivariant']), 
            str(self.codon_type_counter['multidegenerate']), 
            str(self.codon_type_counter['valid_invariant']), 
            str(self.codon_type_counter['valid_variant'])
            ])

    def get_sites_summary(self, degeneracy):
        degeneracy_sites = ['%s_%s' % (degeneracy, site) for site in SITES]
        return "\t".join([
            self.dataset_label, 
            str(sum([self.sites_counter[degeneracy_site] for degeneracy_site in degeneracy_sites])), 
            str(self.sites_counter['%s_invariant' % degeneracy]), 
            str(self.sites_counter['%s_hetA' % degeneracy]),
            str(self.sites_counter['%s_hetB' % degeneracy]),
            str(self.sites_counter['%s_fixed' % degeneracy]), 
            str(self.sites_counter['%s_hetAB' % degeneracy])
            ])

    def get_pi_summary(self, degeneracy):
        return "\t".join([
            self.dataset_label, 
            str(self.pi_metrics['%s_pi_A' % degeneracy]),
            str(self.pi_metrics['%s_pi_B' % degeneracy]),
            str(self.pi_metrics['%s_pi_between' % degeneracy]),
            str(self.pi_metrics['%s_pi_total' % degeneracy])
            ])

    def parse_alignments(self):
        fasta_fs = []
        logging.info("[#] Parsing directory %s ..." % os.path.abspath(self.fasta_dir))
        for f in os.listdir(self.fasta_dir):
            if f.split(".")[-1] in ALN_EXTENSIONS:
                fasta_fs.append(os.path.join(self.fasta_dir, f))
        logging.info("[+]\tFound %s FASTA alignment files ..." % len(fasta_fs))
        logging.info("[#] Parsing FASTA alignment files ...")
        for fasta_f in tqdm(fasta_fs, total=len(fasta_fs), desc="[%] ", ncols=80):
            locus_id = os.path.basename(fasta_f).split(".")[1]
            alnObj = AlnObj(locus_id)
            with open(fasta_f) as fasta_fh:
                data = fasta_fh.read().rstrip("\n").split("\n")
                for idx in range(0, len(data), 2):
                    species_id, sample_id, hap_id, transcript_id = data[idx][1:].split(".")
                    taxon = "%s.%s" % (species_id, sample_id)
                    if not taxon in self.label_by_taxon:
                        sys.exit("[X] Unknown species %s in file %s" % (taxon, fasta_f))
                    label = self.label_by_taxon[taxon]
                    sequence = data[idx + 1]
                    seqObj = SeqObj(label, taxon, hap_id, transcript_id, sequence)
                    alnObj.add_seqObj(seqObj)
            self.alnObjs.append(alnObj)
            for label, codon_counter in alnObj.yield_label_codon_counter():
                self.codon_counter_by_label[label] += codon_counter

    def analyse(self):
        logging.info("[#] Analysing alignments ...")
        for alnObj in tqdm(self.alnObjs, total=len(self.alnObjs), desc="[%] ", ncols=80):
            alnObj.compute_locus_pi()
            self.sites_counter += alnObj.sites_counter
            self.codon_type_counter += alnObj.codon_type_counter
        self.pi_metrics = compute_pi(self.sites_counter)

    def write_results(self):
        # codons_summary
        logging.info("[#] Writing results ...")
        codons_summary_f = '%s.codons_summary.txt' % self.outprefix
        codons_summary_rows = [get_codons_header()]
        codons_summary_rows.append(self.get_codon_summary())
        with open(codons_summary_f, 'w') as codons_summary_fh:
            codons_summary_fh.write("\n".join(codons_summary_rows) + "\n")
        logging.info("[+]\tWritten %s ..." % codons_summary_f)
        # codons_by_locus
        codons_by_locus_f = '%s.codons_by_locus.txt' % self.outprefix
        codons_by_locus_rows = [get_codons_header(locus=True)]
        for alnObj in self.alnObjs:
            codons_by_locus_rows.append(alnObj.get_codon_line())
        with open(codons_by_locus_f, 'w') as codons_by_locus_fh:
            codons_by_locus_fh.write("\n".join(codons_by_locus_rows) + "\n")
        logging.info("[+]\tWritten %s ..." % codons_by_locus_f)
        for degeneracy in DEGENERACIES:
            # sites_summary
            sites_summary_f = '%s.%sd_sites_summary.txt' % (self.outprefix, degeneracy)
            sites_summary_rows = [get_sites_header(degeneracy, self.taxon_by_label)]
            sites_summary_rows.append(self.get_sites_summary(degeneracy))
            with open(sites_summary_f, 'w') as sites_summary_fh:
                sites_summary_fh.write("\n".join(sites_summary_rows) + "\n")
            logging.info("[+]\tWritten %s ..." % sites_summary_f) 
            # pi_summary
            pi_summary_f = '%s.%sd_pi_summary.txt' % (self.outprefix, degeneracy)
            pi_summary_rows = [get_pi_header(degeneracy, self.taxon_by_label)]
            pi_summary_rows.append(self.get_pi_summary(degeneracy))
            with open(pi_summary_f, 'w') as pi_summary_fh:
                pi_summary_fh.write("\n".join(pi_summary_rows) + "\n")
            logging.info("[+]\tWritten %s ..." % pi_summary_f) 
            # pi_by_locus / sites_by_locus
            pi_by_locus_f = '%s.%sd_pi_by_locus.txt' % (self.outprefix, degeneracy)
            pi_by_locus_rows = [get_pi_header(degeneracy, self.taxon_by_label, locus=True)]
            sites_by_locus_f = '%s.%sd_sites_by_locus.txt' % (self.outprefix, degeneracy)
            sites_by_locus_rows = [get_sites_header(degeneracy, self.taxon_by_label, locus=True)]
            for alnObj in self.alnObjs:
                sites_by_locus_rows.append(alnObj.get_sites_line(degeneracy))
                pi_by_locus_rows.append(alnObj.get_pi_line(degeneracy))
            with open(pi_by_locus_f, 'w') as pi_by_locus_fh:
                pi_by_locus_fh.write("\n".join(pi_by_locus_rows) + "\n")
            logging.info("[+]\tWritten %s ..." % pi_by_locus_f)                
            with open(sites_by_locus_f, 'w') as sites_by_locus_fh:
                sites_by_locus_fh.write("\n".join(sites_by_locus_rows) + "\n")
            logging.info("[+]\tWritten %s ..." % sites_by_locus_f)                

if __name__ == "__main__":
    __version__ = 0.3
    args = docopt(__doc__)
    if args['--version']:
        sys.exit("[V] orthodiver.py v%s" % __version__)
    try:
        start_time = timer()
        # logging
        if args['--debug']:
            logging.basicConfig(
                level=logging.DEBUG,
                format='%(asctime)s [%(levelname)-6s] %(message)s',
                datefmt='%m/%d/%Y %H:%M:%S',
                filename="%s.orthodiver.log" % args['--outprefix'],
                filemode='w')
        else:
            logging.basicConfig(
                level=logging.INFO,
                format='%(message)s',
                filename="%s.orthodiver.log" % args['--outprefix'],
                filemode='w')
        console = logging.StreamHandler() # define a Handler which writes INFO messages or higher to the sys.stderr
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter('%(message)s')) # set a format which is simpler for console use
        logging.getLogger('').addHandler(console) # add the handler to the root logger
        
        args = docopt(__doc__)
        dataObj = DataObj(args)
        dataObj.parse_alignments()
        dataObj.analyse()
        dataObj.write_results()
        logging.info("[+] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)

