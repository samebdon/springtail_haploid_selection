#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: minimap2synteny.py -a <STR> -b <STR> -i <STR> [-l <STR> -m <INT> -n <INT> -c <INT> -t <FLT> -h]

  [Options]
    -a, --genomefileA <STR>                     Genomefile for taxon A
    -b, --genomefileB <STR>                     Genomefile for taxon B
    -i, --liftover <STR>                        Liftover file (seq_A, start_A, end_A, seq_B, star_B, end_B, seq_anc)
    -l, --labels <STR>                          Whether to plot labels, choose from True or False [default: True]
    -m, --gapA <INT>                            Gap between genome A chromosomes [default: 10_000_000]
    -n, --gapB <INT>                            Gap between genome B chromosomes [default: 10_000_000]
    -c, --chromosome_width <INT>                Chromosome width [default: 6]
    -t, --alpha <FLT>                           Alpha of alignments [default: 0.1]
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import collections
from scipy.interpolate import make_interp_spline, BSpline
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def generate_genomefile_dict(genomefile, offset, colour):
	genomefile_dict = {}
	cumulative_genome = offset
	with open(genomefile, "r") as fin:
		# for each chromosome, record cumulative coordinates, orientation, and label
		for i, line in enumerate(fin):
			line = line.rstrip()
			chromosome, chromosome_length, orientation, label = line.split("\t")
			chromosome_length = int(chromosome_length)
			genomefile_dict[chromosome] = [cumulative_genome, cumulative_genome + chromosome_length, orientation, label]
			cumulative_genome += chromosome_length
			cumulative_genome += offset
			if colour:
				genomefile_dict[chromosome].append(i)
	return genomefile_dict

# plots chromosomes as lines and adds labels if arg is True
def plot_chromosomes(genomefile_dict, y_coord, labels, chromosome_width):
	for chromosome in genomefile_dict:
		plt.plot([genomefile_dict[chromosome][0], genomefile_dict[chromosome][1]], 
			[y_coord, y_coord], color='slategrey', alpha=1, linewidth=chromosome_width)
		middle_of_chromosome = genomefile_dict[chromosome][0] + \
		((genomefile_dict[chromosome][1] - genomefile_dict[chromosome][0]) / 2)
		plt.plot([genomefile_dict[chromosome][0], genomefile_dict[chromosome][1]], 
			[y_coord*1.02, y_coord*1.02], color='slategrey', alpha=1, linewidth=chromosome_width)
		middle_of_chromosome = genomefile_dict[chromosome][0] + \
		((genomefile_dict[chromosome][1] - genomefile_dict[chromosome][0]) / 2)
		if labels == "True":
			plt.text(middle_of_chromosome, y_coord*1.2, genomefile_dict[chromosome][3], ha='center', va='center', 
				wrap=True, fontsize=20)

def generate_alignment_dicts(liftover, genomefile_A_dict, genomefile_B_dict, alignment_coords):
	with open(liftover, "r") as fin:
		# for each alignment, record coordinates in either genome
		for line in fin:
			alignment = []
			line = line.rstrip()
			seqA, startA, endA, seqB, startB, endB, seqanc = line.split("\t")
			midpointA = int(startA) + ((int(endA) - int(startA)) / 2)
			midpointB = int(startB) + ((int(endB) - int(startB)) / 2)
			average_width = ((int(endA) - int(startA)) + (int(endB) - int(startB))) / 2

			if seqA in genomefile_A_dict.keys():
				# flip coords if orientation is -
				if genomefile_A_dict[seqA][2] == '-':
					midpointA = genomefile_A_dict[seqA][1] - midpointA
				if genomefile_A_dict[seqA][2] == '+':
					midpointA += genomefile_A_dict[seqA][0]
				alignment.append(midpointA)
				
			if seqB in genomefile_B_dict.keys():
				# flip coords if orientation is -
				if genomefile_B_dict[seqB][2] == '-':
					midpointB = genomefile_B_dict[seqB][1] - midpointB
				if genomefile_B_dict[seqB][2] == '+':
					midpointB += genomefile_B_dict[seqB][0]
				alignment.append(midpointB)
				if int(seqanc) == 0:
					line_colour = "lightgrey"
				else:
					line_colour = cm.tab20((int(seqanc) -1) / 7)

			# only interested in alignments on sequences in both genomefiles
			if len(alignment) == 2:
				alignment.append(average_width)
				alignment.append(line_colour)
				alignment_coords.append(alignment)


	return alignment_coords

# code from https://stackoverflow.com/questions/19394505/expand-the-line-with-specified-width-in-data-unit/42972469#42972469
def linewidth_from_data_units(linewidth, axis, reference='x'):
    """
    Convert a linewidth in data units to linewidth in points.

    Parameters
    ----------
    linewidth: float
        Linewidth in data units of the respective reference-axis
    axis: matplotlib axis
        The axis which is used to extract the relevant transformation
        data (data limits and size must not change afterwards)
    reference: string
        The axis that is taken as a reference for the data width.
        Possible values: 'x' and 'y'. Defaults to 'y'.

    Returns
    -------
    linewidth: float
        Linewidth in points
    """
    fig = axis.get_figure()
    if reference == 'x':
        length = fig.bbox_inches.width * axis.get_position().width
        value_range = np.diff(axis.get_xlim())
    elif reference == 'y':
        length = fig.bbox_inches.height * axis.get_position().height
        value_range = np.diff(axis.get_ylim())
    # Convert length to points
    length *= 72
    # Scale linewidth to value range
    return linewidth * (length / value_range)


args = docopt(__doc__)

# generate dicts for each genome with cumulative coordinates
genomefile_A_dict = generate_genomefile_dict(args['--genomefileA'], int(args['--gapA']), colour=False)
genomefile_B_dict = generate_genomefile_dict(args['--genomefileB'], int(args['--gapB']), colour=True)

# set up plot
fig = plt.figure(figsize=(28,4), frameon=False)
ax = fig.add_subplot(111)
ax.axis('off')
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)

# each alignment has coordinates to be recorded
alignment_coords = []
alignment_coords = generate_alignment_dicts(args['--liftover'], genomefile_A_dict, genomefile_B_dict, alignment_coords)

# for each alignment, if the two entries are from the two genomes, plot a curve between them
for alignment in alignment_coords:
	# code is adapted from https://gist.github.com/andrewgiessel/5684769
	upper_x = alignment[0]
	lower_x = alignment[1]
	average_width = alignment[2]
	alignment_colour = alignment[3]
	Y = np.array([-1, -0.568, -0.32, -0.16, -0.056, 0, 0.056, 0.16, 0.32, 0.568, 1])
	X = np.array([lower_x, lower_x + (0.1*(upper_x - lower_x)), 
	lower_x + (0.2*(upper_x - lower_x)), lower_x + (0.3*(upper_x - lower_x)), 
	lower_x + (0.4*(upper_x - lower_x)), lower_x + (0.5*(upper_x - lower_x)), 
	lower_x + (0.6*(upper_x - lower_x)), lower_x + (0.7*(upper_x - lower_x)), 
	lower_x + (0.8*(upper_x - lower_x)), lower_x + (0.9*(upper_x - lower_x)), upper_x])
	# requires sorted arrays, so flip if in the wrong order
	if lower_x > upper_x:
		X = np.flip(X)
		Y = np.flip(Y)
	xnew = np.linspace(X.min(), X.max(), 300) 
	spl = make_interp_spline(X, Y, k=3)
	power_smooth = spl(xnew)
	plt.plot(xnew, power_smooth, color=alignment_colour, alpha=float(args['--alpha']), 
		linewidth=linewidth_from_data_units(average_width, ax, reference='x'))

# plot the chromosomes
plot_chromosomes(genomefile_A_dict, 1, args["--labels"], int(args['--chromosome_width']))
plot_chromosomes(genomefile_B_dict, -1, args["--labels"], int(args['--chromosome_width']))

plt.savefig("minimap2synteny.pdf", format="pdf", bbox_inches='tight')

"""
# example command
~/software/minimap2synteny.py -a brenthis_ino.SP_BI_364.v2_0.sequences.genomefile -b brenthis_daphne.ES_BD_1141.v2_0.sequences.genomefile -i brenthis_daphne.ES_BD_1141.v2_0.sequences.vs.brenthis_ino.SP_BI_364.v2_0.sequences.liftover.txt -l True -m 10_000_000 -n 10_600_000
"""
