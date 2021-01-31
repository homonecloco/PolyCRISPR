import os
import sys
import getopt
import csv
import gzip
import itertools
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from .amplicon import Amplicon
from .amplicon_reads import AmpliconReads

def write_raw_counts(fq, output_prefix, forward_dict, reverse_dict, guides):
	fw_f = fq
	distances_f = output_prefix + "_amplicons.txt"
	with gzip.open(fw_f, "rt") as r1,open(distances_f,"w") as f3:
		#f3.write("fw_name,fw_pos,rev_name,rev_pos,pl_name,pl_pos,orientation\n")
		#seqs = 
		for fw in SeqIO.parse(r1, "fastq") :
			str_seq = str(fw.seq)
			# values = find_plasmid_positions(str_seq, forward_dict, reverse_dict, plasmids, "+")	 
			# if values[6] == ".":
			# 	str_seq = str(fw.seq.reverse_complement())
			# 	values = find_plasmid_positions(str_seq, forward_dict, reverse_dict, plasmids, "-")
			# f3.write(",".join(values)) 
			# f3.write("\n")


def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--guides", default=None,
		help="Guides to search in the amplicons")
	parser.add_argument("-f", "--forward", default=None,
		help="Forward primers in the amplicons")
	parser.add_argument("-r", "--reverse", default=None,
		help="Reverse primers in the amplicons")
	parser.add_argument("-s", "--sequences", default=None,
		help="FastQ file with the amplicons")
	parser.add_argument("-m", "--minimum_coverage", default=200,
		help="FastQ file with the amplicons", type=int)
	parser.add_argument("-o", "--ouptut", default=None, 
		help="Output file. If missing, the ouptut is sent to stdout")
	parser.add_argument("-t", "--target_reference", default=None,
		help="Full sequence of the regions to amplify")
	args = parser.parse_args()
	#print(args)
	return(args)

def main():
	args = parse_arguments()
	ar = AmpliconReads(args.forward, args.reverse, args.guides, args.target_reference)
	counts = ar.amplicon_counts(args.sequences)
	basename = os.path.basename(args.sequences)
	for(k, v ) in counts.items():
		if v.count > args.minimum_coverage:
			#print(k)
			v.find_guides()
			print(basename + "\t" + str(v) )

			#v.best_reference()
	

if __name__ == '__main__':
	main()