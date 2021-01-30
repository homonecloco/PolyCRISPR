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


def get_sequence_counts(fq, forward_dict, reverse_dict):
	fw_f = fq
	ret = defaultdict(int)
	with gzip.open(fw_f, "rt") as r1:
		for fw in SeqIO.parse(r1, "fastq") :
			str_seq = str(fw.seq)
			fw_pos, fw_name, rv_pos, rv_pos = find_primer_positions(str_seq, forward_dict, reverse_dict)
			orientation = "+"
			if(fw_pos == -1 or rv_pos == -1):
				str_seq = str(fw.seq.reverse_complement())
				fw_pos, fw_name, rv_pos, rv_name = find_primer_positions(str_seq, forward_dict, reverse_dict)
				orientation = "-"
			if(fw_pos == -1 and rv_pos == -1):
				orientation = "-"
			else:
				ret[str_seq] += 1
	return ret

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

def write_summary(output_prefix):
	distances_f = output_prefix + "_distances_merge.txt"
	df = pd.read_csv(distances_f)
	summ=df[['fw_name', 'rev_name', 'pl_name', 'pl_pos']].groupby(['fw_name', 'rev_name', 'pl_name']).agg( ['count','mean'])
	summ.to_csv(output_prefix + "summary_merge.csv", sep=',')

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--guides", default=None,
		help="Guides to search in the amplicons")
	parser.add_argument("-f", "--forward", default=None,
		help="Forward primers in the amplicone")
	parser.add_argument("-r", "--reverse", default=None,
		help="Reverse primers in the amplicone")
	parser.add_argument("-s", "--sequences", default=None,
		help="FastQ file with the amplicons")
	parser.add_argument("-m", "--minimum_coverage", default=200,
		help="FastQ file with the amplicons", type=int)
	parser.add_argument("-o", "--ouptut", default=None, 
		help="Output file. If missing, the ouptut is sent to stdout")
	args = parser.parse_args()
	#print(args)
	return(args)

def main():
	args = parse_arguments()
	print("IN MAIN")
	ar = AmpliconReads(args.forward, args.reverse, args.guides)

	# sequences = get_sequence_counts(sequences_path, forward_dict, reverse_dict)
	# for k, v in sequences.items():
	# 	if v > minimum_coverage:
	# 		fw_pos, fw_name, rv_pos, rv_name = find_primer_positions(k, forward_dict, reverse_dict)
	# 		print(fw_pos, fw_name, rv_pos, rv_name )
	# 		print(v)
	# 		print(k)
		#print("\n")

	#write_raw_counts(sequences_path, output_prefix, forward_dict, reverse_dict, plamids_dict)
	#write_summary(output_prefix)

if __name__ == '__main__':
	main()