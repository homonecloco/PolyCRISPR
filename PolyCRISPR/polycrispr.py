import sys
import getopt
import csv
import gzip
import itertools
import pandas as pd
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_seq(seq):
	desc = seq.description
	arr = desc.split(" ")
	strand = arr[1]
	if(strand == "-"):
		seq.seq = seq.seq.reverse_complement()
	return seq

def fasta_to_dict(path):
	ret_dict = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
	ret_dict = {k:reverse_seq(v) for k, v in ret_dict.items()}
	ret_dict = {k:str(v.seq) for k, v in ret_dict.items()}
	return(ret_dict)

def find_hash_position(sequence, tags):
	for key, value in tags.items():
		position = sequence.find(value)
		if position != -1:
			return (key, position)
	return("none", -1)


def find_primer_positions(sequence, forward_dict, reverse_dict):
	fw_name, fw_pos = find_hash_position(sequence, forward_dict)
	rv_name, rv_pos = find_hash_position(sequence, reverse_dict)
	return([fw_pos, fw_name, rv_pos, rv_name])

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

def usage():
	print ("Usage: " + sys.argv[0]  + " --guides=<guides.fa> --primers=<primers.fa> --sequences=<sequences.fq.gz> --output_prefix=<file_prefix> --folder=<folder_with_fastqs>")


def main(argv):
	minimum_coverage = 200
	try: 
		opts, args = getopt.getopt(argv,'g:p:s:o:f:r:m:h', ['guides=','primers=','sequences=', 'output_prefix=','forward=','reverse=','minimum_coverage=','help'])
	except getopt.GetOptError:
		usage()
		sys.exit(2)

	if not opts: 
			print ('No options supplied')
			usage()

	for opt, arg in opts: 
		if opt in ('h', '--help'):
			usage()
			sys.exit(2)
		elif opt in ('g', '--guides'):
			guides_path = arg 
			guides_dict = fasta_to_dict(guides_path)
			print(guides_dict)
		elif opt in ('f', '--forward'):
			primers_path = arg 
			forward_dict = fasta_to_dict(primers_path)
			print(forward_dict)
		elif opt in ('r', '--reverse'):
			primers_path = arg 
			reverse_dict = fasta_to_dict(primers_path)
			print(reverse_dict)
		elif opt in ('s', '--sequences'):
			sequences_path = arg
			print(arg)
		elif opt in ('o', '--output_prefix'):
			output_prefix = arg
		elif opt in ('m', '--minimum_coverage'):
			minimum_coverage = imt(arg)

	sequences = get_sequence_counts(sequences_path, forward_dict, reverse_dict)
	for k, v in sequences.items():
		if v > minimum_coverage:
			fw_pos, fw_name, rv_pos, rv_name = find_primer_positions(k, forward_dict, reverse_dict)
			print(fw_pos, fw_name, rv_pos, rv_name )
			print(v)
			print(k)
		#print("\n")

	#write_raw_counts(sequences_path, output_prefix, forward_dict, reverse_dict, plamids_dict)
	#write_summary(output_prefix)

if __name__ == '__main__':
	main(sys.argv[1:])