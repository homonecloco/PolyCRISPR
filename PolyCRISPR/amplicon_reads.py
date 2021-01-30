from Bio import SeqIO
from Bio.Seq import Seq
from .amplicon import Amplicon

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

class AmpliconReads():
	def __init__(self, forward, reverse, guides):
		self.guides = fasta_to_dict(guides)
		self.reverse = fasta_to_dict(reverse)
		self.forward = fasta_to_dict(forward)

		print(self.guides)
		print(self.reverse)
		print(self.forward)


