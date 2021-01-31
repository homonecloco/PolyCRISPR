from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

import gzip
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
	def __init__(self, forward, reverse, guides, reference):
		self.guides = fasta_to_dict(guides)
		self.reverse = fasta_to_dict(reverse)
		self.forward = fasta_to_dict(forward)
		self.references = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
		self.aligner = Align.PairwiseAligner()

	def setup_amplicons(self):
		pass

	def amplicon_counts(self,fq):
		ret = dict()
		with gzip.open(fq, "rt") as r1:
			for fw in SeqIO.parse(r1, "fastq"):
				amplicon = Amplicon(fw.seq, self)
				if amplicon.orientation == ".":
					next
				str_seq  = str(amplicon.oriented_sequence())
				if str_seq not in ret:
					ret[str_seq] = amplicon
				ret[str_seq].count += 1
#		print(ret)		
		self._amplicon=ret
		return ret



