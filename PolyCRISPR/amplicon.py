from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

class Amplicon():

	def __init__(self, sequence,amplicon_reads):
		self.count = 0
		self.amplicon_reads = amplicon_reads
		self.sequence = sequence
		self.forward_dict = amplicon_reads.forward
		self.reverse_dict = amplicon_reads.reverse
		self.guides       = amplicon_reads.guides
		self.orientation  = "."
		self.reorient()

	def oriented_sequence(self):
		seq = self.sequence
		if(self.orientation == "-"):
			seq = self.sequence.reverse_complement()
		return(seq)

	def find_hash_position(self, tags):
		str_seq = str(self.oriented_sequence())
		for key, value in tags.items():
			position = self.sequence.find(value)
			if position != -1:
				return (key, position)
		return(None, -1)

	def find_primer_positions(self):
		self.fw_name, self.fw_pos = self.find_hash_position(self.forward_dict)
		self.rv_name, self.rv_pos = self.find_hash_position(self.reverse_dict)
		return(self.fw_pos != -1 or self.rv_pos != -1 )

	def reorient(self):
		#This methods sets the orientation to where both primers are found.
		found = self.find_primer_positions()
		if not found:
			self.orientation = "-"
			found = self.find_primer_positions()
		if not found:
			self.orientation = "."
			self.find_primer_positions()

	def __str__(self):
		guides = ""
		if self.found_guides:
			guides = ",".join(self.found_guides)
		return ("\t".join([ 
			str(self.count),
			str(self.fw_name), 
			str(self.rv_name),
			guides,
			str(len(str(self.oriented_sequence()))),
			str(self.oriented_sequence())], 
			)
		)

	def find_guides(self):
		str_seq = str(self.oriented_sequence)
		self.found_guides = []
		for key, value in self.guides.items():
			position = self.sequence.find(value)
			if position != -1:
				self.found_guides.append(key)

	def best_reference(self):
		aligner = self.amplicon_reads.aligner
		print("...")
		print(str(self.oriented_sequence()))
		best = 0
		best_aln = None
		for ref, seq in self.amplicon_reads.references.items():
			alignments = aligner.align(str(seq.seq),
				str( self.oriented_sequence()))
			print(ref)
			i = 0
			for aln in alignments:
				i += 1
				if(best < aln.score):
					best_aln = aln 
					best = aln.score
				if(i > 10):
					break

		print(best)
		if(best > 500):
			print(best_aln)




