from Bio import SeqIO
from Bio.Seq import Seq

class Amplicon():

	def __init__(self, sequence, forward_dict, reverse_dict, guides):
		self.count = 0
		self.sequence = sequence
		self.forward_dict = forward_dict
		self.reverse_dict = reverse_dict
		self.guides       = guides
		self.orientation  = "."
		self.reorient()

	def oriented_sequence(self):
		seq = self.sequence
		if(self.orientation == "-"):
			seq = self.sequence.reverse_complement()
		return(seq)

	def find_hash_position(self, tags):
		for key, value in tags.items():
			str_seq = str(self.oriented_sequence())
			position = self.sequence.find(value)
			if position != -1:
				return (key, position)
		return(None, -1)

	def find_primer_positions(self):
		self.fw_name, self.fw_pos = self.find_hash_position(self.forward_dict)
		self.rv_name, self.rv_pos = self.find_hash_position(self.reverse_dict)

	def reorient(self):
		#This methods sets the orientation to where both primers are found.
		self.find_primer_positions()
		if(self.fw_pos == -1 and self.rv_pos == -1):
			self.orientation = "-"
			self.find_primer_positions()
		if(self.fw_pos == -1 and self.rv_pos == -1):
			self.orientation = "."
			self.find_primer_positions()

	def __str__(self):
		return ("\t".join([
			str(self.oriented_sequence()), 
			str(self.count),
			str(self.fw_name), 
			str(self.rv_name)])
		)



