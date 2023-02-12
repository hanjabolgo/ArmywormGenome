#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
from Bio import SeqIO
dicttest = {}
bed_file = "genome.windows"
dictbed = {}
fasta_file = "Chr.fa"
def get_te_length(fasta_file):
	for key in dictbed.keys():
		Chr_id = dictbed[key][0]
		Chr_sta = int(dictbed[key][1])
		Chr_end = int(dictbed[key][2])
		for seq_record in SeqIO.parse(fasta_file, "fasta"):
			if seq_record.id == Chr_id:
				seq_slice = seq_record.seq[Chr_sta+1 : Chr_end+1]
				seq_sequence = ''.join(seq_slice)
				pattern = re.compile(r'[a-z]')
				seq_lowlettr = pattern.findall(seq_sequence)
				te_length = len(seq_lowlettr)
				te_percent = te_length/(Chr_end-Chr_sta)
				print("{}\t{}\t{}\t{}".format(Chr_id,Chr_sta,Chr_end,te_percent))
get_te_length(fasta_file)
