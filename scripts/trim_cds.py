#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_new_name(old_file_name):
	head, tail = os.path.split(old_file_name)
	new_file_name1 = tail.split(".")[0] + "_reformed_cds.fa"
	return new_file_name1

def change_gen_len(fasta_file):
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		new_file_name = get_new_name(fasta_file)
		if len(seq_record.seq)%3 == 0:
			new_seq = seq_record.seq
			reform_rec = SeqRecord(Seq(new_seq),id = seq_record.id, description = "")
			with open(new_file_name, "a") as handle:
				SeqIO.write(reform_rec, handle, "fasta")
		elif len(seq_record.seq)%3 == 1:
			new_seq = seq_record.seq[1:]
			reform_rec = SeqRecord(Seq(new_seq),id = seq_record.id, description = "")
			with open(new_file_name, "a") as handle:
				SeqIO.write(reform_rec, handle, "fasta")
		else:
			new_seq = seq_record.seq[2:]
			reform_rec = SeqRecord(Seq(new_seq),id = seq_record.id, description = "")
			with open(new_file_name, "a") as handle:
				SeqIO.write(reform_rec, handle, "fasta")

folder_list = os.listdir("./")
for fasta_file in folder_list:
	change_gen_len(fasta_file)
