#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

##get new file name 
def get_new_name(old_file_name):
	head, tail = os.path.split(old_file_name)
	new_file_name1 = tail.split(".")[0] + "_cds.fa"
	return new_file_name1

##extract cds from fasta file 
def extr_sequence(fast_file_name,fasta_file):
	new_file_name = get_new_name(fast_file_name)
	print(new_file_name)
	id_list= []
	for seq_record in SeqIO.parse(fast_file_name, "fasta"):
		sequence_id = str(seq_record.id)
		id_list.append(sequence_id.rstrip().strip('\t'))
	print(id_list)
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		if str(seq_record.id) in id_list:
      seq_remove_stop_coden = seq_record.seq[:-3] #remove stop coden
			rename_rec = SeqRecord(Seq(seq_remove_stop_coden),id = seq_record.id, description = "")
			with open(new_file_name, "a") as handle:
				SeqIO.write(rename_rec, handle, "fasta")

folder_list = os.listdir("./")
folder_list1 = folder_list[:-1]
for fasta_file in folder_list1:
	extr_sequence(fasta_file,"../all.fa")
