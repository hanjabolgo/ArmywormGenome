import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

names = globals() 
folder_list = os.listdir("./")
folder_list1 = folder_list[:-1]
prefix_list= ["Aips", "Apla", "Ayam", "Bmor", "Cexi", "Cpom", "Csup", "Dmel", "Dple", "Harm", "Hvir", "Mlor", "Msep", "Msex", "Obru", "Ofur", "Pbia", "Pxyl", "Sexi", "Sfru", "Slit", "Stie", "Tni"]

for i in prefix_list:
	names[ str(i) ] = ''

for fasta_file in folder_list1:
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		if len(seq_record.seq) > 99:
			seq_record_id = seq_record.id
			seq_id = seq_record_id.split("_", 1)[0]
			if seq_id == "Aips":
				seq_id_seq = seq_record.seq
				Aips = Aips + seq_id_seq
			elif seq_id == "Apla":
				seq_id_seq = seq_record.seq
				Apla = Apla + seq_id_seq
			elif seq_id == "Ayam":
				seq_id_seq = seq_record.seq
				Ayam = Ayam + seq_id_seq
			elif seq_id == "Bmor":
				seq_id_seq = seq_record.seq
				Bmor = Bmor + seq_id_seq
			elif seq_id == "Cexi":
				seq_id_seq = seq_record.seq
				Cexi = Cexi + seq_id_seq
			elif seq_id == "Cpom":
				seq_id_seq = seq_record.seq
				Cpom = Cpom + seq_id_seq
			elif seq_id == "Csup":
				seq_id_seq = seq_record.seq
				Csup = Csup + seq_id_seq
			elif seq_id == "Dmel":
				seq_id_seq = seq_record.seq
				Dmel = Dmel + seq_id_seq
			elif seq_id == "Dple":
				seq_id_seq = seq_record.seq
				Dple = Dple + seq_id_seq
			elif seq_id == "Harm":
				seq_id_seq = seq_record.seq
				Harm = Harm + seq_id_seq
			elif seq_id == "Hvir":
				seq_id_seq = seq_record.seq
				Hvir = Hvir + seq_id_seq
			elif seq_id == "Mlor":
				seq_id_seq = seq_record.seq
				Mlor = Mlor + seq_id_seq
			elif seq_id == "Msep":
				seq_id_seq = seq_record.seq
				Msep = Msep + seq_id_seq
			elif seq_id == "Msex":
				seq_id_seq = seq_record.seq
				Msex = Msex + seq_id_seq
			elif seq_id == "Obru":
				seq_id_seq = seq_record.seq
				Obru = Obru + seq_id_seq
			elif seq_id == "Ofur":
				seq_id_seq = seq_record.seq
				Ofur = Ofur + seq_id_seq
			elif seq_id == "Pbia":
				seq_id_seq = seq_record.seq
				Pbia = Pbia + seq_id_seq
			elif seq_id == "Pxyl":
				seq_id_seq = seq_record.seq
				Pxyl = Pxyl + seq_id_seq
			elif seq_id == "Sexi":
				seq_id_seq = seq_record.seq
				Sexi = Sexi + seq_id_seq
			elif seq_id == "Sfru":
				seq_id_seq = seq_record.seq
				Sfru = Sfru + seq_id_seq
			elif seq_id == "Slit":
				seq_id_seq = seq_record.seq
				Slit = Slit + seq_id_seq
			elif seq_id == "Stie":
				seq_id_seq = seq_record.seq
				Stie = Stie + seq_id_seq
			elif seq_id == "Tni":
				seq_id_seq = seq_record.seq
				Tni = Tni + seq_id_seq

Aips_rec = SeqRecord(Seq(Aips),id = "Aips", description = "")
Apla_rec = SeqRecord(Seq(Apla),id = "Apla", description = "")
Ayam_rec = SeqRecord(Seq(Ayam),id = "Ayam", description = "")
Bmor_rec = SeqRecord(Seq(Bmor),id = "Bmor", description = "")
Cexi_rec = SeqRecord(Seq(Cexi),id = "Cexi", description = "")
Cpom_rec = SeqRecord(Seq(Cpom),id = "Cpom", description = "")
Csup_rec = SeqRecord(Seq(Csup),id = "Csup", description = "")
Dmel_rec = SeqRecord(Seq(Dmel),id = "Dmel", description = "")
Dple_rec = SeqRecord(Seq(Dple),id = "Dple", description = "")
Harm_rec = SeqRecord(Seq(Harm),id = "Harm", description = "")
Hvir_rec = SeqRecord(Seq(Hvir),id = "Hvir", description = "")
Mlor_rec = SeqRecord(Seq(Mlor),id = "Mlor", description = "")
Msep_rec = SeqRecord(Seq(Msep),id = "Msep", description = "")
Msex_rec = SeqRecord(Seq(Msex),id = "Msex", description = "")
Obru_rec = SeqRecord(Seq(Obru),id = "Obru", description = "")
Ofur_rec = SeqRecord(Seq(Ofur),id = "Ofur", description = "")
Pbia_rec = SeqRecord(Seq(Pbia),id = "Pbia", description = "")
Pxyl_rec = SeqRecord(Seq(Pxyl),id = "Pxyl", description = "")
Sexi_rec = SeqRecord(Seq(Sexi),id = "Sexi", description = "")
Sfru_rec = SeqRecord(Seq(Sfru),id = "Sfru", description = "")
Slit_rec = SeqRecord(Seq(Slit),id = "Slit", description = "")
Stie_rec = SeqRecord(Seq(Stie),id = "Stie", description = "")
Tni_rec = SeqRecord(Seq(Tni),id = "Tni", description = "")

with open("test.fa", "a") as handle:
	SeqIO.write(Aips_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Apla_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Ayam_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Bmor_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Cexi_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Cpom_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Csup_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Dmel_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Dple_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Harm_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Hvir_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Mlor_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Msep_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Msex_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Obru_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Ofur_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Pbia_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Pxyl_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Sexi_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Sfru_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Slit_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Stie_rec, handle, "fasta")
with open("test.fa", "a") as handle:
	SeqIO.write(Tni_rec, handle, "fasta") 
