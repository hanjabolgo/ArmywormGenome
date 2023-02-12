################################################################################
# Genome Assembly
################################################################################
base=~/Project/Armyworms/
rawdata=~/Project/Armyworms/00_RawData/
software=~/biosoft
scripts=~/scripts/

cd $base
mkdir 01_Assembly/ 02_Anntation/ 
cd $base/01_Assembly/
mkdir -p GenomeAssembly/MsepF GenomeAssembly/MsepM GenomeAssembly/MlorF GenomeAssembly/MlorM

################################################################################
# 2 Genome assembly 
################################################################################

################################################################################
# 2.1 Using canu to assemble female M. separata genomes
################################################################################
cd $base/01_Assembly/GenomeAssembly/MsepF
$software/canu-2.1.1/bin/canu -p hicanuMsepF -d hicanuMsepF genomeSize=660m maxThreads=40 useGrid=false "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" -pacbio-hifi $rawdata/HiFiReads/MsepF_ccs.fasta.gz &> MsepF.log

################################################################################
# 2.2 Using canu to assemble male M. separata genomes
################################################################################
cd $base/01_Assembly/GenomeAssembly/MsepM
$software/canu-2.1.1/bin/canu -p hicanuMsepM -d hicanuMsepM genomeSize=660m maxThreads=40 useGrid=false "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" -pacbio-hifi $rawdata/HiFiReads/MsepM_ccs.fasta.gz &> MsepM.log

################################################################################
# 2.3 Using canu to assemble female M. loreyi genomes
################################################################################
cd $base/01_Assembly/GenomeAssembly/MlorF
$software/canu-2.1.1/bin/canu -p hicanuMlorF -d hicanuMlorF genomeSize=660m maxThreads=40 useGrid=false "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" -pacbio-hifi $rawdata/HiFiReads/MlorF_ccs.fasta.gz &> MlorF.log

################################################################################
# 2.4 Using canu to assemble male M. loreyi genomes
################################################################################
cd $base/01_Assembly/GenomeAssembly/MlorM
$software/canu-2.1.1/bin/canu -p hicanuMlorM -d hicanuMlorM genomeSize=660m maxThreads=40 useGrid=false "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" -pacbio-hifi $rawdata/HiFiReads/MlorM_ccs.fasta.gz &> MlorM.log

################################################################################
# 3 Purge haplotigs and overlaps
################################################################################
conda activate purge_dups

################################################################################
# 3.1 Using Purge_Dups to purging haplotigs and overlaps in the assembly of female M. separata
################################################################################
cd $base/01_Assembly/GenomeAssembly/MsepF
$software/minimap2.1/minimap2 -t 20 -x map-pb hicanuMsepF.contigs.fasta $rawdata/HiFiReads/MsepF_ccs.fasta.gz | gzip -c - > pb_aln.paf.gz
pbcstat pb_aln.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log
split_fa hicanuMsepF.contigs.fasta > asm.split
$software/minimap2.1/minimap2 -t 20 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed hicanuMsepF.contigs.fasta
mv purged.fa hicanuMsepF.purged.fa

################################################################################
# 3.2 Using Purge_Dups to purging haplotigs and overlaps in the assembly of male M. separata
################################################################################
cd $base/01_Assembly/GenomeAssembly/hicanuMsepM
$software/minimap2.1/minimap2 -t 20 -x map-pb hicanuMsepM.contigs.fasta $rawdata/HiFiReads/MsepM_ccs.fasta.gz | gzip -c - > pb_aln.paf.gz
pbcstat pb_aln.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log
split_fa hicanuMsepM.contigs.fasta > asm.split
minimap2 -t 20 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed hicanuMsepM.contigs.fasta
mv purged.fa hicanuMsepM.purged.fa

################################################################################
# 3.3 Using Purge_Dups to purging haplotigs and overlaps in the assembly of female M. loreyi
################################################################################
cd $base/01_Assembly/GenomeAssembly/hicanuMlorF
$software/minimap2.1/minimap2 -t 20 -x map-pb hicanuMlorF.contigs.fasta $rawdata/HiFiReads/MlorF_ccs.fasta.gz | gzip -c - > pb_aln.paf.gz
pbcstat pb_aln.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log
split_fa hicanuMlorF.contigs.fasta > asm.split
$software/minimap2.1/minimap2 -t 20 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed hicanuMlorF.contigs.fasta
mv purged.fa hicanuMlorF.purged.fa

################################################################################
# 3.4 Using Purge_Dups to purging haplotigs and overlaps in the assembly of male M. loreyi
################################################################################
cd $base/01_Assembly/GenomeAssembly/hicanuMlorM
$software/minimap2.1/minimap2 -t 20 -x map-pb hicanuMlorM.contigs.fasta $rawdata/HiFiReads/MlorM_ccs.fasta.gz | gzip -c - > pb_aln.paf.gz
pbcstat pb_aln.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log
split_fa hicanuMlorM.contigs.fasta > asm.split
$software/minimap2.1/minimap2 -t 20 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed hicanuMlorM.contigs.fasta
mv purged.fa hicanuMlorM.purged.fa

################################################################################
# 4 Remove the contaminating bacterial contigs
################################################################################
cd $base/01_Assembly/
mkdir -p FilterLGT/Animal_DB/ FilterLGT/Bacteria_DB/ FilterLGT/MsepF FilterLGT/MsepM FilterLGT/MlorF FilterLGT/MlorM
conda activate edirect 

################################################################################
# 4.1 Building the blast database
# The GenBank Accession Number is collected from "Olafson, P. U., Aksoy, S., Attardo, G. M., Buckmeier, G., Chen, X., Coates, C. J., . . . Friedrich, M. (2021). The genome of the stable fly, Stomoxys calcitrans, reveals potential mechanisms underlying reproduction, host interactions, and novel targets for pest control. BMC biology, 19(1), 1-31." 
################################################################################
cd $base/01_Assembly/FilterLGT/Bacteria_DB/
epost -input Bacteria_LGT_Pipeline.list -db nucleotide | efetch -format fasta > Bacteria_LGT_Pipeline.fa
$software/ncbi-blast-2.12.0+/bin/makeblastdb -in Bacteria_LGT_Pipeline.fa -dbtype nucl -out Bacteria_LGT_Pipeline -parse_seqids

cd $base/01_Assembly/FilterLGT/Animal_DB/
$software/ncbi-blast-2.12.0+/bin/makeblastdb -in Animal_LGT_Pipeline.fa -dbtype nucl -out Animal_LGT_Pipeline -parse_seqids

################################################################################
# 4.2 Remove the contaminating bacterial contigs in the genome assembly of female M. separata
################################################################################ 
cd $base/01_Assembly/FilterLGT/MsepF
ln -s $base/01_Assembly/GenomeAssembly/hicanuMsepF/hicanuMsepF.purged.fa ./
awk '{if(NR%2==1)a=$1;else{len=length($0); i=1; while(i<=len-1000){print a"_"i;print (substr($0,i,1000));i=i+1000}}}' hicanuMsepF.purged.fa > hicanuMsepF.purged.split.fa
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMsepF.purged.split.fa -db $base/01_Assembly/FilterLGT/Bacteria_DB/Bacteria_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Bacteria_hicanuMsepF
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMsepF.purged.split.fa -db $base/01_Assembly/FilterLGT/Animal_DB/Animal_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Animal_hicanuMsepF

# We removed contigs with strong bacterial bitscores but weak vertebrate scores manually, and got the new file ‘hicanuMsepF.purged.FINAL.corrected.fa’.  
mv hicanuMsepF.purged.FINAL.corrected.fa MsepF.genome.fa

################################################################################
# 4.3 Remove the contaminating bacterial contigs in the genome assembly of male M. separata
################################################################################
cd $base/01_Assembly/FilterLGT/MsepM
ln -s $base/01_Assembly/GenomeAssembly/hicanuMsepM/hicanuMsepM.purged.fa ./
awk '{if(NR%2==1)a=$1;else{len=length($0); i=1; while(i<=len-1000){print a"_"i;print (substr($0,i,1000));i=i+1000}}}' hicanuMsepM.purged.fa > hicanuMsepM.purged.split.fa
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMsepM.purged.split.fa -db $base/01_Assembly/FilterLGT/Bacteria_DB/Bacteria_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Bacteria_hicanuMsepM
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMsepM.purged.split.fa -db $base/01_Assembly/FilterLGT/Animal_DB/Animal_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Animal_hicanuMsepM

# We removed contigs with strong bacterial bitscores but weak vertebrate scores manually, and got the new file ‘hicanuMsepM.purged.FINAL.corrected.fa’.
mv hicanuMsepM.purged.FINAL.corrected.fa MsepM.genome.fa

################################################################################
# 4.4 Remove the contaminating bacterial contigs in the genome assembly of female M. loreyi
################################################################################
cd $base/01_Assembly/FilterLGT/MlorF
ln -s $base/01_Assembly/GenomeAssembly/hicanuMlorF/hicanuMlorF.purged.fa ./
awk '{if(NR%2==1)a=$1;else{len=length($0); i=1; while(i<=len-1000){print a"_"i;print (substr($0,i,1000));i=i+1000}}}' hicanuMlorF.purged.fa > hicanuMlorF.purged.split.fa
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMlorF.purged.split.fa -db $base/01_Assembly/FilterLGT/Bacteria_DB/Bacteria_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Bacteria_hicanuMlorF
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMlorF.purged.split.fa -db $base/01_Assembly/FilterLGT/Animal_DB/Animal_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Animal_hicanuMlorF

# We removed contigs with strong bacterial bitscores but weak vertebrate scores manually, and got the new file ‘hicanuMlorF.purged.FINAL.corrected.fa’.
mv hicanuMlorF.purged.FINAL.corrected.fa MlorF.genome.fa

################################################################################
# 4.5 Remove the contaminating bacterial contigs in the genome assembly of male M. loreyi
################################################################################
cd $base/01_Assembly/FilterLGT/MlorM
ln -s $base/GenomeAssembly/hicanuMlorM/hicanuMlorM.purged.fa ./
awk '{if(NR%2==1)a=$1;else{len=length($0); i=1; while(i<=len-1000){print a"_"i;print (substr($0,i,1000));i=i+1000}}}' hicanuMlorM.purged.fa > hicanuMlorM.purged.split.fa
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMlorM.purged.split.fa -db $base/01_Assembly/FilterLGT/Bacteria_DB/Bacteria_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Bacteria_hicanuMlorM
$software/ncbi-blast-2.12.0+/bin/blastn -query hicanuMlorM.purged.split.fa -db $base/01_Assembly/FilterLGT/Animal_DB/Animal_LGT_Pipeline -evalue 1e-5 -outfmt 6 -num_threads 16 -out Animal_hicanuMlorM

# We removed contigs with strong bacterial bitscores but weak vertebrate scores manually, and got the new file ‘hicanuMlorM.purged.FINAL.corrected.fa’.
mv hicanuMlorM.purged.FINAL.corrected.fa MlorM.genome.fa


################################################################################
# 5 Genome scaffolding 
################################################################################
conda activate 3ddna

################################################################################
# 5.1  Test quality of Hi-C data 
################################################################################
# Check the quality of Hi-C reads of female M. separata
mkdir -p $base/01_Assembly/Hi_C/HiCPro/MsepF/Msep.ref
cd $base/01_Assembly/Hi_C/HiCPro/MsepF/Msep.ref
ln -s $base/01_Assembly/FilterLGT/MsepF/MsepF.genome.fa ./
cd ..
$software/HiC-Pro_3.0.0/bin/utils/digest_genome.py -r A^AGCTT -o MseqF.bed Msep.ref/MsepF.genome.fa
samtools faidx Msep.ref/MsepF.genome.fa && awk '{print $1"\t"$2}' Msep.ref/MsepF.genome.fa.fai > MsepF.genome.fa.size
find $rawdata/Hi-C/MsepF | xargs ls -d > MsepF.reads
$software/HiC-Pro_3.0.0/bin/HiC-Pro --input Msep.reads --output hicpro_output --conf config-hicpro.txt

# Check the quality of Hi-C reads of male M. separata
mkdir -p $base/01_Assembly/Hi_C/HiCPro/MsepM/Msep.ref
cd $base/01_Assembly/Hi_C/HiCPro/MsepM/Msep.ref
ln -s $base/01_Assembly/FilterLGT/MsepM/MsepM.genome.fa ./
cd ..
$software/HiC-Pro_3.0.0/bin/utils/digest_genome.py -r A^AGCTT -o MseqM.bed Msep.ref/MsepM.genome.fa
samtools faidx Msep.ref/MsepM.genome.fa && awk '{print $1"\t"$2}' Msep.ref/MsepM.genome.fa.fai > MsepM.genome.fa.size
find $rawdata/Hi-C/MsepM | xargs ls -d > MsepM.reads
$software/HiC-Pro_3.0.0/bin/HiC-Pro --input Msep.reads --output hicpro_output --conf config-hicpro.txt

# Check the quality of Hi-C reads of female M. loreyi
mkdir -p $base/01_Assembly/Hi_C/HiCPro/MlorF/Mlor.ref
cd $base/01_Assembly/Hi_C/HiCPro/MlorF/Mlor.ref
ln -s $base/01_Assembly/FilterLGT/MlorF/MlorF.genome.fa ./
cd ..
$software/HiC-Pro_3.0.0/bin/utils/digest_genome.py -r A^AGCTT -o MlorF.bed Mlor.ref/MlorF.genome.fa
samtools faidx Msep.ref/MlorF.genome.fa && awk '{print $1"\t"$2}' Mlor.ref/MsepF.genome.fa.fai > MlorF.genome.fa.size
find $rawdata/Hi-C/MlorF | xargs ls -d > MlorF.reads
$software/HiC-Pro_3.0.0/bin/HiC-Pro --input Mlor.reads --output hicpro_output --conf config-hicpro.txt

# Check the quality of Hi-C reads of male M. loreyi
mkdir -p $base/01_Assembly/Hi_C/HiCPro/MlorM/Mlor.ref
cd $base/01_Assembly/Hi_C/HiCPro/MlorM/Mlor.ref
ln -s $base/01_Assembly/FilterLGT/MlorM/MlorM.genome.fa ./
cd ..
$software/HiC-Pro_3.0.0/bin/utils/digest_genome.py -r A^AGCTT -o MlorM.bed Mlor.ref/MlorM.genome.fa
samtools faidx Mlor.ref/MlorM.genome.fa && awk '{print $1"\t"$2}' Mlor.ref/MlorM.genome.fa.fai > MlorM.genome.fa.size
find $rawdata/Hi-C/MlorM | xargs ls -d > MlorM.reads
$software/HiC-Pro_3.0.0/bin/HiC-Pro --input Mlor.reads --output hicpro_output --conf config-hicpro.txt

################################################################################
# 5.2 Genome scaffolding with the Juicer, JuiceBox, 3D-DNA pipeline
################################################################################
# Scaffolding the genome of female M. separata with the Juicer, JuiceBox, 3D-DNA pipeline
cd $base/01_Assembly/Hi_C
mkdir -p MsepF/reference MsepF/restriction_sites MsepF/fastq MsepF/scripts
cd $base/01_Assembly/Hi_C/MsepF/reference
ln -s $base/01_Assembly/FilterLGT/MsepF/MsepF.genome.fa ./
bwa index MsepF.genome.fa
generate_site_positions.py HindIII MsepF MsepF.genome.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' MsepF_HindIII.txt > MsepF.chrom.sizes
mv MsepF* ../restriction_sites
cd ../fastq
ln -s $rawdata/Hi-C/MsepF/Orientalarmyworm_AE742-H02-1_good_1.fq.gz MsepF_R1.fastq.gz
ln -s $rawdata/Hi-C/MsepF/Orientalarmyworm_AE742-H02-1_good_2.fq.gz MsepF_R2.fastq.gz
cd ../scripts
cp $software/juicer/scripts/* ./
cd ..
./scripts/juicer.sh -d ./ -D ./ -z ./references/MsepF.genome.fa -y ./restriction_sites/MsepF_HindIII.txt -p ./restriction_sites/MsepF.chrom.sizes -s HindIII -t 10
$software/3d-dna/run-asm-pipeline.sh -r 0 references/MsepF.genome.fa aligned/merged_nodups.txt

# Editing with Juicebox and running 3D-DNA again
$software/3d-dna/run-asm-pipeline-post-review.sh -r MsepF.genome.final.review.assembly references/MsepF.genome.fa aligned/merged_nodups.txt
mv MsepF.genome.FINAL.corrected.fa MsepF_genome_assembly.fa

# Scaffolding the genome of male M. separata with the Juicer, JuiceBox, 3D-DNA pipeline
cd $base/01_Assembly/Hi_C
mkdir -p MsepM/reference MsepM/restriction_sites MsepM/fastq MsepM/scripts
cd $base/01_Assembly/Hi_C/MsepM/reference
ln -s $base/01_Assembly/FilterLGT/MsepM/MsepM.genome.fa ./
bwa index MsepM.genome.fa
generate_site_positions.py HindIII MsepM MsepM.genome.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' MsepM_HindIII.txt > MsepM.chrom.sizes
mv MsepM* ../restriction_sites
cd ../fastq
ln -s $rawdata/Hi-C/MsepM/Orientalarmyworm_AE742-H01-1_good_1.fq.gz MsepM_R1.fastq.gz
ln -s $rawdata/Hi-C/MsepM/Orientalarmyworm_AE742-H01-1_good_2.fq.gz MsepM_R2.fastq.gz
cd ../scripts
cp $software/juicer/scripts/* ./
cd ..
./scripts/juicer.sh -d ./ -D ./ -z ./references/MsepM.genome.fa -y ./restriction_sites/MsepM_HindIII.txt -p ./restriction_sites/MsepM.chrom.sizes -s HindIII -t 10
$software/3d-dna/run-asm-pipeline.sh -r 0 references/MsepM.genome.fa aligned/merged_nodups.txt

# Editing with Juicebox and running 3D-DNA again
$software/3d-dna/run-asm-pipeline-post-review.sh -r MsepM.genome.final.review.assembly references/MsepM.genome.fa aligned/merged_nodups.txt
mv MsepM.genome.FINAL.corrected.fa MsepM_genome_assembly.fa

# Scaffolding the genome of female M. loreyi with the Juicer, JuiceBox, 3D-DNA pipeline
cd $base/01_Assembly/Hi_C
mkdir -p MlorF/reference MlorF/restriction_sites MlorF/fastq MlorF/scripts
cd $base/01_Assembly/Hi_C/MlorF/reference
ln -s $base/01_Assembly/FilterLGT/MlorF/MlorF.genome.fa ./
bwa index MlorF.genome.fa
generate_site_positions.py HindIII MlorF MlorF.genome.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' MlorF_HindIII.txt > MlorF.chrom.sizes
mv MlorF* ../restriction_sites
cd ../fastq
ln -s $rawdata/Hi-C/MlorF/Leucanialoreyiduponchel_AE742-H04-1_good_1.fq.gz MlorF_R1.fastq.gz
ln -s $rawdata/Hi-C/MlorF/Leucanialoreyiduponchel_AE742-H04-1_good_1.fq.gz MlorF_R2.fastq.gz
cd ../scripts
cp $software/juicer/scripts/* ./
cd ..
./scripts/juicer.sh -d ./ -D ./ -z ./references/MlorF.genome.fa -y ./restriction_sites/MlorF_HindIII.txt -p ./restriction_sites/MlorF.chrom.sizes -s HindIII -t 10
$software/3d-dna/run-asm-pipeline.sh -r 0 references/MlorF.genome.fa aligned/merged_nodups.txt

# Editing with Juicebox and running 3D-DNA again
$software/3d-dna/run-asm-pipeline-post-review.sh -r MlorF.genome.final.review.assembly references/MlorF.genome.fa aligned/merged_nodups.txt
mv MlorF.genome.FINAL.corrected.fa MlorF_genome_assembly.fa
# Then we rename chromosome IDs manuelly

# Scaffolding the genome of female M. loreyi with the Juicer, JuiceBox, 3D-DNA pipeline
cd $base/01_Assembly/Hi_C/
mkdir -p MlorM/reference MlorM/restriction_sites MlorM/fastq MlorM/scripts
cd $base/01_Assembly/Hi_C/MlorM/reference
ln -s $base/01_Assembly/FilterLGT/MlorM/MlorM.genome.fa ./
bwa index MlorM.genome.fa
generate_site_positions.py HindIII MlorM MlorM.genome.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' MlorM_HindIII.txt > MlorM.chrom.sizes
mv MlorM* ../restriction_sites
cd ../fastq
ln -s $rawdata/Hi-C/MlorM/Leucanialoreyiduponchel_AE742-H03-1_good_1.fq.gz MlorM_R1.fastq.gz
ln -s $rawdata/Hi-C/MlorM/Leucanialoreyiduponchel_AE742-H03-1_good_1.fq.gz MlorM_R2.fastq.gz
cd ../scripts
cp $software/juicer/scripts/* ./
cd ..
./scripts/juicer.sh -d ./ -D ./ -z ./references/MlorM.genome.fa -y ./restriction_sites/MlorM_HindIII.txt -p ./restriction_sites/MlorM.chrom.sizes -s HindIII -t 10
$software/3d-dna/run-asm-pipeline.sh -r 0 references/MlorM.genome.fa aligned/merged_nodups.txt

# Editing with Juicebox and running 3D-DNA again
$software/3d-dna/run-asm-pipeline-post-review.sh -r MlorM.genome.final.review.assembly references/MlorM.genome.fa aligned/merged_nodups.txt
mv MlorM.genome.FINAL.corrected.fa MlorM_genome_assembly.fa
# Then we rename chromosome IDs manuelly
