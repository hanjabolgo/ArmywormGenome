################################################################################
# Sex Chromosome Detection
################################################################################
base=~/Project/Armyworms/
rawdata=~/Project/Armyworms/00_RawData/
software=~/biosoft
scripts=~/scripts/

################################################################################
# 9 Sex Chromosome Detection
################################################################################
# Calculate Sequencing depth of male and female M. separata
cd $base
mkdir 03_SexChromosomeDetect && cd 03_SexChromosomeDetect
mkdir -p ZWDetect/Msep ZWDetect/Mlor
cd ZWDetect/Msep
ln -s $base/01_Assembly/Hi_C/MsepF/MsepF_genome_assembly.fa Msep.genome
ln -s $rawdata/resequence/Msep/*gz ./
ln -s $base/02_Annotation/Genepredict/EVM/Msep/Msep.gff3
grep -o -E "Chr[0-9]+" Msep.genome > Chr.id
grep -f Chr.id Msep.genome -A 1 > Chr.fa
bwa mem -t 24 -M -R '@RG\tID:MSEPM2F\tSM:M2F\tLB:MSEP\tPL:Illumina' Chr.fa Orientalarmyworm_AE875-03-R02_good_1.fq.gz Orientalarmyworm_AE875-03-R02_good_2.fq.gz > M2F.sam 2> ./mem-peM.log
samtools view -@ 24 -bS M2F.sam -o M2F.bam
samtools sort -o M2F_sorted.bam M2F.bam
samtools index -@ 24 M2F_sorted.bam
bwa mem -t 6 -M -R '@RG\tID:MSEPF2F\tSM:F2F\tLB:MSEP\tPL:Illumina' Chr.fa Orientalarmyworm_AE875-03-R01_good_1.fq.gz Orientalarmyworm_AE875-03-R01_good_2.fq.gz > F2F.sam 2> ./mem-peF.log
samtools view -@ 6 -bS F2F.sam -o F2F.bam
samtools sort -o F2F_sorted.bam F2F.bam
samtools index -@ 16 F2F_sorted.bam
bedtools  makewindows -g sizes.genome  -w  2000000 -s 500000  > 2M_500k.bed
samtools bedcov 2M_500k.bed F2F_sorted.bam > F2F.counts.txt.tmp
samtools bedcov 2M_500k.bed M2F_sorted.bam > M2F.counts.txt.tmp
samtools faidx Chr.fa
grep "Chr" Msep.gff3 > Msep_chr.gff

# Generate chromosome files 
cat Chr.fa.fai |cut -f 1-2| grep "Chr"|sort -V > Msep.genome.txt
bedtools makewindows -g Msep.genome.txt  -w 2000000 -s 500000  > genome.windows
grep '[[:blank:]]gene[[:blank:]]' Msep_chr.gff | cut -f 1,4,5 | awk '{print $1"\t"$2"\t"$3}'|grep Chr > genes.bed
cat Msep.genome.fai|awk '{print "chr - "$1" "$1" 0 "$2" "$1}' > karyotype.Msep.tmp
sed -i "s/\s/\t/g" karyotype.Msep.tmp
cat karyotype.Msep.tmp | sort -t "r" -n -k 4 > karyotype.Msep.txt

# Calculate GC density 
bedtools nuc -fi Chr.fa -bed Msep.windows |cut -f 1-3,5 | grep -v "#" | awk -vFS="\t" -vOFS="\t" '{print $1,$2,$3,$4}' > circos.gc.density.txt

#  Calculate Gene density 
grep '[[:blank:]]gene[[:blank:]]' Msep_chr.gff | awk '{print $1"\t"$4"\t"$5}'| bedtools coverage -a genome.windows -b - |cut -f 1-4 > circos.gene.density.txt

# Calculate log2 M2F cover 
paste genome.windows M2F.counts.txt.tmp F2F.counts.txt.tmp > M2F_cov_raw.txt
python $script/calculatelog2.py
mv M2F_cov_nor_log2.txt circos.m2fcover.density.txt
sed- i "/^\d+/d" circos.m2fcover.density.txt

#Calculate TE density 
python $script/calculatetedensity.py > circos.TE.density.txt

# plot circos
circos -conf circos.conf

# Calculate Sequencing depth of male and female M. loreyi
cd $base/01_Assembly
mkdir -p ZWDetect/Msep ZWDetect/Mlor
cd ZWDetect/Mlor
ln -s $rawdata/resequence/Mlor/*gz ./
ln -s $base/01_Assembly/Hi_C/MlorF/MlorF_genome_assembly.fa Mlor.genome
ln -s $base/02_Annotation/Genepredict/EVM/Mlor/Mlor.gff3
grep -o -E "Chr[0-9]+" Mlor.genome > Chr.id
grep -f Chr.id Mlor.genome -A 1 > Chr.fa
bwa mem -t 24 -M -R '@RG\tID:MLORM2F\tSM:M2F\tLB:MLOR\tPL:Illumina' Chr.fa Leucanialoreyiduponchel_AE875-03-R03_good_1.fq.gz  Leucanialoreyiduponchel_AE875-03-R03_good_2.fq.gz > M2F.sam 2> ./mem-peM.log
samtools view -@ 24 -bS M2F.sam -o M2F.bam
samtools sort -o M2F_sorted.bam M2F.bam
samtools index -@ 24 M2F_sorted.bam
bwa mem -t 6 -M -R '@RG\tID:MLORF2F\tSM:F2F\tLB:MLOR\tPL:Illumina' Chr.fa Leucanialoreyiduponchel_AE875-03-R04_good_1.fq.gz  Leucanialoreyiduponchel_AE875-03-R04_good_2.fq.gz > F2F.sam 2> ./mem-peF.log
samtools view -@ 6 -bS F2F.sam -o F2F.bam
samtools sort -o F2F_sorted.bam F2F.bam
samtools index -@ 16 F2F_sorted.bam
bedtools  makewindows -g sizes.genome  -w  2000000 -s 500000  > 2M_500k.bed
samtools bedcov 2M_500k.bed F2F_sorted.bam > F2F.counts.txt.tmp
samtools bedcov 2M_500k.bed M2F_sorted.bam > M2F.counts.txt.tmp
samtools faidx Chr.fa
grep "Chr" Mlor.gff3 > Mlor_chr.gff

# Generate chromosome files
cat Chr.fa.fai |cut -f 1-2| grep "Chr"|sort -V > Mlor.genome.txt
bedtools makewindows -g Mlor.genome.txt  -w 2000000 -s 500000  > genome.windows
grep '[[:blank:]]gene[[:blank:]]' Mlor_chr.gff | cut -f 1,4,5 | awk '{print $1"\t"$2"\t"$3}'|grep Chr > genes.bed
cat Mlor.genome.fai|awk '{print "chr - "$1" "$1" 0 "$2" "$1}' > karyotype.Mlor.tmp
sed -i "s/\s/\t/g" karyotype.Mlor.tmp
cat karyotype.Mlor.tmp | sort -t "r" -n -k 4 > karyotype.Mlor.txt

# Calculate GC density
bedtools nuc -fi Mlor.genome -bed Mlor.windows |cut -f 1-3,5 | grep -v "#" | awk -vFS="\t" -vOFS="\t" '{print $1,$2,$3,$4}' > circos.gc.density.txt

#  Calculate Gene density
grep '[[:blank:]]gene[[:blank:]]' Mlor_chr.gff | awk '{print $1"\t"$4"\t"$5}'| bedtools coverage -a Mlor.windows -b - |cut -f 1-4 > circos.gene.density.txt

# Calculate log2 M2F cover
paste genome.windows M2F.counts.txt.tmp F2F.counts.txt.tmp > M2F_cov_raw.txt
python $script/calculatelog2.py
mv M2F_cov_nor_log2.txt circos.m2fcover.density.txt
sed- i "/^\d+/d" circos.m2fcover.density.txt

# Calculate TE density
python $script/calculatetedensity.py > circos.TE.density.txt


# plot circos
circos -conf circos.conf
