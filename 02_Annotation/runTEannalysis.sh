################################################################################
# Genome Assembly
################################################################################
base=~/Project/Armyworms/
rawdata=~/Project/Armyworms/00_RawData/
software=~/biosoft
scripts=~/scripts/

cd $base/
mkdir -p 02_Annotation/TE

################################################################################
# 6 TE analysis
################################################################################

cd $base/02_Annotation/TE_analysis/
ln -s $base/01_Assembly/Hi-C/MsepF/MsepF_genome_assembly.fa ./
ln -s $base/01_Assembly/Hi-C/MlorF/MlorF_genome_assembly.fa ./
ls *.fa|cut -d"_" -f 1| while read id;do $software/RepeatMasker/BuildDatabase -name ${id} ${id}_genome_assembly.fa;$software/RepeatMasker/RepeatModeler -database ${id} -pa 5 -LTRStruct geno -i ${id}_genome_assembly.fa;RepeatMasker -pa 20 -html -a -gff -xsmall -lib ${id}-families.fa ${id}_genome_assembly.fa;$software/RepeatMasker/util/calcDivergenceFromAlign.pl -s ${id}.divsum -noCpGMod -a ${id}.align ${id}_genome_assembly.fa.align;done
