###########################################
# QC of genome assemblies
###########################################
# We ran BUSCO on the genome assembly for each of the ten Noctuoidea moths genomes.

#######################
# BUSCO
#######################

# prepare environment
base=~/Project/Armyworms/01_Assembly/
BUSCOBase=~/DataBase/BUSCO/insecta_odb10

#######################
# run for each species
#######################

mkdir -p $base/Asembly_BUSCO/
cd $base/Asembly_BUSCO/
ln -s ~/DataBase/Insects_Reference/Noctuoidea/Genomes/*.fa ./
ln -s $base/Hi_C/MsepF/MsepF_genome_assembly.fa .
ln -s $base/Hi_C/MsepM/MsepM_genome_assembly.fa .
ln -s $base/Hi_C/MlorF/MlorF_genome_assembly.fa .
ln -s $base/Hi_C/MlorM/MlorM_genome_assembly.fa .

# genome mode
ls *.fa|cut -d"_" -f 1|while read id;do busco -m geno -i ${id}_genome_assembly.fa -l $BUSCOBase -o busco${id} -c 40 -f --offline; done 
