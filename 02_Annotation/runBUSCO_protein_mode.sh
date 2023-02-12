###########################################
# QC of genome annotations
###########################################
# We ran BUSCO on the set of predicted proteins for each of the ten Noctuoidea moths genomes.

#######################
# BUSCO
#######################
base=~/Project/Armyworms/01_Assembly/
export BUSCOBase="/usr/local/home/zhaohanbo/DataBase/BUSCO/insecta_odb10"

#######################
# run for each species
#######################

cd $base/02_Annotation
mkdir Annotation_BUSCO
cd Annotation_BUSCO
ln -s $base/02_Annotation/Functionannotation/Msep/protein.fa ./Msep_pep.fa
ln -s $base/02_Annotation/Functionannotation/Mlor/protein.fa ./Mlor_pep.fa
ln -s ~/DataBase/Insects_Reference/Noctuoidea/Protein/*.fa ./

# protein mode
ls *.fa|cut -d"_" -f 1|while read id;do busco -m prot -i ${id}_pep.fa -l $BUSCOBase -o busco${id} -c 40 -f --offline; done
