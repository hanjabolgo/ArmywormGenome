################################################################################
# Evolution Analysis
################################################################################
base=~/Project/Armyworms/
rawdata=~/Project/Armyworms/00_RawData/
software=~/biosoft
scripts=~/scripts/

cd $base
mkdir 04_Evolution
cd $base/04_Evolution
mkdir Phylogenetic Cafe OrthoFinder Syntenicanalysis

################################################################################
# 10 Phylogenetic analysis
################################################################################

################################################################################
# Run OrthoFinder
################################################################################
mkdir -p OrthoFinder/insect OrthoFinder/Noctuoidea
cd OrthoFinder/insect
ln -s ~/DataBase/Insects_Reference/Noctuoidea/Protein/*.fa ./
ln -s ~/DataBase/Insects_Reference/Insect/Protein/*.fa ./
$software/OrthoFinder-2.5.4/orthofinder -t 24 -a 24 -S diamond_ultra_sens -y -insect -og -o insect_a
a
################################################################################
# Phylogenetic analysis
################################################################################
# Extract cds 
cd  $base/04_Evolution/Phylogenetic
ln -s $base/04_Evolution/OrthoFinder/insect/insect_aa/Orthogroups/ ./
ln -s ~/DataBase/Insects_Reference/Noctuoidea/CDS/*_cds.fa ./
ln -s ~/DataBase/Insects_Reference/Insect/CDS/*_cds.fa ./
cat *_cds.fa > all.fa
mkdir extrasinglecopycds && cd extrasinglecopycds
ln -s $base/04_Evolution/OrthoFinder/insect/insect_aa/Orthogroups/Single_Copy_Orthologue_Sequences/* python $script/extract_cds.py 
cd ..
mkdir reformcds && reformcds
ln -s ../extrasinglecopycds/*_cds.fa ./
python $script/trimc_ds.py
cd ..

# Alignment with prank
mkdir PrankAlignment && cd PrankAlignment
ln -s ../*_reformed_cds.fa ./
ls *cds*|while read id ;do (echo prank -codon -d=${id} -o=${id}.aln -F -t=noc.tree -once);done > com.sh
ParaFly -c com.sh -CPU 40
cd ..

# Filter sequences using Gblocks
mkdir Gblocks && cd Gblocks
ln -s ../PrankAlignment/*aln.best.fas ./
ls *aln.best.fas | while read id ;do (Gblocks ${id} -t=c -b5=n -b3=5);done > gblocks.log
# Then we remove alignments file less than 120 bp manully 
cd ../

# Merge sequences
mkdir Mergesequences && Mergesequences
ln -s ../Gblocks/*-gb ./
python $script/Mergesequences.py
mv test.fasta 23insect.fa
cd ../

# selected the substitution model using ModelTest-NG
mkdir Modeltest && Modeltest
ln -s ../Mergesequences/23insect.fa ./
$software/modeltest-ng/modeltest-ng -i 23insect.fa -d nt
cd ../

# constructed the maximum likelihood phylogenetic tree by RAxML-NG
mkdir Raxml && Raxml
ln -s ../Mergesequences/23insect.fa ./
$software/raxml-ng/raxml-ng --model GTR+I+G4 --all --msa 23insect.fa --bs-trees 1000 --threads 36 --seed 2 --tree test.tree --prefix insect23
mv RAxML_bestTree.insect23 RAxML_bestTree.tre

################################################################################
# Run cafe
################################################################################
cd  $base/04_Evolution/Cafe/
ln -s $base/04_Evolution/OrthoFinder/insect/insect_aa/Orthogroups/Orthogroups.GeneCount.tsv ./
ln -s $base/04_Evolutio/Phylogenetic/raxml/RAxML_bestTree.insect23 ./
awk -v  OFS="\t" '{if($1=="Orthogroup"){print"Descript",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}else{print"Null",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}}' Orthogroups.GeneCount.tsv > GeneCounts.tsv
python $script/filterorthofinder.py > GeneCounts_ex0.tsv
$software/CAFE5/bin/cafe5 -i GeneCounts_ex0.tsv -t RAxML_bestTree.tre --cores 36 -e -o 20211015_error_distribution
ln -s 20211015_error_distribution/Base_error_model.txt
$software/CAFE5/bin/cafe5 -i GeneCounts_ex0.tsv -t FigTree.tre -eBase_error_model.txt -o uni_corr_k2 -k 2 --cores 36 
