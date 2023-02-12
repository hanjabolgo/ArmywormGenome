################################################################################
# RNA-Seq Analysis
################################################################################
base=~/Project/Armyworms/
rawdata=~/Project/Armyworms/00_RawData/
software=~/biosoft
scripts=~/scripts/

cd $base
mkdir 05_RNA-Seq_Analysis
################################################################################
# 12 Gene expression 
################################################################################

# Get gene expression matrix 
cd $base/05_RNA-Seq_Analysis
mkdir -p Expression/Msep/reference Expression/Mlor/reference Expression/Sfru/reference

cd Expression/Msep/reference
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked Msep.fa
hisat2-build Msep.fa Msep
cd ../ 
ln -s $rawdata/Illunima_Seq/Msep/*fastq.gz ./
ln -s $base/02_Annotation/Genepredict/EVM/Msep/Msep.gff3
ls *gz|cut -d"d" -f 1 | sort -u |while read id;do
ls -lh ${id}d_1.fq.gz ${id}d_2.fq.gz
hisat2 -p 36 -x ./reference/Msep -1 ${id}d_1.fq.gz -2 ${id}d_2.fq.gz -S ${id}d.hisat.sam 1>./${id}d.log \
done

ls *.sam|while read id ;do (samtools sort -O bam -@ 24 -o $(basename ${id} ".sam").bam ${id});done
ls *.bam|while read id ;do ( samtools flagstat -@ 24 $id >  $(basename ${id} ".bam").flagstat  );done
ls *.bam|cut -d"d" -f 1 | while read id;do stringtie ${id}d.hisat.bam -eB -G ./Msep.gff3 -o ${id}d.final.stringTie_asm.gtf -A ${id}d.final.gene_abundence.tab -C ${id}d.final.known.cov_refs.gtf; done
rm *.sam *.bam

prepDE.py -i Msep_sample_list.txt -g gene_count_matrix.csv -l 150
ls ./*.tab|while read id;do sed -i "s/Gene Name/Gene_Name/g" ${id}; done
ls ./*.tab|while read id;do sed -i "s/Gene ID/Gene_ID/g" ${id}; done
ls ./*.tab|while read id;do awk 'NR > 1 {print $1 "\t" $8 "\t" $9}' ${id} | sort > ${id}.sort; done
awk '{print "mv" " ./" $2 "_good.gene_abundence.tab.sort " $1 "_good.gene_abundence.tab.sort"}' sampleName_clientId.txt > ChangeSampleName.sh
chmod 755 ChangeSampleName.sh
sh ./ChangeSampleName.sh
ls ./*tab.sort |cut -d"_" -f 1 | while read id; do awk -v simplename=$id '{print $1 "\t" $2 "\t" simplename}' ${id}_good.gene_abundence.tab.sort > ${id}_good.gene_abundence.tab.sort_FPKM.csv; done
ls ./*tab.sort |cut -d"_" -f 1 | while read id; do awk -v simplename=$id '{print $1 "\t" $2 "\t" simplename}' ${id}_good.gene_abundence.tab.sort > ${id}_good.gene_abundence.tab.sort_TPM.csv; done
cat *gene_abundence.tab.sort_FPKM.csv > FPKM.txt
cat *gene_abundence.tab.sort_TPM.csv > TPM.txt
Rscript $software/ReshapeMatrix.R

cd Expression/Mlor/reference
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked Mlor.fa
hisat2-build Mlor.fa Mlor
cd ../
ln -s $rawdata/Illunima_Seq/Mlor/*fastq.gz ./
ln -s $base/02_Annotation/Genepredict/EVM/Mlor/Mlor.gff3
ls *gz|cut -d"d" -f 1 | sort -u |while read id;do
ls -lh ${id}d_1.fq.gz ${id}d_2.fq.gz
hisat2 -p 36 -x ./reference/Mlor -1 ${id}d_1.fq.gz -2 ${id}d_2.fq.gz -S ${id}d.hisat.sam 1>./${id}d.log \
done

ls *.sam|while read id ;do (samtools sort -O bam -@ 24 -o $(basename ${id} ".sam").bam ${id});done
ls *.bam|while read id ;do ( samtools flagstat -@ 24 $id >  $(basename ${id} ".bam").flagstat  );done
ls *.bam|cut -d"d" -f 1 | while read id;do stringtie ${id}d.hisat.bam -eB -G ./Mlor.gff3 -o ${id}d.final.stringTie_asm.gtf -A ${id}d.final.gene_abundence.tab -C ${id}d.final.known.cov_refs.gtf; done
rm *.sam *.bam

prepDE.py -i Mlor_sample_list.txt -g gene_count_matrix.csv -l 150
ls ./*.tab|while read id;do sed -i "s/Gene Name/Gene_Name/g" ${id}; done
ls ./*.tab|while read id;do sed -i "s/Gene ID/Gene_ID/g" ${id}; done
ls ./*.tab|while read id;do awk 'NR > 1 {print $1 "\t" $8 "\t" $9}' ${id} | sort > ${id}.sort; done
awk '{print "mv" " ./" $2 "_good.gene_abundence.tab.sort " $1 "_good.gene_abundence.tab.sort"}' sampleName_clientId.txt > ChangeSampleName.sh
chmod 755 ChangeSampleName.sh
sh ./ChangeSampleName.sh
ls ./*tab.sort |cut -d"_" -f 1 | while read id; do awk -v simplename=$id '{print $1 "\t" $2 "\t" simplename}' ${id}_good.gene_abundence.tab.sort > ${id}_good.gene_abundence.tab.sort_FPKM.csv; done
ls ./*tab.sort |cut -d"_" -f 1 | while read id; do awk -v simplename=$id '{print $1 "\t" $2 "\t" simplename}' ${id}_good.gene_abundence.tab.sort > ${id}_good.gene_abundence.tab.sort_TPM.csv; done
cat *gene_abundence.tab.sort_FPKM.csv > FPKM.txt
cat *gene_abundence.tab.sort_TPM.csv > TPM.txt
Rscript $software/ReshapeMatrix.R

cd Expression/Sfru/reference
ln -s ~/DataBase/Insects_Reference/Noctuoidea/Genomes/Sfru_genome_assembly.fa.masked Sfru.fa
hisat2-build Sfru.fa Sfru
cd ../
ln -s $rawdata/Illunima_Seq/Sfru/*fastq.gz ./
ln -s ~/DataBase/SRA/Sfru/*gz 
ln -s ~/DataBase/Insects_Reference/Noctuoidea/Genomes/Sfru.gff3
ls *gz|cut -d"d" -f 1 | sort -u |while read id;do
ls -lh ${id}d_1.fq.gz ${id}d_2.fq.gz
hisat2 -p 36 -x ./reference/Sfru -1 ${id}d_1.fq.gz -2 ${id}d_2.fq.gz -S ${id}d.hisat.sam 1>./${id}d.log \
done

ls *.sam|while read id ;do (samtools sort -O bam -@ 24 -o $(basename ${id} ".sam").bam ${id});done
ls *.bam|while read id ;do ( samtools flagstat -@ 24 $id >  $(basename ${id} ".bam").flagstat  );done
ls *.bam|cut -d"d" -f 1 | while read id;do stringtie ${id}d.hisat.bam -eB -G ./Sfru.gff3 -o ${id}d.final.stringTie_asm.gtf -A ${id}d.final.gene_abundence.tab -C ${id}d.final.known.cov_refs.gtf; done
rm *.sam *.bam

prepDE.py -i Sfru_sample_list.txt -g gene_count_matrix.csv -l 150
ls ./*.tab|while read id;do sed -i "s/Gene Name/Gene_Name/g" ${id}; done
ls ./*.tab|while read id;do sed -i "s/Gene ID/Gene_ID/g" ${id}; done
ls ./*.tab|while read id;do awk 'NR > 1 {print $1 "\t" $8 "\t" $9}' ${id} | sort > ${id}.sort; done
awk '{print "mv" " ./" $2 "_good.gene_abundence.tab.sort " $1 "_good.gene_abundence.tab.sort"}' sampleName_clientId.txt > ChangeSampleName.sh
chmod 755 ChangeSampleName.sh
sh ./ChangeSampleName.sh
ls ./*tab.sort |cut -d"_" -f 1 | while read id; do awk -v simplename=$id '{print $1 "\t" $2 "\t" simplename}' ${id}_good.gene_abundence.tab.sort > ${id}_good.gene_abundence.tab.sort_FPKM.csv; done
ls ./*tab.sort |cut -d"_" -f 1 | while read id; do awk -v simplename=$id '{print $1 "\t" $2 "\t" simplename}' ${id}_good.gene_abundence.tab.sort > ${id}_good.gene_abundence.tab.sort_TPM.csv; done
cat *gene_abundence.tab.sort_FPKM.csv > FPKM.txt
cat *gene_abundence.tab.sort_TPM.csv > TPM.txt
Rscript $software/ReshapeMatrix.R

################################################################################
#  13 Differential Expression with DESeq2
################################################################################
cd $base/05_RNA-Seq_Analysis/Expression/Msep/
Rscript $software/DESeq2.R

cd $base/05_RNA-Seq_Analysis/Expression/Mlor/
Rscript $software/DESeq2.R

cd $base/05_RNA-Seq_Analysis/Expression/Sfru/
Rscript $software/DESeq2.R
