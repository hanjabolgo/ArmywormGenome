################################################################################
# Gene Annotation
################################################################################
base=~/Project/Armyworms/
rawdata=~/Project/Armyworms/00_RawData/
software=~/biosoft
scripts=~/scripts/

cd $base/02_Annotation
mkdir -p Genepredict/BRAKER Genepredict/Exonerate Genepredict/RNA-Seq  Genepredict/EVM

################################################################################
# 7 Gene prediction
################################################################################

################################################################################
# 7.1 BRAKER2
################################################################################
cd $base/02_Annotation/Genepredict/BRAKER
mkdir -p Msep/AlignIlluminaSeq Msep/AlignIsoseq Mlor/AlignIlluminaSeq Mlor/Align/AlignIsoseq

cd $base/02_Annotation/Genepredict/BRAKER/Msep/AlignIlluminaSeq
ln -s $rawdata/Illunima_Seq/Msep/*fastq.gz ./
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked 
$software/hisat2/bin/hisat2-build Msep MsepF_genome_assembly.fa.masked
paste <(ls -1 *1.fastq.gz ) <(ls -1 *2.fastq.gz) |awk '{print "runAlignConvertSortStats.sh",$1,$2,"~/Project/Armyworms/02_Annotation/Genepredict/BRAKER/Msep/AlignIlluminaSeq Msep"}' >align.sh
chmod 755 align.sh
./align.sh

# Check bam files and remove the bam file with a "PCT_PF_READS_IMPROPER_PAIRS" less than 20.
cat *.bam_alignment.stats > all.bam_alignment.stats
samtools short_reads_all.bam *.bam
stringtie -o transcript.gtf short_reads_all.bam
$software/transdecoder/util/gtf_to_alignment_gff3.pl transcript.gtf > transcript.gff3
$software/gffread/gffread Illumina.gff3 -g MsepF_genome_assembly.fa.masked -M -x transcript_cds.fa -y transcript_protein.fa

cd $base/02_Annotation/Genepredict/BRAKER/Msep/AlignIsoseq
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked
ln -s $rawdata/Iso-Seq/Msep/polished.hq.fasta.gz
gmap_build -d Msep -D . MsepF_genome_assembly.fa.masked
runGmap.sh Msep $base/02_Annotation/Genepredict/BRAKER/Msep/AlignIsoseqTranscripts MsepF_genome_assembly.fa.masked polished.hq.fasta

cd ../
ln -s AlignIsoseqTranscripts/*.bam
ln -s AlignIlluminaSeq/*.bam
ln -s ~/DataBase/Insects_Reference/Noctuoidea/NCBI_transcripts/Noctuoidea_trans.fa
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked
samtools merge all.bam *.bam
braker.pl --cores 40 --species=Msep --genome=MsepF_genome_assembly.fa.masked \
     --softmasking --bam=all.bam \
     --gff3 \
     --prot_seq=Noctuoidea_trans.fa --prg=exonerate \ 

cd $base/02_Annotation/Genepredict/BRAKER/Mlor/AlignIlluminaSeq
ln -s $rawdata/Illunima_Seq/Msep/*fastq.gz ./
ls *gz|while read gz;do gzip -t $gz;done
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked
$software/hisat2/bin/hisat2-build Mlor MlorF_genome_assembly.fa.masked
paste <(ls -1 *1.fastq.gz ) <(ls -1 *2.fastq.gz) |awk '{print "runAlignConvertSortStats.sh",$1,$2,"~/Project/Armyworms/02_Annotation/Genepredict/BRAKER/Mlor/AlignIlluminaSeq Mlor"}' >align.sh
chmod 755 align.sh
./align.sh
# Check bam files and remove the bam file with a "PCT_PF_READS_IMPROPER_PAIRS" less than 20.
cat *.bam_alignment.stats > all.bam_alignment.stats
samtools short_reads_all.bam *.bam
stringtie -o transcript.gtf short_reads_all.bam
$software/transdecoder/util/gtf_to_alignment_gff3.pl transcript.gtf > transcript.gff3
$software/gffread/gffread Illumina.gff3 -g MlorF_genome_assembly.fa.masked -M -x transcript_cds.fa -y transcript_protein.fa

cd $base/02_Annotation/Genepredict/BRAKER/Mlor/AlignIsoseq
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked
ln -s $rawdata/Iso-Seq/Mlor/polished.hq.fasta.gz
gmap_build -d Msep -D . MlorF_genome_assembly.fa.masked
runGmap.sh Msep $base/02_Annotation/Genepredict/BRAKER/Mlor/AlignIsoseqTranscripts MlorF_genome_assembly.fa.masked polished.hq.fasta

cd ../
ln -s AlignIsoseqTranscripts/*.bam
ln -s AlignIlluminaSeq/*.bam
ln -s ~/DataBase/Insects_Reference/Noctuoidea/NCBI_transcripts/Noctuoidea_trans.fa
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked
samtools merge all.bam *.bam
braker.pl --cores 40 --species=Mlor --genome=MlorF_genome_assembly.fa.masked \
     --softmasking --bam=all.bam \
     --gff3 \
     --prot_seq=Noctuoidea_trans.fa --prg=exonerate \

################################################################################
# 7.2 Exonerate
################################################################################
cd $base/02_Annotation/Genepredict/Exonerate
mkdir Msep Mlor

cd Msep
ln -s ~/DataBase/Insects_Reference/Noctuoidea/NCBI_transcripts/Noctuoidea_trans.fa
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked
$software/exonerate-2.4.0/bin/exonerate -t MsepF_genome_assembly.fa.masked -q Noctuoidea_trans.fa --model protein2genome --querytype protein --targettype dna --showvulgar no --softmaskquery yes --softmasktarget yes --minintron 20 --maxintron 100000 --showalignment no --showtargetgff yes --showcigar no --geneseed 250 --score 250 --verbose 0 --gff3 no > exonerate.out
perl $software/EVidenceModeler/EvmUtils/misc/Exonerate_to_evm_gff3.pl exonerate.out > exonerate.gff3

cd ../Mlor
ln -s ~/DataBase/Insects_Reference/Noctuoidea/NCBI_transcripts/Noctuoidea_trans.fa
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked
$software/exonerate-2.4.0/bin/exonerate -t MlorF_genome_assembly.fa.masked -q Noctuoidea_trans.fa --model protein2genome --querytype protein --targettype dna --showvulgar no --softmaskquery yes --softmasktarget yes --minintron 20 --maxintron 100000 --showalignment no --showtargetgff yes --showcigar no --geneseed 250 --score 250 --verbose 0 --gff3 no > exonerate.out
perl $software/EVidenceModeler/EvmUtils/misc/Exonerate_to_evm_gff3.pl exonerate.out > exonerate.gff3

################################################################################
# 7.3 Transcript
################################################################################
cd $base/02_Annotation/Genepredict/RNA-Seq
mkdir Msep Mlor

cd Msep
ln -s $base/02_Annotation/Genepredict/BRAKER/Msep/AlignIlluminaSeq/transcript_cds.fa
ln -s $rawdata/Iso-Seq/Msep/polished.hq.fasta.gz
ln -s $base/02_Annotation/Genepredict/BRAKER/Msep/AlignIlluminaSeq/transcripts.gff3
gzip -d polished.hq.fasta.gz
cat transcript_cds.fa polished.hq.fasta > transcript.fa
$software/TransDecoder/TransDecoder.LongOrfs -t transcript.fa -m 100
$software/TransDecoder/TransDecoder.Predict -t transcript.fa
$software/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

cd ../Mlor
ln -s $base/02_Annotation/Genepredict/BRAKER/Mlor/AlignIlluminaSeq/transcript_cds.fa
ln -s $rawdata/Iso-Seq/Mlor/polished.hq.fasta.gz
ln -s $base/02_Annotation/Genepredict/BRAKER/Mlor/AlignIlluminaSeq/transcripts.gff3
gzip -d polished.hq.fasta.gz
cat transcript_cds.fa polished.hq.fasta > transcript.fa
$software/TransDecoder/TransDecoder.LongOrfs -t transcript.fa -m 100
$software/TransDecoder/TransDecoder.Predict -t transcript.fa
$software/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

################################################################################
# 7.4 EVM
################################################################################
cd $base/02_Annotation/Genepredict/EVM
mkdir Msep Mlor

cd Msep 
ln -s $base/02_Annotation/Genepredict/BRAKER/Msep/braker.gff3 
ln -s $base/02_Annotation/Genepredict/Exonerate/Msep/exonerate.gff3
ln -s $base/02_Annotation/Genepredict/RNA-Seq/Msep/transcripts.fasta.transdecoder.genome.gff3 transcript.gff3
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked
$software/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome MsepF_genome_assembly.fa.masked \
        --gene_predictions braker.gff3 \
        --protein_alignments exonerate.gff3 \
        --transcript_alignments transcript.gff3 \
        --segmentSize 100000 --overlapSize 10000 \
        --partition_listing partitions_list.out

$software/EVidenceModeler/EvmUtils/write_EVM_commands.pl --genome MsepF_genome_assembly.fa.masked --weights `pwd`/weights.txt \
      --gene_predictions braker.gff3 \
      --protein_alignments exonerate.gff3 \
      --transcript_alignments Transcript.gff3 \
      --output_file_name evm.out  \
      --partitions partitions_list.out >  commands.list

parallel --jobs 10 < commands.list

$software/EVidenceModeler/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out \
         --output_file_name evm.out
$software/EVidenceModeler/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out
         --output evm.out  --genome MsepF_genome_assembly.fa.masked
find . -regex ".*evm.out.gff3" -exec cat {} \; | bedtools sort -i - > EVM.all.gff
$software/gffread/gffread EVM.all.gff -g MsepF_genome_assembly.fa.masked -y tr_cds.fa
bioawk -c fastx '$seq < 50 {print $comment}' tr_cds.fa | cut -d '=' -f 2 > short_aa_gene_list.txt
grep -v -w -f short_aa_gene_list.txt EVM.all.gff > filter.gff
grep -v -w -f short_aa_gene_list.txt EVM.all.gff > filter.gff
python $script/sort_EVM.py filter.gff filter.gff EVM_sort.gff
python $script/rename_gff.py -g EVM_sort.gff -c Mlor_key.list
mv result.rename.gff3 Msep.gff3
# Then we corrected "Msep.gff3" by apollo manully

cd ../Mlor
ln -s $base/02_Annotation/Genepredict/BRAKER/Mlor/braker.gff3
ln -s $base/02_Annotation/Genepredict/Exonerate/Mlor/exonerate.gff3
ln -s $base/02_Annotation/Genepredict/RNA-Seq/Mlor/transcripts.fasta.transdecoder.genome.gff3 transcript.gff3
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked
$software/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome MlorF_genome_assembly.fa.masked \
        --gene_predictions braker.gff3 \
        --protein_alignments exonerate.gff3 \
        --transcript_alignments transcript.gff3 \
        --segmentSize 100000 --overlapSize 10000 \
        --partition_listing partitions_list.out

$software/EVidenceModeler/EvmUtils/write_EVM_commands.pl --genome MlorF_genome_assembly.fa.masked --weights `pwd`/weights.txt \
      --gene_predictions braker.gff3 \
      --protein_alignments exonerate.gff3 \
      --transcript_alignments Transcript.gff3 \
      --output_file_name evm.out  \
      --partitions partitions_list.out >  commands.list

parallel --jobs 10 < commands.list

$software/EVidenceModeler/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out \
         --output_file_name evm.out
$software/EVidenceModeler/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out
         --output evm.out  --genome MlorF_genome_assembly.fa.masked
find . -regex ".*evm.out.gff3" -exec cat {} \; | bedtools sort -i - > EVM.all.gff
$software/gffread/gffread EVM.all.gff -g MlorF_genome_assembly.fa.masked -y tr_cds.fa
bioawk -c fastx '$seq < 50 {print $comment}' tr_cds.fa | cut -d '=' -f 2 > short_aa_gene_list.txt
grep -v -w -f short_aa_gene_list.txt EVM.all.gff > filter.gff
python $script/sort_EVM.py filter.gff filter.gff EVM_sort.gff
python $script/rename_gff.py -g EVM_sort.gff -c Mlor_key.list
mv result.rename.gff3 Mlor.gff3
# Then we corrected "Mlor.gff3" by apollo manully

################################################################################
# 7.5 ncRNA prediction
################################################################################
cd $base/02_Annotation/ncRNA/
mkdir -p RNAmmer/Msep RNAmmer/Mlor tRNAscan-SE/Msep tRNAscan-SE/Mlor

cd RNAmmer/Msep
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked
$software/RNAmmer-1.2/bin/rnammer -S euk -m lsu,ssu,tsu -h Msep_hmmreport -gff Msep_rRNA.gff -f Msep_rRNA.frn MsepF_genome_assembly.fa.masked
$base/tRNAscan-SE-1.3.1/tRNAscan-SE --thread 20 -o Msep_tRNA.out -f DFNC_F_tRNA.ss -m Msep_tRNA.stats MsepF_genome_assembly.fa.masked

cd ../Mlor
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked
$software/RNAmmer-1.2/bin/rnammer -S euk -m lsu,ssu,tsu -h Mlor_hmmreport -gff Mlor_rRNA.gff -f Mlor_rRNA.frn MlorF_genome_assembly.fa.masked
$base/tRNAscan-SE-1.3.1/tRNAscan-SE --thread 20 -o Msep_tRNA.out -f DFNC_F_tRNA.ss -m Msep_tRNA.stats MlorF_genome_assembly.fa.masked

################################################################################
# 8 Functional annotation
################################################################################
cd ~/Database/NR/
$software/Diamond-2.0.6/bin/diamond makedb --in nr -d Ref_nr --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp
cd ~/DataBasw/Uniprot/
$software/Diamond-2.0.6/bin/diamond -in uniprot_sprot.fasta -d swissport

cd $base/02_Annotation/
mkdir -p Functionannotation/Msep Functionannotation/Mlor
cd Functionannotation/Msep
ln -s $base/02_Annotation/TE_analysis/MsepF_genome_assembly.fa.masked
ln -s $base/02_Annotation/Genepredict/EVM/Msep/Msep.gff3
$software/gffread/gffread Msep.gff3 -g MsepF_genome_assembly.fa.masked -M -x cds.fa -y protein.fa
$software/Diamond-2.0.6/bin/diamond blastp --db ~/DataBase/NR/Ref_nr.dmnd --out Msep2nr.out --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles sphylums --query protein.fa --evalue 1e-10 --threads 40 --ultra-sensitive
python -m jcvi.formats.blast best -n 1 Msep2nr.out
$software/Diamond-2.0.6/bin/diamond blastp --db ~/DataBasw/Uniprot/swissport.dmnd --out Msep2sw.out --outfmt 6 --query protein.fa --evalue 1e-10 --threads 40 --ultra-sensitive
python -m jcvi.formats.blast best -n 1 Msep2sw.out

cd ../Mlor
ln -s $base/02_Annotation/TE_analysis/MlorF_genome_assembly.fa.masked
ln -s $base/02_Annotation/Genepredict/EVM/Mlor/Mlor.gff3
$software/gffread/gffread Mlor.gff3 -g MlorF_genome_assembly.fa.masked -M -x cds.fa -y protein.fa
$software/Diamond-2.0.6/bin/diamond blastp --db ~/DataBase/NR/Ref_nr.dmnd --out Mlor2nr.out --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles sphylums --query protein.fa --evalue 1e-10 --threads 40 --ultra-sensitive
python -m jcvi.formats.blast best -n 1 Mlor2nr.out
$software/Diamond-2.0.6/bin/diamond blastp --db ~/DataBasw/Uniprot/swissport.dmnd --out Mlor2sw.out --outfmt 6 --query protein.fa --evalue 1e-10 --threads 40 --ultra-sensitive
python -m jcvi.formats.blast best -n 1 Mlor2sw.out
