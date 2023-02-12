# ArmywormGenome
Code for the armymorm genome project

# Introduction

This repository holds code that was used in for the genome project of *Mythimna separata* and *Mythimna loreyi*. The details were in following docs and the reference.



# Genome assembly

We documented the pipeline used in the paper to scaffold contigs assembled by Pacbio HiFi reads. 



## Genome survey

Estimating genome size and heterozygosity.

```
01_Assembly/runGCE.sh
```



## Genome assembly

Assemly and  clean the genome. 

```
01_Assembly/runGenomeAssembly.sh
```



## Assembly assessment

Code for quality assessment of assemblies.

```
01_Assembly/runBUSCO_genome_mode.sh
01_Assembly/plotBUSCO_genome_mode.sh
```



# Genome Annotation

## Transposable elements analysis

```
02_Annotation/runTEannalysis.sh
```



## Gene prediction

```
02_Annotation/runGenomeAnnotation.sh
```



## Annotation assessment

```
02_Annotation/runBUSCO_protein_mode.sh
02_Annotation/plotBUSCO_protein_mode.sh
```



# Sex Chromosome Detect

```
03_SexChromosomeDetect/runSexChromosomeDetection.sh
```



# Evolution Analysis

## Phylogenetic analysis

Identify orthologs and orthogroups between different species. Alignment of single copy orthologs. Phylogenetic inference with single copy orthologs. Tree divergence dating. Gene Family Evolution. 

```
04_Evolution/runEvolution.sh
```

## Syntenic analysis

```
04_Evolution/runSyntenicanalysis.Rmd
```



# Transcriptome analysis

```
05_RNA-Seq_Analysis/runRNASeqAnalysis.sh
```




