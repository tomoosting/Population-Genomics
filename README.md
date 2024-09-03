# Population-Genomics

This repository contains a wide-range of bash en Rscripts to perform population genomic analyses.
I have split this repository up input two sections, genotyping and population analyses.

## 1. Genotyping
Genotyping will start with taking raw Illumina read data (fastq.gz) from individual samples and ultimately produce varient call format (VCF) files which will all the relevent SNP data for population analyses.
This pipeline assumed you have access to reference genome. If you want to assembly one you can have a look at [this repository](https://github.com/tomoosting/ONT_Genome_Assembly)

steps include
1. Quality control
2. read alignment
3. genotyping
4. varient filtering

Following all steps will produce multiple VCF files that can be used for different kinds of analyses.
1. all_sites (both variant and invariant sites)
2. qc (all high quality SNPs)
3. outlier (SNPs under putative selection)
4. neutral (selectively neutral and independantly segregating SNPs)
Each population analyses will mention which dataset to use and why.

## 2. population analyses
Here I will outline a range of different analyses you perform.

analyses include
1. geographic distance (requires GPS coordinates)
2. Heterozygosity estimates
3. Analyse of molecular varience (AMOVA)
4. Principal component analyses (PCA)
5. ADMIXTURE
6. Pairwise genetic differentiation (FST)
7. Isolation-by-distance (IBD)
8. (Relative) migration estmates
9. Identification of hybrids
10. Genomescan

If you end up using any of this code, please cite ...

