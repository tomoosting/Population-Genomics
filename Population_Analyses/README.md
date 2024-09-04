# Performing population genomic analyses

This directory lists all population genomic analyses available in this repository.
Detailed instructions on how to run each analyses are located in the various subdirectories.
Example scripts that can be modified to fit the users needs are located in the scripts folder. NOTE: most Rmarkdown scripts require the file 'custom_script.R' to be present in the same directory.

## file types
Besides a VCF file, we will encounter the following file types
1. gds
2. genlight
3. plink
4. genepop

These can all be created form a VCF file as they contain the same information (genotypes). Please see the convert scipts to generate the various file formats.
I highly recommand taking some time to understand how a VCF file is constructed because it will allow you start playing with the data yourself much faster. NOTE: all of the other file formats only contain the FORMAT/GT field stored in a VCF file, representing the individual genotypes for all samples and loci. 

These analyses generally assume you have have a file named 'sample_information.tsv' that stores metadata for individual samples. Most common fields include IND,SITE,POP,LONG,LATH reprsenting sample name, sample site, population, longitude and latitude. To present mistakes it is best to sort this file in the same order as samples are listed in your VCF file (genrally alphabetical).

Below is a short summary for each analyses including which datasets, file types, and metadata you need to perform the analyses.

## Analyses included in the repository 
### 1. geographic distance (requires GPS coordinates)
metadata:IND,SITE,LATH,LONG

dataset:(all_sites,qc,outlier,neutral)

filetype:(VCF,gds,genlight,plink,genepop)

An Rmarkdown script that allows you to estimate the distannce (km) between provided GSP locations (does not require genetic data). This can be either from sample locations, or for indidual samples. You will need this to estimate isolation-by-distance and GEA analysis. The output will be a distance matrix.

### 2. Estimate sequencing depth statistics
filetype: bam files
Calculate depth statistics to, this is actually run on bam files individual samples. 
This can be useful to check whether some samples have a much higher coverage, which can create biases when estimating genetic diveristy metrics such as heterozygosity. 

### 3. Heterozygosity estimates
Estimate levels of heterozygosity for each sample in your dataset and test for significnat differences between sample locations or populations.

### 4. Analyse of molecular varience (AMOVA)

### 5. Principal component analyses (PCA)

### 6. ADMIXTURE

![](figures/snapper_norm_q2_ggplot.png)

### 7. Pairwise genetic differentiation (FST)

### 8. Isolation-by-distance (IBD)

### 9. (Relative) migration estmates

### 10. Identification of hybrids

### 11. Genomescan

### 12. Genotype-by-environment analysis (GEA)