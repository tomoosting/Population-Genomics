# Performing population genomic analyses

This directory lists all population genomic analyses available in this repository.
Detailed instructions on how to run each analyses are located in the various subdirectories.
Example scripts that can be modified to fit the users needs are located in the scripts folder. NOTE: most Rmarkdown scripts require the file 'custom_Rfunctions.R' to be present in the same directory.

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
A PCA allows you to test for genetic structure without making an a priori assumption about population structure 
(in contrast, pairwise-FST requires you to group samples). This makes it an easy first step if you don't have any information on genetic structure.

![](Figures/5.PCA.png)

### 6. ADMIXTURE
ADMIXTURE estimates the most likely number op K population present in your dataset. Similar to a PCA, this approach also doesn't make any assumptions about population structure. Ones you have determined the most likely number of populations you can plot the estimated levels of ancestry to each of K populations. Output is visualised using what is commenly referred to as a structure plot. 

![](Figures/6.ADMIXTURE.png)

### 7. Pairwise genetic differentiation (FST)
Estimate the level of genetic differentiation between groups of individuals, usually sample locations or populations. Output is visualised using the heatmap. 

![](Figures/7.PWFST_heatmap.png)

### 8. Isolation-by-distance (IBD)
Test if genetic distance is positively correlated with geographic distance. Here we plot genetic distance (liniarised FST) against geogrphic distance and plot a regresion line, and perform a mantal test to test for IBD.

![](Figures/8.IBD_plot.png)

### 9. Migration Rates
Estimate the rate of migration between populations. I'll present two methods that you can apply using SNP data. 1) BayesAss3 estimates the fraction of a population that comes from another population per generation, 2) divmigrate (Rpackage diveRsity) estimates relative migration rates. 

### 10. Identification of hybrids
NewHybrids estimates the probability of individuals being assigned to a genotype frequency class, including first (F1) and second (F2) generation hybrids. 

### 11. Genomescan (Selectionscan)
Test for genome regions under selecion by estimating a range of summary statistics (fst, dxy, π) along the genome using a sliding window. Output is visuailsed using a Manhattan plot which 

### 12. Genotype-by-environment analysis (GEA)
perform (partial) RDA analyses to test for correlations between genetic variation and environmental factors (e.g. temperature or pH). 