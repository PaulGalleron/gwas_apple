# Required packages

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)

if(!require(gplots)) install.packages("gplots")
if(!require(LDheatmap)) install.packages("LDheatmap")
if(!require(genetics)) install.packages("genetics")
if(!require(ape)) install.packages("ape")
if(!require(compiler)) install.packages("compiler")
if(!require(grid)) install.packages("grid")
if(!require(bigmemory)) install.packages("bigmemory")
if(!require(EMMREML)) install.packages("EMMREML")
if(!require(scatterplot3d)) install.packages("scatterplot3d")
if(!require(lme4)) install.packages("lme4")

if(!'multtest'%in% installed.packages()[,"Package"]){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("multtest")
  BiocManager::install("snpStats")
}

source("http://zzlab.net/GAPIT/emma.txt")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/GAPIT.library.R")


######## To use if the LDheatmap package was not properly downloaded

#install.packages("LDheatmap")
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("snpStats","rtracklayer","GenomicRanges","GenomInfoDb","IRanges"))

#library("devtools")
#library("pkgload")
# Install the latest development version from GitHub with
#devtools::install_github("SFUStatgen/LDheatmap")


################################################################################
## Read phenotypic and genotypic data and group information                   ##
################################################################################

# Phenotypic data
pheno <- "../data/BLUP_tavelure_multi_2023_automne.csv"
pheno <- read.table(pheno, header = TRUE,
                    sep = ",", stringsAsFactors = TRUE,
                    check.names = FALSE)

# Select the phenotypic variable to study
pheno <- pheno[, 1:2]

# Genotypic data
geno <-"../data/geno_320k_GDDH13.txt"
geno <- fread(geno, sep = "\t")

# Group information
panel <- "../data/panel.csv"
panel <- read.table(panel, header = TRUE,
                    sep = ",", stringsAsFactors = TRUE,
                    check.names = FALSE)


################################################################################
################## Obtain prepared data for GAPIT ##############################
################################################################################

# Call the function used to prepare data for GWAS analyses with GAPIT

source("../functions/functions_before_gwas.R")

result <- process_data(
  pheno = pheno, # Phenotypic dataset
  geno = geno, # Genotypic dataset
  panel = panel, # Panel with groups
  filter_groups = c("cidre", "dessert_moderne", "dessert_ancien"),# Selected groups
  pheno_percent = 100, # Percentage representing the number of studied genotypes
  # (50% corresponds to half of the genotypes, those with the highest values)
  maf_threshold = 0.05, # maf
  pheno_column = "BLUP_AUDPC_reel" # Name of the phenotypic variable
)

# Extract obtained data
pheno_exit <- result$myY
geno_exit <- result$myG

# Ensure the geno_exit data is in the correct format
geno_exit <- as.data.frame(geno_exit)
colnames(geno_exit) <- NULL


### Save obtained data
file_path_geno <- "../data/geno_480k.txt"
fwrite(geno_exit, file = file_path_geno, sep = "\t", quote = FALSE,
       na = "NA", col.names = FALSE)

file_path_pheno <- "../data/pheno_480k.txt"
fwrite(pheno_exit, file = file_path_pheno, sep = "\t", quote = FALSE,
       na = "NA")


######### Alternative if you don't want to save myG and myY for space reasons
#myY <- result$myY
#myG <- result$myG
#myY <- as.data.frame(myY)
#myG <- as.data.frame(myG)
# Convert columns from the 12th onward to integer
#myG[ , 12:ncol(myG)] <- lapply(myG[ , 12:ncol(myG)], as.integer)

##########################################################

# Read the data

myG <- fread(file_path_geno , sep = "\t", header = FALSE)
myG <- as.data.frame(myG)

myY <- fread(file_path_pheno, sep = "\t")
myY <- as.data.frame(myY)

# Replace values in the second column starting from the second row with NA
# It appears that sometimes when GAPIT doesn't work, the allele coding column needs
# to be replaced by NA, as this column is not necessary
#myG_gapit <- myG
#myG_gapit[-1, 2] <- NA


################################################################################
############################# GWAS with GAPIT ##################################
################################################################################

# Use GAPIT with BLINK to perform GWAS
myGAPIT <- GAPIT(
  Y = myY, # Phenotypic data
  G = myG, # Genotypic data (or myG_gapit)
  model = "BLINK", # Model choice (MLMM, FarmCPU, BLINK)
  PCA.total = 3 # Structure, first PCA components
)


################################################################################
################### New Figures Derived from GAPIT #############################
################################################################################

source("../functions/functions_after_gwas.R")

#to use if GAPIT detects significants snp

# General function that returns multiple figures based on GAPIT outputs
run_all_plots(
  myG = myG, # Genotypic dataset
  myY = myY, # Phenotypic dataset
  model = "BLINK", # Model choice
  trait_pheno = "BLUP_AUDPC_reel", # Studied phenotypic variable
  y_axis_title = "BLUP of AUDPC at Day 28" # Y-axis title
)
# Returns Manhattan plot, box plot, violin plot, a spreadsheet with significant SNPs and genotypes, 
# the chromosome map with marker density, and the chromosome map with significant SNPs

# PCA with groups
plot_pca_scatter("GAPIT.Genotype.PCA.csv", 
                 panel,
                 "Principal Component 1 (2.96%)",
                 "Principal Component 2 (2.07%)")
# Percentages in titles are from the GAPIT-generated figure

# Box plot of significant SNPs with groups
generate_box_plots_with_groups(myG = myG,
                               myY = myY, 
                               panel = panel,
                               model = "BLINK",
                               trait_pheno = "BLUP_AUDPC_reel",
                               y_axis_title = "BLUP of AUDPC")
