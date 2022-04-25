#folder containing RNAseq data saved from GTEx database
RAW_file_source <- "~/Adipose_files/GTEx_v8_tissue_subsets/"
#output folder
output_overall <- "~/Adipose_files/RESULTS/"
#Tissue files, Either the Adipose Subcutaneous or Adipose Visceral data, saved from the GTEx database as .csv files.
tissue_file <- "Adipose - Subcutaneous.csv"
# tissue_file <- "Adipose - Visceral (Omentum).csv"
#A list of Gene names to run correlations for.
seeds_file <- c("ADIPOQ","LIPE","PLIN1",
                "FKBP10","COL6A1","COL6A2",
                "UPK3B","MSLN","KRT19",
                "MMRN2","ESAM","CDH5",
                "KCNMB1","CNN1","MYH11",
                "CD68","C1QC","FCER1G",
                "CSF3R","FCGR3B","CXCR2",
                "CPA3","TPSB2","TPSAB1",
                "TRBC2","CD6","CD3E",
                "IGKC","JCHAIN","MZB1")
#Reference file containing ENSEMBL IDs and their equivilent Genenames
ENSEMBL <- "ENSEMBL_Genes_V8"
#Correlation method
CorrType <- "spearman"


run_all_corr <- function(tissue_file, seeds_file, RAW_file_source,  ENSEMBL, output_overall, CorrType){
# 0. Load in libraries
library(readr)
library(dplyr)
library(psych)
library(tools)
# 1. Reading in files --------------------------------------------------------
#Read datafile from "Raw untouched files" folder
#make seeds file in excel, just a list of the gene names you want to calculate correlations againast.
print("Reading dataset file")
dataset <- read_csv(paste0(RAW_file_source, tissue_file))
print("Reading seeds")
seeds <- seeds_file
print("reading ENSEMBL_lookup file")
ENSEMBL_key <- read.csv(paste0(RAW_file_source,ENSEMBL, sep=";"))


# 2. Applying functions to data ----------------------------------------------
print("making output folder")
output_directory <- paste0(output_overall, "/", "Correlation_", file_path_sans_ext(basename(seeds_file)), "/")
dir.create(path = output_directory)

#Remove names and labels from samples
print("remove phenotype data")
remove_phen <- function(x){
  subset(x, select = -c(X1, Name, Tissue.x))
}
dataset22 <- remove_phen(dataset)
#Creating full list of ENSEMBL IDs and gene names for recombining with data later
print("make genelist")
make_genelist <- function(d1, ENSEMBL_key){
  genesX <- data.frame(colnames(d1))
  genesX <- data.frame(genesX[-c(1:3),])
  colnames(genesX) <- "Name"
  left_join(genesX, ENSEMBL_key)
}
genelist <- make_genelist(dataset, ENSEMBL_key)
#Convert seeds from gene names into ENSIGID
print("seeds2_FPKM")
genes_to_ENSEMBL <- function(d, ensm, s){
  colnames(s) <- "GENENAME"
  seeds2 <- left_join(s, ensm)
  seeds2$GENENAME <- NULL
  seeds2$ENSIG_ID <- NULL
  d[,colnames(d) %in% seeds2$Name]
}
seeds2_FPKM <- genes_to_ENSEMBL(dataset22, ENSEMBL_key, seeds)
#### calculating correlation values ####
#Define the correlation functions ####
correlation_frame <- function(d, s){
  corr_data_fdr <- corr.test(d, s, method=CorrType, adjust="fdr", ci=FALSE)
}
make_corr_table <- function(d){
  corr_data_corrValues <- data.frame(d$r)
  setDT(corr_data_corrValues, keep.rownames = TRUE)[]
  setnames(corr_data_corrValues, 1, "Name")
}
make_FDR_table <- function(d){
  corr_data_corrValues <- data.frame(d$p)
  setDT(corr_data_corrValues, keep.rownames = TRUE)[]
  setnames(corr_data_corrValues, 1, "Name")
}
#Running the functions ####
print("calculating correlations")
Merged_Corr <- correlation_frame(dataset22, seeds2_FPKM)
print("writing data to file")
write_delim(right_join(genelist, make_corr_table(Merged_Corr)), file =  paste0(output_directory, "FULL_CORR_", file_path_sans_ext(basename(tissue_file)), ".csv"), delim = ",")
write_delim(right_join(genelist, make_FDR_table(Merged_Corr)), file =  paste0(output_directory, "FULL_FDR_", file_path_sans_ext(basename(tissue_file)), ".csv"), delim = ",")
# Create metadata files for full dataset---------------------------------------------------
print("calculating mean TPMs")
metadata_tissue <- data_frame(
  "Name" = c(colnames(dataset22)), 
  "TPM_Mean" = c(colMeans(dataset22)), 
  "TPM_SD" = c(sapply(dataset22, sd)), 
  "TPM_MAX" = c(sapply(dataset22, max)), 
  "TPM_MIN" = c(sapply(dataset22, min)), 
  "TPM=0" = c(sapply(dataset22, function(x) sum(x == 0))), 
  "TPM<0.1" = c(sapply(dataset22, function(x) sum(x < 0.1))),
  "n_samples" = c(nrow(dataset22))
  )
#add percentages for 0 values
metadata_tissue$"TPM=0_percent" <- 100*metadata_tissue$"TPM=0"/metadata_tissue$n_samples
metadata_tissue$"TPM<0.1_percent" <- 100*metadata_tissue$"TPM<0.1"/metadata_tissue$n_samples
metadata_tissue$"Coeff_Variation" <- metadata_tissue$TPM_SD/metadata_tissue$TPM_Mean
metadata_tissue <- metadata_tissue[c(1:5,11,6:10)]
metadata_tissue <- right_join(genelist, metadata_tissue)
write_delim(metadata_tissue, file = paste0(output_directory, "FULL_METADATA_", file_path_sans_ext(basename(tissue_file)), ".csv"), delim = ",")

}

run_all_corr(tissue_file, seeds_file, RAW_file_source, ENSEMBL, output_overall, CorrType)
