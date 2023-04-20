###Load libraries needed
library(optparse)
library(data.table)
library(Rfast)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--test_ids"), type="character")
outparse <- parse_args(parser)

input <- outparse$input
data_id <- outparse$data_id
test_ids <- outparse$test_ids

###Read in data
options(datatable.fread.datatable=FALSE)
pheno <- readRDS(paste0("../data/phenotypes/simulated/", input, "_pheno_", data_id, ".rds"))$Y
test_ids <- fread(paste0("../data/phenotypes/", test_ids, "_", data_id, ".txt"), showProgress=FALSE)

###Keep only training individuals 
training_inds_pheno <- which(!(rownames(pheno) %in% test_ids[,2])) ##Get only training individuals
pheno <- pheno[training_inds_pheno, ]

###Compute phenotypic covariance
covY <- cova(pheno, center=TRUE, large = TRUE)

###Save file
saveRDS(covY, file=paste0("../output/misc/", input, "_phenotypic_cov_", data_id, ".rds"))
