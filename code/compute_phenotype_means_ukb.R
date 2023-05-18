###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)
library(bigstatsr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--outprefix"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

pheno_dat <- outparse$pheno
test_ids <- outparse$test_ids
output <- outparse$output
seed <- outparse$seed

###Set seed
set.seed(seed)

###Read in data
options(datatable.fread.datatable=FALSE)
test_ids <- fread(test_ids, showProgress=FALSE)
pheno <- readRDS(pheno_dat,)$Y
training_inds_pheno <- which(!(rownames(pheno) %in% test_ids[,2])) ##Get only training individuals
pheno <- pheno[training_inds_pheno, ]

###Compute phenotype means
pheno_means <- colMeans(pheno)

saveRDS(pheno_means, file=output)

