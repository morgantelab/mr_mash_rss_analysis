###Load libraries needed
library(optparse)
library(mice)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--prop_missing"), type="numeric")
parser <- add_option(parser, c("--prop_missing_by_cases"), type="logical")
parser <- add_option(parser, c("--mechanism"), type="character")
parser <- add_option(parser, c("--pattern"), type="character")
outparse <- parse_args(parser)

pheno_dat <- outparse$pheno
test_ids <- outparse$test_ids
output <- outparse$output
seed <- outparse$seed
prop_missing_pheno <- outparse$prop_missing
prop_missing_by_cases <- outparse$prop_missing_by_cases
missing_mechanism <- outparse$mechanism
pattern <- outparse$pattern

###Set seed
set.seed(seed)

###Get phenotype for training individuals only
options(datatable.fread.datatable=FALSE)
test_ids <- fread(test_ids, showProgress=FALSE, header=FALSE)
pheno <- readRDS(pheno_dat)$Y
training_inds_pheno <- which(!(rownames(pheno) %in% test_ids[,2])) ##Get only training individuals
pheno_train <- pheno[training_inds_pheno, ]

###Assign missing values to the phenotype matrix for the training individuals
if(pattern == "full"){
  missingness_pattern <- as.matrix(do.call(CJ, replicate(ncol(pheno_train), 0:1, FALSE)))
  missingness_pattern <- missingness_pattern[-c(1, nrow(missingness_pattern)),]
} else if(pattern == "simple") {
  missingness_pattern <- NULL
}

pheno_miss <- ampute(data=pheno_train, patterns=missingness_pattern, prop=prop_missing_pheno, 
                     bycases=prop_missing_by_cases, mech=missing_mechanism)$amp

###Add missing pheno matrix for training individuals to original matrix
pheno[training_inds_pheno, ] <- as.matrix(pheno_miss)

###Write out results
saveRDS(list(Y=pheno), file=output)
