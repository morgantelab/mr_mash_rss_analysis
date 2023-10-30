###Load libraries needed
library(optparse)
library(data.table)
library(Rfast)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--normalize"), type="logical")
outparse <- parse_args(parser)

input <- outparse$pheno
output <- outparse$output
seed <- outparse$seed
test_ids <- outparse$test_ids
normalize <- outparse$normalize

set.seed(seed)

###Read in data
options(datatable.fread.datatable=FALSE)
pheno <- readRDS(input)$Y
test_ids <- fread(test_ids, showProgress=FALSE, header=FALSE)

###Keep only training individuals 
training_inds_pheno <- which(!(rownames(pheno) %in% test_ids[,2])) ##Get only training individuals
pheno <- pheno[training_inds_pheno, ]

if(normalize){
  ###Function for quantile normalization
  inv_normalise <- function(x) { #this would also tolerate NAs
    return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
  }
  
  pheno <- apply(pheno, 2, inv_normalise)
}

###Compute phenotypic covariance
covY <- cova(pheno, center=TRUE, large = TRUE)

###Save file
saveRDS(covY, file=output)
