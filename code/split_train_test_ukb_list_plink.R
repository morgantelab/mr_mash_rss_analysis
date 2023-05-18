###Load libraries need
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--n_train"), type="integer")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

input <- outparse$input
output <- outparse$output
n_train <- outparse$n_train
seed <- outparse$seed

###Set seed
set.seed(seed)

###Read in data
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(input, showProgress=FALSE)

###Spilt train-test
idx_train <- sample(x=1:nrow(geno_fam), size=n_train)

test_ids <- geno_fam[-idx_train, 1:2]

###Write out the data
fwrite(test_ids, file=output, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
