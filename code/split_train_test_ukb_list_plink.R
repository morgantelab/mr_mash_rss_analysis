###Load libraries need
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--n_train"), type="integer")
parser <- add_option(parser, c("--data_id"), type="integer")
outparse <- parse_args(parser)

input <- outparse$input
n_train <- outparse$n_train
data_id <- outparse$data_id

###Set seed
set.seed(data_id)

###Read in data
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(paste0("../data/genotypes/", input, ".fam"), showProgress=FALSE)

###Spilt train-test
idx_train <- sample(x=1:nrow(geno_fam), size=n_train)

test_ids <- geno_fam[-idx_train, 1:2]

###Write out the data
fwrite(test_ids, file=paste0("../data/phenotypes/simulated/", input, "_test_ids_", data_id, ".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
