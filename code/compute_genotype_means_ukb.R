###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)
library(bigstatsr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--chr"), type="integer")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical")
parser <- add_option(parser, c("--ncores"), type="integer")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
test_ids <- outparse$test_ids
output <- outparse$output
seed <- outparse$seed
chr <- outparse$chr
impute_missing <- outparse$impute_missing
ncores <- outparse$ncores

###Set seed
set.seed(seed)

###Read in data
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(paste0("..", unlist(strsplit(geno_dat, ".", fixed=TRUE))[3], ".fam"), showProgress=FALSE)
test_ids <- fread(test_ids, showProgress=FALSE)
training_inds_geno <- which(!(geno_fam[,2] %in% test_ids[,2])) ##Get only training individuals

tmp <- tempfile(tmpdir="/data2/morgante_lab/fabiom/tmp")
rds <- snp_readBed2(geno_dat, ind.row=training_inds_geno, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

rm(list=c("geno_fam", "test_ids", "training_inds_geno"))

###Impute missing genotypes --> needed because this is the input to mr.mash.rss (via summary stats)
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

p <- geno$genotypes$ncol

if(chr==0){
  inds <- seq_len(p)
  chr <- "All"
} else {
  inds <- which(geno$map$chromosome==chr)
}
  
geno_means <- big_colstats(X, ind.col = inds, ncores = ncores)$sum / geno$genotypes$nrow
  
saveRDS(geno_means, file=output)

###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))
