###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)
library(bigstatsr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--outprefix"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical")
parser <- add_option(parser, c("--ncores"), type="integer")
outparse <- parse_args(parser)

pheno_dat <- outparse$pheno
geno_dat <- outparse$geno
test_ids <- outparse$test_ids
outprefix <- outparse$outprefix
data_id <- outparse$data_id
chr <- outparse$chr
impute_missing <- outparse$impute_missing
ncores <- outparse$ncores

###Set seed
set.seed(data_id)

###Read in data
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(paste0("../data/genotypes/", geno_dat, ".fam"), showProgress=FALSE)
test_ids <- fread(paste0("../data/phenotypes/", test_ids, "_", data_id, ".txt"), showProgress=FALSE)
training_inds_geno <- which(!(geno_fam[,2] %in% test_ids[,2])) ##Get only training individuals

pheno <- readRDS(paste0("../data/phenotypes/", pheno_dat, "_", data_id, ".rds"))$Y
training_inds_pheno <- which(!(rownames(pheno) %in% test_ids[,2])) ##Get only training individuals
pheno <- pheno[training_inds_pheno, ]

tmp <- tempfile(tmpdir="/data2/morgante_lab/fabiom/tmp")
rds <- snp_readBed2(paste0("../data/genotypes/", geno_dat, ".bed"), ind.row=training_inds_geno, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

rm(list=c("geno_fam", "test_ids", "training_inds_geno", "training_inds_pheno"))

###Compute phenotype means
pheno_means <- colMeans(pheno)

saveRDS(pheno_means, file=paste0("../output/misc/", outprefix, "_pheno_means_", data_id, ".rds"))

###Compute genotype means
chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

###Impute missing genotypes --> needed because this is the input to mr.mash.rss (via summary stats)
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

p <- geno$genotypes$ncol

for(i in chrs){ 

  if(chr==0){
    inds <- seq_len(p)
  } else {
    inds <- which(geno$map$chromosome==i)
  }
  
  geno_means <- big_colstats(X, ind.col = inds, ncores = ncores)$sum / geno$genotypes$nrow
  
  if(length(chrs)==1 && chrs==0){
    i <- "All"
  }
  
  saveRDS(geno_means, file=paste0("../output/misc/", outprefix, "_chr", i, "_geno_means_", data_id, ".rds"))
}

