###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--window_size"), type="integer", default=0)
parser <- add_option(parser, c("--thr_r2"), type="numeric", default=0)
parser <- add_option(parser, c("--impute_missing"), type="logical")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
test_ids <- outparse$test_ids
data_id <- outparse$data_id
chr <- outparse$chr
ncores <- outparse$ncores
window_size <- outparse$window_size
thr_r2 <- outparse$thr_r2
impute_missing <- outparse$impute_missing


###Set seed
set.seed(data_id)

###Read in data
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(paste0("../data/genotypes/", geno_dat, ".fam"), showProgress=FALSE)
test_ids <- fread(paste0("../data/phenotypes/", test_ids, "_", data_id, ".txt"), showProgress=FALSE)
training_inds <- which(!(geno_fam[,2] %in% test_ids[,2])) ##Get only training individuals
rm(list=c("geno_fam", "test_ids"))

rds <- snp_readBed2(paste0("../data/genotypes/", geno_dat, ".bed"), ind.row=training_inds, ncores=ncores)
geno <- snp_attach(rds)

p <- geno$genotypes$ncol

###Impute missing genotypes
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

###Set window size
if(window_size==0){
  win <- p
} else {
  win <- window_size
}

###Compute LD
chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

##Loop over chromosomes
for(i in chrs){ 
  
  #Set chromosome
  if(i==0){
    inds <- seq_len(p)
  } else {
    inds <- which(geno$map$chromosome==i)
  }
  
  #Compute LD
  LD <- snp_cor(Gna=X, ind.col=inds, size=win, thr_r2=thr_r2, ncores=ncores)

  #Save to binary file    
  if(length(chrs)==1 && chrs==0){
    i <- "All"
  }
  conn <- file(paste0("../data/LD_matrices/", geno_dat, "_chr", i, "_LD_", data_id, ".ld.bin"), "wb")
  writeBin(as.numeric(LD), conn)
  close(conn)
}
