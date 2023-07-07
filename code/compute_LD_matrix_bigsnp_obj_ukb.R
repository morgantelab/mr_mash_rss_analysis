###Load libraries needed
library(optparse)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--genetic_map"), type="character", default=NULL)
parser <- add_option(parser, c("--chr"), type="integer")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--window_size"), type="integer", default=0)
parser <- add_option(parser, c("--thr_r2"), type="numeric", default=0)
parser <- add_option(parser, c("--impute_missing"), type="logical", default=FALSE)
parser <- add_option(parser, c("--sparse"), type="logical")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
genetic_map <- outparse$genetic_map
seed <- outparse$seed
chr <- outparse$chr
ncores <- outparse$ncores
window_size <- outparse$window_size
thr_r2 <- outparse$thr_r2
impute_missing <- outparse$impute_missing
sparse <- outparse$sparse
output <- outparse$output

###Set seed
set.seed(seed)

###
if(chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9")){
  chr <- paste0("0", chr)
}

###Read in data
geno <- snp_attach(geno_dat)

p <- geno$genotypes$ncol

###Impute missing geno1types
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

###Set chromosome
if(chr==0){
  inds <- seq_len(p)
} else {
  inds <- which(geno$map$chromosome==chr)
}

###Compute LD
if(sparse){
  if(!is.null(genetic_map)){
    POS <- snp_asGeneticPos(as.integer(geno$map$chromosome), geno$map$physical.pos, 
                            dir=genetic_map, 
                            ncores=ncores)
  } else {
    POS <- geno$map$physical.pos
  }
  
  LD <- snp_cor(Gna=X, ind.col=inds, size=win, infos.pos=POS[inds], thr_r2=thr_r2, ncores=ncores)
} else {
  LD <- snp_cor(Gna=X, ind.col=inds, size=win, thr_r2=thr_r2, ncores=ncores)
}
  
###Save to binary file    
conn <- file(output, "wb")
writeBin(as.numeric(LD), conn)
close(conn)

