###Load libraries needed
library(optparse)
library(Matrix)
library(bigsnpr)

###Function to make low rank approximation
make_lowrank_R <- function(R, tol=0.8){
  tol <- 1 - sqrt(tol)
  
  eig <- eigen(R, symmetric = TRUE)
  d <- eig$values
  to_remove <- which(abs(d) < tol)
  
  if(length(to_remove)>0){
    d[to_remove] <- 0
    
    R <- eig$vectors %*% (t(eig$vectors)*d)
    
  } else {
    warning("No eigenvalue is smaller than the threshold set. The original matrix is returned.")
  }    
  
  return(list(R=R, d=d))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--blocks"), type="character", default=NULL)
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--thr_maf"), type="numeric", default=0)
parser <- add_option(parser, c("--thr_r2_eigen"), type="numeric", default=-1)
parser <- add_option(parser, c("--impute_missing"), type="logical", default=FALSE)
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--temp_dir"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
blocks <- outparse$blocks
seed <- outparse$seed
ncores <- outparse$ncores
thr_maf <- outparse$thr_maf
thr_r2_eigen <- outparse$thr_r2_eigen
impute_missing <- outparse$impute_missing
output <- outparse$output
temp_dir <- outparse$temp_dir

###Set seed
set.seed(seed)

###Read in data
tmp <- tempfile(tmpdir=temp_dir)
rds <- snp_readBed2(geno_dat, backingfile=tmp, ncores=ncores)
block_data <- read.table(blocks, header=TRUE, sep="\t")
geno <- snp_attach(rds)

p <- geno$genotypes$ncol

mapp <- geno$map

###Impute missing genotypes
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

###Compute maf and filter out variants by it
if(thr_maf>0){
  maf <- snp_MAF(G=X, ncores=ncores)
  inds <- which(maf>=thr_maf)
} else {
  inds <- seq_len(p)
}

LD_denoised <- vector("list", nrow(block_data))

###Loop over blocks
for(i in 1:nrow(block_data)){

  ###Select variants in the current block and intersect them with those surviving MAF threshold
  block_ind <- which(mapp$physical.pos >= block_data[i, "start"] & mapp$physical.pos < block_data[i, "stop"])
  block_ind <- block_ind[which(block_ind %in% inds)]
  
  if(length(block_ind)==1){
    LD_denoised[[i]] <- matrix(1,1,1)
  } else {
    ###Compute LD
    tmp1 <- tempfile(tmpdir=temp_dir)
    LD <- big_cor(X=X, ind.col=block_ind, backingfile=tmp1)
    
    ###Denoise via eigen decomposition
    if(thr_r2_eigen>0){
      LD1 <- make_lowrank_R(LD[], tol=thr_r2_eigen)$R
      diag(LD1) <- 1
      LD_denoised[[i]] <- LD1
    } else {
      LD_denoised[[i]] <- LD
    }
    
    rm(LD, LD1)
  }
}

LD <- bdiag(LD_denoised)

###Filter map, add allele frequency and save it to a file
map_filtered <- mapp[inds, ]
statz <- bigstatsr::big_colstats(X=X, ind.col=inds, ncores=ncores)
map_filtered$freq <- (statz$sum/nrow(X))/2
write.table(map_filtered, paste0(output, ".map"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  
###Save to binary file    
conn <- file(paste0(output, ".ld.bin"), "wb")
writeBin(as.numeric(LD), conn)
close(conn)

###Save to rds
saveRDS(LD, paste0(output, ".rds"))

###Remove temporary files
file.remove(paste0(tmp, c(".bk", ".rds")))

