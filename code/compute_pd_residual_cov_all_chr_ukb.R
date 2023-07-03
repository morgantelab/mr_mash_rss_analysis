###Load libraries needed
library(optparse)
library(Matrix)
library(Rfast)

###Functions needed
is_strong <- function(x){
  strong <- any(abs(x) > 2)
  return(strong)
}

makepd <- function(mat=NULL, tol=0.0001) {
  rn <- rownames(mat)
  e <- eigen(mat)
  U <- e$vectors
  e <- e$values
  e[e<tol] <- tol #set eigenvalues smaller than tol to tol
  D <- diag(e)
  mat <- U%*%D%*%t(U)
  colnames(mat) <- rownames(mat) <- rn
  
  return(list(mat=mat, evals_adj=e))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats_prefix"), type="character")
parser <- add_option(parser, c("--sumstats_suffix"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

sumstats_prefix <- outparse$sumstats_prefix
sumstats_suffix <- outparse$sumstats_suffix
chr <- outparse$chr
output <- outparse$output
seed <- outparse$seed

set.seed(seed)

chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

###Extract weak effects
for(i in chrs){
  dat <- readRDS(paste0(sumstats_prefix, i, sumstats_suffix))
  Z <- dat$Bhat/dat$Shat
  strong <- which(apply(Z, 1, is_strong))
  
  if(i==1){
    Z_weak <- Z[-strong, ]
  } else {
    Z_weak <- rbind(Z_weak, Z[-strong, ])
  }
  
}

###Estimate V
Vhat <- matrix(0, nrow=ncol(Z_weak), ncol=ncol(Z_weak))

for(j in 1:nrow(Z_weak)){
  Vhat <- Vhat + tcrossprod(Z_weak[j,])
}

Vhat <- Vhat/nrow(Z_weak)

Vhat <- makepd(mat=Vhat, tol=1e-6)$mat

if(!is.symmetric(Vhat)){
  Vhat <- as.matrix(forceSymmetric(Vhat))
}

###Save file
saveRDS(Vhat, file=output)
