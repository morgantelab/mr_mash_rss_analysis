###Load libraries needed
library(optparse)
library(qgg)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_scores"), type="character")
parser <- add_option(parser, c("--sblub_fit_prefix"), type="character", default=NULL)
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--method"), type="character")
parser <- add_option(parser, c("--meff"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--traits"), type="character", default="-1")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_scores <- outparse$LD_scores
sblub_fit_prefix <- outparse$sblub_fit_prefix
n <- outparse$n
method <- outparse$method
meff <- outparse$meff
output <- outparse$output
seed <- outparse$seed
data_id <- outparse$data_id
traits <- eval(parse(text=outparse$traits))

###Set seed
set.seed(seed)

###Read in data
univ_sumstats <- readRDS(sumstats)
ld_scores <- read.table(gzfile(LD_scores), header=TRUE, sep="\t")
ldscores <- ld_scores$L2
names(ldscores) <- ld_scores$SNP

###Select traits
if(length(traits)>1){
  univ_sumstats$Bhat <- univ_sumstats$Bhat[, traits]
  univ_sumstats$Shat <- univ_sumstats$Shat[, traits]
  univ_sumstats$n <- univ_sumstats$n[traits]
}

###Compute Z scores
Bhat <- univ_sumstats$Bhat
Z <- univ_sumstats$Bhat / univ_sumstats$Shat
r <- ncol(Z)
p <- nrow(Z)

###Compute genetic parameters
genpar <- ldsc(ldscores=ldscores, z=Z, n=rep(n, r), what="rg")

###If h2 is <= 0, set it to 1e-8
genpar$h2[which(genpar$h2<=0)] <- 1e-8

###Load SBLUP results
if(method=="blup"){
  Bhat_sblup <- matrix(as.numeric(NA), nrow=p, ncol=r)
  
  it <- 0  
  for(i in traits){
    it <- it+1
    fit_sblup <- readRDS(paste0(sblub_fit_prefix, "_trait", i, "_", data_id, ".rds"))
    
    Bhat_sblup[, it] <- fit_sblup$res$Effect
  }
  
  Bhat <- Bhat_sblup
}

###Fit wMT-BLUP
tic <- proc.time()[[3]]

fit_wMT_BLUP <- mtadj(h2=genpar$h2, rg=genpar$rg, b=Bhat, n=rep(n, ncol(Z)), meff=meff, method=method, statistics="b")

toc <- proc.time()[[3]]

fit_wMT_BLUP$elapsed_time <- toc-tic

saveRDS(fit_wMT_BLUP, file=output)
