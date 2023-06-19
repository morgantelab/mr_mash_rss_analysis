###Load libraries needed
library(optparse)
library(qgg)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--method"), type="character")
parser <- add_option(parser, c("--vg"), type="numeric", default=NULL)
parser <- add_option(parser, c("--vb"), type="numeric", default=NULL)
parser <- add_option(parser, c("--ve"), type="numeric", default=NULL)
parser <- add_option(parser, c("--ssb_prior"), type="numeric", default=NULL)
parser <- add_option(parser, c("--sse_prior"), type="numeric", default=NULL)
parser <- add_option(parser, c("--lambda"), type="numeric", default=NULL)
parser <- add_option(parser, c("--h2"), type="numeric", default=NULL)
parser <- add_option(parser, c("--pi"), type="numeric", default=NULL)
parser <- add_option(parser, c("--updateB"), type="logical")
parser <- add_option(parser, c("--updateE"), type="logical")
parser <- add_option(parser, c("--updateG"), type="logical")
parser <- add_option(parser, c("--updatePi"), type="logical")
parser <- add_option(parser, c("--adjustE"), type="logical")
parser <- add_option(parser, c("--nub"), type="integer", default=4)
parser <- add_option(parser, c("--nue"), type="integer", default=4)
parser <- add_option(parser, c("--nit"), type="integer", default=1000)
parser <- add_option(parser, c("--nburn"), type="integer", default=200)
parser <- add_option(parser, c("--algorithm"), type="integer", default=1)
parser <- add_option(parser, c("--trait"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_matrix <- outparse$LD_matrix
n <- outparse$n
method <- outparse$method
vg <- outparse$vg
vb <- outparse$vb
ve <- outparse$ve
ssb_prior <- outparse$ssb_prior
sse_prior <- outparse$sse_prior
lambda <- outparse$lambda
h2 <- outparse$h2
Pi <- outparse$pi
updateB <- outparse$updateB
updateE <- outparse$updateE
updatePi <- outparse$updatePi 
updateG <- outparse$updateG
adjustE <- outparse$adjustE
nub <- outparse$nub
nue <- outparse$nue
nit <- outparse$nit
nburn <- outparse$nburn
algorithm <- outparse$algorithm
trait <- outparse$trait
output <- outparse$output
seed <- outparse$seed

###Set seed
set.seed(seed)

###Read in data
univ_sumstats <- readRDS(sumstats)
p <- nrow(univ_sumstats$Bhat)
r <- ncol(univ_sumstats$Bhat)
LD <- matrix(readBin(LD_matrix, what="numeric", n=p^2), nrow=p, ncol=p, byrow=TRUE)

###Organize summary stats in the format required by qgg
stat <- data.frame(b=univ_sumstats$Bhat[, trait], seb=univ_sumstats$Shat[, trait], 
                   n=n)

###Fit sBayesX
fit_sbayes <- qgg:::sbayes(stat=stat, LD=LD, method=method, lambda=lambda, vg=vg, vb=vb, ve=ve,  
                          h2=h2, pi=Pi, ssb_prior=ssb_prior, sse_prior=sse_prior, nub=nub, nue=nue,
                          updateB=updateB, updateE=updateE, updatePi=updatePi, updateG=updateG, 
                          adjustE=adjustE, nit=nit, nburn=nburn, algorithm=algorithm) 

saveRDS(fit_sbayes, file=output)


