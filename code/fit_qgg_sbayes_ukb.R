###Load libraries needed
library(optparse)
library(qgg)

###Function needed
gbayes_wrap <- function(bhat, shat, n, R, method, vg=NULL, vb=NULL, ve=NULL, 
                        ssb_prior=NULL, sse_prior=NULL, lambda=NULL, h2=NULL, 
                        pi=NULL, updateB=TRUE, updateE=TRUE, updatePi=TRUE, 
                        updateG=TRUE, adjustE=FALSE, nub=4, nue=4, 
                        nit=1000, nburn=200, algorithm=1){
  
  # Process LD
  LDvalues <- split(R, rep(1:ncol(R), each = nrow(R)))
  LDindices <- lapply(1:ncol(R),function(x) { (1:ncol(R))-1 } )
  
  # Prepare parameter input
  methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
  method <- match(method, methods) - 1
  if( !sum(method%in%c(0:5))== 1 ) stop("Method specified not valid")
  
  m <- length(bhat)
  b <- rep(0,m)
  mask <- rep(FALSE, m)
  
  # Compute XtX, Xty, yty 
  ww <- 1/(shat^2 + bhat/n)
  wy <- bhat*ww
  b2 <- bhat^2
  seb2 <- shat^2
  yy <- (b2 + (stat$n-2)*seb2)*ww
  yy <- median(yy)
  
  n <- as.integer(median(n))
  
  # Fit BLR model
  out <- qgg:::sbayes_sparse(yy=yy,
                             wy=wy,
                             ww=ww,
                             b=b,
                             LDvalues=LDvalues,
                             LDindices=LDindices,
                             n=n,
                             m=m,
                             mask=mask,
                             method=method,
                             vg=vg,
                             vb=vb,
                             ve=ve, 
                             ssb_prior=ssb_prior,
                             sse_prior=sse_prior,
                             lambda=lambda,
                             h2=h2,
                             pi=pi,
                             updateB=updateB,
                             updateE=updateE,
                             updatePi=updatePi, 
                             updateG=updateG,
                             adjustE=adjustE,
                             nub=nub,
                             nue=nue,
                             nit=nit,
                             nburn=nburn,
                             algorithm=algorithm)
  
  return(out)
}


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

###Fit sBayesX
fit_sbayes <- gbayes_wrap(bhat=univ_sumstats$Bhat[, trait], shat=univ_sumstats$Shat[, trait], n=n, R=LD, 
                          method=method, vg=vg, vb=vb, ve=ve, ssb_prior=ssb_prior, sse_prior=sse_prior, 
                          lambda=lambda, h2=h2, pi=Pi, updateB=updateB, updateE=updateE, updatePi=updatePi, 
                          updateG=updateG, adjustE=adjustE, nub=nub, nue=nue, nit=nit, nburn=nburn, 
                          algorithm=algorithm) 

saveRDS(fit_sbayes, file=output)


