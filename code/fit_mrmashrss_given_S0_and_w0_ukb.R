###Load libraries needed
library(optparse)
library(mr.mash.alpha)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--pheno_cov"), type="character")
parser <- add_option(parser, c("--residual_cov"), type="character")
parser <- add_option(parser, c("--prior"), type="character", default=NULL)
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--w0_null_scale"), type="numeric", default=0)
parser <- add_option(parser, c("--update_w0"), type="logical", default=FALSE)
parser <- add_option(parser, c("--update_w0_max_iter"), type="numeric", default=Inf)
parser <- add_option(parser, c("--tol"), type="numeric", default=1e-2)
parser <- add_option(parser, c("--check_R"), type="logical", default=FALSE)
parser <- add_option(parser, c("--max_iter"), type="integer", default=5000)
parser <- add_option(parser, c("--standardize"), type="logical", default=TRUE)
parser <- add_option(parser, c("--verbose"), type="logical", default=TRUE)
parser <- add_option(parser, c("--update_V"), type="logical", default=TRUE)
parser <- add_option(parser, c("--update_V_method"), type="character", default="full")
parser <- add_option(parser, c("--w0_threshold"), type="numeric", default=0)
parser <- add_option(parser, c("--ca_update_order"), type="character", default="consecutive")
parser <- add_option(parser, c("--X_colmeans"), type="character", default=NULL)
parser <- add_option(parser, c("--Y_colmeans"), type="character", default=NULL)
parser <- add_option(parser, c("--mu1_init"), type="character", default=NULL)
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--traits"), type="character", default="-1")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_matrix <- outparse$LD_matrix
pheno_cov <- outparse$pheno_cov
residual_cov <- outparse$residual_cov
prior <- outparse$prior
n <- outparse$n
update_w0 <- outparse$update_w0
update_w0_max_iter <- outparse$update_w0_max_iter
w0_null_scale <- outparse$w0_null_scale
tol <- outparse$tol
check_R <- outparse$check_R
max_iter <- outparse$max_iter
standardize <- outparse$standardize
verbose <- outparse$verbose
update_V <- outparse$update_V
update_V_method <- outparse$update_V_method
w0_threshold <- outparse$w0_threshold
ca_update_order <- outparse$ca_update_order
X_colmeans <- outparse$X_colmeans
Y_colmeans <- outparse$Y_colmeans
mu1_init <- outparse$mu1_init
ncores <- outparse$ncores
output <- outparse$output
traits <- eval(parse(text=outparse$traits))
seed <- outparse$seed

###Set seed
set.seed(seed)

###Read in data
univ_sumstats <- readRDS(sumstats)
p <- nrow(univ_sumstats$Bhat)
r <- ncol(univ_sumstats$Bhat)
if(tail(unlist(strsplit(LD_matrix, ".", fixed = TRUE)), 1) == "bin"){
  LD <- matrix(readBin(LD_matrix, what="numeric", n=p^2), nrow=p, ncol=p, byrow=TRUE)
} else if(tail(unlist(strsplit(LD_matrix, ".", fixed = TRUE)), 1) == "rds"){
  LD <- readRDS(LD_matrix)
}
covY <- readRDS(pheno_cov)
V <- readRDS(residual_cov)
prior <- tryCatch(readRDS(prior), 
                     error = function(e) {
                       return(NULL)
                     },
                     warning = function(w) {
                       return(NULL)
                     }
)
mu1_init <- tryCatch(readRDS(mu1_init), 
                      error = function(e) {
                        return(NULL)
                      },
                      warning = function(w) {
                        return(NULL)
                      }
)
X_colmeans <- tryCatch(readRDS(X_colmeans), 
                     error = function(e) {
                       return(NULL)
                     },
                     warning = function(w) {
                       return(NULL)
                     }
)
Y_colmeans <- tryCatch(readRDS(Y_colmeans), 
                     error = function(e) {
                       return(NULL)
                     },
                     warning = function(w) {
                       return(NULL)
                     }
)

###Get the prior covariances and weights
S0 <- prior$S0
w0 <- prior$w0

if(w0_null_scale > 0){
  w0[1] <- sum(w0[-1])*w0_null_scale
  w0 <- w0/sum(w0)
}

###Subsets traits
if(length(traits)>1){
  univ_sumstats$Bhat <- univ_sumstats$Bhat[, traits]
  univ_sumstats$Shat <- univ_sumstats$Shat[, traits]
  univ_sumstats$n <- univ_sumstats$n[traits]
  S0 <- lapply(S0, function(x, sel){x[sel, sel]}, traits)
  covY <- covY[traits, traits]
  V <- V[traits, traits]
  if(!is.null(mu1_init)){
    mu1_init <- mu1_init[, traits]
  }
  if(!is.null(Y_colmeans)){
    Y_colmeans <- Y_colmeans[traits]
  }
  
}

###Fit mr.mash.rss
tic <- proc.time()[[3]]
fit_mrmash_rss <- mr.mash.rss(Bhat=univ_sumstats$Bhat, Shat=univ_sumstats$Shat, covY=covY, R=LD, n=n, S0=S0, 
                              w0=w0, update_w0=update_w0, tol=tol, V=V, check_R=check_R, max_iter=max_iter,
                              convergence_criterion="ELBO", compute_ELBO=TRUE, standardize=standardize,
                              verbose=verbose, update_V=update_V, update_V_method=update_V_method, update_w0_max_iter=update_w0_max_iter,
                              w0_threshold=w0_threshold, ca_update_order=ca_update_order, mu1_init=mu1_init,
                              X_colmeans=X_colmeans, Y_colmeans=Y_colmeans, nthreads=ncores)
toc <- proc.time()[[3]]

fit_mrmash_rss$elapsed_time <- toc-tic


saveRDS(fit_mrmash_rss, file=output)


