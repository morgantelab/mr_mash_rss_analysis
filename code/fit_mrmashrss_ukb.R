###Load libraries needed
library(optparse)
library(mr.mash.alpha)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--pheno_cov"), type="character")
parser <- add_option(parser, c("--residual_cov"), type="character")
parser <- add_option(parser, c("--canonical_cov"), type="logical")
parser <- add_option(parser, c("--data_driven_cov"), type="character", default=NULL)
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--prop_nonzero"), type="numeric", default=0.05)
parser <- add_option(parser, c("--update_w0"), type="logical", default=TRUE)
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
data_driven_cov <- outparse$data_driven_cov
canonical_cov <- outparse$canonical_cov
n <- outparse$n
prop_nonzero <- outparse$prop_nonzero
update_w0 <- outparse$update_w0
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
S0_data <- tryCatch(readRDS(data_driven_cov), 
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

###Compute prior covariance matrices (S0)
if(canonical_cov){
  S0_can <- compute_canonical_covs(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
}

if(!is.null(S0_data) && canonical_cov){
  S0 <- c(S0_data, S0_can)
} else if(is.null(S0_data) && canonical_cov){
  S0 <- S0_can
} else if(!is.null(S0_data) && !canonical_cov){
  S0 <- S0_data
}

grid <- autoselect.mixsd(univ_sumstats, mult=sqrt(2))^2
S0 <- expand_covs(S0, grid, zeromat=TRUE)

###Compute initial estimates of prior weights (w0)
w0 <- c((1-prop_nonzero), rep(prop_nonzero/(length(S0)-1), (length(S0)-1)))

###Subsets traits
if(length(traits)>1){
  univ_sumstats <- lapply(univ_sumstats, function(x, sel){x[, sel]}, traits)
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
                              verbose=verbose, update_V=update_V, update_V_method=update_V_method,
                              w0_threshold=w0_threshold, ca_update_order=ca_update_order, mu1_init=mu1_init,
                              X_colmeans=X_colmeans, Y_colmeans=Y_colmeans, nthreads=ncores)
toc <- proc.time()[[3]]

fit_mrmash_rss$elapsed_time <- toc-tic


saveRDS(fit_mrmash_rss, file=output)


