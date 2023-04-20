###Load libraries needed
library(optparse)
library(mr.mash.alpha)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--pheno_cov"), type="character")
parser <- add_option(parser, c("--residual_cov"), type="character")
parser <- add_option(parser, c("--n"), type="integer")
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
parser <- add_option(parser, c("--outprefix"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_matrix <- outparse$LD_matrix
pheno_cov <- outparse$pheno_cov
residual_cov <- outparse$residual_cov
n <- outparse$n
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
outprefix <- outparse$outprefix
data_id <- outparse$data_id

###Set seed
set.seed(data_id)

###Read in data
univ_sumstats <- readRDS(paste0("../output/summary_statistics/", sumstats, "_", data_id, ".rds"))
p <- nrow(univ_sumstats$Bhat)
r <- ncol(univ_sumstats$Bhat)
LD <- matrix(readBin(paste0("../data/LD_matrices/", LD_matrix, "_", data_id, ".ld.bin"), what="numeric", n=p^2), nrow=p, ncol=p, byrow=TRUE)
covY <- readRDS(paste0("../output/misc/", pheno_cov, "_", data_id, ".rds"))
V <- readRDS(paste0("../output/misc/", residual_cov, "_", data_id, ".rds"))
mu1_init <- tryCatch(readRDS(paste0("../output/misc/", mu1_init, "_", data_id, ".rds")), 
                      error = function(e) {
                        return(NULL)
                      },
                      warning = function(w) {
                        return(NULL)
                      }
)
X_colmeans <- tryCatch(readRDS(paste0("../output/misc/", X_colmeans, "_", data_id, ".rds")), 
                     error = function(e) {
                       return(NULL)
                     },
                     warning = function(w) {
                       return(NULL)
                     }
)
Y_colmeans <- tryCatch(readRDS(paste0("../output/misc/", Y_colmeans, "_", data_id, ".rds")), 
                     error = function(e) {
                       return(NULL)
                     },
                     warning = function(w) {
                       return(NULL)
                     }
)

###Compute prior covariance matrices (S0) with canonical matrices
S0_can <- compute_canonical_covs(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
grid <- autoselect.mixsd(univ_sumstats, mult=sqrt(2))^2
S0 <- expand_covs(S0_can, grid, zeromat=TRUE)

###Compute initial estimates of prior weights (w0)
prop_nonzero <- 0.05
w0 <- c((1-prop_nonzero), rep(prop_nonzero/(length(S0)-1), (length(S0)-1)))

###Fit mr.mash.rss
fit_mrmash_rss <- mr.mash.rss(Bhat=univ_sumstats$Bhat, Shat=univ_sumstats$Shat, covY=covY, R=LD, n=n, S0=S0, 
                              w0=w0, update_w0=update_w0, tol=tol, V=V, check_R=check_R, max_iter=max_iter,
                              convergence_criterion="ELBO", compute_ELBO=TRUE, standardize=standardize,
                              verbose=verbose, update_V=update_V, update_V_method=update_V_method,
                              w0_threshold=w0_threshold, ca_update_order=ca_update_order, mu1_init=mu1_init,
                              X_colmeans=X_colmeans, Y_colmeans=Y_colmeans, nthreads=ncores)

saveRDS(fit_mrmash_rss, file=paste0("../output/mr_mash_rss_fit/", outprefix, "_mr_mash_rss_fit_", data_id, ".rds"))
