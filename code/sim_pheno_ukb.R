###Load libraries
library(Matrix)
library(bigsnpr)
library(optparse)
# library(data.table)

###Simulate data from given big X
simulate_mr_mash_data_from_given_big_X <- function(X, p_causal, r, r_causal, intercepts,
                                               pve, B_cor, B_scale, w, V_cor, seed, ncores=ncores){
  ##Check that the inputs are correct
  if(!inherits(X, 'FBM')){
    stop("X must be an object of class FBM")
  }
  if(length(intercepts)!=r)
    stop("intercepts must be of length equal to r.")
  if(any(sapply(r_causal, length)>r))
    stop("r_causal cannot be greater than r.")
  if(!(length(B_cor)==length(B_scale) & length(w)==length(B_cor) & length(B_scale)==length(w)))
    stop("B_cor, B_scale, and w must have the same length.")
  if(abs(sum(w) - 1) > 1e-10)
    stop("Elements of w must sum to 1.")
  if(length(pve)!=1 & length(pve)!=r)
    stop("pve must be of length equal to 1 or r.")
  
  ##Get number of mixture components, samples, and variables
  n <- nrow(X)
  p <- ncol(X)
  K <- length(w)
  
  ##Simulate true effects from N_r(0, Sigma) or \sum_K w_k N_r(0, Sigma_k) where Sigma and Sigma_k are given 
  ##covariance matrices across traits and w_k is the mixture proportion associated to Sigma_k
  Sigma <- vector("list", K)
  for(i in 1:K){
    r_mix_length <- length(r_causal[[i]])
    Sigma_offdiag <- B_scale[i]*B_cor[i]
    Sigma[[i]] <- matrix(Sigma_offdiag, nrow=r_mix_length, ncol=r_mix_length)
    diag(Sigma[[i]]) <- B_scale[i]
  }
  #Sample effects from a mixture of MVN distributions or a single MVN distribution
  B_causal <- matrix(0, nrow=p_causal, ncol=r)
  if(K>1){
    mixcomps <- sample(x=1:K, size=p_causal, prob=w, replace=TRUE)
    for(j in 1:p_causal){
      comp_to_use <- mixcomps[j]
      r_causal_mix <- r_causal[[comp_to_use]]
      B_causal[j, r_causal_mix] <- mvtnorm::rmvnorm(n=1, mean=rep(0, length(r_causal_mix)), sigma=Sigma[[comp_to_use]])
    }
  } else {
    r_causal_length <- length(r_causal[[1]])
    r_causal_index <- r_causal[[1]]
    B_causal[, r_causal_index] <- mvtnorm::rmvnorm(n=p_causal, mean=rep(0, r_causal_length), sigma=Sigma[[1]])
  }
  B <- matrix(0, ncol=r, nrow=p)
  causal_variables <- sample(x=(1:p), size=p_causal)
  B[causal_variables, ] <- B_causal
  
  ##Center X
  X_colmeans <- bigstatsr::big_colstats(X, ncores=ncores)$sum/nrow(X)
  
  ##Compute G and its variance
  G <- bigstatsr::big_prodMat(X, B, center=X_colmeans, ncores=ncores)
  B <- as(B, "CsparseMatrix")
  Var_G <- matrixStats::colVars(G)
  
  ##Compute residual covariance
  Var_E <- ((1/pve)-1)*Var_G
  Var_E[which(Var_E<=.Machine$double.eps)] <- 1
  D <- diag(x=sqrt(Var_E))
  V_cor_mat <- matrix(V_cor, nrow=r, ncol=r)
  diag(V_cor_mat) <- 1
  V <- D %*% V_cor_mat %*% D
  
  ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
  Y <- mr.mash.alpha:::matrix_normal_indep_rows(G + matrix(intercepts, n, r, byrow=TRUE), V, seed)
  
  ##Compile output
  causal_responses <- r_causal
  names(causal_responses) <- paste0("Component", 1:K)
  names(Sigma) <- paste0("Component", 1:K)
  out <- list(Y=Y, B=B, V=V, Sigma=Sigma, intercepts=intercepts, causal_responses=causal_responses)
  if(K>1){
    if(p_causal>1){
      causal_variables_mixcomps <- cbind(causal_variables, mixcomps)
      causal_variables_mixcomps <- causal_variables_mixcomps[order(causal_variables_mixcomps[, 1]), ]
      out$causal_variables <- causal_variables_mixcomps[, 1]
      out$causal_vars_to_mixture_comps <- causal_variables_mixcomps[, 2]
    } else {
      out$causal_variables <- causal_variables
      out$causal_vars_to_mixture_comps <- mixcomps
    }
  } else {
    out$causal_variables <- sort(causal_variables)
    out$causal_vars_to_mixture_comps <- rep(1, p_causal)
  }
  
  return(out)
}

###Parse option
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--p_causal"), type="integer")
parser <- add_option(parser, c("--r"), type="integer")
parser <- add_option(parser, c("--r_causal"), type="character")
parser <- add_option(parser, c("--pve"), type="character")
parser <- add_option(parser, c("--B_scale"), type="character")
parser <- add_option(parser, c("--B_cor"), type="character")
parser <- add_option(parser, c("--V_cor"), type="numeric")
parser <- add_option(parser, c("--w"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--temp_dir"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
outparse <- parse_args(parser)

geno <- outparse$geno
output <- outparse$output
p_causal <- outparse$p_causal
r <- outparse$r
r_causal <- eval(parse(text=outparse$r_causal))
pve <- eval(parse(text=outparse$pve))
B_scale <- eval(parse(text=outparse$B_scale))
B_cor <- eval(parse(text=outparse$B_cor))
w <- eval(parse(text=outparse$w))
V_cor <- outparse$V_cor
seed <- outparse$seed
temp_dir <- outparse$temp_dir
ncores <- outparse$ncores

###Set seed
set.seed(seed)

###Read in genotype data
tmp <- tempfile(tmpdir=temp_dir)
rds <- snp_readBed2(geno, backingfile=tmp, ncores=ncores)
dat <- snp_attach(rds)

###Impute missing values
X <- snp_fastImputeSimple(dat$genotypes, method="mean2", ncores=ncores)

###Convert from FBM to matrix
FID <- dat$fam$family.ID
IID <- dat$fam$sample.ID

###Simulate phenotype
out_sim <- simulate_mr_mash_data_from_given_big_X(X=X, p_causal=p_causal, r=r, r_causal=r_causal, intercepts=rep(1, r),
                                                  pve=pve, B_cor=B_cor, B_scale=B_scale, w=w, V_cor=V_cor, seed=seed,
                                                  ncores=ncores)
colnames(out_sim$Y) <- paste0("y", 1:r)
rownames(out_sim$Y) <- IID

pheno_file <- data.frame(FID=FID, IID=IID, out_sim$Y)

###Write out phenotypes
saveRDS(out_sim, file=output)
# fwrite(x=pheno_file, file=paste0((output, ".txt"), sep = "\t", 
#        row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)

###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))


