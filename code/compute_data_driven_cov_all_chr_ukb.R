###Load libraries needed
library(optparse)
library(mashr)

###Function needed
is_strong <- function(x, thresh){
  strong <- any(abs(x) > thresh)
  return(strong)
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats_prefix"), type="character")
parser <- add_option(parser, c("--sumstats_suffix"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--residual_cov"), type="character")
parser <- add_option(parser, c("--strong_Z_thresh"), type="integer")
parser <- add_option(parser, c("--n_PCs"), type="integer")
parser <- add_option(parser, c("--flash_remove_singleton"), type="logical")
parser <- add_option(parser, c("--ED_algorithm"), type="character")
parser <- add_option(parser, c("--ted_zero_thresh"), type="numeric", default=1e-10)
parser <- add_option(parser, c("--canonical_cov"), type="logical")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

sumstats_prefix <- outparse$sumstats_prefix
sumstats_suffix <- outparse$sumstats_suffix
output <- outparse$output
chr <- outparse$chr
residual_cov <- outparse$residual_cov
strong_thresh <- outparse$strong_Z_thresh
npcs <- outparse$n_PCs
flash_remove_singleton <- outparse$flash_remove_singleton
ED_algorithm <- outparse$ED_algorithm
ted_zero_thresh <- outparse$ted_zero_thresh
canonical_cov <- outparse$canonical_cov
seed <- outparse$seed

###Set seed
set.seed(seed)

chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

###Extract strong effects
for(i in chrs){
  dat <- readRDS(paste0(sumstats_prefix, i, sumstats_suffix))
  Z <- dat$Bhat/dat$Shat
  strong <- which(apply(Z, 1, is_strong, strong_thresh))
  
  if(i==1){
    Z_strong <- Z[strong, ]
  } else {
    Z_strong <- rbind(Z_strong, Z[strong, ])
  }
}

if(nrow(Bhat_strong) < 20){
  stop("Too few strong effects with current strong_Z_thresh value.")
}

###Estimate data-driven prior
if(residual_cov=="diagonal"){
  V <- diag(ncol(Bhat_strong))
} else {
  V <- readRDS(residual_cov)
}

dat_mash <- mash_set_data(Z_strong, V=cov2cor(V))

if(ED_algorithm=="bovy"){
  ##Compute factor analyses
  U_pca <- cov_pca(data=dat_mash, npc=npcs)
  U_flash_default <- cov_flash(data=dat_mash, factors="default", tag="default",
                               remove_singleton=flash_remove_singleton, output_model=NULL)
  # U_flash_nonneg <- cov_flash(data=dat_mash, factors="nonneg", tag="nonneg",
  #                             remove_singleton=flash_remove_singleton, output_model=NULL)
  U_emp <- crossprod(Z_strong) / nrow(Z_strong)

  ##De-noise data-driven matrices via extreme deconvolution
  U_datadriven <- c(list(BB=U_emp), U_pca)

  if(!is.null(U_flash_default)){
    U_datadriven <- c(U_datadriven, U_flash_default)
  }

  # if(!is.null(U_flash_nonneg)){
  #   U_datadriven <- c(U_datadriven, U_flash_nonneg)
  # }
  
  if(canonical_cov){
    U_can <- cov_canonical(dat_mash)
    U_datadriven <- c(U_datadriven, U_can)
  }
  
  U_ed <- cov_ed(dat_mash, U_datadriven)
} else if(ED_algorithm=="ted"){
  library(udr)
  # 1. Initialize 50 unconstrained covariance matrices for udr.
  R <- nrow(V)
  K <- 50

  # 2. Add small amount of penalty(inverse-Wishart), strength = R. When sample size is large, 
  # the penalty amount R is small.
  fit0 <- ud_init(dat_mash, n_unconstrained = K)
  fit1 <- ud_fit(fit0, control = list(unconstrained.update = "ted", lambda = R, penalty.type = "iw",                                      maxiter=5000, tol = 1e-2, tol.lik = 1e-3), verbose=TRUE,
                                      zero.threshold=ted_zero_thresh)

  U_ed <- fit1$U
}


###Save file
saveRDS(U_ed, file=output)


