###Load libraries needed
library(optparse)
library(mashr)
library(data.table)

###Function needed
is_weak <- function(x, thresh){
  weak <- all(abs(x) < thresh)
  return(weak)
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--regions_dir"), type="character")
parser <- add_option(parser, c("--n_strong"), type="integer")
parser <- add_option(parser, c("--n_weak"), type="integer")
parser <- add_option(parser, c("--residual_cov"), type="character", default="full")
parser <- add_option(parser, c("--n_PCs"), type="integer")
parser <- add_option(parser, c("--flash_remove_singleton"), type="logical")
parser <- add_option(parser, c("--ED_algorithm"), type="character")
parser <- add_option(parser, c("--ted_zero_thresh"), type="numeric", default=1e-10)
parser <- add_option(parser, c("--canonical_cov"), type="logical")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
output <- outparse$output
regions_dir <- outparse$regions_dir
n_strong <- outparse$n_strong
n_weak <- outparse$n_weak
residual_cov <- outparse$residual_cov
npcs <- outparse$n_PCs
flash_remove_singleton <- outparse$flash_remove_singleton
ED_algorithm <- outparse$ED_algorithm
ted_zero_thresh <- outparse$ted_zero_thresh
canonical_cov <- outparse$canonical_cov
seed <- outparse$seed

###Set seed
set.seed(seed)

filenames <- list.files(regions_dir, pattern="*.txt")

dat <- readRDS(sumstats)
Z <- dat$Bhat/dat$Shat

###Extract effects
it <- 0
for(nam in filenames){
  it <- it+1
  
  ##Select only variants in that region
  rsids <- fread(paste(regions_dir, nam, sep="/"), showProgress=FALSE, 
               colClasses = "character", header=TRUE, data.table = FALSE)
  Z_sel <- Z[rsids$ID, ]
  
  ##Obtain strong signals
  strong_id <- names(sort(apply(abs(Z_sel), 1, max), decreasing = TRUE)[1:n_strong])
  
  ##Obtain weak signals
  if(residual_cov=="full"){
    weak_id <- names(sample(x=which(apply(Z_sel, 1, is_weak, 2)), size=n_weak))
  }
  
  if(it==1){
    Z_strong <- Z_sel[strong_id, ]
    
    if(residual_cov=="full"){
      Z_weak <- Z_sel[weak_id, ]
    }
  } else {
    Z_strong <- rbind(Z_strong, Z_sel[strong_id, ])
    
    if(residual_cov=="full"){
      Z_weak <- rbind(Z_weak, Z_sel[weak_id, ])
    }
  }
}

###Estimate data-driven prior
##Estimate residual covariance
if(residual_cov=="diagonal"){
  V <- diag(ncol(Z_strong))
} else if(residual_cov=="full"){
  V <- matrix(0, nrow=ncol(Z_weak), ncol=ncol(Z_weak))
  
  for(j in 1:nrow(Z_weak)){
    V <- V + tcrossprod(Z_weak[j,])
  }
  
  V <- V/nrow(Z_weak)
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


