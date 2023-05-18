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
parser <- add_option(parser, c("--ED_algorithm"), type="character", )
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
    Bhat_strong <- dat$Bhat[strong, ]
    Shat_strong <- dat$Shat[strong, ]
  } else {
    Bhat_strong <- rbind(Bhat_strong, dat$Bhat[strong, ])
    Shat_strong <- rbind(Shat_strong, dat$Shat[strong, ])
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

dat_mash <- mash_set_data(Bhat_strong, Shat_strong, V=cov2cor(V))

##Compute data-driven matrices
U_pca <- cov_pca(data=dat_mash, npc=npcs)
U_flash_default <- cov_flash(data=dat_mash, factors="default", tag="default",
                              remove_singleton=flash_remove_singleton, output_model=NULL)
U_flash_nonneg <- cov_flash(data=dat_mash, factors="nonneg", tag="nonneg",
                             remove_singleton=flash_remove_singleton, output_model=NULL)
U_emp <- mr.mash.alpha:::cov_empirical(data=dat_mash)

##De-noise data-driven matrices via extreme deconvolution
U_datadriven <- c(U_pca, U_flash_default, U_flash_nonneg, list(BB=U_emp))
if(ED_algorithm=="bovy"){
  U_ed <- cov_ed(dat_mash, U_datadriven)
} else if(ED_algorithm=="ted" || ED_algorithm=="ed"){
  library(udr)
  f0 <- ud_init(X = as.matrix(Bhat_strong/Shat_strong), V = V, U_scaled = list(), 
                U_unconstrained = U_datadriven, n_rank1=0)
  res <-ud_fit(f0, X = na.omit(f0$X), control = list(unconstrained.update = ED_algorithm, 
                                                     resid.update = "none",
                                                     maxiter=5000, tol.lik = 1e-3), 
                                                     verbose=TRUE)

U_ed <- lapply(res$U,"[[",2)
}


###Save file
saveRDS(U_ed, file=output)


