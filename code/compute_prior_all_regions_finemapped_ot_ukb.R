###Load libraries needed
library(optparse)
library(mashr)
library(data.table)
library(bigsnpr)
library(bigstatsr)
library(susieR)

options(warn=1)

###Function needed
is_weak <- function(x, thresh){
  weak <- all(abs(x) < thresh)
  return(weak)
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--fold_file"), type="character")
parser <- add_option(parser, c("--sample_file"), type="character")
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--fold"), type="integer")
parser <- add_option(parser, c("--regions_dir"), type="character")
parser <- add_option(parser, c("--n_weak"), type="integer")
parser <- add_option(parser, c("--residual_cov"), type="character", default="full")
parser <- add_option(parser, c("--n_PCs"), type="integer")
parser <- add_option(parser, c("--flash_remove_singleton"), type="logical")
parser <- add_option(parser, c("--ED_algorithm"), type="character")
parser <- add_option(parser, c("--ted_zero_thresh"), type="numeric", default=1e-10)
parser <- add_option(parser, c("--canonical_cov"), type="logical")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--temp_dir"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
fold_file <- outparse$fold_file
sample_file <- outparse$sample_file
sumstats <- outparse$sumstats
output <- outparse$output
fold <- outparse$fold
regions_dir <- outparse$regions_dir
n_weak <- outparse$n_weak
residual_cov <- outparse$residual_cov
npcs <- outparse$n_PCs
flash_remove_singleton <- outparse$flash_remove_singleton
ED_algorithm <- outparse$ED_algorithm
ted_zero_thresh <- outparse$ted_zero_thresh
canonical_cov <- outparse$canonical_cov
seed <- outparse$seed
temp_dir <- outparse$temp_dir

###Set seed
set.seed(seed)

filenames <- list.files(regions_dir, pattern="*.txt")

###Read in data
geno <- snp_attach(geno_dat)
dat <- readRDS(sumstats)
Z <- dat$Bhat/dat$Shat

##This is only needed to get the individuals in the training set for LD computation 
fold_info <- readRDS(fold_file)
rownames(fold_info) <- fold_info$id
fold_info$id <- NULL
ids <- bigreadr::fread2(sample_file)
ids <- ids[-1, ]
ids <- ids[which(as.character(ids$ID_2) %in% rownames(fold_info)), "ID_2"]
fold_info <- fold_info[as.character(ids), ]
train_inds <- which(fold_info$fold != fold)

rm(fold_info, ids)

###Extract effects
it <- 0
for(nam in filenames){
  it <- it+1
  
  ##Select only variants in that region
  rsids <- fread(paste(regions_dir, nam, sep="/"), showProgress=FALSE, 
               colClasses = "character", header=TRUE, data.table = FALSE)

  rsids_common <- intersect(rownames(Z), rsids$ID)
  
  Z_sel <- Z[rsids_common, ]
  Bhat_sel <- dat$Bhat[rsids_common, ]
  Shat_sel <- dat$Shat[rsids_common, ]
  snp_sel_geno <- which(geno$map$rsid %in% rsids_common)
  
  ##Compute LD
  tmp <- tempfile(tmpdir=temp_dir)
  LD <- big_cor(X=geno$genotypes, backingfile=tmp, ind.col=snp_sel_geno, ind.row=train_inds)
  LD <- LD[]

  ##Fine mapping with susie_rss
  high_pip_snps_cs_idx <- c()
  snps_cs_idx <- c()
    
  for(i in 1:ncol(Z_sel)){
    fit_susie_rss <- suppressMessages(susie_rss(bhat=Bhat_sel[, i], shat=Shat_sel[, i], 
                                                var_y=1, R=LD, n=length(train_inds), L=10, 
                                                estimate_residual_variance=TRUE))
    if(!is.null(fit_susie_rss$sets$cs)){
      cs <- fit_susie_rss$sets$cs
      for(j in 1:length(cs)){
        which_max_pip <- which.max(fit_susie_rss$pip[cs[[j]]])
        high_pip_snps_cs_idx <- c(high_pip_snps_cs_idx, cs[[j]][which_max_pip])
      }
      
      snps_cs_idx <- c(snps_cs_idx, unlist(cs, use.names=FALSE))
    } else {
      next
    }
  }
  
  if(length(high_pip_snps_cs_idx)>0){
    high_pip_snps_idx <- high_pip_snps_cs_idx
  } else {
    file.remove(paste0(tmp, ".bk"))
    warning(paste0("Region ", nam, " does not have any finemapped SNP."))
    cat(sprintf("Finished analyzing region %d.\n", it))
    next
  }
  
  ##Obtain strong signals
  strong_idx <- sort(unique(high_pip_snps_idx))
  
  ##Obtain weak signals
  if(residual_cov=="full"){
    Z_sel_no_strong <- Z_sel[-strong_idx, ]
      
    weak_id <- names(sample(x=which(apply(Z_sel_no_strong, 1, is_weak, 2)), size=n_weak))
  }
  
  ###Extract strong and weak signals
  if(it==1){
    Z_strong <- Z_sel[strong_idx, , drop=FALSE]
    
    if(residual_cov=="full"){
      Z_weak <- Z_sel[weak_id, , drop=FALSE]
    }

  } else {
    Z_strong <- rbind(Z_strong, Z_sel[strong_idx, , drop=FALSE])

    if(residual_cov=="full"){
      Z_weak <- rbind(Z_weak, Z_sel[weak_id, , drop=FALSE])
    }
  }
  
  file.remove(paste0(tmp, ".bk"))
  cat(sprintf("Finished analyzing region %d.\n", it))
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
  fit1 <- ud_fit(fit0, control = list(unconstrained.update = "ted", lambda = R, penalty.type = "iw",
                                      maxiter=5000, tol = 1e-2, tol.lik = 1e-3), verbose=TRUE,
                                      zero.threshold=ted_zero_thresh)

  U_ed <- fit1$U
}


###Save file
saveRDS(U_ed, file=output)


