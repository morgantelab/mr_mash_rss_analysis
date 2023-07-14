###Load libraries needed
library(optparse)
library(bigsnpr)
library(RhpcBLASctl)

###Function to compute summary stats
compute_univariate_sumstats_bigsnp <- function(geno_obj, Y, Z=NULL, chr=0, standardize=FALSE, 
                                               impute_missing=TRUE, train_inds=NULL,
                                               normalize=FALSE, standardize.response=FALSE, 
                                               mc.cores=1){
  
  if(impute_missing){
    X <- bigsnpr::snp_fastImputeSimple(geno_obj$genotypes, method = "mean2", ncores = mc.cores)
  } else {
    X <- geno_obj$genotypes
  }

  r <- ncol(Y)
  n <- nrow(Y)
  p <- geno_obj$genotypes$ncol
  
  variable_names <- geno_obj$map$marker.ID
  response_names <- colnames(Y)
  
  if(is.null(train_inds)){
    row_inds <- seq_len(n)
  } else {
    row_inds <- train_inds
    Y <- Y[row_inds, ]
    if(!is.null(Z)){
      Z <- Z[row_inds, ]
    }
  }
  
  if(standardize.response){
    Y <- mr.mash.alpa::scale_fast2(Y, scale=TRUE)$M
  }
  
  if(chr=="0"){
    col_inds <- seq_len(p)
  } else {
    col_inds <- which(geno_obj$map$chromosome==chr)
  }
  
  if(normalize){
    
    inv_normalise <- function(x) { #this would also tolerate NAs
      return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
    }
    
    res_mat_train <- matrix(as.numeric(NA), nrow=length(row_inds), ncol=r)
    colnames(res_mat_train) <- colnames(Y)
    
    for(s in 1:r){
      fit_norm <- lm(Y[, s] ~ as.matrix(Z))
      res <- resid(fit_norm) + coef(fit_norm)[1]
      res_mat_train[, s] <- inv_normalise(res)
    }
    
    Z <- NULL
    Y <- res_mat_train
  }
  
  linreg <- function(i, X, Y, Z, cols, rows){
    
    fit <- bigstatsr::big_univLinReg(X=X, y=Y[, i], covar.train=Z,
                                     ind.col=cols, ind.train=rows,  ncores = 1)
    bhat <- fit$estim
    shat <- fit$std.err
    
    return(list(bhat=bhat, shat=shat))
  }
  
  if(!is.null(Z)){
    Z <- as.matrix(Z)
  }
  
  if(mc.cores>1){
    cl <- parallel::makeCluster(mc.cores)
    out <- parallel::parLapply(cl, 1:r, linreg, X, Y, Z, col_inds, row_inds)
    parallel::stopCluster(cl)
  } else {
    out <- lapply(1:r, linreg, X, Y, Z, col_inds, row_inds)
  }
  
  Bhat <- sapply(out,"[[","bhat")
  Shat <- sapply(out,"[[","shat")
  colnames(Bhat) <- colnames(Shat) <- response_names
  rownames(Bhat) <- rownames(Shat) <- variable_names[col_inds]
  
  ###Put coefficients and ses on the standardized X scale   
  if(standardize){
    sds <- sqrt(bigstatsr::big_colstats(X, ind.col=col_inds, ind.row=row_inds)$var)
  	
    Bhat <- Bhat*sds
    Shat <- Shat*sds
  }
  
  return(list(Bhat=Bhat, Shat=Shat))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--fold"), type="integer")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--standardize"), type="logical")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical")
parser <- add_option(parser, c("--normalize"), type="logical")
outparse <- parse_args(parser)

chr <- outparse$chr
fold <- outparse$fold
standardize <- outparse$standardize
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing
normalize <- outparse$normalize

###Set MKL threads
RhpcBLASctl::blas_set_num_threads(1)

###Set seed
set.seed(1)

###
if(chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9")){
  chr <- paste0("0", chr)
}

###Read in data
pheno_dat <- readRDS("../data/phenotypes/ukb_tiezzi_cleaned_BP_covar_pheno.rds")
rownames(pheno_dat) <- pheno_dat$ID
pheno_dat$ID <- NULL
ukb_data_location <- "/data2/morgante_lab/data/ukbiobank/genotypes/imputed/"
ids <- bigreadr::fread2(paste0(ukb_data_location, "ukb22828_c1_b0_v3_s487271.sample"))
ids <- ids[-1, ]

###Order phenotype and covariate data according to genotype data
ids <- ids[which(as.character(ids$ID_2) %in% rownames(pheno_dat)), "ID_2"]
pheno_dat <- pheno_dat[as.character(ids), ]

###Get training individuals indexes
train_inds <- which(pheno_dat$fold != fold)

###Load genotype data
geno <- snp_attach("/scratch1/fabiom/ukb_geno_imp_HM3_tiezzi.rds")

###Compute summary stats
out <- compute_univariate_sumstats_bigsnp(geno_obj=geno, Y=pheno_dat[, 1:6], Z=pheno_dat[, 7:28],
                                          train_inds=train_inds, chr=chr, standardize=standardize,
                                          impute_missing=impute_missing, standardize.response=FALSE, 
                                          normalize=normalize, mc.cores=ncores)

if(chr=="0"){
  chr <- "All"
} else if(chr %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09")){
  chr <- unlist(strsplit(chr, split=""))[2]
}

saveRDS(out, file=paste0("../output/summary_statistics/ukb_tiezzi_BP_chr", chr, "_sumstats_", fold, ".rds"))

