###Load libraries needed
library(optparse)
library(bigsnpr)
library(bigstatsr)

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
  
  variable_names <- geno_obj$map$rsid
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
      dat <- data.frame(y=Y[, s], Z)
      fit_norm <- lm(y ~ ., dat)
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
  
  if(!is.null(Z) && !is.matrix(Z)){
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
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--sample_file"), type="character")
parser <- add_option(parser, c("--fold"), type="integer")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--standardize"), type="logical")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical")
parser <- add_option(parser, c("--normalize"), type="logical")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

geno <- outparse$geno
pheno <- outparse$pheno
sample_file <- outparse$sample_file
chr <- outparse$chr
fold <- outparse$fold
standardize <- outparse$standardize
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing
normalize <- outparse$normalize
output <- outparse$output

###Set seed
set.seed(1)

###
if(chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9")){
  chr <- paste0("0", chr)
}

###Read in data
pheno_dat <- readRDS(pheno)
rownames(pheno_dat) <- pheno_dat$id
pheno_dat$id <- NULL
ids <- bigreadr::fread2(sample_file)
ids <- ids[-1, ]

###Order phenotype and covariate data according to genotype data
ids <- ids[which(as.character(ids$ID_2) %in% rownames(pheno_dat)), "ID_2"]
pheno_dat <- pheno_dat[as.character(ids), ]

###Get training individuals indexes
train_inds <- which(pheno_dat$fold != fold)

###Load genotype data
geno <- snp_attach(geno)

###Compute summary stats
#Prepare phenotypes
pheno <- pheno_dat[, c("DPa", "SPa", "height", "hip", "BMI", "BMR", "WHR")]

#Prepare covariates
covar <- pheno_dat[, c("sex", "assessment_centre", "age", "genotype_measurement_batch",
                       paste0("pc_genetic", 1:10))]
covar$age2 <- covar$age^2
covar <- covar[, c("sex", "assessment_centre", "age", "age2", "genotype_measurement_batch",
                       paste0("pc_genetic", 1:10))]
covar$sex <- as.factor(covar$sex)
covar$assessment_centre <- as.factor(covar$assessment_centre)
covar$genotype_measurement_batch <- as.factor(covar$genotype_measurement_batch)

if(!normalize){
  covar <- covar_from_df(covar)
}

#Fit linear model
out <- compute_univariate_sumstats_bigsnp(geno_obj=geno, Y=pheno, Z=covar,
                                          train_inds=train_inds, chr=chr, standardize=standardize,
                                          impute_missing=impute_missing, standardize.response=FALSE, 
                                          normalize=normalize, mc.cores=ncores)

saveRDS(out, file=output)

