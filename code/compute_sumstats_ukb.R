###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)

###Function to compute summary stats
compute_univariate_sumstats_bigsnp <- function(geno_obj, Y, chr=0, standardize=FALSE, 
                                               impute_missing=TRUE,
                                               normalize=FALSE, mc.cores=1){
  
  if(impute_missing){
    X <- bigsnpr::snp_fastImputeSimple(geno_obj$genotypes, method = "mean2", ncores = mc.cores)
  } else {
    X <- geno_obj$genotypes
  }

  r <- ncol(Y)
  p <- geno_obj$genotypes$ncol
  
  variable_names <- geno_obj$map$marker.ID
  response_names <- colnames(Y)
  
  if(normalize){
    ###Function for quantile normalization
    inv_normalise <- function(x) { #this would also tolerate NAs
      return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
    }

    Y <- apply(Y, 2, inv_normalise)
  }
  
  if(chr==0){
    inds <- seq_len(p)
  } else {
    inds <- which(geno_obj$map$chromosome==chr)
  }
  
  linreg <- function(i, X, Y, cols){
    
    fit <- bigstatsr::big_univLinReg(X=X, y=Y[, i], ind.col=cols, ncores = 1)
    bhat <- fit$estim
    shat <- fit$std.err
    
    return(list(bhat=bhat, shat=shat))
  }
  
  if(mc.cores>1){
    cl <- parallel::makeCluster(mc.cores)
    out <- parallel::parLapply(cl, 1:r, linreg, X, Y, inds)
    parallel::stopCluster(cl)
  } else {
    out <- lapply(1:r, linreg, X, Y, inds)
  }
  
  Bhat <- sapply(out,"[[","bhat")
  Shat <- sapply(out,"[[","shat")
  colnames(Bhat) <- colnames(Shat) <- response_names
  rownames(Bhat) <- rownames(Shat) <- variable_names[inds]
  
  ###Add some information
  statz <- bigstatsr::big_colstats(X=X, ind.col=inds, ncores=mc.cores)
  a1f <- (statz$sum/nrow(X))/2
  alleles_df <- data.frame(chr=chr, bp=geno_obj$map$physical.pos[inds], a1=geno_obj$map$allele1[inds], a2=geno_obj$map$allele2[inds], freq=a1f)
  rownames(alleles_df) <- variable_names[inds]
  
  ###Put coefficients and ses on the standardized X scale   
  if(standardize){
    sds <- sqrt(statz$var)
    
    Bhat <- Bhat*sds
    Shat <- Shat*sds
  }
  
  return(list(Bhat=Bhat, Shat=Shat))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--chr"), type="integer")
parser <- add_option(parser, c("--standardize"), type="logical")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical")
parser <- add_option(parser, c("--temp_dir"), type="character")
parser <- add_option(parser, c("--normalize"), type="logical")
outparse <- parse_args(parser)

pheno_dat <- outparse$pheno
geno_dat <- outparse$geno
test_ids <- outparse$test_ids
output <- outparse$output
seed <- outparse$seed
chr <- outparse$chr
standardize <- outparse$standardize
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing
temp_dir <- outparse$temp_dir
normalize <- outparse$normalize

###Set seed
set.seed(seed)

###Read in data and filter them
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(paste0(unlist(strsplit(geno_dat, ".", fixed=TRUE))[1], ".fam"), showProgress=FALSE, header=FALSE)
test_ids <- fread(test_ids, showProgress=FALSE, header=FALSE)
training_inds_geno <- which(!(geno_fam[,2] %in% test_ids[,2])) ##Get only training individuals

pheno <- readRDS(pheno_dat)$Y
training_inds_pheno <- which(!(rownames(pheno) %in% test_ids[,2])) ##Get only training individuals
pheno <- pheno[training_inds_pheno, ]

tmp <- tempfile(tmpdir=temp_dir)
rds <- snp_readBed2(geno_dat, ind.row=training_inds_geno, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

rm(list=c("geno_fam", "test_ids", "training_inds_geno", "training_inds_pheno"))

###Compute summary stats
out <- compute_univariate_sumstats_bigsnp(geno_obj=geno, Y=pheno, chr=chr, standardize=standardize,
                                          impute_missing=impute_missing,
                                          normalize=normalize, mc.cores=ncores)

if(chr==0){
  chr <- "All"
}

saveRDS(out, file=output)


###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))
