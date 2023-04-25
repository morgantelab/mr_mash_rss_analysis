###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)

###Function to compute accuracy
compute_accuracy <- function(Y, Yhat) {
  bias <- rep(as.numeric(NA), ncol(Y))
  r2 <- rep(as.numeric(NA), ncol(Y))
  rmse <- rep(as.numeric(NA), ncol(Y))
  
  for(i in 1:ncol(Y)){
    if(var(Yhat[, i])>0){
      fit  <- lm(Y[, i] ~ Yhat[, i])
      bias[i] <- coef(fit)[2] 
      r2[i] <- summary(fit)$r.squared
      rmse[i] <- sqrt(mean((Y[, i] - Yhat[, i])^2))
    } else {
      bias[i] <- NA
      r2[i] <- NA
      rmse[i] <- sqrt(mean((Y[, i] - Yhat[, i])^2))
    }
    
  }
  
  return(list(bias=bias, r2=r2, rmse=rmse, scaled_rmse=rmse/matrixStats::colSds(Y)))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--model"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--pheno_means"), type="character")
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--prefix"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical", default=TRUE)
outparse <- parse_args(parser)

model <- outparse$model
chr <- outparse$chr
pheno_means <- outparse$pheno_means
pheno_dat <- outparse$pheno
geno_dat <- outparse$geno
test_ids <- outparse$test_ids
prefix <- outparse$prefix
data_id <- outparse$data_id
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing


###Set seed
set.seed(data_id)

###Read in data
options(datatable.fread.datatable=FALSE)
test_ids <- fread(paste0("../data/phenotypes/", test_ids, "_", data_id, ".txt"), showProgress=FALSE)

pheno <- readRDS(paste0("../data/phenotypes/", pheno_dat, "_", data_id, ".rds"))$Y
test_inds_pheno <- which(rownames(pheno) %in% test_ids[,2]) ##Get only test individuals
pheno_test <- pheno[test_inds_pheno, ]

geno_fam <- fread(paste0("../data/genotypes/", geno_dat, ".fam"), showProgress=FALSE)
test_inds_geno <- which(geno_fam[,2] %in% test_ids[,2]) ##Get only test individuals
tmp <- tempfile()
rds <- snp_readBed2(paste0("../data/genotypes/", geno_dat, ".bed"), ind.row=test_inds_geno, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

pheno_means <- readRDS(paste0("../output/misc/", pheno_means, "_", data_id, ".rds"))

rm(list=c("geno_fam", "test_ids", "test_inds_geno", "test_inds_pheno", "pheno"))

###Impute missing genotypes
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

p <- geno$genotypes$ncol

###Compute predictions
chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

pheno_pred <- vector("list", length=length(chrs))

it <- 0

for(i in chrs){ 
  
  it <- it+1
  
  ##Chromosome SNP indexes
  if(chr==0){
    inds <- seq_len(p)
  } else {
    inds <- which(geno$map$chromosome==i)
  }
  
  ##Read in model fit
  model_fit <- readRDS(paste0("../output/", model, "_fit/", prefix, "_chr", i, "_", model, "_fit_", data_id, ".rds"))

  ##Compute predictions
  pheno_pred[[it]] <- t(t(big_prodMat(X, model_fit$mu1, ind.col=inds)) + model_fit$intercept)
}

###Compute predictions and accuracy
pheno_pred <- Reduce("+", pheno_pred)
pheno_pred <- t(t(pheno_pred) - (pheno_means*(length(chrs)-1)))

accuracy <- compute_accuracy(pheno_test, pheno_pred)

###Save results to file
saveRDS(accuracy, file=paste0("../output/accuracy/", prefix, "_", model, "_pred_acc_", data_id, ".rds"))


