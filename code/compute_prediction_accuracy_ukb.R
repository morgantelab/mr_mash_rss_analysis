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
parser <- add_option(parser, c("--trait"), type="character")
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
trait <- outparse$trait
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
tmp <- tempfile(tmpdir="/data2/morgante_lab/fabiom/tmp")
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
##Get chromosomes to analyze
chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

##Get traits to analyze
traitscar <- unlist(strsplit(trait, ":"))

if(length(traitscar)==1){
  traits <- as.integer(traitscar)
} else {
  traits <- as.integer(traitscar[1]):as.integer(traitscar[2])
}

r <- length(traits)

pheno_pred <- vector("list", length=length(chrs))
Bhat_all <- vector("list", length=length(chrs))

##Loop over chromosomes
it <- 0

for(i in chrs){ 
  
  it <- it+1
  
  ##Chromosome SNP indexes
  if(chr==0){
    inds <- seq_len(p)
  } else {
    inds <- which(geno$map$chromosome==i)
  }
  
  if(model=="mr_mash_rss"){
    ##Read in model fit
    model_fit <- readRDS(paste0("../output/", model, "_fit/", prefix, "_chr", i, "_", model, "_fit_", data_id, ".rds"))
    
    ##Compute predictions
    pheno_pred[[it]] <- t(t(big_prodMat(X, model_fit$mu1[, traits], ind.col=inds)) + model_fit$intercept[traits])

    ##Store effects
    Bhat_all[[it]] <- model_fit$mu1[, traits]
    
    } else if(model=="ldpred2_auto"){
    
    it2 <- 0
    
    for(j in traits){
      
      it2 <- it2+1
      
      ##Read in model fit
      model_fit <- readRDS(paste0("../output/", model, "_fit/", prefix, "_chr", i, "_", model, "_fit_trait", j, "_", data_id, ".rds"))
      ##Quality control over chains
      range_corr <- sapply(model_fit, function(auto) diff(range(auto$corr_est)))
      to_keep <- (range_corr > (0.95 * quantile(range_corr, 0.95)))
      ##Compute posterior mean after QC
      if(it2==1){
        Bhat <- rowMeans(sapply(model_fit[to_keep], function(auto) auto$beta_est))
      } else {
        Bhat <- cbind(Bhat, rowMeans(sapply(model_fit[to_keep], function(auto) auto$beta_est)))
      }
    }
    
    ##Compute predictions
    pheno_pred[[it]] <- big_prodMat(X, Bhat, ind.col=inds)
    
    ##Store effects
    Bhat_all[[it]] <- Bhat
  }
}

###Compute predictions and accuracy
pheno_pred <- Reduce("+", pheno_pred)

##Remove mean of Y_training which was over counted in the per chromosome estimate of the intercept
if(model=="mr_mash_rss"){
  pheno_pred <- t(t(pheno_pred) - (pheno_means[traits]*(length(chrs)-1)))
}

accuracy <- compute_accuracy(pheno_test, pheno_pred)

###Save results to file
saveRDS(accuracy, file=paste0("../output/prediction_accuracy/", prefix, "_", model, "_pred_acc_", data_id, ".rds"))

###Save effects to file
effects <- do.call("rbind", Bhat_all)
saveRDS(effects, file=paste0("../output/estimated_effects/", prefix, "_", model, "_effects_", data_id, ".rds"))

###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))
