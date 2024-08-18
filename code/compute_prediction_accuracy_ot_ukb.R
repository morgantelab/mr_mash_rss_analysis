###Load libraries needed
library(optparse)
library(bigreadr)
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
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--sample_file"), type="character")
parser <- add_option(parser, c("--samples_to_keep"), type="character")
parser <- add_option(parser, c("--model_fit_dir"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical", default=TRUE)
parser <- add_option(parser, c("--output_eff"), type="character")
parser <- add_option(parser, c("--output_pred_acc"), type="character")
parser <- add_option(parser, c("--prefix"), type="character")
parser <- add_option(parser, c("--fold"), type="integer")
outparse <- parse_args(parser)

model <- outparse$model
chrs <- eval(parse(text=outparse$chr))
traits <- eval(parse(text=outparse$trait))
pheno_dat <- outparse$pheno
geno_dat <- outparse$geno
sample_file <- outparse$sample_file
samples_to_keep <- outparse$samples_to_keep
model_fit_dir <- outparse$model_fit_dir
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing
output_eff <- outparse$output_eff
output_pred_acc <- outparse$output_pred_acc
prefix <- outparse$prefix
fold <- outparse$fold

###Set seed
set.seed(1)

###Read in genotypes
geno <- snp_attach(geno_dat)

###Impute missing genotypes
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

p <- geno$genotypes$ncol

###Read in sample file
ids <- bigreadr::fread2(sample_file)
ids <- ids[-1, ]

###Read in list of individuals selected for the analysis
ids1 <- data.table::fread(samples_to_keep, data.table=FALSE, 
                          colClasses = "character", header=FALSE)

###Read in phenotype
pheno_test <- readRDS(pheno_dat)

###Get only bc selected individuals
ids_sel <- as.character(ids[as.character(ids$ID_2) %in% ids1[,2], 2])

###Get only test individuals indeces and ids
test_inds_geno <- which(ids_sel %in% as.character(rownames(pheno_test))) ##Get only test individuals
test_ids_geno <- ids_sel[test_inds_geno]

###Order pheno data according to geno data
pheno_test <-  pheno_test[as.character(test_ids_geno), ]

###Compute predictions
##Get traits to analyze
r <- length(traits)

pheno_pred <- vector("list", length=length(chrs))
Bhat_all <- vector("list", length=length(chrs))

##Loop over chromosomes
it <- 0

for(i in chrs){
  
  it <- it+1
  
  ##Chromosome SNP indexes
  if(i==0){
    inds <- seq_len(p)
    i <- "All"
  } else {
    if(i %in% 1:9){
      ii <- paste0("0", i)
    } else {
      ii <- as.character(i)
    }
    inds <- which(geno$map$chromosome==ii)
  }
  
  if(model %in% c("mr_mash_rss_sparse_LD_V_all_chr_ldpred2_auto_init_prior_finemapped",
                  "mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped")){
    ##Read in model fit
    model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_", fold, ".rds"))
    
    ##Store effects
    Bhat <- model_fit$mu1[, traits]
    
  } else if(model %in% c("ldpred2_auto","bayesN","bayesA","bayesL","bayesC","bayesR")){
    
    Bhat <- matrix(as.numeric(NA), nrow=length(inds), ncol=length(traits))
    
    for(j in traits){
      
      ##Read in model fit
      model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_trait", j, "_", fold, ".rds"))
      
      if(model=="ldpred2_auto"){
        
        model_fit$elapsed_time <- NULL
        
        ##Quality control over chains
        range_corr <- sapply(model_fit, function(auto) diff(range(auto$corr_est)))
        to_keep <- (range_corr > (0.95 * quantile(range_corr, 0.95, na.rm = TRUE)))
        
        if(all(is.na(to_keep)) || all(to_keep==FALSE)){
          warning(paste0("All chains dropped because of divergence! Skipping trait ", j, "for chr", i))
          next
        } else if(any(is.na(to_keep))){
          to_keep[which(is.na(to_keep))] <- FALSE
          warning("Some chains dropped because of divergence.")
        } else if(any(to_keep==FALSE)){
          warning("Some chains dropped because of divergence.")
        }
        
        ##Compute posterior mean after QC
        if(sum(to_keep)==1){
          Bhat[, j] <- sapply(model_fit[to_keep], function(auto) auto$beta_est)
        } else {
          Bhat[, j] <- rowMeans(sapply(model_fit[to_keep], function(auto) auto$beta_est))
        }
      } else if(model %in% c("bayesN","bayesA","bayesL","bayesC","bayesR")){
        Bhat[, j] <- model_fit$bm
      }    
    }
  }
  
  ##Store effects
  Bhat_all[[it]] <- Bhat
  
  if(model=="ldpred2_auto"){
    saveRDS(Bhat, file=paste0("../output/estimated_effects/", prefix, "_chr", i, "_ldpred2_auto_effects_", fold, ".rds"))
  }
  
  if(model=="bayesR"){
    saveRDS(Bhat, file=paste0("../output/estimated_effects/", prefix, "_chr", i, "_bayesR_effects_", fold, ".rds"))
  }

  ##Compute predictions
  pheno_pred[[it]] <- big_prodMat(X, Bhat, ind.row=test_inds_geno, ind.col=inds, ncores=ncores)
}

###Compute predictions and accuracy
pheno_pred <- Reduce("+", pheno_pred)
accuracy <- compute_accuracy(pheno_test[, traits], pheno_pred)

###Save accuracy to file
saveRDS(accuracy, file=output_pred_acc)

###Save effects to file
effects <- do.call("rbind", Bhat_all)
saveRDS(effects, file=output_eff)
