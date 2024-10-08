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
    if(all(is.na(Yhat[, i]))){
      bias[i] <- NA
      r2[i] <- NA
      rmse[i] <- NA
    } else if(var(Yhat[, i])>0){
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
parser <- add_option(parser, c("--model_fit_dir"), type="character")
parser <- add_option(parser, c("--sumstats_dir"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical", default=TRUE)
parser <- add_option(parser, c("--output_eff"), type="character")
parser <- add_option(parser, c("--output_pred_acc"), type="character")
parser <- add_option(parser, c("--temp_dir"), type="character")
parser <- add_option(parser, c("--prefix"), type="character")
parser <- add_option(parser, c("--normalize"), type="logical")
outparse <- parse_args(parser)

model <- outparse$model
chr <- eval(parse(text=outparse$chr))
traits <- eval(parse(text=outparse$trait))
pheno_means <- outparse$pheno_means
pheno_dat <- outparse$pheno
geno_dat <- outparse$geno
test_ids <- outparse$test_ids
model_fit_dir <- outparse$model_fit_dir
sumstats_dir <- outparse$sumstats_dir
data_id <- outparse$data_id
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing
output_eff <- outparse$output_eff
output_pred_acc <- outparse$output_pred_acc
temp_dir <- outparse$temp_dir
prefix <- outparse$prefix
normalize <- outparse$normalize

###Set seed
set.seed(data_id)

###Read in data
options(datatable.fread.datatable=FALSE)
test_ids <- fread(test_ids, showProgress=FALSE, header=FALSE)

pheno <- readRDS(pheno_dat)$Y
test_inds_pheno <- which(rownames(pheno) %in% test_ids[,2]) ##Get only test individuals
pheno_test <- pheno[test_inds_pheno, ]

if(normalize){
  ###Function for quantile normalization
  inv_normalise <- function(x) { #this would also tolerate NAs
    return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
  }
  
  pheno_test <- apply(pheno_test, 2, inv_normalise)
}

geno_fam <- fread(paste0(unlist(strsplit(geno_dat, ".", fixed=TRUE))[1], ".fam"), showProgress=FALSE, header=FALSE)
test_inds_geno <- which(geno_fam[,2] %in% test_ids[,2]) ##Get only test individuals
tmp <- tempfile(tmpdir=temp_dir)
rds <- snp_readBed2(geno_dat, ind.row=test_inds_geno, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

pheno_means <- tryCatch(readRDS(pheno_means), 
                       error = function(e) {
                         return(NULL)
                       },
                       warning = function(w) {
                         return(NULL)
                       }
)


rm(list=c("geno_fam", "test_ids", "test_inds_geno", "test_inds_pheno", "pheno"))

###Impute missing genotypes
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

p <- geno$genotypes$ncol

###Compute predictions
##Get traits to analyze
r <- length(traits)

pheno_pred <- vector("list", length=length(chr))
Bhat_all <- vector("list", length=length(chr))

##Loop over chromosomes
it <- 0

for(i in chr){
  
  it <- it+1
  
  ###Load summary stats to get the SNP names
  univ_sumstats <- readRDS(paste0(sumstats_dir, prefix, "_chr", i, "_sumstats_matched_", data_id, ".rds"))
  
  inds <- which(geno$map$marker.ID %in% rownames(univ_sumstats$Bhat))
  
  if(model=="mr_mash_rss"){
    ##Read in model fit
    model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_", data_id, ".rds"))
    
    ##Compute predictions
    if(!is.null(pheno_means)){
      pheno_pred[[it]] <- t(t(big_prodMat(X, model_fit$mu1[, traits], ind.col=inds, ncores=ncores)) + model_fit$intercept[traits])
    } else {
      pheno_pred[[it]] <- big_prodMat(X, model_fit$mu1[, traits], ind.col=inds, ncores=ncores)
    }

    ##Store effects
    Bhat_all[[it]] <- model_fit$mu1[, traits]
    
  } else if(model %in% c("mvbayesN","mvbayesA","mvbayesL","mvbayesC","mvbayesC_rest","mvbayesR")){
    ##Read in model fit
    model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_", data_id, ".rds"))
    
    ##Compute predictions
    pheno_pred[[it]] <- big_prodMat(X, model_fit$bm[, traits], ind.col=inds, ncores=ncores)
    
    ##Store effects
    Bhat_all[[it]] <- model_fit$bm[, traits]
  
  } else if(model=="wmt_sblup"){
    model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_", data_id, ".rds"))
    
    ##Compute predictions
    pheno_pred[[it]] <- big_prodMat(X, model_fit$b[, traits], ind.col=inds, ncores=ncores)
    
    ##Store effects
    Bhat_all[[it]] <- model_fit$b[, traits]
  } else if(model %in% c("ldpred2_auto","mtag_ldpred2_auto","bayesN","bayesA","bayesL","bayesC","bayesR")){
    
    Bhat <- matrix(as.numeric(NA), nrow=length(inds), ncol=length(traits))
    
    for(j in traits){
      
      ##Read in model fit
      model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_trait", j, "_", data_id, ".rds"))
      
      if(model %in% c("ldpred2_auto", "mtag_ldpred2_auto")){
        
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
    
    ##Compute predictions
    pheno_pred[[it]] <- big_prodMat(X, Bhat, ind.col=inds, ncores=ncores)
    
    ##Store effects
    Bhat_all[[it]] <- Bhat
  }
}

###Compute predictions and accuracy
pheno_pred <- Reduce("+", pheno_pred)

##Remove mean of Y_training which was over counted in the per chromosome estimate of the intercept
if(!is.null(pheno_means)){
  pheno_pred <- t(t(pheno_pred) - (pheno_means[traits]*(length(chrs)-1)))
}

accuracy <- compute_accuracy(pheno_test, pheno_pred)

###Save accuracy to file
saveRDS(accuracy, file=output_pred_acc)

###Save effects to file
effects <- do.call("rbind", Bhat_all)
saveRDS(effects, file=output_eff)

###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))
