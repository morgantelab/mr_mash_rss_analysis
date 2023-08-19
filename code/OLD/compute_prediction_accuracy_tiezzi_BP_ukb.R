###Load libraries needed
library(optparse)
library(bigreadr)
library(dplyr)
library(vctrs)
library(glue)
library(parallel)
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
parser <- add_option(parser, c("--geno_dir"), type="character")
parser <- add_option(parser, c("--sample_file"), type="character")
parser <- add_option(parser, c("--model_fit_dir"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical", default=TRUE)
parser <- add_option(parser, c("--output_eff"), type="character")
parser <- add_option(parser, c("--output_pred_acc"), type="character")
parser <- add_option(parser, c("--prefix"), type="character")
parser <- add_option(parser, c("--temp_dir"), type="character")
parser <- add_option(parser, c("--fold"), type="integer")
outparse <- parse_args(parser)

model <- outparse$model
chr <- outparse$chr
trait <- outparse$trait
pheno_dat <- outparse$pheno
geno_dir <- outparse$geno_dir
sample_file <- outparse$sample_file
model_fit_dir <- outparse$model_fit_dir
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing
output_eff <- outparse$output_eff
output_pred_acc <- outparse$output_pred_acc
temp_dir <- outparse$temp_dir
prefix <- outparse$prefix
fold <- outparse$fold

###Set seed
set.seed(1)

###Read in sample file
ids <- bigreadr::fread2(sample_file)
ids <- ids[-1, ]

###Read in phenotype
pheno_test <- readRDS(pheno_dat)

##Get chromosomes to analyze
chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

###Get list of HM3 SNPs to use
map_hapmap3 <- as.data.frame(readRDS("../data/misc/map_hm3.rds"))
ukb_data_location <- geno_dir

cl <- makeCluster(ncores)
clusterExport(cl, c("map_hapmap3", "ukb_data_location"))
list_snp_id <- parLapply(cl, chrs, function(chr) {
  mfi <- paste0(ukb_data_location, "ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
  joined <- dplyr::inner_join(cbind(chr = chr, infos_chr), map_hapmap3[, c("chr", "rsid")],
                              by = c("chr" = "chr", "V2" = "rsid"))
  with(joined[!vctrs::vec_duplicate_detect(joined$V2), ],
       paste(chr, V3, V4, V5, sep = "_"))
})
stopCluster(cl)

###Get only test individuals indeces and ids
test_inds_geno <- which(as.character(ids[,2]) %in% as.character(rownames(pheno_test))) ##Get only test individuals
test_ids_geno <- ids[which(as.character(ids[,2]) %in% as.character(rownames(pheno_test))), 2]

###Order pheno data according to geno data
pheno_test <-  pheno_test[as.character(test_ids_geno), ]

###Read in genotypes
tmp <- tempfile(tmpdir=temp_dir)
geno <- snp_readBGEN(
  bgenfiles   = glue::glue(paste0(ukb_data_location, "ukb22828_c{chr}_b0_v3.bgen"), chr = chrs), 
  backingfile = tmp,
  list_snp_id = list_snp_id,
  ind_row     = test_inds_geno,
  ncores      = ncores
)  

geno <- snp_attach(geno)

###Impute missing genotypes
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

p <- geno$genotypes$ncol

###Compute predictions
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
  if(i==0){
    inds <- seq_len(p)
    i <- "All"
  } else {
    if(i %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9")){
      ii <- paste0("0", i)
    } else {
      ii <- i
    }
    inds <- which(geno$map$chromosome==ii)
  }
  
  if(model %in% c("mr_mash_rss", "mr_mash_rss_init")){
    ##Read in model fit
    model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_", fold, ".rds"))
    
    ##Store effects
    Bhat <- model_fit$mu1[, traits]
    
    } else if(model %in% c("ldpred2_auto","bayesN","bayesA","bayesL","bayesC","bayesR")){
    
    it2 <- 0
    
    for(j in traits){
      
      it2 <- it2+1
      
      ##Read in model fit
      model_fit <- readRDS(paste0(model_fit_dir, prefix, "_chr", i, "_", model, "_fit_trait", j, "_", fold, ".rds"))
      
      if(model=="ldpred2_auto"){
        ##Quality control over chains
        range_corr <- sapply(model_fit, function(auto) diff(range(auto$corr_est)))
        to_keep <- (range_corr > (0.95 * quantile(range_corr, 0.95)))
        ##Compute posterior mean after QC
        if(it2==1){
          Bhat <- rowMeans(sapply(model_fit[to_keep], function(auto) auto$beta_est))
        } else {
          Bhat <- cbind(Bhat, rowMeans(sapply(model_fit[to_keep], function(auto) auto$beta_est)))
        }
      } else if(model %in% c("bayesN","bayesA","bayesL","bayesC","bayesR")){
        if(it2==1){
          Bhat <- model_fit$bm
        } else {
          Bhat <- cbind(Bhat, model_fit$bm)
        }
      }
    }
  }
  
  ##Store effects
  Bhat_all[[it]] <- Bhat
  
  ##Compute predictions
  pheno_pred[[it]] <- big_prodMat(X, Bhat, ind.col=inds)
}

###Compute predictions and accuracy
pheno_pred <- Reduce("+", pheno_pred)
accuracy <- compute_accuracy(pheno_test[, traits], pheno_pred)

###Save accuracy to file
saveRDS(accuracy, file=output_pred_acc)

###Save effects to file
effects <- do.call("rbind", Bhat_all)
saveRDS(effects, file=output_eff)

###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))
