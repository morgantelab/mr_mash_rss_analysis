###Load libraries needed
library(optparse)
library(RhpcBLASctl)

###Function needed
inv_normalise <- function(x) { #this would also tolerate NAs
  return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--normalize"), type="logical")
outparse <- parse_args(parser)

normalize <- outparse$normalize

###Set MKL threads
RhpcBLASctl::blas_set_num_threads(1)

###Set seed
set.seed(1)

###Read in data
pheno_dat <- readRDS("../data/phenotypes/ukb_tiezzi_cleaned_BP_covar_pheno.rds")
rownames(pheno_dat) <- pheno_dat$ID
pheno_dat$ID <- NULL

###Set up 
nfolds <- 5
traits <- c("Basophill_perc", "Eosinophill_perc", "Haemoglobin_conc", "HLR_perc", 
            "Lymphocyte_perc", "MCV", "Monocyte_perc", "MSCV", "Neutrophill_perc", "Platelet_dw", 
            "Platelet_count", "Platelet_crit", "RBC_count", "RBC_dw", "Reticulocyte_perc", 
            "WBC_count")

for(fold in 1:nfolds){
  ###Get training and test individuals indexes
  train_inds <- which(pheno_dat$fold != fold)
  test_inds <- which(pheno_dat$fold == fold)
  
  res_mat_train <- matrix(as.numeric(NA), nrow=length(train_inds), ncol=length(traits))
  rownames(res_mat_train) <- rownames(pheno_dat)[train_inds]
  res_mat_test <- matrix(as.numeric(NA), nrow=length(test_inds), ncol=length(traits))
  rownames(res_mat_test) <- rownames(pheno_dat)[test_inds]
  colnames(res_mat_train) <- colnames(res_mat_test) <- traits
 
  ###Loop over traits 
  i <- 0
  for(trait in traits){
    i <- i+1
    form <- as.formula(paste0(trait, "~ AOP + Sex_SI + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +   
                      PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
                      PC16 + PC17 + PC18 + PC19 + PC20"))
    
    ###Fit model to the training set and get residuals
    fit_train <- lm(form, data=pheno_dat, subset=train_inds)
    res_train <- resid(fit_train) + coef(fit_train)[1]

    ###Fit model to the test set and get residuals
    fit_test <- lm(form, data=pheno_dat, subset=test_inds)
    res_test <- resid(fit_test) + coef(fit_test)[1]

    if(normalize){
      res_mat_train[, i] <- inv_normalise(res_train)  
      res_mat_test[, i] <- inv_normalise(res_test)
    } else {
      res_mat_train[, i] <- res_train
      res_mat_test[, i] <- res_test
    }
  }
  
  saveRDS(cov(res_mat_train), file=paste0("../output/misc/ukb_tiezzi_BP_phenotypic_cov_", fold, ".rds"))
  saveRDS(res_mat_test, file=paste0("../data/phenotypes/ukb_tiezzi_cleaned_BP_adjusted_pheno_test_", fold, ".rds"))
}


