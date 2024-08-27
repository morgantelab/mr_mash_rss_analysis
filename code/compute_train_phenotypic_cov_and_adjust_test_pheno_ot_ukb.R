###Load libraries needed
library(optparse)

###Function needed
inv_normalise <- function(x) { #this would also tolerate NAs
  return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--normalize"), type="logical")
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--traits"), type="character")
parser <- add_option(parser, c("--output_pheno"), type="character")
parser <- add_option(parser, c("--output_cov"), type="character")
parser <- add_option(parser, c("--fold"), type="integer")
outparse <- parse_args(parser)

normalize <- outparse$normalize
pheno <- outparse$pheno
traits <- eval(parse(text=outparse$traits))
output_pheno <- outparse$output_pheno
output_cov <- outparse$output_cov
fold <- outparse$fold

###Set seed
set.seed(1)

###Read in data
pheno_dat <- readRDS(pheno)
rownames(pheno_dat) <- pheno_dat$id
pheno_dat$id <- NULL

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
  form <- as.formula(paste0(trait, "~ as.factor(sex) + as.factor(assessment_centre) + age + 
                            I(age^2) + as.factor(genotype_measurement_batch) + pc_genetic1 +
                            pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 + pc_genetic6 +
                            pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10"))
    
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
  
saveRDS(cov(res_mat_train), file=output_cov)
saveRDS(res_mat_test, file=output_pheno)



