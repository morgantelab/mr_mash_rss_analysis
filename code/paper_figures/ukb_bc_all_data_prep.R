set.seed(1)

prefix <- "../output/prediction_accuracy/ukb_bc"
prefix_2 <- "../output/prediction_accuracy/ukb_bc_sampled"

repz <- 1:5
methodz <- c("mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped", "ldpred2_auto", "bayesR")
metric <- "r2"

n_col <- 5
n_row <- length(repz) * length(methodz) * 16

###Full data
pheno_full <- readRDS("../data/phenotypes/ukb_cleaned_bc_adjusted_pheno_test_1.rds")

res_full <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res_full) <- c("rep", "method", "trait", "metric", "score")

it <- 0

for(repp in repz){
  for(met in methodz){
    dat <- readRDS(paste0(prefix, "_", met, "_pred_acc_", repp, ".rds"))
    
    for(trait in 1:length(dat$r2)){
      it <- it + 1
      
      res_full[it, 1] <- repp
      res_full[it, 2] <- met
      res_full[it, 3] <- colnames(pheno_full)[trait]
      res_full[it, 4] <- metric
      res_full[it, 5] <- dat$r2[trait]
    }
  }
}

saveRDS(res_full, "../output/paper_figures_intermediate_data/ukb_bc_full_data.rds")


###Sampled data
pheno_sampled <- readRDS("../data/phenotypes/ukb_cleaned_bc_sampled_adjusted_pheno_test_1.rds")

res_sampled <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res_sampled) <- c("rep", "method", "trait", "metric", "score")

it <- 0

for(repp in repz){
  for(met in methodz){
    dat <- readRDS(paste0(prefix_2, "_", met, "_pred_acc_", repp, ".rds"))
    
    for(trait in 1:length(dat$r2)){
      it <- it + 1
      
      res_sampled[it, 1] <- repp
      res_sampled[it, 2] <- met
      res_sampled[it, 3] <- colnames(pheno_sampled)[trait]
      res_sampled[it, 4] <- metric
      res_sampled[it, 5] <- dat$r2[trait]
    }
  }
}

saveRDS(res_sampled, "../output/paper_figures_intermediate_data/ukb_bc_sampled_data.rds")
