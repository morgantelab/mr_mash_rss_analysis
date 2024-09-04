set.seed(1)


prefix <- "../output/prediction_accuracy/ukb_caucasian_white_british_unrel_100000"
metric <- "r2"

scenarioz <- c("external_LD_equal_effects_indep_resid", "missing_pheno_02_full_equal_effects_indep_resid",
               "missing_pheno_08_full_equal_effects_indep_resid")
methodz <- c("mr_mash_rss", "mvbayesC", "mvbayesC_rest", "wmt_sblup", "mtag_ldpred2_auto", 
             "ldpred2_auto", "bayesR")

nreps <- 20

n_col <- 6
n_row <- nreps * length(methodz) * length(scenarioz) * 5

res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "scenario", "method", "trait", "metric", "score")


###External LD
sce <- scenarioz[1]
repz <- c(1:7, 9:20, 22)

i <- 0

for(met in methodz){
  for(repp in repz){
    dat <- readRDS(paste0(prefix, "_", sce, "_", met, "_pred_acc_", repp, ".rds"))
      
    for(trait in 1:length(dat$r2)){
      i <- i + 1
        
      res[i, 1] <- repp
      res[i, 2] <- sce
      res[i, 3] <- met
      res[i, 4] <- trait
      res[i, 5] <- metric
      res[i, 6] <- dat$r2[trait]
    }
  }
}

###External LD
repz <- 1:21

for(sce in scenarioz[2:3]){
  for(met in methodz){
    for(repp in repz){
      if(sce=="missing_pheno_02_full_equal_effects_indep_resid" && repp==18){
        next
      }

      if(sce=="missing_pheno_08_full_equal_effects_indep_resid" && repp==21){
        next
      }
      
      dat <- readRDS(paste0(prefix, "_", sce, "_", met, "_pred_acc_", repp, ".rds"))
      
      for(trait in 1:length(dat$r2)){
        i <- i + 1
        
        res[i, 1] <- repp
        res[i, 2] <- sce
        res[i, 3] <- met
        res[i, 4] <- trait
        res[i, 5] <- metric
        res[i, 6] <- dat$r2[trait]
      }
    }
  }
}


saveRDS(res, "../output/paper_figures_intermediate_data/ukb_caucasian_white_british_unrel_100000_robustness_sims_scenarios.rds")


