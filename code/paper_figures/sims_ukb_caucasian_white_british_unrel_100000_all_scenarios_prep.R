set.seed(1)

repz <- 1:20
prefix <- "../output/prediction_accuracy/ukb_caucasian_white_british_unrel_100000"
metric <- "r2"

scenarioz <- c("equal_effects_indep_resid", "trait_1_only_effects_indep_resid",
               "blocks_shared_effects_indep_resid", "equal_effects_low_pve_indep_resid",
               "equal_effects_50000causal_indep_resid", "equal_effects_10traits_indep_resid")
methodz <- c("mr_mash_rss", "mvbayesC", "mvbayesC_rest", "wmt_sblup", "mtag_ldpred2_auto", 
             "ldpred2_auto", "bayesR")

n_col <- 6
n_row <- length(repz)
n_row <- (n_row * length(methodz) * (length(scenarioz) - 1) * 5) + n_row * (length(methodz) - 1) * 10


res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "scenario", "method", "trait", "metric", "score")

i <- 0

for(sce in scenarioz){
  for(met in methodz){
    for(repp in repz){
      if(sce=="equal_effects_10traits_indep_resid" && met=="mvbayesC"){
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

saveRDS(res, "../output/paper_figures_intermediate_data/ukb_caucasian_white_british_unrel_100000_all_sims_scenarios.rds")


