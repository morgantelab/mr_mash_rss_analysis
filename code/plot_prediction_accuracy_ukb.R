###Load libraries
library(ggplot2)
library(cowplot)

repz <- 1:20
prefix <- "../output/prediction_accuracy/ukb_caucasian_white_british_unrel_100000"
scenarioz <- "equal_effects_indep_resid"
methodz <- c("mr_mash_rss", "ldpred2_auto", "ldpred2_auto_gwide")
metric <- "r2"
traitz <- 1:5

i <- 0

n_col <- 6
n_row <- length(repz) * length(scenarioz) * length(methodz) * length(traitz)
res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "scenario", "method", "trait", "metric", "score")



for(sce in scenarioz){
  for(met in methodz){
    for(repp in repz){
      dat <- readRDS(paste0(prefix, "_", sce, "_", met, "_pred_acc_", repp, ".rds"))
      
      for(trait in traitz){
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

res <- transform(res, scenario=as.factor(scenario),
                      method=as.factor(method),
                      trait=as.factor(trait))

p_methods_shared <- ggplot(res, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  stat_summary(fun=mean, geom="point", shape=23,
               position = position_dodge2(width = 0.87,   
                                          preserve = "single")) +
  ylim(0.3, 0.51) +
  scale_fill_manual(values = c("pink", "red", "green"))+ #, labels = c("g-lasso", "smt-lasso", "e-net")) +
  labs(x = "Trait", y = expression(italic(R)^2), title = "Equal effects", fill="Method") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

print(p_methods_shared)
