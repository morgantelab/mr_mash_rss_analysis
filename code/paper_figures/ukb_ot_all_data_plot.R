###Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

set.seed(1)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###Full data
##Prediction accuracy
res_full1 <- readRDS("../output/paper_figures_intermediate_data/ukb_ot_full_data.rds")

res_full <- transform(res_full1,
                 method=factor(method, levels=c("bayesR", "ldpred2_auto", 
                                                "mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped",
                                                "mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped_mash_weights"),
                               labels=c("SBayesR", "LDpred2", "mr.mash-rss", "mr.mash-rss mash")),
                 trait=factor(trait, levels=c("DPa", "SPa", "BMI", "weight", "hip", "waist", "BFP", "TFM"),
                              labels=c("DP", "SP", "BMI", "weight", "hip", "waist", "BFP", "TFM")))

p <- ggplot(res_full, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.75) +
  scale_fill_manual(values = cbPalette) +
  ylim(0.05, 0.145) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.80,0.9))

#print(p)

ggsave("../analysis/paper_figures/FigS3.pdf", plot=p, device="pdf", units="in", height=9, width=11)

