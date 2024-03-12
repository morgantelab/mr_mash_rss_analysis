###Load libraries
library(ggplot2)
library(cowplot)

set.seed(1)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###Full data
res_full <- readRDS("../output/paper_figures_intermediate_data/ukb_bc_full_data.rds")

res_full <- transform(res_full,
                 method=factor(method, levels=c("bayesR", "ldpred2_auto", 
                                                "mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped"),
                               labels=c("SBayesR", "LDpred2", "mr.mash-rss")),
                 trait=factor(trait, levels=c("RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV", 
                                                "Reticulocyte_perc", "HLR_perc",
                                                "Platelet_count", "Plateletcrit", "PDW", 
                                                "WBC_count", "Lymphocyte_perc", "Monocyte_perc", "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc"),
                              labels=c("RBC#", "HGB", "MCV", "RDW", "MSCV", 
                                       "RET%", "HLR%",
                                       "PLT#", "PCT", "PDW", 
                                       "WBC#", "LYMPH%", "MONO%", "NEUT%", "EO%", "BASO%")))

p_3 <- ggplot(res_full, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.75) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,
  #                                         preserve = "single")) +
  scale_fill_manual(values = cbPalette) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.80,0.9))

#print(p_3)

ggsave("../analysis/paper_figures/Fig3.pdf", plot=p_3, device="pdf", units="in", height=9, width=11)


###Sampled data
res_sampled <- readRDS("../output/paper_figures_intermediate_data/ukb_bc_sampled_data.rds")

res_sampled <- transform(res_sampled,
                 method=factor(method, levels=c("bayesR", "ldpred2_auto", 
                                                "mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped"),
                               labels=c("SBayesR", "LDpred2", "mr.mash-rss")),
                 trait=factor(trait, levels=c("RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV", 
                                              "Reticulocyte_perc", "HLR_perc",
                                              "Platelet_count", "Plateletcrit", "PDW", 
                                              "WBC_count", "Lymphocyte_perc", "Monocyte_perc", "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc"),
                              labels=c("RBC#", "HGB", "MCV", "RDW", "MSCV", 
                                       "RET%", "HLR%",
                                       "PLT#", "PCT", "PDW", 
                                       "WBC#", "LYMPH%", "MONO%", "NEUT%", "EO%", "BASO%")))

p_4 <- ggplot(res_sampled, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.75) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,
  #                                         preserve = "single")) +
  scale_fill_manual(values = cbPalette) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.8,0.9))

#print(p_4)

ggsave("../analysis/paper_figures/Fig4.pdf", plot=p_4, device="pdf", units="in", height=9, width=11)

