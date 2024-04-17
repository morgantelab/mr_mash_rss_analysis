###Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

set.seed(1)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###Full data
##Prediction accuracy
res_full1 <- readRDS("../output/paper_figures_intermediate_data/ukb_bc_full_data.rds")

res_full <- transform(res_full1,
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

###Percentage increase in R2 vs h2g
h2_dat <- readRDS("../output/paper_figures_intermediate_data/ukb_bc_full_data_h2.rds")
means_full <- res_full1 %>% group_by(method, trait) %>% summarise(avg = mean(score)) %>% as.data.frame(.)
means_full_bayesR <- means_full[which(means_full$method=="bayesR"), ]
means_full_LDpred2 <- means_full[which(means_full$method=="ldpred2_auto"), ]
means_full_mrmashrss <- means_full[which(means_full$method=="mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped"), ]

change_bayesR_mrmashrss_full <- ((means_full_mrmashrss$avg - means_full_bayesR$avg)/means_full_bayesR$avg)*100
names(change_bayesR_mrmashrss_full) <- means_full_mrmashrss$trait 
change_LDpred2_mrmashrss_full <- ((means_full_mrmashrss$avg - means_full_LDpred2$avg)/means_full_LDpred2$avg)*100
names(change_LDpred2_mrmashrss_full) <- means_full_mrmashrss$trait

dat_change_LDpred2_mrmashrss_full <- data.frame(trait=names(change_LDpred2_mrmashrss_full), perc_change=change_LDpred2_mrmashrss_full)
dat_LDpred2_full <- inner_join(dat_change_LDpred2_mrmashrss_full, h2_dat, by=join_by(trait))
dat_change_bayesR_mrmashrss_full <- data.frame(trait=names(change_bayesR_mrmashrss_full), perc_change=change_bayesR_mrmashrss_full)
dat_bayesR_full <- inner_join(dat_change_bayesR_mrmashrss_full, h2_dat, by=join_by(trait))

###Plot of relative improvement of mr.mash over the enet based on RMSE vs sample size by tissue
p_bayesR <- ggplot(dat_bayesR_full, aes(x=h2, y=perc_change)) + 
  geom_point(shape=16) +
  geom_smooth(method='lm', formula= y~x) +
  labs(x = expression(italic(h["g"]^2)), y = expression("Relative difference in" ~ italic(R)^2), title="mr.mash-rss vs SBayesR") + 
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

p_LDpred2 <- ggplot(dat_LDpred2_full, aes(x=h2, y=perc_change)) + 
  geom_point(shape=16) +
  geom_smooth(method='lm', formula= y~x) +
  labs(x = expression(italic(h["g"]^2)), y = expression("Relative difference in" ~ italic(R)^2), title="mr.mash-rss vs LDpred2") + 
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

p_4 <- plot_grid(p_LDpred2, p_bayesR,
                 align = 'vh',
                 labels = c("A", "B"),
                 hjust = -1,
                 nrow = 1
)

#print(p_4)

ggsave("../analysis/paper_figures/Fig4.pdf", plot=p_4, device="pdf", units="in", height=5, width=11)

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

p_5 <- ggplot(res_sampled, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.75) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,
  #                                         preserve = "single")) +
  scale_fill_manual(values = cbPalette) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position=c(0.8,0.9))

#print(p_5)

ggsave("../analysis/paper_figures/Fig5.pdf", plot=p_5, device="pdf", units="in", height=9, width=11)
