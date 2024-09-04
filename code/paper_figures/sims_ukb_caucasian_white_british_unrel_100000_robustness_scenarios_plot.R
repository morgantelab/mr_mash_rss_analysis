###Load libraries
library(ggplot2)
library(cowplot)

set.seed(1)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dat <- readRDS("../output/paper_figures_intermediate_data/ukb_caucasian_white_british_unrel_100000_robustness_sims_scenarios.rds")


dat_external_LD <- dat[which(dat$scenario=="external_LD_equal_effects_indep_resid"), ]
dat_missing_pheno_02 <- dat[which(dat$scenario=="missing_pheno_02_full_equal_effects_indep_resid"), ]
dat_missing_pheno_08 <- dat[which(dat$scenario=="missing_pheno_08_full_equal_effects_indep_resid"), ]

###External LD
res_external_LD <- transform(dat_external_LD, scenario=as.factor(scenario),
                                    method=factor(method, 
                                                  levels = c("bayesR", "ldpred2_auto", 
                                                             "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                             "mtag_ldpred2_auto", "wmt_sblup")),
                                    trait=as.factor(trait))

p_external_LD <- ggplot(res_external_LD, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="External LD") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))


#print(p_external_LD)

###Fig S1
ggsave("../analysis/paper_figures/FigS1.pdf", plot=p_external_LD, device="pdf", units="in", height=5, width=8)



###Missing pheno 0.2
res_missing_pheno_02 <- transform(dat_missing_pheno_02, scenario=as.factor(scenario),
                             method=factor(method, 
                                           levels = c("bayesR", "ldpred2_auto", 
                                                      "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                      "mtag_ldpred2_auto", "wmt_sblup")),
                             trait=as.factor(trait))

p_missing_pheno_02 <- ggplot(res_missing_pheno_02, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="20% individuals with missing values") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_missing_pheno_02)

###Missing pheno 0.8
res_missing_pheno_08 <- transform(dat_missing_pheno_08, scenario=as.factor(scenario),
                                  method=factor(method, 
                                                levels = c("bayesR", "ldpred2_auto", 
                                                           "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                           "mtag_ldpred2_auto", "wmt_sblup")),
                                  trait=as.factor(trait))

p_missing_pheno_08 <- ggplot(res_missing_pheno_08, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="80% individuals with missing values") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_missing_pheno_08)


###Figure S2
#Extract legend
legend_1 <- get_legend(
  p_missing_pheno_02 + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position=c(0.04,0.8),legend.direction = "horizontal")
)

#Make the multi panel plot 
p_1 <- plot_grid(p_missing_pheno_02 + theme(legend.position="none"),
                 p_missing_pheno_08 + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("A", "B"),
                 hjust = -1,
                 nrow = 1
)

p_1 <- plot_grid(p_1, legend_1, ncol = 1, rel_heights = c(1, .1))

ggsave("../analysis/paper_figures/FigS2.pdf", plot=p_1, device="pdf", units="in", height=5, width=11)



