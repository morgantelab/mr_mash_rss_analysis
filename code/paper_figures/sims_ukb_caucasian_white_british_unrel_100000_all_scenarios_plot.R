###Load libraries
library(ggplot2)
library(cowplot)

set.seed(1)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dat <- readRDS("../output/paper_figures_intermediate_data/ukb_caucasian_white_british_unrel_100000_all_sims_scenarios.rds")


dat_equal_effects <- dat[which(dat$scenario=="equal_effects_indep_resid"), ]
dat_mostly_null <- dat[which(dat$scenario=="trait_1_only_effects_indep_resid"), ]
dat_shared_subgroups <- dat[which(dat$scenario=="blocks_shared_effects_indep_resid"), ]
dat_low_h2g <- dat[which(dat$scenario=="equal_effects_low_pve_indep_resid"), ]
dat_highly_polygenic <- dat[which(dat$scenario=="equal_effects_50000causal_indep_resid"), ]
dat_more_traits <- dat[which(dat$scenario=="equal_effects_10traits_indep_resid"), ]


###Equal Effects
res_equal_effects <- transform(dat_equal_effects, scenario=as.factor(scenario),
                                    method=factor(method, 
                                                  levels = c("bayesR", "ldpred2_auto", 
                                                             "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                             "mtag_ldpred2_auto", "wmt_sblup")),
                                    trait=as.factor(trait))

p_equal_effects <- ggplot(res_equal_effects, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="Equal Effects") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_equal_effects)

###Mostly null
res_mostly_null <- transform(dat_mostly_null, scenario=as.factor(scenario),
                             method=factor(method, 
                                           levels = c("bayesR", "ldpred2_auto", 
                                                      "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                      "mtag_ldpred2_auto", "wmt_sblup")),
                             trait=as.factor(trait))

p_mostly_null <- ggplot(res_mostly_null, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="Mostly Null") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_mostly_null)

###Shared in Subgroups
res_shared_subgroups <- transform(dat_shared_subgroups, scenario=as.factor(scenario),
                                  method=factor(method, 
                                                levels = c("bayesR", "ldpred2_auto", 
                                                           "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                           "mtag_ldpred2_auto", "wmt_sblup")),
                                  trait=as.factor(trait))

p_shared_subgroups <- ggplot(res_shared_subgroups, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="Shared Effects in Subgroups") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  geom_hline(yintercept=0.3, linetype="dashed", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_shared_subgroups)

###Low h2g
res_low_h2g <- transform(dat_low_h2g, scenario=as.factor(scenario),
                         method=factor(method, 
                                       levels = c("bayesR", "ldpred2_auto", 
                                                  "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                  "mtag_ldpred2_auto", "wmt_sblup")),
                         trait=as.factor(trait))

p_low_h2g <- ggplot(res_low_h2g, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title=expression(bold("Low") ~ bolditalic(h["g"]^2))) +
  
  geom_hline(yintercept=0.3, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_low_h2g)

###High Polygenicity
res_highly_polygenic <- transform(dat_highly_polygenic, scenario=as.factor(scenario),
                                  method=factor(method, 
                                                levels = c("bayesR", "ldpred2_auto", 
                                                           "mr_mash_rss", "mvbayesC", "mvbayesC_rest", 
                                                           "mtag_ldpred2_auto", "wmt_sblup")),
                                  trait=as.factor(trait))

p_highly_polygenic <- ggplot(res_highly_polygenic, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette, labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC", "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="High Polygenicity") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_highly_polygenic)

###More Phenotypes
res_more_traits <- transform(dat_more_traits, scenario=as.factor(scenario),
                             method=factor(method, 
                                           levels = c("bayesR", "ldpred2_auto", 
                                                      "mr_mash_rss", "mvbayesC_rest", 
                                                      "mtag_ldpred2_auto", "wmt_sblup")),
                             trait=as.factor(trait))

p_more_traits <- ggplot(res_more_traits, aes(x = trait, y = score, fill = method)) +
  geom_boxplot(color = "black", outlier.size = 0.5, width = 0.85) +
  # stat_summary(fun=mean, geom="point", shape=23,
  #              position = position_dodge2(width = 0.87,   
  #                                         preserve = "single")) +
  ylim(0, 0.51) +
  scale_fill_manual(values = cbPalette[-4], labels = c("SBayesR", "LDpred2", "mr.mash-rss", 
                                                   "SmvBayesC-rest", 
                                                   "MTAG+LDpred2", "wMT-SBLUP")) +
  labs(x = "Phenotype", y = expression(italic(R)^2), fill="Method", title="More Phenotypes") +
  geom_hline(yintercept=0.5, linetype="dotted", linewidth=1, color = "black") +
  theme_cowplot(font_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=14))

#print(p_more_traits)

###Figure 1
#Extract legend
legend_1 <- get_legend(
  p_equal_effects + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position=c(0.04,0.8),legend.direction = "horizontal")
)


#Make the multi panel plot 
p_1 <- plot_grid(p_equal_effects + theme(legend.position="none"),
                 p_mostly_null + theme(legend.position="none"),
                 p_shared_subgroups + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("A", "B", "C"),
                 hjust = -1,
                 nrow = 1
)

p_1 <- plot_grid(p_1, legend_1, ncol = 1, rel_heights = c(1, .1))

ggsave("../analysis/paper_figures/Fig1.eps", plot=p_1, device="eps", units="in", height=5, width=11)

###Figure 2
#Extract legend
legend_2 <- get_legend(
  p_low_h2g + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position=c(0.04,0.8),legend.direction = "horizontal")
)


#Make the multi panel plot 
p_2 <- plot_grid(p_low_h2g + theme(legend.position="none"),
                 p_highly_polygenic + theme(legend.position="none"),
                 p_more_traits + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("A", "B", "C"),
                 hjust = -1,
                 nrow = 1
)

p_2 <- plot_grid(p_2, legend_2, ncol = 1, rel_heights = c(1, .1))
ggsave("../analysis/paper_figures/Fig2.eps", plot=p_2, device="eps", units="in", height=5, width=11)


