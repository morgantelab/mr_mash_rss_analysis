library(dplyr)

set.seed(1)

###Load data
res_full <- readRDS("../output/paper_figures_intermediate_data/ukb_bc_full_data.rds")
res_sampled <- readRDS("../output/paper_figures_intermediate_data/ukb_bc_sampled_data.rds")
pheno <- readRDS("../data/phenotypes/ukb_cleaned_bc_adjusted_pheno_test_1.rds")

options(pillar.sigfig = 4)

###Compute mean R2 by methods and traits
means_full <- res_full %>% group_by(method, trait) %>% summarise(avg = mean(score)) %>% as.data.frame(.)
means_sampled <- res_sampled %>% group_by(method, trait) %>% summarise(avg = mean(score)) %>% as.data.frame(.)

###Compute relative change in R2 
means_full_bayesR <- means_full[which(means_full$method=="bayesR"), ]
means_full_LDpred2 <- means_full[which(means_full$method=="ldpred2_auto"), ]
means_full_mrmashrss <- means_full[which(means_full$method=="mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped"), ]

change_bayesR_mrmashrss_full <- ((means_full_mrmashrss$avg - means_full_bayesR$avg)/means_full_bayesR$avg)*100
names(change_bayesR_mrmashrss_full) <- means_full_mrmashrss$trait 

summary(change_bayesR_mrmashrss_full)

change_LDpred2_mrmashrss_full <- ((means_full_mrmashrss$avg - means_full_LDpred2$avg)/means_full_LDpred2$avg)*100
names(change_LDpred2_mrmashrss_full) <- means_full_mrmashrss$trait

summary(change_LDpred2_mrmashrss_full)

means_sampled_bayesR <- means_sampled[which(means_sampled$method=="bayesR"), ]
means_sampled_LDpred2 <- means_sampled[which(means_sampled$method=="ldpred2_auto"), ]
means_sampled_mrmashrss <- means_sampled[which(means_sampled$method=="mr_mash_rss_sparse_LD_V_all_chr_bayesR_init_prior_finemapped"), ]

change_bayesR_mrmashrss_sampled <- ((means_sampled_mrmashrss$avg - means_sampled_bayesR$avg)/means_sampled_bayesR$avg)*100
names(change_bayesR_mrmashrss_sampled) <- means_sampled_mrmashrss$trait

summary(change_bayesR_mrmashrss_sampled)

change_LDpred2_mrmashrss_sampled <- ((means_sampled_mrmashrss$avg - means_sampled_LDpred2$avg)/means_sampled_LDpred2$avg)*100
names(change_LDpred2_mrmashrss_sampled) <- means_sampled_mrmashrss$trait 

summary(change_LDpred2_mrmashrss_sampled)

