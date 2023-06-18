###Commands to run these scripts (from /data2/morgante_lab/ukbiobank_projects/mr_mash_rss/run)

sbatch ../scripts/04_run_split_train_test_sample_filter_ukb_caucasian_white_british_unrel_list_plink.sbatch 11
sbatch --array=1-22 ../scripts/05_run_compute_geno_means_ukb_caucasian_white_british_unrel_by_chr.sbatch 11
sbatch ../scripts/05_run_sim_pheno_sample_filter_ukb_caucasian_white_british_unrel.sbatch 11
sbatch --array=1-22 ../scripts/06_run_compute_sumstats_sample_filter_ukb_caucasian_white_british_unrel_by_chr.sbatch 11
sbatch ../scripts/07_run_compute_pheno_means_ukb_caucasian_white_british_unrel.sbatch 11
sbatch ../scripts/07_run_compute_phenotypic_cov_ukb_caucasian_white_british_unrel.sbatch 11
sbatch ../scripts/07_run_compute_residual_cov_all_chr_ukb_caucasian_white_british_unrel.sbatch 11
sbatch --array=1-22 ../scripts/08_run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.sbatch 11 1 #trait
sbatch --array=1-22 ../scripts/08_run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.sbatch 11 2 #trait
sbatch --array=1-22 ../scripts/08_run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.sbatch 11 3 #trait
sbatch --array=1-22 ../scripts/08_run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.sbatch 11 4 #trait
sbatch --array=1-22 ../scripts/08_run_fit_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel_by_chr_and_trait.sbatch 11 5 #trait
sbatch ../scripts/08_run_compute_data_driven_cov_all_chr_ukb_caucasian_white_british_unrel.sbatch 11
sbatch --array=1-22 ../scripts/09_run_fit_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel_by_chr.sbatch 11
sbatch ../scripts/09_run_compute_pred_accuracy_ldpred2auto_sample_filter_ukb_caucasian_white_british_unrel.sbatch 11
sbatch ../scripts/10_run_compute_pred_accuracy_mrmashrss_sample_filter_ukb_caucasian_white_british_unrel.sbatch 11

