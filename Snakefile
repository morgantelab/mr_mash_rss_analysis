rule all:
  input:
    expand(config["DEST"]+"output/estimated_effects/"+config["DATA_PREFIX"]+"_mr_mash_rss_effects_{DATASET}.rds",DATASET=config["DATASET"]),
    expand(config["DEST"]+"output/prediction_accuracy/"+config["DATA_PREFIX"]+"_mr_mash_rss_pred_acc_{DATASET}.rds",DATASET=config["DATASET"]),
    expand(config["DEST"]+"output/estimated_effects/"+config["DATA_PREFIX"]+"_ldpred2_auto_effects_{DATASET}.rds",DATASET=config["DATASET"]),
    expand(config["DEST"]+"output/prediction_accuracy/"+config["DATA_PREFIX"]+"_ldpred2_auto_pred_acc_{DATASET}.rds",DATASET=config["DATASET"]),
    expand(config["DEST"]+"output/estimated_effects/"+config["DATA_PREFIX"]+"_ldpred2_auto_gwide_effects_{DATASET}.rds",DATASET=config["DATASET"]),
    expand(config["DEST"]+"output/prediction_accuracy/"+config["DATA_PREFIX"]+"_ldpred2_auto_gwide_pred_acc_{DATASET}.rds",DATASET=config["DATASET"])

rule split_train_test:
  input:
    config["FAM"]
  output:
    config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt"
  params:
    script=config["SCRIPT"]+"split_train_test_ukb_list_plink.R",
    n=config["N_TRAIN"],
    dataset="{DATASET}"
  resources: cpus=1, mem_mb=2048, time_min=10
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        Rscript {params.script} \
                                  --input {input} \
                                  --n_train {params.n} \
                                  --seed {params.dataset} \
                                  --output {output}
        module unload R/4.2.3
        """

rule compute_genotype_means:
  input:
    geno=config["BED"],
    tids=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt"
  output:
    config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_chr{CHR}_geno_means_{DATASET}.rds"
  params:
    script=config["SCRIPT"]+"compute_genotype_means_ukb.R",
    dataset="{DATASET}",
    imp_flag=config["IMPUTE_FLAG"],
    tmp=config["TMP"],
    chr="{CHR}"
  resources: cpus=4, mem_mb=100000, time_min=1440
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        echo "Start analyzing chromosome {params.chr} ..."
        Rscript {params.script} \
                                  --geno {input.geno} \
                                  --test_ids {input.tids} \
                                  --output {output} \
                                  --seed {params.dataset} \
                                  --chr {params.chr} \
                                  --impute_missing {params.imp_flag} \
                                  --temp_dir {params.tmp} \
                                  --ncores 4
        module unload R/4.2.3
        """

rule simulate_phenotype:
  input:
    config["BED"]
  output:
    config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_pheno_{DATASET}.rds"
  params:
    p_causal=config["P_CAUSAL"],
    r=config["R"],
    pve=config["PVE"],
    B_cor=config["B_COR"],
    V_cor=config["V_COR"],
    dataset="{DATASET}",
    tmp=config["TMP"],
    script=config["SCRIPT"]+"sim_pheno_ukb.R"
  resources: cpus=16, mem_mb=200000, time_min=240
  shell:
        """
        module load gsl/2.7
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        Rscript {params.script} \
                                    --geno {input} \
                                    --p_causal {params.p_causal} \
                                    --r {params.r} \
                                    --pve {params.pve} \
                                    --B_cor {params.B_cor} \
                                    --V_cor {params.V_cor} \
                                    --seed {params.dataset} \
                                    --temp_dir {params.tmp} \
                                    --output {output}
        module unload R/4.2.3
        module unload gsl/2.7
        """

rule compute_summary_stats:
  input:
    pheno=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_pheno_{DATASET}.rds",
    geno=config["BED"],
    tids=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt"
  output:
    config["DEST"]+"output/summary_statistics/"+config["DATA_PREFIX"]+"_chr{CHR}_sumstats_{DATASET}.rds"
  params:
    dataset="{DATASET}",
    chr="{CHR}",
    stand=config["STANDARDIZE"],
    imp_flag=config["IMPUTE_FLAG"],
    tmp=config["TMP"],
    script=config["SCRIPT"]+"compute_sumstats_ukb.R"
  resources: cpus=2, mem_mb=100000, time_min=1440
  shell:
        """
        module load gsl/2.7
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        echo "Start analyzing chromosome {params.chr} ..."
        Rscript {params.script} \
                                    --pheno {input.pheno} \
                                    --geno {input.geno} \
                                    --test_ids {input.tids} \
                                    --output {output} \
                                    --seed {params.dataset} \
                                    --chr {params.chr} \
                                    --standardize {params.stand} \
                                    --impute_missing {params.imp_flag} \
                                    --temp_dir {params.tmp} \
                                    --ncores 2
        echo "Finished analyzing chromosome {params.chr}."
        module unload R/4.2.3
        module unload gsl/2.7
        """

rule compute_phenotype_means:
  input:
    pheno=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_pheno_{DATASET}.rds",
    tids=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt"
  output:
    config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_pheno_means_{DATASET}.rds"
  params:
    dataset="{DATASET}",
    script=config["SCRIPT"]+"compute_phenotype_means_ukb.R"
  resources: cpus=1, mem_mb=2000, time_min=30
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        Rscript {params.script} \
                                    --pheno {input.pheno} \
                                    --test_ids {input.tids} \
                                    --output {output} \
                                    --seed {params.dataset}
        module unload R/4.2.3
        """

rule compute_phenotype_covariance:
  input:
    pheno=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_pheno_{DATASET}.rds",
    tids=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt"
  output:
    config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_phenotypic_cov_{DATASET}.rds"
  params:
    dataset="{DATASET}",
    script=config["SCRIPT"]+"compute_phenotypic_cov_ukb_sim.R"
  resources: cpus=8, mem_mb=6000, time_min=120
  shell:
        """
        module load gsl/2.7
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        Rscript {params.script} \
                                    --pheno {input.pheno} \
                                    --test_ids {input.tids} \
                                    --seed {params.dataset} \
                                    --output {output}
        module unload R/4.2.3
        module unload gsl/2.7
        """

def get_sumstats_input_files(wildcards):
    DATASET = wildcards.DATASET
    return expand(config["DEST"]+"output/summary_statistics/"+config["DATA_PREFIX"]+"_chr{CHR}_sumstats_{DATASET}.rds", DATASET=DATASET, CHR=config["CHR"])

rule compute_residual_covariance_all_chr:
  input:
    lambda wildcards: get_sumstats_input_files(wildcards)
  output:
    config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_residual_cov_{DATASET}.rds"
  params:
    prefix="output/summary_statistics/chr",
    suffix="_sumstats_{DATASET}.rds",
    chr_r=config["CHR_RANGE"],
    dataset="{DATASET}",
    script=config["SCRIPT"]+"compute_residual_cov_all_chr_ukb.R"
  resources: cpus=1, mem_mb=4096, time_min=60
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        Rscript {params.script} \
                                    --sumstats_prefix {params.prefix} \
                                    --sumstats_suffix {params.suffix} \
                                    --chr {params.chr_r} \
                                    --output {output} \
                                    --seed {params.dataset}
        module unload R/4.2.3
        """

rule compute_data_driven_covariance_all_chr:
  input:
    sum_stat=lambda wildcards: get_sumstats_input_files(wildcards),
    res_cov=config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_residual_cov_{DATASET}.rds"
  output:
    config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_"+config["DD_ED"]+"_data_driven_cov_{DATASET}.rds"
  params:
    prefix="output/summary_statistics/chr",
    suffix="_sumstats_{DATASET}.rds",
    chr_r=config["CHR_RANGE"],
    str_z_thr=config["DD_STRONG_Z_THR"],
    n_pcs=config["DD_N_PCS"],
    frs=config["DD_FLASH_REMOVE_SINGLE"],
    ed=config["DD_ED"],
    dataset="{DATASET}",
    script=config["SCRIPT"]+"compute_data_driven_cov_all_chr_ukb.R"
  resources: cpus=8, mem_mb=7000, time_min=1440
  shell:
        """
        module load gsl/2.7
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        Rscript {params.script} \
                                    --sumstats_prefix {params.prefix} \
                                    --sumstats_suffix {params.suffix} \
                                    --chr {params.chr_r} \
                                    --residual_cov {input.res_cov} \
                                    --strong_Z_thresh {params.str_z_thr} \
                                    --n_PCs {params.n_pcs} \
                                    --flash_remove_singleton {params.frs} \
                                    --ED_algorithm {params.ed} \
                                    --output {output} \
                                    --seed {params.dataset}
        module unload R/4.2.3
        module unload gsl/2.7
        """

rule fit_ldpred2auto_by_chr_trait:
  input:
    sumstats=config["DEST"]+"output/summary_statistics/"+config["DATA_PREFIX"]+"_chr{CHR}_sumstats_{DATASET}.rds",
    ldmat=config["LDMAT"]+config["LDMAT_PREFIX"]+"{CHR}_LD.ld.bin"
  output:
    config["DEST"]+"output/ldpred2_auto_fit/"+config["DATA_PREFIX"]+"_chr{CHR}_ldpred2_auto_fit_trait{TRAIT}_{DATASET}.rds"
  params:
    n=config["N_TRAIN"],
    h2_init=config["H2_INIT"],
    burn_in=config["LDPRED_BURN_IN"],
    num_iter=config["LDPRED_ITER"],
    trait="{TRAIT}",
    dataset="{DATASET}",
    chr="{CHR}",
    script=config["SCRIPT"]+"fit_ldpred2auto_ukb.R"
  resources: cpus=4, mem_mb=100000, time_min=720
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        echo "Start analyzing trait {params.trait} and chromosome {params.chr} ..."
        Rscript {params.script} \
                                    --sumstats {input.sumstats} \
                                    --LD_matrix {input.ldmat} \
                                    --n {params.n} \
                                    --h2_init {params.h2_init} \
                                    --burn_in {params.burn_in} \
                                    --num_iter {params.num_iter} \
                                    --trait {params.trait} \
                                    --ncores 4 \
                                    --seed {params.dataset} \
                                    --output {output}
        module unload R/4.2.3
        module unload gsl/2.7
        """

rule fit_mrmash_rss:
  input:
    sumstats=config["DEST"]+"output/summary_statistics/"+config["DATA_PREFIX"]+"_chr{CHR}_sumstats_{DATASET}.rds",
    ldmat=config["LDMAT"]+config["LDMAT_PREFIX"]+"{CHR}_LD.ld.bin",
    pheno_cov=config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_phenotypic_cov_{DATASET}.rds",
    res_cov=config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_residual_cov_{DATASET}.rds",
    dd_cov=config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_"+config["DD_ED"]+"_data_driven_cov_{DATASET}.rds",
    x_colm=config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_chr{CHR}_geno_means_{DATASET}.rds",
    y_colm=config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_pheno_means_{DATASET}.rds"
  output:
    config["DEST"]+"output/mr_mash_rss_fit/"+config["DATA_PREFIX"]+"_chr{CHR}_mr_mash_rss_fit_{DATASET}.rds"
  params:
    can_cov=config["MRMASH_CAN_COV"],
    n=config["N_TRAIN"],
    prop_nz=config["MRMASH_PROP_NONZERO"],
    stand=config["STANDARDIZE"],
    upd_v_met=config["MRMASH_UPDATE_V_METHOD"],
    w0_thresh=config["MRMASH_W0_THRESH"],
    script=config["SCRIPT"]+"fit_mrmashrss_ukb.R",
    chr="{CHR}",
    dataset="{DATASET}"
  resources: cpus=4, mem_mb=100000, time_min=720
  shell:
        """
        module load gsl/2.7
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        echo "Start analyzing chromosome {params.chr} ..."
        Rscript {params.script} \
                                    --sumstats {input.sumstats} \
                                    --LD_matrix {input.ldmat} \
                                    --pheno_cov {input.pheno_cov} \
                                    --residual_cov {input.res_cov} \
                                    --data_driven_cov {input.dd_cov} \
                                    --canonical_cov {params.can_cov} \
                                    --n {params.n} \
                                    --prop_nonzero {params.prop_nz} \
                                    --standardize {params.stand} \
                                    --update_V_method {params.upd_v_met} \
                                    --w0_threshold {params.w0_thresh} \
                                    --X_colmeans {input.x_colm} \
                                    --Y_colmeans {input.y_colm} \
                                    --ncores 4 \
                                    --seed {params.dataset} \
                                    --output {output}
        echo "Finished analyzing chromosome {params.chr}."
        module unload R/4.2.3
        module unload gsl/2.7
        """

rule fit_ldpred2auto_gwide_by_trait:
  input:
    sumstats=lambda wildcards: get_sumstats_input_files(wildcards),
    ldmat=expand(config["LDMAT"]+config["LDMAT_PREFIX"]+"{CHR}_LD.ld.bin", CHR=config["CHR"])
  output:
    config["DEST"]+"output/ldpred2_auto_gwide_fit/"+config["DATA_PREFIX"]+"_chrAll_ldpred2_auto_fit_trait{TRAIT}_{DATASET}.rds"
  params:
    ss_prefix="output/summary_statistics/chr",
    ss_suffix="_sumstats_{DATASET}.rds",
    ld_prefix="data/LD_matrices/"+config["LDMAT_PREFIX"],
    ld_suffix="_LD.ld.bin",
    n=config["N_TRAIN"],
    h2_init=config["H2_INIT"],
    burn_in=config["LDPRED_BURN_IN"],
    num_iter=config["LDPRED_ITER"],
    trait="{TRAIT}",
    dataset="{DATASET}",
    chr_r=config["CHR_RANGE"],
    tmp=config["TMP"],
    script=config["SCRIPT"]+"fit_ldpred2auto_gwide_ukb.R"
  resources: cpus=4, mem_mb=150000, time_min=1440
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        echo "Start analyzing trait {params.trait} ..."
        Rscript {params.script} \
                                    --sumstats_prefix {params.ss_prefix} \
                                    --sumstats_suffix {params.ss_suffix} \
                                    --LD_matrix_prefix {params.ld_prefix} \
                                    --LD_matrix_suffix {params.ld_suffix} \
                                    --chr {params.chr_r} \
                                    --n {params.n} \
                                    --h2_init {params.h2_init} \
                                    --burn_in {params.burn_in} \
                                    --num_iter {params.num_iter} \
                                    --trait {params.trait} \
                                    --ncores 4 \
                                    --seed {params.dataset} \
                                    --temp_dir {params.tmp} \
                                    --output {output}
        module unload R/4.2.3
        module unload gsl/2.7
        """

def get_pred_input_files(wildcards):
    DATASET = wildcards.DATASET
    return expand(config["DEST"]+"output/ldpred2_auto_fit/"+config["DATA_PREFIX"]+"_chr{CHR}_ldpred2_auto_fit_trait{TRAIT}_{DATASET}.rds", DATASET=DATASET, CHR=config["CHR"],TRAIT=config["TRAIT"])

rule compute_prediction_accuracy_ldpred:
  input:
    pheno=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_pheno_{DATASET}.rds",
    geno=config["BED"],
    tids=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt",
    list=lambda wildcards: get_pred_input_files(wildcards)
  output:
    eff=config["DEST"]+"output/estimated_effects/"+config["DATA_PREFIX"]+"_ldpred2_auto_effects_{DATASET}.rds",
    pred_acc=config["DEST"]+"output/prediction_accuracy/"+config["DATA_PREFIX"]+"_ldpred2_auto_pred_acc_{DATASET}.rds"
  params:
    model=config["LDPRED_MODEL"],
    model_fit_dir=config["LDPRED_MODEL_FIT_DIR"],
    dataset="{DATASET}",
    chr_r=config["CHR_RANGE"],
    trait_r=config["TRAIT_RANGE"],
    imp_flag=config["IMPUTE_FLAG"],
    script=config["SCRIPT"]+"compute_prediction_accuracy_ukb.R",
    tmp=config["TMP"],
    prefix=config["DATA_PREFIX"]
  resources: cpus=4, mem_mb=100000, time_min=120
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        Rscript {params.script} \
                                    --model {params.model} \
                                    --model_fit_dir {params.model_fit_dir} \
                                    --pheno {input.pheno} \
                                    --geno {input.geno} \
                                    --test_ids {input.tids} \
                                    --data_id {params.dataset} \
                                    --chr {params.chr_r} \
                                    --trait {params.trait_r} \
                                    --impute_missing {params.imp_flag} \
                                    --output_eff {output.eff} \
                                    --output_pred_acc {output.pred_acc} \
                                    --temp_dir {params.tmp} \
                                    --ncores 4 \
                                    --prefix {params.prefix}
        module unload R/4.2.3
        """
       
def get_mrmashrss_input_files(wildcards):
    DATASET = wildcards.DATASET
    return expand(config["DEST"]+"output/mr_mash_rss_fit/"+config["DATA_PREFIX"]+"_chr{CHR}_mr_mash_rss_fit_{DATASET}.rds", DATASET=DATASET, CHR=config["CHR"])

rule compute_prediction_accuracy_mrmashrss:
  input:
    pheno=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_pheno_{DATASET}.rds",
    geno=config["BED"],
    pheno_means=config["DEST"]+"output/misc/"+config["DATA_PREFIX"]+"_pheno_means_{DATASET}.rds",
    tids=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt",
    list=lambda wildcards: get_mrmashrss_input_files(wildcards)
  output:
    eff=config["DEST"]+"output/estimated_effects/"+config["DATA_PREFIX"]+"_mr_mash_rss_effects_{DATASET}.rds",
    pred_acc=config["DEST"]+"output/prediction_accuracy/"+config["DATA_PREFIX"]+"_mr_mash_rss_pred_acc_{DATASET}.rds"
  params:
    model=config["MRMASH_MODEL"],
    model_fit_dir=config["MRMASH_MODEL_FIT_DIR"],
    dataset="{DATASET}",
    chr_r=config["CHR_RANGE"],
    trait_r=config["TRAIT_RANGE"],
    imp_flag=config["IMPUTE_FLAG"],
    prefix=config["DATA_PREFIX"],
    tmp=config["TMP"],
    script=config["SCRIPT"]+"compute_prediction_accuracy_ukb.R"
  resources: cpus=4, mem_mb=100000, time_min=120
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        Rscript {params.script} \
                                    --model {params.model} \
                                    --model_fit_dir {params.model_fit_dir} \
                                    --pheno {input.pheno} \
                                    --geno {input.geno} \
                                    --pheno_means {input.pheno_means} \
                                    --test_ids {input.tids} \
                                    --data_id {params.dataset} \
                                    --chr {params.chr_r} \
                                    --trait {params.trait_r} \
                                    --impute_missing {params.imp_flag} \
                                    --output_eff {output.eff} \
                                    --output_pred_acc {output.pred_acc} \
                                    --prefix {params.prefix} \
                                    --temp_dir {params.tmp} \
                                    --ncores 4
        module unload R/4.2.3
        """
        
def get_ldpred_gwide_input_files(wildcards):
    DATASET = wildcards.DATASET
    return expand(config["DEST"]+"output/ldpred2_auto_gwide_fit/"+config["DATA_PREFIX"]+"_chrAll_ldpred2_auto_fit_trait{TRAIT}_{DATASET}.rds", DATASET=DATASET, TRAIT=config["TRAIT"])

rule compute_prediction_accuracy_ldpred_gwide:
  input:
    pheno=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_pheno_{DATASET}.rds",
    geno=config["BED"],
    tids=config["DEST"]+"data/phenotypes/simulated/"+config["DATA_PREFIX"]+"_test_ids_{DATASET}.txt",
    list=lambda wildcards: get_ldpred_gwide_input_files(wildcards)
  output:
    eff=config["DEST"]+"output/estimated_effects/"+config["DATA_PREFIX"]+"_ldpred2_auto_gwide_effects_{DATASET}.rds",
    pred_acc=config["DEST"]+"output/prediction_accuracy/"+config["DATA_PREFIX"]+"_ldpred2_auto_gwide_pred_acc_{DATASET}.rds"
  params:
    model=config["LDPRED_MODEL"],
    model_fit_dir=config["LDPRED_MODEL_FIT_GWIDE_DIR"],
    dataset="{DATASET}",
    chr="0",
    trait_r=config["TRAIT_RANGE"],
    imp_flag=config["IMPUTE_FLAG"],
    script=config["SCRIPT"]+"compute_prediction_accuracy_ukb.R",
    tmp=config["TMP"],
    prefix=config["DATA_PREFIX"]
  resources: cpus=4, mem_mb=100000, time_min=120
  shell:
        """
        source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
        module load R/4.2.3
        export MKL_NUM_THREADS=1
        export OMP_NUM_THREADS=1
        Rscript {params.script} \
                                    --model {params.model} \
                                    --model_fit_dir {params.model_fit_dir} \
                                    --pheno {input.pheno} \
                                    --geno {input.geno} \
                                    --test_ids {input.tids} \
                                    --data_id {params.dataset} \
                                    --chr {params.chr} \
                                    --trait {params.trait_r} \
                                    --impute_missing {params.imp_flag} \
                                    --output_eff {output.eff} \
                                    --output_pred_acc {output.pred_acc} \
                                    --temp_dir {params.tmp} \
                                    --ncores 4 \
                                    --prefix {params.prefix}
        module unload R/4.2.3
        """
