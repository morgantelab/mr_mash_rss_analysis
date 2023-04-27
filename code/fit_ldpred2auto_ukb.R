###Load libraries needed
library(optparse)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--sparse"), type="logical", default=FALSE)
parser <- add_option(parser, c("--h2_init"), type="numeric", default=0.5)
parser <- add_option(parser, c("--shrink_corr"), type="numeric", default=1)
parser <- add_option(parser, c("--allow_jump_sign"), type="logical", default=TRUE)
parser <- add_option(parser, c("--verbose"), type="logical", default=FALSE)
parser <- add_option(parser, c("--burn_in"), type="integer", default=500)
parser <- add_option(parser, c("--num_iter"), type="integer", default=500)
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--outprefix"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--trait"), type="integer")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_matrix <- outparse$LD_matrix
n <- outparse$n
sparse <- outparse$sparse
h2_init <- outparse$h2_init
shrink_corr <- outparse$shrink_corr
allow_jump_sign <- outparse$allow_jump_sign
verbose <- outparse$verbose
burn_in <- outparse$burn_in
num_iter <- outparse$num_iter
ncores <- outparse$ncores
outprefix <- outparse$outprefix
data_id <- outparse$data_id
trait <- outparse$trait

vec_p_init <- seq_log(1e-4, 0.2, length.out = 30)

###Set seed
set.seed(data_id)

###Read in data
univ_sumstats <- readRDS(paste0("../output/summary_statistics/", sumstats, "_", data_id, ".rds"))
p <- nrow(univ_sumstats$Bhat)
LD <- matrix(readBin(paste0("../data/LD_matrices/", LD_matrix, "_", data_id, ".ld.bin"), what="numeric", n=p^2), nrow=p, ncol=p, byrow=TRUE)

###Prepare the sumstats and LD matrix
df_beta <- data.frame(beta=univ_sumstats$Bhat[, trait],
                      beta_se=univ_sumstats$Shat[, trait],
                      n_eff=rep(n, times=p))
corr <- as_SFBM(as(LD, "CsparseMatrix"))
rm(LD)

###Fit mr.mash.rss
fit_ldpred2_auto <- snp_ldpred2_auto(corr=corr, df_beta=df_beta, h2_init=h2_init, vec_p_init=vec_p_init,
                                     burn_in=burn_in, num_iter=num_iter, sparse=sparse, verbose=verbose,
                                     allow_jump_sign=allow_jump_sign, shrink_corr=shrink_corr, ncores=ncores)

saveRDS(fit_ldpred2_auto, file=paste0("../output/ldpred2_auto_fit/", outprefix, "_ldpred2_auto_fit_trait", trait, "_", data_id, ".rds"))
