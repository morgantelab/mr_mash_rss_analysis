###Load libraries needed
library(optparse)
library(bigsnpr)
library(mashr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--sumstats_prefix"), type="character")
parser <- add_option(parser, c("--sumstats_suffix"), type="character")
parser <- add_option(parser, c("--prior_matrices"), type="character")
parser <- add_option(parser, c("--residual_cov"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--rnd_size"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--LD_clump_r2_thresh"), type="numeric")
parser <- add_option(parser, c("--temp_dir"), type="character")
parser <- add_option(parser, c("--fold_file"), type="character")
parser <- add_option(parser, c("--sample_file"), type="character")
parser <- add_option(parser, c("--fold"), type="integer")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
sumstats_prefix <- outparse$sumstats_prefix
sumstats_suffix <- outparse$sumstats_suffix
prior_matrices <- outparse$prior_matrices
residual_cov <- outparse$residual_cov
chrs <- eval(parse(text=outparse$chr))
rnd_size <- outparse$rnd_size
ncores <- outparse$ncores
output <- outparse$output
seed <- outparse$seed
LD_clump_r2_thresh <- outparse$LD_clump_r2_thresh
temp_dir <- outparse$temp_dir
fold_file <- outparse$fold_file
sample_file <- outparse$sample_file
fold <- outparse$fold

###Set seed
set.seed(seed)

###Read in data
geno <- snp_attach(geno_dat)
V <- readRDS(residual_cov)
prior_mat <- readRDS(prior_matrices)
fold_info <- readRDS(fold_file)

##This is only needed to get the individuals in the training set for LD computation 
rownames(fold_info) <- fold_info$id
fold_info$id <- NULL
ids <- bigreadr::fread2(sample_file)
ids <- ids[-1, ]
ids <- ids[which(as.character(ids$ID_2) %in% rownames(fold_info)), "ID_2"]
fold_info <- fold_info[as.character(ids), ]
train_inds <- which(fold_info$fold != fold)
rm(ids, fold_info)

Bhat_all <- vector("list", length(chrs))
Shat_all <- vector("list", length(chrs))

for(i in chrs){
  
  ###Read in summary stats
  sumstats <- paste0(sumstats_prefix, i, sumstats_suffix)
  univ_sumstats <- readRDS(sumstats)
  
  if(i %in% c(1:9)){
    chr <- paste0("0", i)
  } else {
    chr <- as.character(i)
  }
  
  
  ###Filter geno data by chromosome
  to_keep <- which(geno$map$chromosome==chr)
  tmp <- tempfile(tmpdir=temp_dir)
  geno_filt_path <- snp_subset(geno, ind.row = train_inds, ind.col = to_keep, backingfile = tmp)
  geno_filt <- snp_attach(geno_filt_path)
  
  ###LD clumping on MAF
  snps_keep <- snp_clumping(geno_filt$genotypes, infos.chr = as.integer(geno_filt$map$chromosome),
                            infos.pos = geno_filt$map$physical.pos,
                            thr.r2 = LD_clump_r2_thresh, ncores = ncores)
  
  ###Extract summary stats for selected SNPs
  Bhat_sel <- univ_sumstats$Bhat[snps_keep, ]
  Shat_sel <- univ_sumstats$Shat[snps_keep, ]
  
  ###Assign selected summary stats to list
  Bhat_all[[i]] <- Bhat_sel
  Shat_all[[i]] <- Shat_sel
  
  file.remove(paste0(tmp, c(".bk",".rds")))
  print(paste("Finished chromosome ", i))
}

Bhat_all_sel <- do.call("rbind", Bhat_all)
Shat_all_sel <- do.call("rbind", Shat_all)

###Randomly sample desired number of SNPs
if(nrow(Bhat_all_sel) > rnd_size){
  sel_keep <- sort(x=sample(1:nrow(Bhat_all_sel), size=rnd_size))
  
  Bhat_all_sel_rnd <- Bhat_all_sel[sel_keep, ]
  Shat_all_sel_rnd <- Shat_all_sel[sel_keep, ]
} else {
  Bhat_all_sel_rnd <- Bhat_all_sel
  Shat_all_sel_rnd <- Shat_all_sel
}

###Compute grid and set data for mash
grid <- mr.mash.alpha::autoselect.mixsd(list(Bhat=Bhat_all_sel_rnd, Shat=Shat_all_sel_rnd), mult=sqrt(2)) ##do not square it here! mash does that internally.
dat_mash <- mash_set_data(Bhat_all_sel_rnd, Shat_all_sel_rnd, V=cov2cor(V), alpha = 0)

###Fit mash
fit <- mash(dat_mash, Ulist=prior_mat, grid=grid, usepointmass=TRUE)

###Remodel mash output needed
fitted_pi <- fit$fitted_g$pi
names(fitted_pi) <- gsub(".", "_grid", names(fitted_pi), fixed=TRUE)

###Prepare prior to save
##Obtain expanded prior matrices
S0 <- mr.mash.alpha::expand_covs(mats=prior_mat, grid=grid^2, zeromat=TRUE)
fitted_pi <- fitted_pi[names(S0)]

##Keep only mixture components with weight greater than 0
comps_to_keep <- which(fitted_pi > 0)
out <- list(S0=S0[comps_to_keep], w0=fitted_pi[comps_to_keep])

saveRDS(out, file=output)

