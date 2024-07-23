###Load libraries needed
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--sumstats_mtag_prefix"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--traits"), type="character")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
sumstats_mtag_prefix <- outparse$sumstats_mtag_prefix
output <- outparse$output
traits <- eval(parse(text=outparse$traits))

###Read sumstats
univ_sumstats <- readRDS(sumstats)

###Compiled MTAG results
mtag_ss_Bhat <- matrix(as.numeric(NA), nrow=nrow(univ_sumstats$Bhat), ncol=length(traits))
colnames(mtag_ss_Bhat) <- colnames(univ_sumstats$Bhat)[traits]
rownames(mtag_ss_Bhat) <- rownames(univ_sumstats$Bhat)
mtag_ss_Shat <- matrix(as.numeric(NA), nrow=nrow(univ_sumstats$Shat), ncol=length(traits))
colnames(mtag_ss_Shat) <- colnames(univ_sumstats$Shat)[traits]
rownames(mtag_ss_Shat) <- rownames(univ_sumstats$Shat)

for(i in traits){
  ##Read MTAG results
  mtag_ss <- read.table(paste0(sumstats_mtag_prefix, i, ".txt"), header = TRUE, sep="\t")
  mtag_ss_bhat <- mtag_ss$mtag_beta
  names(mtag_ss_bhat) <- mtag_ss$SNP
  mtag_ss_shat <- mtag_ss$mtag_se
  names(mtag_ss_shat) <- mtag_ss$SNP
  
  ##Intersect the MTAG results with the univariate results
  in_common <- which(rownames(mtag_ss_Bhat) %in% names(mtag_ss_bhat))
  mtag_ss_Bhat[in_common, i] <- mtag_ss_bhat
  mtag_ss_Shat[in_common, i] <- mtag_ss_shat
  
  ##Assign univariate effect sizes and se to variants missing in the MTAG results (MTAG does not deal with indels)
  missing_mtag <- which(is.na(mtag_ss_Bhat[, i]))
  mtag_ss_Bhat[missing_mtag, i] <- univ_sumstats$Bhat[missing_mtag, i]
  mtag_ss_Shat[missing_mtag, i] <- univ_sumstats$Shat[missing_mtag, i]
}

###Save results
saveRDS(list(Bhat=mtag_ss_Bhat, Shat=mtag_ss_Shat), file=output)
