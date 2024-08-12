###Load libraries needed
library(optparse)
library(data.table)
library(dplyr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--mtag_sumstats"), type="character")
parser <- add_option(parser, c("--plink_sumstats"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

mtag_sumstats <- outparse$mtag_sumstats
plink_sumstats <- outparse$plink_sumstats
output <- outparse$output
ncores <- outparse$ncores

###Set seed
set.seed(1)

setDTthreads(ncores)

##Read MTAG and plink results
dat_mtag <- fread(mtag_sumstats, data.table=FALSE, showProgress = FALSE)
dat_plink <- fread(plink_sumstats, data.table=FALSE, showProgress = FALSE)
  
##Create unique id for MTAG results
dat_mtag$unique_id <- with(dat_mtag, paste(CHR, paste(BP, A2, A1, sep = "_"), sep = ":"))
  
##Merge data using mtag results as backbone, since they do not contain ambiguous SNPs and duplicated rsid
dat <- left_join(dat_mtag, dat_plink, by = join_by(unique_id==ID))
  
dat <- dat[, c("SNP", "unique_id", "CHR", "BP", "REF", "ALT", "A1.y", "A1_FREQ", "N", "BETA", "SE", "T_STAT", "P",
               "mtag_beta", "mtag_se", "mtag_z", "mtag_pval")]
colnames(dat) <- c("rsid", "unique_id", "chr", "bp", "ref", "alt", "a1", "a1_freq", "n", "plink_beta", "plink_se", "plink_z", "plink_p",
                   "mtag_beta", "mtag_se", "mtag_z", "mtag_p")
  
##Write out the results
fwrite(dat, file=output, quote=FALSE, na="NA", row.names=FALSE, col.names=TRUE, showProgress=FALSE)

