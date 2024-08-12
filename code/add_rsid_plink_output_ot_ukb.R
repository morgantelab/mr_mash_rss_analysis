###Load libraries needed
library(optparse)
library(data.table)
library(dplyr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--mfi"), type="character")
parser <- add_option(parser, c("--plink_sumstats"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

mfi <- outparse$mfi
plink_sumstats <- outparse$plink_sumstats
output <- outparse$output
ncores <- outparse$ncores
chr <- outparse$chr

###Set seed
set.seed(1)

###Set data.table threads
setDTthreads(4)

###Read in the data
dat <- fread(mfi, data.table = FALSE, showProgress=FALSE)
dat_ss <- fread(plink_sumstats, data.table = FALSE, showProgress=FALSE)

###Create unique ID
dat$unique_id <- paste(chr, with(dat, paste(V3, V4, V5, sep = "_")), sep=":")

dat_joined <- inner_join(dat_ss, dat, by = join_by(ID == unique_id))
dat_joined$V1 <- dat_joined$V3 <- dat_joined$V4 <- dat_joined$V5 <- dat_joined$V6 <- dat_joined$V7 <- dat_joined$V8 <- NULL
colnames(dat_joined)[17] <- "rsID"

###Write output
write.table(dat_joined, file=output, col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
