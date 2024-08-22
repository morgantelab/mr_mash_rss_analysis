###Load libraries needed
library(optparse)
library(data.table)
library(vctrs)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--mfi"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--info_thresh"), type="numeric")
parser <- add_option(parser, c("--only_biallelic"), type="logical")
parser <- add_option(parser, c("--snps_only"), type="logical")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

mfi <- outparse$mfi
output <- outparse$output
ncores <- outparse$ncores
chr <- outparse$chr
info_thresh <- outparse$info_thresh
only_biallelic <- outparse$only_biallelic
snps_only <- outparse$snps_only

###Set seed
set.seed(1)

###Set data.table threads
setDTthreads(ncores)

###Read in the data
dat <- fread(mfi, data.table = FALSE, showProgress=FALSE)

###Filter for INFO score > info_thresh
dat_filt <- dat[which(dat$V8>=info_thresh), ]

###Filter out multiallelic variants
if(only_biallelic){
  multi_all <- which(vec_duplicate_detect(dat_filt$V3))
  dat_filt <- dat_filt[-multi_all, ]
}

###Filter out indels
if(snps_only){
  indels_ref <- which(!(dat_filt$V4 %in% c("A", "G", "T", "C")))
  dat_filt <- dat_filt[-indels_ref, ]
  indels_alt <- which(!(dat_filt$V5 %in% c("A", "G", "T", "C")))
  dat_filt <- dat_filt[-indels_alt, ]
}

###Create unique ID###
unique_id <- paste(chr, with(dat_filt, paste(V3, V4, V5, sep = "_")), sep=":")

###Write output
write.table(matrix(unique_id, ncol=1), file=output, col.names = FALSE, row.names = FALSE, quote=FALSE, sep="\t")

