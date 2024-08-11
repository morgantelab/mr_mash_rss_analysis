###Load libraries needed
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--mfi"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--info_thresh"), type="numeric")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

mfi <- outparse$mfi
output <- outparse$output
ncores <- outparse$ncores
chr <- outparse$chr
info_thresh <- outparse$info_thresh

###Set seed
set.seed(1)

###Set data.table threads
setDTthreads(4)

###Read in the data
dat <- fread(mfi, data.table = FALSE, showProgress=FALSE)

###Filter for INFO score > 0.6
dat_filt <- dat[which(dat$V8>=info_thresh), ]

###Create unique ID###
unique_id <- paste(chr, with(dat_filt, paste(V3, V4, V5, sep = "_")), sep=":")

###Write output
write.table(matrix(unique_id, ncol=1), file=output, col.names = FALSE, row.names = FALSE, quote=FALSE, sep="\t")

