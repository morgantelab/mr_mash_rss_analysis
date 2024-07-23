###Load libraries needed
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

input <- outparse$input
output <- outparse$output

###Write variants info files
dat <- read.table(input, sep="\t", header=FALSE)
colnames(dat) <- c("CHROM", "ID", "GEN_POS", "POS", "REF", "ALT")
dat1 <- dat[, c(1, 4, 5, 6, 2)]
  
write.table(dat1, output, sep="\t", col.names=FALSE,
            row.names=FALSE, quote=FALSE)
