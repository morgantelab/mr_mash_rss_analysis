###Load libraries needed
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--prefix"), type="character")
parser <- add_option(parser, c("--suffix"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

prefix <- outparse$prefix
suffix <- outparse$suffix
output <- outparse$output
ncores <- outparse$ncores
chr <- eval(parse(text=outparse$chr))

###Set seed
set.seed(1)

###Set data.table threads
setDTthreads(8)

###Merge cromosome files
ss_all <- vector("list", length=length(chr))

for(i in chr){
  ss_all[[i]] <- fread(paste0(prefix, i, suffix), data.table = FALSE, showProgress=FALSE)
}

ss_all_df <- do.call("rbind", ss_all)

###Write output
write.table(ss_all_df, file=output, col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
