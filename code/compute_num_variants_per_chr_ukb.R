###Load libraries needed
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
outparse <- parse_args(parser)

input <- outparse$input

###Read in data
options(datatable.fread.datatable=FALSE)

nsnps_per_chr <- data.frame(chr=1:22, n_variants=vector("numeric", 22))

for(i in 1:22){
  freq <- fread(paste0("../data/LD_matrices/number_of_variants/", input, "_chr", i, ".frq"), header=TRUE, showProgress=FALSE)
  nsnps_per_chr[i, 2] <- nrow(freq)
}

fwrite(nsnps_per_chr, file=paste0("../data/LD_matrices/number_of_variants/", input, "_num_variants_per_chr.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")