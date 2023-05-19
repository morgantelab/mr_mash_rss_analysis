###Load libraries needed
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno_fam"), type="character")
parser <- add_option(parser, c("--withdrawn_inds"), type="character")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno_fam
withdrawn_inds <- outparse$withdrawn_inds
output <- outparse$output

###Read in the data
options(datatable.fread.datatable=FALSE)

geno_fam <- fread(geno_dat, showProgress=FALSE, header=FALSE)[, 1:2]
withdr <- fread(withdrawn_inds, showProgress=FALSE, header=FALSE)

###Merge the data
total <- rbind(geno_fam, withdr)
total <- total[!duplicated(total), ]

fwrite(total, file=output, row.names=FALSE, col.names=FALSE, sep=" ", quote=FALSE, showProgress=FALSE)
