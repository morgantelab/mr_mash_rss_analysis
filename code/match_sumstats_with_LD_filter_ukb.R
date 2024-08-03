###Load libraries needed
library(optparse)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--map"), type="character")
parser <- add_option(parser, c("--strand_flip"), type="logical")
parser <- add_option(parser, c("--remove_dups"), type="logical")
parser <- add_option(parser, c("--join_by_pos"), type="logical")
parser <- add_option(parser, c("--match_min_prop"), type="numeric", default=0.2)
parser <- add_option(parser, c("--maf_sumstats"), type="numeric")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--traits"), type="character")
parser <- add_option(parser, c("--output_LD_prefix"), type="character")
parser <- add_option(parser, c("--output_sumstats"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

sumstats_file <- outparse$sumstats
map_file <- outparse$map
LD_matrix <- outparse$LD_matrix
strand_flip <- outparse$strand_flip
remove_dups <- outparse$remove_dups
join_by_pos <- outparse$join_by_pos
match_min_prop <- outparse$match_min_prop
frq <- outparse$maf_sumstats
n <- outparse$n
traits <- eval(parse(text=outparse$traits))
output_LD_prefix <- outparse$output_LD_prefix
output_sumstats <- outparse$output_sumstats
seed <- outparse$seed

###Set seed
set.seed(seed)

###Read in the data
sumstats <- readRDS(sumstats_file)
mapp <- read.table(map_file, header=TRUE, sep="\t")
mapp <- mapp[, -3]
colnames(mapp) <- c("chr", "rsid", "pos", "a1", "a0", "freq")
LD <- readRDS(LD_matrix)

###Compute sd of allele freq in the LD ref panel
sd_ref <- sqrt(2 * mapp$freq * (1 - mapp$freq))

is_bad <- c()

###Loop over traits
for(trait in traits){
  ##Create sumstats data to match
  sumstats_df <- data.frame(rsid=rownames(sumstats$Bhat), chr=sumstats$alleles$chr, 
                            pos=sumstats$alleles$bp, a1=sumstats$alleles$a1, a0=sumstats$alleles$a2,
                            beta=sumstats$Bhat[, trait], beta_se=sumstats$Shat[, trait], n_eff=n)
  
  ##Match sumstats with LD panel
  df_beta <- snp_match(sumstats_df, mapp, strand_flip=strand_flip, remove_dups=remove_dups, 
                       return_flip_and_rev=TRUE, join_by_pos=join_by_pos, match.min.prop=match_min_prop)
  
  ###Get list of bad variants (https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html)
  sd_ss <- with(df_beta, 1 / sqrt(n_eff * beta_se^2 + beta^2))
  is_bad <- c(is_bad, which(sd_ss < (0.5 * sd_ref) | sd_ss > (sd_ref + 0.1) | sd_ss < 0.05 | sd_ref < 0.05))
  
  ##If no variant was allele-reversed and are not in the last trait, go to next trait
  if(all(df_beta$`_REV_`==FALSE) & trait != tail(traits, 1)){
    next
  ##If no variant was allele-reversed and are in the last trait, subset the data
  } else if(all(df_beta$`_REV_`==FALSE) & trait == tail(traits, 1)){
    sumstats$Bhat <- sumstats$Bhat[df_beta$`_NUM_ID_.ss`, ]
    sumstats$Shat <- sumstats$Shat[df_beta$`_NUM_ID_.ss`, ]
    sumstats$alleles <- sumstats$alleles[df_beta$`_NUM_ID_.ss`, ]
    mapp <- mapp[df_beta$`_NUM_ID_`, ]
  ##If some variants were allele-reversed and are in the last trait, error out
  } else if(any(df_beta$`_REV_`==TRUE)) {
    stop("Some variants were allele-reversed. Need to account for that")
  }
}

###Additional quality control filters
isBad <- which(sumstats$alleles$freq<frq | sumstats$alleles$freq>(1-frq) | abs(sumstats$alleles$freq - mapp$freq)>0.05)
is_Bad <- sort(unique(c(isBad, is_bad)))

###Filter sumstats and LD matrix
if(length(isBad) > 1 && all.equal(rownames(sumstats$alleles), mapp$rsid)){
  sumstats$Bhat <- sumstats$Bhat[-is_Bad, ]
  sumstats$Shat <- sumstats$Shat[-is_Bad, ]
  sumstats$alleles <- sumstats$alleles[-is_Bad, ]
  
  LD <- LD[-is_Bad, -is_Bad]
}

###Save LD to binary file    
conn <- file(paste0(output_LD_prefix, ".ld.bin"), "wb")
writeBin(as.numeric(LD), conn)
close(conn)

###Save LD to rds
saveRDS(LD, paste0(output_LD_prefix, ".rds"))

###Save sumstats to rds
saveRDS(sumstats, output_sumstats) 


