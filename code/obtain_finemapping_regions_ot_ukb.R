###Load libraries needed
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats_prefix"), type="character")
parser <- add_option(parser, c("--sumstats_suffix"), type="character")
parser <- add_option(parser, c("--traits"), type="character")
parser <- add_option(parser, c("--region_size"), type="integer")
parser <- add_option(parser, c("--mhc"), type="logical")
parser <- add_option(parser, c("--chr"), type="integer")
parser <- add_option(parser, c("--sig_threshold"), type="numeric")
parser <- add_option(parser, c("--gwas_method"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output_prefix_snp"), type="character")
parser <- add_option(parser, c("--output_prefix_boundaries"), type="character")
outparse <- parse_args(parser)

sumstats_prefix <- outparse$sumstats_prefix
sumstats_suffix <- outparse$sumstats_suffix
traits <- eval(parse(text=outparse$trait))
region_size <- outparse$region_size
mhc <- outparse$mhc
chr <- outparse$chr
sig_threshold <- outparse$sig_threshold
gwas_method <- outparse$gwas_method
output_prefix_snp <- outparse$output_prefix_snp
output_prefix_boundaries <- outparse$output_prefix_boundaries
ncores <- outparse$ncores

###Set seed
set.seed(1)

setDTthreads(ncores)

regions_snp_all_traits <- vector("list", length(traits))
regions_boundaries_all_traits <- vector("list", length(traits))
names(regions_snp_all_traits) <- names(regions_boundaries_all_traits) <- traits

for(k in 1:length(traits)){
  
  trait <- traits[k]

  dat <- fread(paste0(sumstats_prefix, trait, sumstats_suffix), data.table=FALSE, showProgress = FALSE)

  ###Add fine-mapping regions start stop
  dat$start <- dat$bp - region_size
  dat$stop <- dat$bp + region_size

  ###Extract chromosome of interest
  dat_chr <- dat[which(dat$chr == chr), ]
  
  rm(dat); gc()
  
  ###Extract significant SNPs for that chr
  if(gwas_method=="mtag"){
    dat_chr_sig <- dat_chr[which(dat_chr$mtag_p < sig_threshold), ]
  } else if(gwas_method=="plink"){
    dat_chr_sig <- dat_chr[which(dat_chr$plink_p < sig_threshold), ]
  }
  
  if(nrow(dat_chr_sig)==0){
    next
  }
  
  ###Assign SNPs to regions identified by an interval around the significant SNPs
  regions_snp <- vector("list", nrow(dat_chr_sig))
  regions_boundaries <- vector("list", nrow(dat_chr_sig))
  
  for(j in 1:nrow(dat_chr_sig)){
    regions_snp[[j]] <- dat_chr[which(dat_chr$bp > dat_chr_sig[j, "start"] & dat_chr$bp <= dat_chr_sig[j, "stop"]), "rsid"]
    regions_boundaries[[j]] <- c(dat_chr_sig[j, "chr"], dat_chr_sig[j, "start"], dat_chr_sig[j, "stop"])
  }
  
  ###Merge regions if they share one or more SNPs
  for(j in 2:length(regions_snp)){
    if(j==2){i <- 1}
    
    if(length(intersect(regions_snp[[i]], regions_snp[[j]])) > 0){
      
      regions_snp[[i]] <- unique(c(regions_snp[[i]], regions_snp[[j]]))
      regions_snp[[j]] <- "merged"
      regions_boundaries[[i]] <- c(regions_boundaries[[i]][1], regions_boundaries[[i]][2], regions_boundaries[[j]][3])
      regions_boundaries[[j]] <- "merged"
      
    } else {
      i <- j
    }
  }
  
  ###Filter out trash
  regions_snp_merged <- list()
  regions_boundaries_merged <- list()
  
  it <- 1
  for(i in 1:length(regions_snp)){
    if(all(regions_snp[[i]] != "merged")){
      regions_snp_merged[[it]] <- regions_snp[[i]]
      regions_boundaries_merged[[it]] <- regions_boundaries[[i]]
      it <- it + 1
    }
  }
  
  regions_snp_all_traits[[k]] <- regions_snp_merged
  regions_boundaries_all_traits[[k]] <- regions_boundaries_merged
}

###Prepare boundaries data by sorting regions across traits
regions_boundaries_all_traits_df <- as.data.frame(do.call("rbind", unlist(regions_boundaries_all_traits, recursive = FALSE)))
regions_boundaries_all_traits_df$pos <- 1:nrow(regions_boundaries_all_traits_df)
regions_boundaries_all_traits_df_ordered <- regions_boundaries_all_traits_df[order(regions_boundaries_all_traits_df[, 2]), ]

###Prepare SNP data by sorting regions across traits
regions_snp_all_traits_unlist <- unlist(regions_snp_all_traits, recursive = FALSE)
regions_snp_all_traits_unlist_ordered <- regions_snp_all_traits_unlist[regions_boundaries_all_traits_df_ordered$pos]


###Merge regions if they share one or more SNPs across traits
regions_boundaries_all_traits_df_ordered$pos <- NULL

for(j in 2:length(regions_snp_all_traits_unlist_ordered)){
  if(j==2){i <- 1}
  
  if(length(intersect(regions_snp_all_traits_unlist_ordered[[i]], regions_snp_all_traits_unlist_ordered[[j]])) > 0){
    
    regions_snp_all_traits_unlist_ordered[[i]] <- unique(c(regions_snp_all_traits_unlist_ordered[[i]], regions_snp_all_traits_unlist_ordered[[j]]))
    regions_snp_all_traits_unlist_ordered[[j]] <- "merged"
    regions_boundaries_all_traits_df_ordered[i, ] <- c(regions_boundaries_all_traits_df_ordered[i, 1], regions_boundaries_all_traits_df_ordered[i, 2], 
                                                       regions_boundaries_all_traits_df_ordered[j, 3])
    regions_boundaries_all_traits_df_ordered[j, ] <- c(0,0,0)
    
  } else {
    i <- j
  }
}

###Clean up output
merged_out_regions <- which(apply(regions_boundaries_all_traits_df_ordered, 1, function(x){all(x==0)}))

regions_boundaries_all_traits_df_ordered <- regions_boundaries_all_traits_df_ordered[-merged_out_regions, ]
colnames(regions_boundaries_all_traits_df_ordered) <- c("chr", "start", "stop")

regions_snp_all_traits_unlist_ordered <- regions_snp_all_traits_unlist_ordered[-merged_out_regions]

regions_boundaries_all_traits_df_ordered$n_snps <- sapply(regions_snp_all_traits_unlist_ordered, length)

###Write out output
for(i in 1:length(regions_snp_all_traits_unlist_ordered)){
  df <- data.frame(ID=regions_snp_all_traits_unlist_ordered[[i]])
  
  ##Skip MHC region if requested
  if(!mhc){
    if(regions_boundaries_all_traits_df_ordered[i,1] == 6 && 
       (regions_boundaries_all_traits_df_ordered[i,2] >= 25000000 & regions_boundaries_all_traits_df_ordered[i,3] <= 36000000)){
      next
    }
  }
  
  write.table(df, file=paste0(output_prefix_snp, "_chr", regions_boundaries_all_traits_df_ordered[i,1],
                              ".", regions_boundaries_all_traits_df_ordered[i,2], ".", regions_boundaries_all_traits_df_ordered[i,3], "_snp_names.txt"),
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

write.table(regions_boundaries_all_traits_df_ordered, file=paste0(output_prefix_boundaries, "_chr", chr, "_regions_info.txt"), 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


