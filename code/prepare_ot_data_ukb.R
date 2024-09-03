###Load libraries needed
library(optparse)
library(data.table)
library(dplyr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--raw_pheno"), type="character")
parser <- add_option(parser, c("--bc_pheno"), type="character")
parser <- add_option(parser, c("--n_fold"), type="integer", default=-1)
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output_inds_list"), type="character")
parser <- add_option(parser, c("--output_pheno_data"), type="character")
outparse <- parse_args(parser)

ncores <- outparse$ncores
raw_pheno <- outparse$raw_pheno
bc_pheno <- outparse$bc_pheno
n_fold <- outparse$n_fold
output_pheno_data <- outparse$output_pheno_data
output_inds_list <- outparse$output_inds_list

###Set seed
set.seed(1)

###Set data.table threads
setDTthreads(ncores)

###Load data
##Raw phenotypes
dat <- fread(raw_pheno, data.table = FALSE, showProgress=FALSE, colClasses = "character")

##Cleaned blood cell traits
dat_bc <- readRDS(bc_pheno)
dat_covar <- dat_bc[, c("id", "sex", "assessment_centre", "age", "genotype_measurement_batch", paste0("pc_genetic", 1:10), "fold")]
rm(dat_bc); gc()

###Select columns to keep and filter the raw data
#weight, waist, hip, BMI, TFM, BFP
cols_to_keep <- c("eid", "6153-0.0", "6153-0.1", "6153-0.2", "6153-0.3", "6177-0.0", "6177-0.1", "6177-0.2", 
                  "4079-0.0", "4080-0.0", "48-0.0", "49-0.0", "23128-0.0")

dat_filt <- dat[, cols_to_keep]
colnames(dat_filt) <- c("id", "med_F1", "med_F2", "med_F3", "med_F4", "med_M1", "med_M2", "med_M3", 
                        "DP", "SP", "waist", "hip", "TFM")

rm(dat); gc()

###Convert numerical columns except the first one (the sample ids) to numeric values, and set all empty strings to NA.
n <- length(dat_filt)
for (i in 2:n) {
  x          <- dat_filt[,i]
  x[x == ""] <- as.character(NA)
  if(grepl('med', colnames(dat_filt)[i])){
    dat_filt[,i]    <- x
  }else{
    dat_filt[,i]    <- as.numeric(x)
  }
}

###Adjust blood pressure for blodd pressure medication (+10 to DP and +15 to SP)
# 6153-0.0	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones (Females)
dat_filt$mediblood1 <- as.numeric(NA)
dat_filt$mediblood1 <- replace(dat_filt$mediblood1, dat_filt$med_F1%in%c("-7", "1", "3", "4", "5"), 0)
dat_filt$mediblood1 <- replace(dat_filt$mediblood1, dat_filt$med_F1%in%c("-1", "-3"), NA)
dat_filt$mediblood1 <- replace(dat_filt$mediblood1, dat_filt$med_F1%in%c("2"), 1)

# 6153-0.1	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones (Females)
dat_filt$mediblood3 <- as.numeric(NA)
dat_filt$mediblood3 <- replace(dat_filt$mediblood3, dat_filt$med_F2%in%c("-7","1","3", "4", "5"), 0)
dat_filt$mediblood3 <- replace(dat_filt$mediblood3, dat_filt$med_F2%in%c("-1", "-3"), NA)
dat_filt$mediblood3 <- replace(dat_filt$mediblood3, dat_filt$med_F2%in%c("2"), 1)

# 6153-0.2	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones (Females)
dat_filt$mediblood5 <- as.numeric(NA)
dat_filt$mediblood5 <- replace(dat_filt$mediblood5, dat_filt$med_F3%in%c("-7","1","3", "4", "5"), 0)
dat_filt$mediblood5 <- replace(dat_filt$mediblood5, dat_filt$med_F3%in%c("-1", "-3"), NA)
dat_filt$mediblood5 <- replace(dat_filt$mediblood5, dat_filt$med_F3%in%c("2"), 1)

# 6153-0.3	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones (Females)
dat_filt$mediblood6 <- as.numeric(NA)
dat_filt$mediblood6 <- replace(dat_filt$mediblood6, dat_filt$med_F4%in%c("-7","1","3", "4", "5"), 0)
dat_filt$mediblood6 <- replace(dat_filt$mediblood6, dat_filt$med_F4%in%c("-1", "-3"), NA)
dat_filt$mediblood6 <- replace(dat_filt$mediblood6, dat_filt$med_F4%in%c("2"), 1)

# 6177-0.0	226928	Categorical (multiple) (Males)
dat_filt$mediblood2 <- as.numeric(NA)
dat_filt$mediblood2 <- replace(dat_filt$mediblood2, dat_filt$med_M1%in%c("-7","1","3"), 0)
dat_filt$mediblood2 <- replace(dat_filt$mediblood2, dat_filt$med_M1%in%c("-1", "-3"), NA)
dat_filt$mediblood2 <- replace(dat_filt$mediblood2, dat_filt$med_M1%in%c(2), 1)

# 6177-0.1	226928	Categorical (multiple) (Males)
dat_filt$mediblood4 <- as.numeric(NA)
dat_filt$mediblood4 <- replace(dat_filt$mediblood4, dat_filt$med_M2%in%c("-7","1","3"), 0)
dat_filt$mediblood4 <- replace(dat_filt$mediblood4, dat_filt$med_M2%in%c("-1", "-3"), NA)
dat_filt$mediblood4 <- replace(dat_filt$mediblood4, dat_filt$med_M2%in%c("2"), 1)

# 6177-0.2	226928	Categorical (multiple) (Males)
dat_filt$mediblood7 <- as.numeric(NA)
dat_filt$mediblood7 <- replace(dat_filt$mediblood7, dat_filt$med_M3%in%c("-7","1","3"), 0)
dat_filt$mediblood7 <- replace(dat_filt$mediblood7, dat_filt$med_M3%in%c("-1", "-3"), NA)
dat_filt$mediblood7 <- replace(dat_filt$mediblood7, dat_filt$med_M3%in%c("2"), 1)

dat_filt$add.med10 <- dat_filt$add.med15 <- as.numeric(NA)
dat_filt$add.med10 <- replace(dat_filt$add.med10, dat_filt$mediblood1==0 | dat_filt$mediblood2==0 |
                                dat_filt$mediblood3==0 | dat_filt$mediblood4==0 | dat_filt$mediblood5==0 |
                                dat_filt$mediblood6==0 | dat_filt$mediblood7==0, 0)
dat_filt$add.med15 <- replace(dat_filt$add.med15, dat_filt$mediblood1==0 | dat_filt$mediblood2==0 |
                                dat_filt$mediblood3==0 | dat_filt$mediblood4==0 | dat_filt$mediblood5==0 |
                                dat_filt$mediblood6==0 | dat_filt$mediblood7==0, 0)
dat_filt$add.med10 <- replace(dat_filt$add.med10, dat_filt$mediblood1==1 | dat_filt$mediblood2==1 |
                                dat_filt$mediblood3==1 | dat_filt$mediblood4==1 | dat_filt$mediblood5==1 |
                                dat_filt$mediblood6==1 | dat_filt$mediblood7==1, 10)
dat_filt$add.med15 <- replace(dat_filt$add.med15, dat_filt$mediblood1==1 | dat_filt$mediblood2==1 |
                                dat_filt$mediblood3==1 | dat_filt$mediblood4==1 | dat_filt$mediblood5==1 |
                                dat_filt$mediblood6==1 | dat_filt$mediblood7==1, 15)

dat_filt$SPa <- dat_filt$SP + dat_filt$add.med15
dat_filt$DPa <- dat_filt$DP + dat_filt$add.med10

###Create pulse pressure
dat_filt$PPa <- dat_filt$SPa - dat_filt$DPa

###Keep only variables of interest
dat_filt <- dat_filt[, c("id", "waist", "hip", "TFM", "DPa", "SPa", "PPa")]


###Join bc covariate data with BP data
dat_filt1 <- inner_join(dat_covar, dat_filt, by = join_by(id))
rm(dat_filt); gc()

# library(reshape2)
# library(ggplot2)
# 
# corr <- cor(dat_filt1[, c("waist", "hip", "TFM", "DPa", "SPa", "PPa")], use="pairwise.complete.obs")
# 
# melted_cormat <- melt(corr)
# 
# ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        midpoint = 0, limit = c(-1,1), space = "Lab",
#                        name="Pearson\nCorrelation") +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1,
#                                    size = 12, hjust = 1))+
#   coord_fixed()

###Select phenotypes for the analysis and keep only complete cases
dat_filt1 <- dat_filt1[complete.cases(dat_filt1), ]

###Remove individuals with "abnormal" measurements
pheno_names <- c("waist", "hip", "TFM", "DPa", "SPa", "PPa")
pheno_names_rint <- paste(pheno_names, "rint", sep="_")

##RINT
dat_filt2 <- matrix(as.numeric(NA), nrow=nrow(dat_filt1), ncol=length(pheno_names_rint))
colnames(dat_filt2) <- pheno_names_rint
dat_filt_final <- cbind(dat_filt1, dat_filt2)

for(name in pheno_names){
  id_rint <- which(colnames(dat_filt_final) == paste(name, "rint", sep="_"))
  id <- which(colnames(dat_filt_final) == name)
  dat_filt_final[, id_rint] = qnorm((rank(dat_filt_final[,id],na.last="keep")-0.5)/sum(!is.na(dat_filt_final[,id])))
}

##Compute empirical covariance matrix
covy <- dat_filt_final %>% select(all_of(pheno_names_rint)) %>% cov
D2 <- stats::mahalanobis(dat_filt_final %>% select(all_of(pheno_names_rint)), center=0, cov=covy) ## mahalanobis distance
dat_filt_final <- dat_filt_final[D2 < qchisq(0.01, df=length(pheno_names), lower.tail = FALSE),]

###Drop unneeded columns
dat_filt_final <- dat_filt_final[, -which(colnames(dat_filt_final) %in% pheno_names_rint)]

###Sample within fold, if requested
if(n_fold > 0){
  ##Write individuals id BEFORE sampling
  indiv_list <- cbind(dat_filt_final$id, dat_filt_final$id)
  
  base_out <- unlist(strsplit(output_inds_list, ".", fixed=TRUE))
  base_out_full <- paste(base_out[1], "full_data", sep="_")
  output_inds_list_full <- paste(base_out_full, base_out[2], sep=".")
  
  fwrite(x=indiv_list, file=output_inds_list_full, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
  
  ###Write full data BEFORE sampling
  base_out_pheno <- unlist(strsplit(output_pheno_data, ".", fixed=TRUE))
  base_out_pheno_full <- paste(base_out_pheno[1], "full_data", sep="_")
  output_pheno_data_full <- paste(base_out_pheno_full, base_out_pheno[2], sep=".")
  
  saveRDS(dat_filt_final, file=output_pheno_data_full)
  
  ##Sample
  dat_filt_final <- dat_filt_final %>% group_by(fold) %>% slice_sample(n=n_fold) %>% ungroup() %>% 
    arrange(as.numeric(id)) %>% as.data.frame()
}

cat("The final data set has", nrow(dat_filt_final), "individuals.\n")

cat("The training set corresponding to each test set has the following numbers of individuals:\n", 
    (nrow(dat_filt_final) - table(dat_filt_final$fold)))

###Prepare output
indiv_list <- cbind(dat_filt_final$id, dat_filt_final$id)

###Write out output
saveRDS(dat_filt_final, file=output_pheno_data)

fwrite(x=indiv_list, file=output_inds_list, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)

