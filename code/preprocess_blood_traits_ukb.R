library(data.table)

options(datatable.fread.datatable=FALSE)

###Load data
dat <- fread("/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb673826.csv", showProgress=FALSE, header=TRUE)

###Select only initial measurement
dat <- dat[, c(1, seq(2, 94, by=3))]
colnames(dat) <- c("ID", "WBC_count", "RBC_count", "Haemoglobin_conc", "Haematocrit_perc", "MCV", "MCH", "MCH_conc",
                   "RBC_dw", "Platelet_count", "Platelet_crit", "MPV", "Platelet_dw", "Lymphocyte_count",
                   "Monocyte_count", "Neutrophill_count", "Eosinophill_count", "Basophill_count", "NRBC_count",
                   "Lymphocyte_perc", "Monocyte_perc", "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
                   "NRBC_perc", "Reticulocyte_perc", "Reticulocyte_count", "MRV", "MSCV", "IRF", "HLR_perc",
                   "HLR_count")

###Write out data
saveRDS(dat, "../data/raw/ukb_blood_traits.rds")
