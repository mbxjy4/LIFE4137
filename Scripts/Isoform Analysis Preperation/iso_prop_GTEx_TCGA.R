Iso_prop <- read.table(gzfile("TcgaTargetGtex_rsem_isopct.gz"), sep = "\t", header = TRUE, check.names = FALSE)
#Read in list of GTEx and TCGA samples
samples <- readLines("GTEx_TCGA_CR_Samples.txt")
#Get list of columns which are in iso_prop and are in samples vector
matched_cols <- sapply(colnames(Iso_prop), function(col) any(startsWith(col, samples)))
#Filter to only keep those columns
Iso_prop_filtered <- Iso_prop[, matched_cols]
#Set row names of Iso_prop to be filtered df row names
rownames(Iso_prop_filtered) <- Iso_prop[, 1]
#Save df
write.table(Iso_prop_filtered, file = "GTEx_TCGA_samples.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
