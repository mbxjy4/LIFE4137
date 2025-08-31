#Load packages
library(biomaRt)
library(tidyr)
library(dplyr)
library(ggplot2)

###########################################################################################
#Isoform analysis between GTEx normal samples and CMS subtypes
###########################################################################################

#Load ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#Get transcript table - ALL TRANSCRIPTS (inc. non-coding)
goi_df <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ccds',
                 'transcript_biotype', 'transcript_version'),
  filters = 'hgnc_symbol',
  values = 'ZEB1',
  mart = ensembl
)


#get list of ENST transcripts
transcripts <-goi_df$ensembl_transcript_id


#Load in iso prop table and cms subtype info for samples
Iso_prop <- read.table("GTEx_TCGA_samples.tsv", sep = "\t", header = TRUE, check.names = FALSE)
CMS_class <- read.table("TCGA_CMSclass.tsv", sep = "\t", header = TRUE, check.names = FALSE)

#transpose table to set samples as rows
Tiso_prop <- as.data.frame(t(Iso_prop))
#filter df to only include transcripts for GOI
Fil_tiso_prop <- Tiso_prop %>%
  dplyr::select(matches(paste0("^(", paste(transcripts, collapse = "|"), ")")))
#Get list of ENST columns
iso_cols <- colnames(Fil_tiso_prop)[startsWith(colnames(Fil_tiso_prop), "ENST")]
#Remove samples where sum of iso proportion = 0
Fil_tiso_prop_sums <- rowSums(Fil_tiso_prop[iso_cols], na.rm = TRUE)
Fil_tiso_prop <- Fil_tiso_prop[Fil_tiso_prop_sums != 0, ]

#Remove outliers from Fil_Tiso_prop
outliers <- function(x, multiplier = 8) {
  av <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  (x < (av - multiplier * sd)) | (x > (av + multiplier * sd))
}

Fil_tiso_prop1 <- Fil_tiso_prop  # Copy original df
# Loop over each isoform column
for (iso in iso_cols) {
  outlier_mask <- outliers(Fil_tiso_prop1[[iso]])
  Fil_tiso_prop1[[iso]][outlier_mask] <- NA
}

#Get number of outliers for each gene - used to decide multiplier
sapply(Fil_tiso_prop1[iso_cols], function(x) sum(is.na(x)))

####################################


#add -01 onto end of CMS class df samples to match fil_Tran_iso_GTEx ready for merging
CMS_class$shared_samples <- paste0(CMS_class$shared_samples, "-01")
#Create column for tumour type (normal or tumour)
Fil_tiso_prop1$Type <- NA
#assign 'normal' or 'tumour' to type depending on last 2 digits of sample name
Fil_tiso_prop1$Type <- ifelse(grepl("-01$", rownames(Fil_tiso_prop1)), "Tumour", "Normal")
#create new shared_samples column to be able to add cms class from CMS_class
Fil_tiso_prop1$shared_samples <- row.names(Fil_tiso_prop1)
#filter df to remove normal TCGA samples
Fil_tiso_prop2 <- Fil_tiso_prop1 %>%
  filter(!grepl("-11$", shared_samples)) %>%
  filter(!grepl("-02$", shared_samples))

#merge CMS class df with iso exp matrix
FinalDF <- full_join(Fil_tiso_prop2, CMS_class, by = "shared_samples")

#Set paper_res to be normal if it is a normal sample
FinalDF$paper_res <- ifelse(FinalDF$Type == "Normal", "Normal", FinalDF$paper_res)
#Set CMS class and type (normal/tumour) to be factors
FinalDF$paper_res <- as.factor(FinalDF$paper_res)
FinalDF$Type <- as.factor(FinalDF$Type)

#Filter DF to create seperate DFs for each CMS subtype and normal samples

CMS1 <- FinalDF %>%
  filter(paper_res == "CMS1" | paper_res == "Normal")
CMS1 <- droplevels(CMS1)

CMS2 <- FinalDF %>%
  filter(paper_res == "CMS2" | paper_res == "Normal")
CMS2 <- droplevels(CMS2)

CMS3 <- FinalDF %>%
  filter(paper_res == "CMS3" | paper_res == "Normal")
CMS3 <- droplevels(CMS3)

CMS4 <- FinalDF %>%
  filter(paper_res == "CMS4" | paper_res == "Normal")
CMS4 <- droplevels(CMS4)


sink("ZEB1_Transcripts.txt")

for (iso in iso_cols) {
  # Run t-test
  print("######################################################################")
  print("T Test between normal GTEx samples and CMS subtype for Iso:")
  print(iso)
  ttestCMS1 <- t.test(CMS1[[iso]] ~ CMS1$paper_res)
  ttestCMS2 <- t.test(CMS2[[iso]] ~ CMS2$paper_res)
  ttestCMS3 <- t.test(CMS3[[iso]] ~ CMS3$paper_res)
  ttestCMS4 <- t.test(CMS4[[iso]] ~ CMS4$paper_res)
  print("CMS1")
  print(ttestCMS1)
  print("CMS2")
  print(ttestCMS2)
  print("CMS3")
  print(ttestCMS3)
  print("CMS4")
  print(ttestCMS4)
}

sink()

###########################################################################################
#Pairwise analysis TCGA tumour vs normal samples
###########################################################################################
#Create new DF for paired samples (Fil_tiso_prop1 contains normal TCGA sampels as well)
PairedDF <- Fil_tiso_prop1
#Filter to only include isoforms for GOI
PairedDF <-  PairedDF%>%
  dplyr::select(matches(paste0("^(", paste(transcripts, collapse = "|"), ")")))
#create new 'type' column
PairedDF$Type <- NA
#assign 'normal' or 'tumour' to type depending on last 2 digits of sample name
PairedDF$Type <- ifelse(grepl("-01$", rownames(PairedDF)), "Tumour", "Normal")
#create new shared_samples column to be able to add cms class from CMS_class
PairedDF$shared_samples <- row.names(PairedDF)
#Reload CMS class so that the shared samples don't include -01 or -11
CMS_class <- read.table("TCGA_CMSclass.tsv", sep = "\t", header = TRUE, check.names = FALSE)
#Drop the last -01 or -11 from share_samples
PairedDF$shared_samples <- substr(PairedDF$shared_samples, 1, nchar(PairedDF$shared_samples) - 3)
#Merging iso prop df and CMS subtype df
PairedCMSDF <- merge(PairedDF, CMS_class, by = "shared_samples")

#Set CMS class and type (normal/tumour) to be factors
PairedCMSDF$paper_res <- as.factor(PairedCMSDF$paper_res)
PairedCMSDF$Type <- as.factor(PairedCMSDF$Type)

#Get samples which have a sample for tumour and normal tissue
paired_ids <- intersect(
  PairedCMSDF$shared_samples[PairedCMSDF$Type == "Tumour"],
  PairedCMSDF$shared_samples[PairedCMSDF$Type == "Normal"]
)

#create df with only paired samples
FinalPairedDF <- PairedCMSDF[PairedCMSDF$shared_samples %in% paired_ids, ]

#pairwise t test for samples with normal and tumour samples

sink("ZEB1_Transcripts.txt", append = TRUE)

for (iso in iso_cols) {
  #pivot wider the df based off of the sample and type
  print("##############################################################################")
  print("Pairwise T Test between tumour and normal TCGA samples for Iso:")
  print(iso)
  
  wide_df <- FinalPairedDF %>%
    dplyr::select(shared_samples, Type, all_of(iso)) %>%
    pivot_wider(names_from = Type, values_from = iso)
  
  #Remove NA rows
  wide_df <- wide_df %>% filter(!is.na(Tumour) & !is.na(Normal))
  
  #Run pairwise t test, printing ENST number and t test result
  result <- t.test(wide_df$Tumour, wide_df$Normal, paired = TRUE)
  print(iso)
  print(result)
}
sink()

#########################################################################################
#Differential isoform expression between CMS subtypes for all TCGA tumour samples
#########################################################################################

#create df with ENST iso props for only TCGA tumour samples
TumourDF <- FinalDF %>%
  filter(Type == "Tumour")
#Drop factor levels
TumourDF <- droplevels(TumourDF)

sink("ZEB1_Transcripts.txt", append = TRUE)
for (iso in iso_cols) {
  print("###########################################################################")
  print("Differential isoform expression between CMS subtypes for all TCGA tumour samples")
  print(iso)
  formula <- reformulate("paper_res", response = iso)
  aov_res <- aov(formula, data = TumourDF)
  print((summary(aov_res)))
  print((TukeyHSD(aov_res)))
}
sink()

#######################################################################################
#Classifying isoforms as protein/non-protein coding
#######################################################################################

#Create new DF
FinalDF1 <- FinalDF
#Filter DF to only include Tumour Samples
FinalDF1 <- FinalDF1 %>%
  filter(Type == "Tumour") 
#Create Columns for different types of trnascripts
FinalDF1$High <- NA
FinalDF1$Low <- NA
FinalDF1$Non <- NA

#Set row names of df goi_df to be ENST number
rownames(goi_df) <- goi_df$ensembl_transcript_id

#High conf. protein coding vector == ENSTs that have a CCDS and are protein coding
high_pc <- rownames(goi_df)[goi_df$transcript_biotype == "protein_coding" & startsWith(goi_df$ccds, "CCDS")]
#Non protein coding vector = ENSTs which aren't classed as protein coding
non_pc <- rownames(goi_df)[goi_df$transcript_biotype != "protein_coding"]
#Low conf. protein coding vector == ENSTs which are protein coding but don't have a CCDS
low_pc <- rownames(goi_df)[goi_df$transcript_biotype == "protein_coding" & goi_df$ccds == ""]

#Remove version number from ENST IDs - so that it ENSTs in vectors match ENSTs in Exp. DF
colnames(FinalDF1) <- sub("\\..*$", "", colnames(FinalDF1))


#Sum the proportions of all low conf. pc ENSTs (if only one then set it directly to that)
if (length(low_pc) == 1) {
  FinalDF1$Low <- FinalDF1[[low_pc]]
} else { 
  FinalDF1$Low <- rowSums(FinalDF1[, colnames(FinalDF1) %in% low_pc], na.rm = TRUE)
}

#Sum the proportions of all non pc ENSTs (if only one then set it directly to that)
if (length(non_pc) == 1) {
  FinalDF1$Non <- FinalDF1[[non_pc]]
} else { 
  FinalDF1$Non <- rowSums(FinalDF1[, colnames(FinalDF1) %in% non_pc], na.rm = TRUE)
}

#Sum the proportions of all high conf. pc ENSTs (if only one then set it directly to that)
if (length(high_pc) == 1) {
  FinalDF1$High <- FinalDF1[[high_pc]]
} else { 
  FinalDF1$High <- rowSums(FinalDF1[, colnames(FinalDF1) %in% high_pc], na.rm = TRUE)
}

#Check sums add up to ~100
FinalDF1$Sum <- rowSums(FinalDF1[, c("High", "Low", "Non")], na.rm = TRUE)
mean(FinalDF1$Sum, na.rm = TRUE)

#######################################################################################
#Run differential expression analysis for changes in type of transcript in CMS subtypes
#######################################################################################

#Create vector of column names for low, high, non
iso_type <- c("High", "Low", "Non")

sink("ZEB1_Iso_OUTrm.txt", append = TRUE)
for (type in iso_type) {
  print("###########################################################################")
  print("Differential isoform expression between CMS subtypes for all TCGA tumour samples")
  print(type)
  formula <- reformulate("paper_res", response = type)
  aov_res <- aov(formula, data = FinalDF1)
  print((summary(aov_res)))
  print((TukeyHSD(aov_res)))
}

sink()

#######################################################################################
#Pairwise analysis for isoform subtypes
#######################################################################################

#Take paired sample DF and Create Columns for different types of trnascripts
FinalPairedDF$High <- NA
FinalPairedDF$Low <- NA
FinalPairedDF$Non <- NA

#Remove version number from ENST IDs - so that it ENSTs in vectors match ENSTs in Exp. DF
colnames(FinalPairedDF) <- sub("\\..*$", "", colnames(FinalPairedDF))

#Sum the proportions of all low conf. pc ENSTs (if only one then set it directly to that)
if (length(low_pc) == 1) {
  FinalPairedDF$Low <- FinalPairedDF[[low_pc]]
} else { 
  FinalPairedDF$Low <- rowSums(FinalPairedDF[, colnames(FinalPairedDF) %in% low_pc], na.rm = TRUE)
}

#Sum the proportions of all non pc ENSTs (if only one then set it directly to that)
if (length(non_pc) == 1) {
  FinalPairedDF$Non <- FinalPairedDF[[non_pc]]
} else { 
  FinalPairedDF$Non <- rowSums(FinalPairedDF[, colnames(FinalPairedDF) %in% non_pc], na.rm = TRUE)
}

#Sum the proportions of all high conf. pc ENSTs (if only one then set it directly to that)
if (length(high_pc) == 1) {
  FinalPairedDF$High <- FinalPairedDF[[high_pc]]
} else { 
  FinalPairedDF$High <- rowSums(FinalPairedDF[, colnames(FinalPairedDF) %in% high_pc], na.rm = TRUE)
}

#Check sums add up to ~100
FinalPairedDF$Sum <- rowSums(FinalPairedDF[, c("High", "Low", "Non")], na.rm = TRUE)
mean(FinalPairedDF$Sum, na.rm = TRUE)

#Run pairwise comparison
sink("ZEB1_Iso_OUTrm.txt", append = TRUE)

for (type in iso_type) {
  #pivot wider the df based off of the sample and type
  print("##############################################################################")
  print("Pairwise T Test between tumour and normal TCGA samples for type of isoform expression")
  
  wide_df <- FinalPairedDF %>%
    dplyr::select(shared_samples, Type, all_of(type)) %>%
    pivot_wider(names_from = Type, values_from = type)
  
  #Remove NA rows
  wide_df <- wide_df %>% filter(!is.na(Tumour) & !is.na(Normal))
  
  #Run pairwise t test, printing ENST number and t test result
  result <- t.test(wide_df$Tumour, wide_df$Normal, paired = TRUE)
  print(type)
  print(result)
}
sink()

#########################################################################################
#CMS subtype vs GTEx normal comparison for differential expression of isoform subtypes. 
#########################################################################################

#Create Columns for different types of transcripts
FinalDF$High <- NA
FinalDF$Low <- NA
FinalDF$Non <- NA

#Remove version number from ENST IDs - so that it ENSTs in vectors match ENSTs in Exp. DF
colnames(FinalDF) <- sub("\\..*$", "", colnames(FinalDF))

#Sum the proportions of all low conf. pc ENSTs (if only one then set it directly to that)
if (length(low_pc) == 1) {
  FinalDF$Low <- FinalDF[[low_pc]]
} else { 
  FinalDF$Low <- rowSums(FinalDF[, colnames(FinalDF) %in% low_pc], na.rm = TRUE)
}

#Sum the proportions of all non pc ENSTs (if only one then set it directly to that)
if (length(non_pc) == 1) {
  FinalDF$Non <- FinalDF[[non_pc]]
} else { 
  FinalDF$Non <- rowSums(FinalDF[, colnames(FinalDF) %in% non_pc], na.rm = TRUE)
}

#Sum the proportions of all high conf. pc ENSTs (if only one then set it directly to that)
if (length(high_pc) == 1) {
  FinalDF$High <- FinalDF[[high_pc]]
} else { 
  FinalDF$High <- rowSums(FinalDF[, colnames(FinalDF) %in% high_pc], na.rm = TRUE)
}

FinalDF$Sum <- rowSums(FinalDF[, c("High", "Low", "Non")], na.rm = TRUE)
mean(FinalDF$Sum, na.rm = TRUE)

#Filter DFs to create seperate DFs for each CMS subtype and normal samples

CMS1 <- FinalDF %>%
  filter(paper_res == "CMS1" | paper_res == "Normal")
CMS1 <- droplevels(CMS1)

CMS2 <- FinalDF %>%
  filter(paper_res == "CMS2" | paper_res == "Normal")
CMS2 <- droplevels(CMS2)

CMS3 <- FinalDF %>%
  filter(paper_res == "CMS3" | paper_res == "Normal")
CMS3 <- droplevels(CMS3)

CMS4 <- FinalDF %>%
  filter(paper_res == "CMS4" | paper_res == "Normal")
CMS4 <- droplevels(CMS4)

sink("ZEB1_Iso_OUTrm.txt", append = TRUE)
for (type in iso_type) {
  # Run t-test
  print("######################################################################")
  print("T Test between normal GTEx samples and CMS subtype for Iso:")
  print(type)
  ttestCMS1_2 <- t.test(CMS1[[type]] ~ CMS1$paper_res)
  ttestCMS2_2 <- t.test(CMS2[[type]] ~ CMS2$paper_res)
  ttestCMS3_2 <- t.test(CMS3[[type]] ~ CMS3$paper_res)
  ttestCMS4_2 <- t.test(CMS4[[type]] ~ CMS4$paper_res)
  print("CMS1")
  print(type)
  print(ttestCMS1_2)
  print("CMS2")
  print(type)
  print(ttestCMS2_2)
  print("CMS3")
  print(type)
  print(ttestCMS3_2)
  print("CMS4")
  print(type)
  print(ttestCMS4_2)
}

sink()

#load package for creating significance bars
library(ggsignif)

#create graph for HCPC
ggplot(FinalDF, aes(x = paper_res, y = High, fill = paper_res)) +
    geom_boxplot(outlier.shape = NA, notch = TRUE) +
    geom_jitter(width = 0.1, alpha = 0.4) +
    labs(
      title = "HCPC ZEB1 Expression in Normal and CRC samples",
      y = "Proportion of MUC4 Expression (%)",
      x = "CMS Subtype",
      fill = "Subtype Group"
    ) +
    coord_cartesian(ylim = c(25, 100)) +
    theme_classic() +
    geom_signif(comparisons = list(c("CMS4", "CMS1"),
                                   c("CMS4", "CMS2"),
                                   c("CMS4", "CMS3")),
                map_signif_level = TRUE, 
                y_position = c(96, 92, 88))

ggsave("HCPC_ZEB1.png")

  geom_point()+
  geom_line(aes(group=shared_samples))

