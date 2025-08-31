#Install and load packages
#install.packages('tidyverse')
library(tidyverse)

#Load in Gene Ex. matrix (previously filtered to only include TCGA primary tumour samples (end in -01))
allgene_ex <- read.table("TCGA_GeneEx.tsv", sep = "\t", header = TRUE, check.names = FALSE)

#Create Vector of ENSG IDs for genes of interest (GOI)
ensg_ids <- c(
  "ENSG00000133703", "ENSG00000198848", "ENSG00000101384",
  "ENSG00000138755", "ENSG00000169245", "ENSG00000169248",
  "ENSG00000100453", "ENSG00000120217", "ENSG00000089692",
  "ENSG00000131203", "ENSG00000197646", "ENSG00000136997",
  "ENSG00000168036", "ENSG00000168646", "ENSG00000139292",
  "ENSG00000141736", "ENSG00000103460", "ENSG00000122359",
  "ENSG00000152256", "ENSG00000117394", "ENSG00000067225",
  "ENSG00000121879", "ENSG00000165556", "ENSG00000109670",
  "ENSG00000198286", "ENSG00000157764", "ENSG00000141510",
  "ENSG00000134982", "ENSG00000105329", "ENSG00000148516",
  "ENSG00000169554", "ENSG00000111371", "ENSG00000122691",
  "ENSG00000100985", "ENSG00000115414", "ENSG00000108821",
  "ENSG00000121594", "ENSG00000114013", "ENSG00000163606",
  "ENSG00000141646", "ENSG00000117713", "ENSG00000213281",
  "ENSG00000105173", "ENSG00000091664", "ENSG00000135914",
  "ENSG00000139926", "ENSG00000039068", "ENSG00000158481",
  "ENSG00000188676", "ENSG00000113520", "ENSG00000112116",
  "ENSG00000115008", "ENSG00000277632", "ENSG00000046774",
  "ENSG00000186081", "ENSG00000124469", "ENSG00000134258",
  "ENSG00000026508", "ENSG00000066468", "ENSG00000077782",
  "ENSG00000110092", "ENSG00000145113", "ENSG00000148848",
  "ENSG00000163359", "ENSG00000156510", "ENSG00000122786",
  "ENSG00000085733")

#Create Gene symbol vector for GOIs
ensg_sym <- c(
  "KRAS", "CES1", "JAG1",
  "CXCL9", "CXCL10", "CXCL11",
  "GZMB", "CD274", "LAG3",
  "IDO1", "PDL2", "MYC",
  "CTNNB1", "AXIN2", "LGR5",
  "ERBB6", "TOX3", "AUXA11",
  "PDK1", "GLUT1", "PKM2",
  "PIK3CA", "CDX2", "FBXW7",
  "CARD11", "BRAF", "TP53",
  "APC", "TGFB1", "ZEB1",
  "ZEB2", "SNAT1", "TWIST1",
  "MMP9", "FN1", "COL1A1",
  "CD80", "CD88",
  "CD200R1", "SMAD4", "ARID1A",
  "NRAS", "CCNE1",
  "SLC17A6", "HTR2B", "FRMD",
  "CDH1", "CD1C", "IDO2",
  "IL4", "IL17F", "IL1A",
  "CCL3", "MAGEC2", "KRT5",
  "CEACAM8", "VTCN1", "CD44",
  "FGFR2", "FGFR1", "CCND1",
  "MUC4", "ADAM12", "COL6A3",
  "HKDC1", "CALD1", "CTTN")

#Create df for ENSG IDs and Gene Symbols
ensg_id_sym <- as.data.frame(ensg_ids, ensg_sym)

#Filter gene ex. matrix to only include GOI
fil_gene_ex <- allgene_ex[grepl("^ENSG00000133703|^ENSG00000198848|^ENSG00000101384|
  |^ENSG00000138755|^ENSG00000169245|^ENSG00000169248|
  |^ENSG00000100453|^ENSG00000120217|^ENSG00000089692|
  |^ENSG00000131203|^ENSG00000197646|^ENSG00000136997|
  |^ENSG00000168036|^ENSG00000168646|^ENSG00000139292|
  |^ENSG00000141736|^ENSG00000103460|^ENSG00000122359|
  |^ENSG00000152256|^ENSG00000117394|^ENSG00000067225|
  |^ENSG00000121879|^ENSG00000165556|^ENSG00000109670|
  |^ENSG00000198286|^ENSG00000157764|^ENSG00000141510|
  |^ENSG00000134982|^ENSG00000105329|^ENSG00000148516|
  |^ENSG00000169554|^ENSG00000111371|^ENSG00000122691|
  |^ENSG00000100985|^ENSG00000115414|^ENSG00000108821|
  |^ENSG00000121594|^ENSG00000114013|^ENSG00000163606|
  |^ENSG00000141646|^ENSG00000117713|^ENSG00000213281|
  |^ENSG00000105173|^ENSG00000091664|^ENSG00000135914|
  |^ENSG00000139926|^ENSG00000039068|^ENSG00000158481|
  |^ENSG00000188676|^ENSG00000113520|^ENSG00000112116|
  |^ENSG00000115008|^ENSG00000277632|^ENSG00000046774|
  |^ENSG00000186081|^ENSG00000124469|^ENSG00000134258|
  |^ENSG00000026508|^ENSG00000066468|^ENSG00000077782|
  |^ENSG00000110092|^ENSG00000145113|^ENSG00000148848|
  |^ENSG00000163359|^ENSG00000156510|^ENSG00000122786|
  |^ENSG00000085733|^ENSG00000085733", rownames(allgene_ex)), ]


###############################

#Load in TCGA CMS subtype DF
TCGA_CMS <- read.table("TCGA_CMSclass.tsv", sep = "\t", header = TRUE, check.names = FALSE)
#Add -01 to the end of every sample to match format of Gene ex. matrix
TCGA_CMS$shared_samples <- paste0(TCGA_CMS$shared_samples, "-01")
#Create vector of only shared sample IDs
TCGA_shared_samples <- colnames(fil_gene_ex)
#Filter CMS classification so that it only includes shared samples
TCGA_CMS_fil <- TCGA_CMS %>% 
  filter(shared_samples %in% TCGA_shared_samples)
#Transposing gene ex df so that sample ID = rows and genes become col headers
Tran_Gene_Ex <- as.data.frame(t(fil_gene_ex))
#Creating TCGA shared samples col ready for merging gene ex and CMS class data
Tran_Gene_Ex$shared_samples <- rownames(Tran_Gene_Ex)
#Merging Gene Ex df and CMS subtype df
CMS_GeneEx <- merge(Tran_Gene_Ex, TCGA_CMS_fil, by = "shared_samples")
#Adding sample IDs back to be row names
rownames(CMS_GeneEx) <- rownames(Tran_Gene_Ex)
#Save as .tsv
write.table(CMS_GeneEx, file = "TCGA_gene_ex.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

#################################################################################

#Get only ENSG IDs (gene col names) except for shared_samples and paper_res
gene_cols <- setdiff(colnames(CMS_GeneEx), c("paper_res", "shared_samples"))

#Remove outliers
outliers <- function(x, multiplier = 5) {
  av <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  (x < (av - multiplier * sd)) | (x > (av + multiplier * sd))
}

#Copy df
CMS_GeneEx1 <- CMS_GeneEx  

#Loop over each gene column removing outliers
for (ensg in gene_cols) {
  outlier_mask <- outliers(CMS_GeneEx1[[ensg]])
  CMS_GeneEx1[[ensg]][outlier_mask] <- NA
}

#Get number of outliers for each gene
sapply(CMS_GeneEx1[gene_cols], function(x) sum(is.na(x)))

#Create output file
output_file <- "anova_tukey_results.txt"
#sink(output_file)

#Loop through each gene
for (gene in gene_cols) {
  cat("Analysing Gene:", gene, "\n")
  
  
  #Define formula as 'ENSG ID, ~paper_res' -> downstream will analyse gene in repsect to differences in CMS class (paper_res)
  formula <- as.formula(paste(gene, "~ paper_res"))
  
  #Run the 1-way ANOVA for that gene and get summary stats (overall p vlue)
  anova_result <- aov(formula, data = CMS_GeneEx1)
  cat("\n--- ANOVA Summary ---\n")
  print(summary(anova_result))
  
  #Run Tukey post hoc test for sig dif. between each CMS subtype (seperate p-values)
  cat("\n--- Tukey HSD ---\n")
  print(TukeyHSD(anova_result))
  
  #Visualising - replace ENSG number & Gene ID to create sep. graph for each protein
  ggplot(CMS_GeneEx1, aes(x = paper_res, y = .data[[gene]], fill = paper_res)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste("Expression of", gene, "by CMS Subtype"),
         y = "Normalised Gene Expression (log2(TPM + 0.001))",
         x = "CMS Subtype",
         fill = "CMS Subtype") +
    geom_jitter(width = 0.1, alpha = 0.4)
  
  #Save each plot (change name)
  ggsave(filename = paste0(gene, "_CMS_ex.png"))
}

#sink()

#Visualising - replace ENSG number & Gene ID to create sep. graph for each protein
ggplot(CMS_GeneEx1, aes(x = paper_res, y = ENSG00000108821.13, fill = paper_res)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  labs(title = paste("COL1A1 CMS Gene Expression"),
       y = "Normalised Gene Expression (log2(TPM + 0.001))",
       x = "CMS Subtype",
       fill = "CMS Subtype") +
  geom_jitter(width = 0.1, alpha = 0.4) + theme_classic()

#Save each plot (change name)
ggsave(filename = "COL1A1_Gene_Ex.png")
