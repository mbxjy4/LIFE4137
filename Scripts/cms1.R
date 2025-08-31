#Install packages and dependencies
install.packages("devtools")
install.packages("shiny")
#Load devtools
library(devtools)
#Load dplyr for any future processing
library(dplyr)
#Install CMSclassifier
install_github("Sage-Bionetworks/CMSclassifier", force = TRUE)
#Load CMSclassifier
library(CMSclassifier)
##############################################################

#Running RF CMSclassifier on Gene Ex matrix from Synapse

#Load Gene Ex data table for all samples as df
df <- read.table("formatted_crc_data.tsv", sep = "\t", header = TRUE, check.names = FALSE)
#Run RF CMS on df (transform as col names = entrez gene IDs - needs to be other way round)
Rfcms <- CMSclassifier::classifyCMS(t(df), method="RF")[[3]]
#Confirm all 273 genes used in RF are in model genes - 273/273(100%) are.
length(intersect(mgenes, colnames(df)))

##############################################################

#Comparing Results of my CMS run with those reported on Synapse

#Load in results from Synapse
paper_res <- read.table("clinical_molecular_public_all.tsv", sep = "\t", header = TRUE, check.names = FALSE)
#Set sample names in paper_res to be the row names
rownames(paper_res) <- paper_res$sample
#Get the shared samples only
shared_samples <- intersect(rownames(Rfcms), rownames(paper_res))
#Create vectors of CMS classifications for shared samples for paper and self
Rfcms1 <- Rfcms[shared_samples, "RF.predictedCMS"]
paper_res1 <- paper_res[shared_samples, "cms_label"]
#Create matrix with results from self RF and paper for each sample
df1 <- cbind(shared_samples, Rfcms1, paper_res1)
#Turn into DF
df1 <- as.data.frame(df1)
#Filter to remove rows where RFcms1 = NA or paper_res = NOLBL (no label)
df2 <- df1 %>%
  filter(!is.na(Rfcms1) & paper_res1 != "NOLBL")

sum(df2$Rfcms1 == df2$paper_res1) #= 1571/1638 match (96%)

#########
#Conduct SS CMS
SScms <- CMSclassifier::classifyCMS(t(df),method="SSP")[[3]]
#Get shared sample SS classification
SScms1 <- SScms[shared_samples, "SSP.predictedCMS"]
#Create matrix with results from self SS and paper for each sample
df3 <- cbind(shared_samples, SScms1, paper_res1)
#Turn into df
df3 <- as.data.frame(df3)
#Filter to remove rows where SScms1 = NA or paper_res = NOLBL (no label)
df4 <- df3 %>%
  filter(!is.na(SScms1) & paper_res1 != "NOLBL")

sum(df4$SScms1 == df4$paper_res1) # =1525/1598 match = 95%

#Combining both models to one df

df5 <- full_join(df2, df4, "shared_samples")
#Populate NAs in paper_res1.x with the value from paper_res1.y
df5$paper_res1.x <- coalesce(df5$paper_res1.x, df5$paper_res1.y)
#Confirm no NAs in paper_resx.1
sum(is.na(df5$paper_res1.x)) #No NAs present
#Rename paper_res1.x to just paper_res
df5 <- df5 %>%
  rename(paper_res = paper_res1.x)
#Remove paper_res1.y col
df5 <- df5 %>% select(-paper_res1.y)
#Create df which keeps samples which match with CMSclassifier paper eiter through RF or SS
df6 <- df5%>%
  filter(Rfcms1 == paper_res | SScms1 == paper_res)
#Save this df
write.table(df6, file = "RForSS_matchedsamples.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
#Create df with just sample ID and CMSclassification cols only
df7 <- df6 %>% select(-Rfcms1, -SScms1)
#Save this df
write.table(df7, file = "Sample_CMSclass.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

#Filter df7 to only include TCGA samples
sum(grepl("^TCGA", df7$shared_samples)) # = 474 TCGA samples
df8 <- df7 %>%
  filter(grepl("^TCGA", shared_samples)) #df8 has 474 TCGA samples
#Save df8 (TCGA sample and CMS subtype)
write.table(df8, file = "TCGA_CMSclass.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

#Create .txt file with all TCGA samples for which I have a CMS classification for
TCGA_samples <- df8$shared_samples
writeLines(TCGA_samples, "TCGA_samples.txt")
