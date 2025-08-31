library(dplyr)
library(ggplot2)
library(survival)
#install.packages("survminer")
library(survminer)

#Load in Gene Ex df for TCGA samples w/CMS class - created in multi gene project
#Already has shared_samples and paper_res columns
GeneExDF <- read.table("TCGA_gene_ex.tsv", sep = "\t", header = TRUE)
#Load in survival data df
SurvivalDF <- read.table("TCGA_survival_data.tsv", sep = "\t", header = TRUE)

#Get vector of Gene IDs
ENSGs <- colnames(GeneExDF)[startsWith(colnames(GeneExDF), "ENSG")]
#Convert values of Gene Ex to be 2^TPM + 1
GeneExDF[, colnames(GeneExDF) %in% ENSGs] <- 
  (2^GeneExDF[, colnames(GeneExDF) %in% ENSGs]) + 1

#Create 'sample' column from shared samples in GeneExDF to match col name of SurvivalDF
GeneExDF$sample <- GeneExDF$shared_samples
#Get character vector of samples which I have survival and CMS data for
shared_samples <- intersect(GeneExDF$sample, SurvivalDF$sample)
#Filter survival Df to only include those samples
SurvivalDF <- SurvivalDF %>%
  filter(SurvivalDF$sample %in% shared_samples)
#Merge dfs keeping only samples which have both CMS and survival data
DF <- inner_join(SurvivalDF, GeneExDF, by = "sample")
#Filter to only include CMS'x' samples
DF <- DF %>% 
  filter(DF$paper_res == "CMS4")

#Create column in DF for GOI and set it to be the respective ENSG ID
DF$CD274 <- DF$ENSG00000120217

#Divide times by 30 - Convert days to months
DF$OS.time <- DF$OS.time / 30
DF$DSS.time <- DF$DSS.time / 30
DF$DFI.time <- DF$DFI.time / 30
DF$PFI.time <- DF$PFI.time / 30

dir.create("out/CD274/CMS1/Final", recursive = TRUE)

#########################################################################################

DF$CD274_median <- NA
DF$CD274_median[DF$CD274 <= quantile(DF$CD274, 0.5, na.rm=TRUE)] <- "LOW"
DF$CD274_median[DF$CD274 > quantile(DF$CD274, 0.5, na.rm=TRUE)] <- "HIGH"

DF$CD274_25p <- NA
DF$CD274_25p[DF$CD274 <= quantile(DF$CD274, 0.25, na.rm=TRUE)] <- "LOW"
DF$CD274_25p[DF$CD274 > quantile(DF$CD274, 0.25, na.rm=TRUE)] <- "HIGH"

DF$CD274_75p <- NA
DF$CD274_75p[DF$CD274 <= quantile(DF$CD274, 0.75, na.rm=TRUE)] <- "LOW"
DF$CD274_75p[DF$CD274 > quantile(DF$CD274, 0.75, na.rm=TRUE)] <- "HIGH"

DF$CD274_66p <- NA
DF$CD274_66p[DF$CD274 <= quantile(DF$CD274, 0.66, na.rm=TRUE)] <- "LOW"
DF$CD274_66p[DF$CD274 > quantile(DF$CD274, 0.66, na.rm=TRUE)] <- "HIGH"

DF$CD274_33p <- NA
DF$CD274_33p[DF$CD274 <= quantile(DF$CD274, 0.33, na.rm=TRUE)] <- "LOW"
DF$CD274_33p[DF$CD274 > quantile(DF$CD274, 0.33, na.rm=TRUE)] <- "HIGH"

DF$CD274_33_66 <- NA
DF$CD274_33_66[DF$CD274 <= quantile(DF$CD274, 0.33, na.rm=TRUE)] <- "LOW"
DF$CD274_33_66[DF$CD274 >= quantile(DF$CD274, 0.66, na.rm=TRUE)] <- "HIGH"

datasetName <- "CD274"
# SURVIVAL CURVES - CD274
for(strata in c("CD274_median", "CD274_25p", "CD274_75p", "CD274_66p", "CD274_33p", "CD274_33_66")){
  
  H_gTxt <- paste("HIGH CD274", sep="")
  L_gTxt <- paste("LOW CD274", sep="")
  suffix <- paste(datasetName, strata, sep="_")
  
  sc = "OS"
  tmp.sub <- subset(DF, !is.na(OS) & !is.na(OS.time))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(OS.time, OS) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF","black"), title=paste(strata,sep=" - "), xlab="Months", ylab="Overall Survival (OS)", 
                   break.time.by = 30,
                   xlim = c(0, 120),
                   font.x = c(18, "black"),
                   font.y = c(18, "black"),
                   font.tickslab = c(18, "plain", "black"),
                   legend = c(0.25, 0.15),
                   legend.title = "",
                   legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                  paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(18, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord=c(70, 80), pval.size = 5,
                   size = 1, censor.size = 2)
  
  tiff(filename = paste("out/CD274/CMS1/Final/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
  
  sc = "DSS"
  tmp.sub <- subset(DF, !is.na(DSS) & !is.na(DSS.time))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit <- survfit(Surv(DSS.time, DSS) ~ usedSubset, data = tmp.sub)
  p <- ggsurvplot(fit, palette = c("#ED0000FF", "black"),title=paste(strata,sep=" - "), xlab="Months", ylab="Disease Specific Survival (DSS)", 
                  break.time.by = 30, 
                  xlim = c(0, 120),
                  font.x = c(18, "black"),
                  font.y = c(18, "black"),
                  font.tickslab = c(18, "plain", "black"),
                  legend = c(0.25, 0.15),
                  legend.title = "",
                  legend.labs =c(paste0 (H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                 paste0 (L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                  font.legend = c(18, "black"),
                  fun = function(y) y*100,
                  pval = TRUE, pval.coord = c(70, 80), pval.size = 5,
                  size = 1, censor.size = 2)
  
  tiff(filename = paste ("out/CD274/CMS1/Final/KM_", sc, "_", suffix, "_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12,compression = c("none"),bg = "white")
  print(p)
  dev.off()
  
  
  sc = "PFI"
  tmp.sub <- subset(DF,!is.na(PFI) & !is.na(PFI.time))
  tmp.sub$usedSubset <- as.factor(tmp.sub[[strata]])
  tmp.sub <- subset(tmp.sub,!is.na(usedSubset))
  
  fit<- survfit(Surv(PFI.time, PFI) ~ usedSubset, data= tmp.sub)
  p <-  ggsurvplot(fit, palette = c("#ED0000FF","black"),title=paste(strata,sep=" - "), xlab="Months", ylab="Progression Free Interval (PFI)",
                   break.time.by = 30, 
                   xlim = c(0, 120),
                   font.x = c(18, "black"),
                   font.y = c(18, "black"),
                   font.tickslab = c(18, "plain", "black"),
                   legend = c(0.25, 0.15),
                   legend.title = "",
                   legend.labs = c(paste0(H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""),
                                   paste0(L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep="")),
                   font.legend = c(18, "black"),
                   fun = function(y) y*100,
                   pval = TRUE, pval.coord = c(70, 80), pval.size = 5,
                   size = 1, censor.size = 2)
  
  tiff(filename =  paste("out/CD274/CMS1/Final/KM_", sc, "_", suffix,"_v1.tiff", sep = ""),   
       width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
  print(p)
  dev.off()
}

######################################################################################
#For cms subtypes - only showing overall survival for each cms subtype 
#nothing to do with gene ex. 

#CMS1_gTxt <- paste("CMS1 GOI", sep="")
#CMS1_gTxt <- paste("CMS1 GOI", sep="")
#CMS1_gTxt <- paste("CMS1 GOI", sep="")
#CMS1_gTxt <- paste("CMS1 GOI", sep="")
#suffix <- paste(datasetName, strata, sep="_")

#fit<- survfit(Surv(PFI.time, PFI) ~ usedSubset, data= tmp.sub)
#p <-  ggsurvplot(fit, palette = c("#ED0000FF","black", "green"),title=paste(strata,sep=" - "), xlab="Months", ylab="Progression Free Interval (PFI)",
#                break.time.by = 30, 
#               xlim = c(0, 2000),
#              font.x = c(18, "black"),
#             font.y = c(18, "black"),
#            font.tickslab = c(18, "plain", "black"),
#           legend = c(0.25, 0.15),
#          legend.title = "",
#         legend.labs = c(paste0(H_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="HIGH"))[[1]],")",sep=""), 
#                         paste0(L_gTxt, " (n = ", dim(subset(tmp.sub,usedSubset=="LOW"))[[1]],")",sep=""),
#                        paste0(CMS1_gTxt, " (n = ", dim(subset(tmp.sub,paper_res=="CMS1"))[[1]],")",sep="")),
#      font.legend = c(18, "black"),
#     fun = function(y) y*100,
#    pval = TRUE, pval.coord = c(70, 80), pval.size = 5,
#   size = 1, censor.size = 2)
#  

#tiff(filename =  paste("out/Final/CMS1/KM_", sc, "_", suffix,"_v2.tiff", sep = ""),   
#    width = 530, height = 500, units = "px", pointsize = 12, compression = c("none"), bg = "white")
#print(p)
#dev.off()

#31
#39
#67
#96
#122
#149