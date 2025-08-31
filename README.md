**LIFE4137**
Materials and Methods used in Completion of LIFE4137 Individual Research Project

This Repo contains the full materials and methods including all code used in completion of the individual research project. Ensure required files for each script are in the CWD for smooth running of code. 

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Running CMSclassifier](#Running-CMSclassifier)
- [Gene Expression Analysis](#gene-expression-analysis)
- [Isoform Proportion Analysis Preperation in Ada](#isoform-proportion-analysis-preperation-in-ada)
- [Isoform Analysis](#isoform-analysis)
- [Survival Analysis](#survival-analysis)

- <!-- TOC end --> 

<!-- TOC --><a name="Running-CMSclassifier"></a>
## Running CMSclassifier

Gene expression data used in the modelling for CMSclassifier was downloaded from [synapse](https://www.synapse.org/Synapse:syn4983432). RF and SS CMSclassifier methods were both run on the complete dataset in R using ____________CMS1.R. Only TCGA samples where RF or SS matched the [described CMS class](https://www.synapse.org/Synapse:syn4978510) assigned in the original paper were reatined and saved as [TCGA_CMSclass.tsv](https://github.com/mbxjy4/LIFE4137/blob/main/Additional%20Files/TCGA_CMSclass.tsv). 

Script Used: [cms1.R](https://github.com/mbxjy4/LIFE4137/blob/main/Scripts/cms1.R)

<!-- TOC --><a name="gene-expression-analysis"></a>
## Gene Expression Analysis

Gene expression data was downloaded from [Xenabroser](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). Data was filtered prior to loading into R to only include TCGA primary tumour samples (ending in -01).

Script used: [Gene_Ex_outlierRM.R](https://github.com/mbxjy4/LIFE4137/blob/main/Scripts/Gene_Ex_outlierRM.R)

67 genes were selected for analysis based on published literature which suggested difference in gene expression in CRC versus normal tissue or between CMS subtypes. Data was filtered to only include expression for these 67 genes and further filtered to only include samples for which I had a CMS classification for (using TCGA_CMSclass.tsv). Outliers were removed and one-way ANOVA performed with additional Tukey multiple comparisons of the means to test for signifcant difference in gene expression between each CMS subtype relationship. 31 genes showed significant difference between at least one subtype and all others and were taken forward for downstream transcriptomic and suvival analysis. Full results available [here](https://github.com/mbxjy4/LIFE4137/blob/main/Gene%20Expression%20Results/anova_tukey_results.txt).

<!-- TOC --><a name="isoform-proportion-analysis-preperation-in-ada"></a>
## Isoform Proportion Analysis Preperation in Ada

Isoform proportion data was downloaded from [Xenabrowser](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_isopct&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). Due to the size of the file, the data was trimmed down using Ada to remove TCGA samples which I didn't have a CMS classification for. 

Scripts Used:
[iso_prop_GTEx_TCGA.R](https://github.com/mbxjy4/LIFE4137/blob/main/Scripts/Isoform%20Analysis%20Preperation/iso_prop_GTEx_TCGA.R)
[iso_prop_GTEx_TCGA.sh](https://github.com/mbxjy4/LIFE4137/blob/main/Scripts/Isoform%20Analysis%20Preperation/iso_prop_GTEx_TCGA.sh)

<!-- TOC --><a name="isoform-analysis"></a>
## Isoform Analysis


The above preparation created GTEx_TCGA_samples.tsv which was used for isoform expression analysis alongside TCGA_CMSclass.tsv. Isoform proportion data was filtered to only include transcripts for each GOI which were then classified as HCPC, LCPC or NPC. Three seperate analyses was then performed:

1) T-tests comparing each CMS subtype expression (TCGA) for HCPC, LCPC and NPC with normal (GTEx) expression.
2) Pairwsie T-tests comparing CRC expression (TCGA) in primary tumour samples with solid normal tissue samples from the same individual.
3) Non-pairwise T-tests comparing expression between each CMS subtype with Tukey multiple comparisons of the means.

Script Used: [all_stats_outliersRM.R](https://github.com/mbxjy4/LIFE4137/blob/main/Scripts/all_stats_outliersRM.R)

To run the code for each specific gene, replace the gene symbol in the 'values =' condition for getBM to create the goi_df. Update file names as appropriate. 

NOTE: The above analysis was also performed on specific transcripts if isoform groups showed particular interest (e.g. MUC4 and ZEB1). This analysis is at the start of the script as it was originally conducted before realising the full scale of transcript analysis for every gene. 

Full results available [here](https://github.com/mbxjy4/LIFE4137/tree/main/Isoform%20Expression%20Results) 

<!-- TOC --><a name="survival-analysis"></a>
## Survival Analysis

Survival data was downloaded from [Xenabrowser](https://xenabrowser.net/datapages/?dataset=TCGA_survival_data&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). Code was provided by Heshmat Borhani with minor edits performed. Due to word limitations, only DSS was evaluated.

Script used: [survival_CMS_specific.R](https://github.com/mbxjy4/LIFE4137/blob/main/Scripts/survival_CMS_specific.R) 

Update file names and directories as appropriate. To run the code for each gene, create a new column in DF with the desired gene's ENSG ID. E.g.

```bash
#Create column in DF for GOI and set it to be the respective ENSG ID
DF$CD274 <- DF$ENSG00000120217
```

Due to the number of results produced these are not available online however can be requested if required. 





