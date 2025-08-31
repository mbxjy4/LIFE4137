**LIFE4137**
Materials and Methods used in Completion of LIFE4137 Individual Research Project

This Repo contains the full materials and methods including all code used in completion of the individual research project. 

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Running CMSclassifier](#Running-CMSclassifier)
- [Gene Expression Analysis](#gene-expression-analysis)
- [Isoform Proportion Analysis Preperation in Ada](#transcript-preperation-analysis-in-ada)
- [Isoform Analysis](#isoform-analysis)
- [Survival Analysis](#survival-analysis)

- <!-- TOC end --> 

<!-- TOC --><a name="Running-CMSclassifier"></a>
## Running CMSclassifier

Gene expression data used in the modelling for CMSclassifier was downloaded from [synapse](https://www.synapse.org/Synapse:syn4983432). RF and SS CMSclassifier methods were both run on the complete dataset in R using ____________CMS1.R. Only TCGA samples where RF or SS matched the [described CMS class](https://www.synapse.org/Synapse:syn4978510) assigned in the original paper were reatined and saved as ____________ TCGA_CMSclass.tsv. 

<!-- TOC --><a name="gene-expression-analysis"></a>
## Gene Expression Analysis

Gene expression data was downloaded from [Xenabroser](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). Data was filtered prior to loading into R to only include TCGA primary tumour samples (ending in -01). 67 genes were selected for analysis based on published literature which suggested difference in gene expression in CRC versus normal tissue or between CMS subtypes. Data was filtered to only include expression for these 67 genes and further filtered to only include samples for which I had a CMS classification for (using TCGA_CMSclass.tsv). Outliers were removed and one-way ANOVA performed with additional Tukey multiple comparisons of the means to test for signifcant difference in gene expression between each CMS subtype relationship. 31 genes showed significant difference between at least one subtype and all others and were taken forward for downstream transcriptomic and suvival analysis. 




