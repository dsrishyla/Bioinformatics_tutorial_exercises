#script to run survival analysis using TCGA data
##I did not independently write this script, I followed this tutorial: https://www.youtube.com/watch?v=DnygUTAZFmM

#setwd("C:/Users/diksh/Downloads")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(survival)
install.packages("survminer")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")

library(DESeq2)
library(survminer)

#Getting clinical data for TCGA-BRCA cohort
clinical_brca <- GDCquery_clinic("TCGA-BRCA")

#select the data required for survival analysis
any(colnames(clinical_brca) %in% c("vital_status","days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_brca) %in% c("vital_status","days_to_last_follow_up", "days_to_death"))
clinical_brca[,c(46,54,55)]

#how many patients died vs alive
table(clinical_brca$vital_status)

#censor out patients who haven't died yet
clinical_brca$deceased <- ifelse(clinical_brca$vital_status == "Alive", FALSE, TRUE)
table(clinical_brca$deceased)

#re-formatting time to event
clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == 'Alive', clinical_brca$days_to_last_follow_up, clinical_brca$days_to_death)
clinical_brca$overall_survival

#obtain gene expression data
query_brca_all = GDCquery(
  project="TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open"
)

output_brca <- getResults(query_brca_all)

#get the ids of 20 tumor samples from the large dataset
tumor <- output_brca[output_brca$sample_type == "Primary Tumor", "cases"][1:20]
tumor

#prepare query obtain gene expression data from 20 primary tumors
query_brca = GDCquery(
  project="TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open",
  barcode=tumor)

#download from query
GDCdownload(query_brca)

#get counts
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment=TRUE)
brca_matrix <- assay(tcga_brca_data, "unstranded")
brca_matrix[1:10, 1:10]

#extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data)) #(gene) #contains gene names
coldata <- as.data.frame(colData(tcga_brca_data)) #(sample)

#vst transform counts to be used in survival analysis -----
# setting up countData object
dds <- DESeqDataSetFromMatrix(countData= brca_matrix,
                              colData = coldata,
                              design = ~1)

#Remove genes with sum total of 10 reads across all samples

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

#Variance stabilization
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10, 1:10]

#Changing the shape of the data
#Get data ONLY  from the TP53 gene and add the metadata info to it
brca_tp53 <- brca_matrix_vst %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key='case_id', value='counts', -gene_id) %>%
  left_join(.,gene_metadata,by='gene_id') %>%
  filter(gene_name == 'TP53')


#Classify into low vs high expression groups based on median values
median_value <- median(brca_tp53$counts) 

brca_tp53$strata <- ifelse(brca_tp53$counts >= median_value, "HIGH", "LOW")

#Replace case_id column with abbreviated gene id
brca_tp53$case_id <- gsub('-01.*','',brca_tp53$case_id)

#Merge id from brca_tp53 with clinical_brca
brca_tp53 <- merge(brca_tp53, clinical_brca, by.x='case_id', by.y='submitter_id')

#fit survival curve
fit <- survfit(Surv(overall_survival, deceased)~strata, data=brca_tp53)
fit
ggsurvplot(fit, data=brca_tp53, pval=T, risk.table=T) 

#to compare survival probabilities between two groups
fit2 <- survdiff(Surv(overall_survival, deceased)~strata, data=brca_tp53)
fit2
