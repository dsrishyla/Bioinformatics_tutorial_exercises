#I followed this tutorial: https://www.youtube.com/watch?v=OzNzO8qwwp0

#######################
#load libraries
library(DESeq2)
library(tidyverse)

#read in counts data
counts_data <- read.csv('counts_data.csv')

#read in sample info - to know which columns correspond to which treatment
coldata <- read.csv('sample_info.csv')

#make sure rownames in coldata data matches colnames in counts data
all(colnames(counts_data) %in% rownames(coldata))

#are they in the same order
all(colnames(counts_data)==rownames(colData))

#construct DeSeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = coldata,
                       design = ~dexamethasone)

dds 


#pre-filtering: removing rows with low gene counts
#keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#set factor level - use untreated as reference level, to compare
relevel(dds$dexamethasone, ref="untreated")

#NOTE: Collapse technical replicates, NOT biological replicaes

#Step 3: Run DeSeq

dds <- DESeq(dds)
res <- results(dds)
res

#Explore results

summary(res)

#changing p threshold
res0.01 <- results(dds, alpha=0.01)
summary(res0.01)

#contrasts
resultsNames(dds)

#MA Plot - log2fold change vs mean of normalized counts; scatterplots; blue = sig. diff. expressed genes - the ones to the right, with high counts + change are of interest
plotMA(res)





