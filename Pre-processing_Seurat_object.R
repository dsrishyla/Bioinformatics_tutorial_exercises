#Below, I replicated the code from another researcher's script while teaching myself how to preprocess a Seurat object
#The code is replicated from the file "Fig_1A-G_and_Fig_S1A-G.R" in this repository: https://github.com/alonmillet/apoe-ad-age-atlas 
#I also used this tutorial for guidance: https://www.youtube.com/watch?v=xbX49h7BiUU 
###########################################################################################

install.packages('Seurat')
install.packages('hdf5r')
install.packages('flexmix')
install.packages('R.utils')
install.packages('devtools')
devtools::install_github('satijalab/seurat-wrappers')

library(Seurat)
library(ggplot2)
library(stringr)
library(flexmix)
library(SeuratWrappers)

#READ IN DATA, CREATE SEURAT OBJECT, FILTER MITOCHONDRIAL DNA, 
mtx_directory = "C:/Users/diksh/Downloads/GSE225503_RAW"
dir_list = list.files(mtx_directory)

adapoe = NULL

for (i in seq_along(dir_list)){
  tmp = Read10X_h5(paste0(mtx_directory,"/",dir_list[1]))
  if(strsplit(dir_list[10],"_")[[1]][2] == "96wk"){
    tmpseu = CreateSeuratObject(counts = tmp)
    tmpseu$age = "2yr"
  }else{
    tmpseu = CreateSeuratObject(counts = tmp$`Gene Expression`)
    tmpseu[["HTO"]] = CreateAssayObject(counts = tmp$`Multiplexing Capture`)
    tmpseu$age = strsplit(dir_list[10],"_")[[1]][2]
  }
  tmpseu$orig.ident = dir_list[i]
  tmpseu = RenameCells(object = tmpseu, new.names = paste0(dir_list[i], "_", rownames(tmpseu[[]])))
  # MiQC filtering
  tmpseu[["percent.mt"]] = PercentageFeatureSet(tmpseu, pattern = "^mt-")
  tmpseu = RunMiQC(tmpseu, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", 
                   posterior.cutoff = 0.75, model.slot = "flexmix_model", model.type = "spline")
  tmpseu = subset(tmpseu, miQC.keep == "keep")
  # Merge as output
  if (i == 1){
    adapoe = tmpseu
  } else {
    adapoe = merge(x = adapoe, y = tmpseu)
  }
}
#NORMALIZE USING SC TRANSFORM




#COMPARE WITH AUTHOR'S PROCESSED RDS FILE
