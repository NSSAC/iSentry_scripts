library("Seurat")
library("cowplot")
library("ggplot2")

#Set up Seurat objects
wild_neg.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_neg_wt_1234/outs/filtered_feature_bc_matrix.h5")
wild_neg <- CreateSeuratObject(counts = wild_neg.data,project = "Wild Negative")
wild_neg$strain = "Wild"
wild_neg$concentration = "Negative"
wild_neg$experiment = "Wild Negative"
#wild_neg <- NormalizeData(wild_neg,verbose=F)
#wild_neg <- FindVariableFeatures(wild_neg,selection.method = "vst")

mutant_neg.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_neg_mt_1234/outs/filtered_feature_bc_matrix.h5")
mutant_neg <- CreateSeuratObject(counts = mutant_neg.data,project = "Mutant Negative")
mutant_neg$strain = "Mutant"
mutant_neg$concentration = "Negative"
mutant_neg$experiment = "Mutant Negative"
#mutant_neg <- NormalizeData(mutant_neg,verbose=F)
#mutant_neg <- FindVariableFeatures(mutant_neg,selection.method = "vst")

wild_low.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_low_wt_1234/outs/filtered_feature_bc_matrix.h5")
wild_low <- CreateSeuratObject(counts = wild_low.data,project = "Wild Low")
wild_low$strain = "Wild"
wild_low$concentration = "Low"
wild_low$experiment = "Wild Low"
#wild_low <- NormalizeData(wild_low,verbose=F)
#wild_low <- FindVariableFeatures(wild_low,selection.method = "vst")

mutant_low.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_low_mt_1234/outs/filtered_feature_bc_matrix.h5")
mutant_low <- CreateSeuratObject(counts = mutant_low.data,project = "Mutant Low")
mutant_low$strain = "Mutant"
mutant_low$concentration = "Low"
mutant_low$experiment = "Mutant Low"
#mutant_low <- NormalizeData(mutant_low,verbose=F)
#mutant_low <- FindVariableFeatures(mutant_low,selection.method = "vst")

wild_med.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_med_wt_1234/outs/filtered_feature_bc_matrix.h5")
wild_med <- CreateSeuratObject(counts = wild_med.data,project = "Wild Medium")
wild_med$strain = "Wild"
wild_med$concentration = "Medium"
wild_med$experiment = "Wild Medium"
#wild_med <- NormalizeData(wild_med,verbose=F)
#wild_med <- FindVariableFeatures(wild_med,selection.method = "vst")

mutant_med.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_med_mt_1234/outs/filtered_feature_bc_matrix.h5")
mutant_med <- CreateSeuratObject(counts = mutant_med.data,project = "Mutant Medium")
mutant_med$strain = "Mutant"
mutant_med$concentration = "Medium"
mutant_med$experiment = "Mutant Medium"
#mutant_med <- NormalizeData(mutant_med,verbose=F)
#mutant_med <- FindVariableFeatures(mutant_med,selection.method = "vst")

wild_high.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_high_wt_1234/outs/filtered_feature_bc_matrix.h5")
wild_high <- CreateSeuratObject(counts = wild_high.data,project = "Wild High")
wild_high$strain = "Wild"
wild_high$concentration = "High"
wild_high$experiment = "Wild High"
#wild_high <- NormalizeData(wild_high,verbose=F)
#wild_high <- FindVariableFeatures(wild_high,selection.method = "vst")

mutant_high.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_high_mt_1234/outs/filtered_feature_bc_matrix.h5")
mutant_high <- CreateSeuratObject(counts = mutant_high.data,project = "Mutant High")
mutant_high$strain = "Mutant"
mutant_high$concentration = "High"
mutant_high$experiment = "Mutant High"
#mutant_high <- NormalizeData(mutant_high,verbose=F)
#mutant_high <- FindVariableFeatures(mutant_high,selection.method = "vst")

#Merge Seurat Objects
merge.object <- merge(wild_neg, y = c(mutant_neg,wild_low,mutant_low,wild_med,mutant_med,wild_high,mutant_high), add.cell.ids = c("WN","MN","WL","ML","WM","MM","WH","MH"), project = "Combined Analysis")

#Split single object into datasets
merge.experiment <- SplitObject(merge.object,split.by = "experiment")

#Normalize and Find Variable Features
for (i in 1:length(merge.experiment)) {
    merge.experiment[[i]] <- NormalizeData(merge.experiment[[i]], verbose = F)
    merge.experiment[[i]] <- FindVariableFeatures(merge.experiment[[i]],selection.method = "vst")
}

#Find Integration Anchors
reference.list <- merge.experiment[c("Wild Negative","Mutant Negative","Wild Low","Mutant Low","Wild Medium","Mutant Medium","Wild High", "Mutant High")]
merge.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) #Will terminate if a data file does not have a minimum number of cells (k.filter) 
merge.integrated <- IntegrateData(anchorset = merge.anchors, dims = 1:30)

#switch to integrated assay. Can switch back to unintigrated using the "RNA" assay
DefaultAssay(merge.integrated) <- "integrated"

#Run standard visualization and clustering workflow
merge.integrated <- ScaleData(merge.integrated, verbose = F)
merge.integrated <- RunPCA(merge.integrated, npcs=30, verbose = F)
merge.integrated <- RunUMAP(merge.integrated, reduction = "pca", dims = 1:30)
merge.integrated <- FindNeighbors(alldata,reduction="pca",dims = 1:20)
merge.integrated <- FindClusters(alldata,resolution=0.5)

#Save data object
saveRDS(merge.integrated,file="All_Datasets.rds")

