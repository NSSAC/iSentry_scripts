logFile <- file("/home/cc8dm/SingleCell/Rlog_Objects.out")
cat("starting Seurat log: Loading Seurat library\n",file=logFile)

library("Seurat")
cat("Loading cowplot library\n",file="/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
library("cowplot")

#Set up Seurat Objects
#cat("Setting up Seurat object: Negative Wild Type\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wild_neg.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_neg_wt_1234/outs/filtered_feature_bc_matrix.h5")
wild_neg <- CreateSeuratObject(counts = wild_neg.data,project = "Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v1[["percent.mt"]] <- PercentageFeatureSet(v1,pattern = "^MT-")
#v1 <- subset(v1,subset = percent.mt < 5)
#insert condition to compare: example: v1$type = "Mutant", v1$concentration = "Low"
wild_neg$strain = "Wild"
wild_neg$concentration = "Negative"
wild_neg <- NormalizeData(wild_neg,verbose=F)
wild_neg <- FindVariableFeatures(wild_neg,selection.method = "vst") #other selection methods exist, help("FindVariableFeatures")

#cat("Setting up Seurat object: Negative Mutant Type\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mutant_neg.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_neg_mt_1234/outs/filtered_feature_bc_matrix.h5")
mutant_neg <- CreateSeuratObject(counts = mutant_neg.data,project = "Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v1[["percent.mt"]] <- PercentageFeatureSet(v1,pattern = "^MT-")
#v1 <- subset(v1,subset = percent.mt < 5)
#insert condition to compare: example: v1$type = "Mutant", v1$concentration = "Low"
mutant_neg$strain = "Mutant"
mutant_neg$concentration = "Negative"
mutant_neg <- NormalizeData(mutant_neg,verbose=F)
mutant_neg <- FindVariableFeatures(mutant_neg,selection.method = "vst") #other selection methods exist, help("FindVariableFeatures")

#cat("Setting up Seurat object: Low Wild Type\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wild_low.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_low_wt_1234/outs/filtered_feature_bc_matrix.h5")
wild_low <- CreateSeuratObject(counts = wild_low.data,project = "Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v1[["percent.mt"]] <- PercentageFeatureSet(v1,pattern = "^MT-")
#v1 <- subset(v1,subset = percent.mt < 5)
#insert condition to compare: example: v1$type = "Mutant", v1$concentration = "Low"
wild_low$strain = "Wild"
wild_low$concentration = "Low"
wild_low <- NormalizeData(wild_low,verbose=F)
wild_low <- FindVariableFeatures(wild_low,selection.method = "vst") #other selection methods exist, help("FindVariableFeatures")

#cat("Setting up Seurat object: Low Mutant Type\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mutant_low.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_low_mt_1234/outs/filtered_feature_bc_matrix.h5")
mutant_low <- CreateSeuratObject(counts = mutant_low.data,project = "Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v1[["percent.mt"]] <- PercentageFeatureSet(v1,pattern = "^MT-")
#v1 <- subset(v1,subset = percent.mt < 5)
#insert condition to compare: example: v1$type = "Mutant", v1$concentration = "Low"
mutant_low$strain = "Mutant"
mutant_low$concentration = "Low"
mutant_low <- NormalizeData(mutant_low,verbose=F)
mutant_low <- FindVariableFeatures(mutant_low,selection.method = "vst") #other selection methods exist, help("FindVariableFeatures")

#cat("Setting up Seurat object: Medium Wild Type\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wild_med.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_med_wt_1234/outs/filtered_feature_bc_matrix.h5")
wild_med <- CreateSeuratObject(counts = wild_med.data,project = "Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v1[["percent.mt"]] <- PercentageFeatureSet(v1,pattern = "^MT-")
#v1 <- subset(v1,subset = percent.mt < 5)
#insert condition to compare: example: v1$type = "Mutant", v1$concentration = "Low"
wild_med$strain = "Wild"
wild_med$concentration = "Medium"
wild_med <- NormalizeData(wild_med,verbose=F)
wild_med <- FindVariableFeatures(wild_med,selection.method = "vst") #other selection methods exist, help("FindVariableFeatures")

#cat("Setting up Seurat object: Medium Mutant Type\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mutant_med.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_med_mt_1234/outs/filtered_feature_bc_matrix.h5")
mutant_med <- CreateSeuratObject(counts = mutant_med.data,project = "Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v1[["percent.mt"]] <- PercentageFeatureSet(v1,pattern = "^MT-")
#v1 <- subset(v1,subset = percent.mt < 5)
#insert condition to compare: example: v1$type = "Mutant", v1$concentration = "Low"
mutant_med$strain = "Mutant"
mutant_med$concentration = "Medium"
mutant_med <- NormalizeData(mutant_med,verbose=F)
mutant_med <- FindVariableFeatures(mutant_med,selection.method = "vst") #other selection methods exist, help("FindVariableFeatures")

cat("Finished Setting up Seurat objects\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)

#Find Integration Anchors, Integrate Data, 
cat("Setting Up Integration Anchors: Negative Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
nwm.anchors <- FindIntegrationAnchors(object.list = list(wild_neg,mutant_neg), dims = 1:20) #Look into dims??
cat("Integrating Data: Negative Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
nwm.combined <- IntegrateData(anchorset = nwm.anchors, dims = 1:20)
DefaultAssay(nwm.combined) <- "integrated"
cat("Scale Data and PCA: Negative Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
nwm.combined <- ScaleData(nwm.combined,verbose=F)
nwm.combined <- RunPCA (nwm.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Negative Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
nwm.combined <- RunUMAP(nwm.combined,reduction = "pca", dims = 1:20)
nwm.combined <- FindNeighbors(nwm.combined,reduction = "pca", dims = 1:20)
nwm.combined <- FindClusters(nwm.combined,resolution = 0.5)
saveRDS(nwm.combined,file="Negative_Wild_Mutant.rds")
rm(nwm.anchors)
rm(nwm.combined)

cat("Setting Up Integration Anchors: Low Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
lwm.anchors <- FindIntegrationAnchors(object.list = list(wild_low,mutant_low), dims = 1:20) 
cat("Integrating Data: Low Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
lwm.combined <- IntegrateData(anchorset = lwm.anchors, dims = 1:20)
DefaultAssay(lwm.combined) <- "integrated"
cat("Scale Data and PCA: Low Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
lwm.combined <- ScaleData(lwm.combined,verbose=F)
lwm.combined <- RunPCA (lwm.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Low Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
lwm.combined <- RunUMAP(lwm.combined,reduction = "pca", dims = 1:20)
lwm.combined <- FindNeighbors(lwm.combined,reduction = "pca", dims = 1:20)
lwm.combined <- FindClusters(lwm.combined,resolution = 0.5)
saveRDS(lwm.combined,file="Low_Wild_Mutant.rds")
rm(lwm.anchors)
rm(lwm.combined)

cat("Setting Up Integration Anchors: Medium Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mwm.anchors <- FindIntegrationAnchors(object.list = list(wild_med,mutant_med), dims = 1:20) 
cat("Integrating Data: Medium Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mwm.combined <- IntegrateData(anchorset = mwm.anchors, dims = 1:20)
DefaultAssay(mwm.combined) <- "integrated"
cat("Scale Data and PCA: Medium Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mwm.combined <- ScaleData(mwm.combined,verbose=F)
mwm.combined <- RunPCA (mwm.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Medium Wild vs Mutant\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mwm.combined <- RunUMAP(mwm.combined,reduction = "pca", dims = 1:20)
mwm.combined <- FindNeighbors(mwm.combined,reduction = "pca", dims = 1:20)
mwm.combined <- FindClusters(mwm.combined,resolution = 0.5)
saveRDS(mwm.combined,file="Medium_Wild_Mutant.rds")
rm(mwm.anchors)
rm(mwm.combined)

cat("Setting Up Integration Anchors: Wild Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wln.anchors <- FindIntegrationAnchors(object.list = list(wild_low,wild_neg), dims = 1:20) 
cat("Integrating Data: Wild Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wln.combined <- IntegrateData(anchorset = wln.anchors, dims = 1:20)
DefaultAssay(wln.combined) <- "integrated"
cat("Scale Data and PCA: Wild Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wln.combined <- ScaleData(wln.combined,verbose=F)
wln.combined <- RunPCA (wln.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Wild Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wln.combined <- RunUMAP(wln.combined,reduction = "pca", dims = 1:20)
wln.combined <- FindNeighbors(wln.combined,reduction = "pca", dims = 1:20)
wln.combined <- FindClusters(wln.combined,resolution = 0.5)
saveRDS(wln.combined,file="Wild_Low_Negative.rds")
rm(wln.anchors)
rm(wln.combined)

cat("Setting Up Integration Anchors: Wild Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wmn.anchors <- FindIntegrationAnchors(object.list = list(wild_med,wild_neg), dims = 1:20)  
cat("Integrating Data: Wild Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wmn.combined <- IntegrateData(anchorset = wmn.anchors, dims = 1:20)
DefaultAssay(wmn.combined) <- "integrated"
cat("Scale Data and PCA: Wild Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wmn.combined <- ScaleData(wmn.combined,verbose=F)
wmn.combined <- RunPCA (wmn.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Wild Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wmn.combined <- RunUMAP(wmn.combined,reduction = "pca", dims = 1:20)
wmn.combined <- FindNeighbors(wmn.combined,reduction = "pca", dims = 1:20)
wmn.combined <- FindClusters(wmn.combined,resolution = 0.5)
saveRDS(wmn.combined,file="Wild_Medium_Negative.rds")
rm(wmn.anchors)
rm(wmn.combined)

cat("Setting Up Integration Anchors: Wild Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wml.anchors <- FindIntegrationAnchors(object.list = list(wild_med,wild_low), dims = 1:20) 
cat("Integrating Data: Wild Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wml.combined <- IntegrateData(anchorset = wml.anchors, dims = 1:20)
DefaultAssay(wml.combined) <- "integrated"
cat("Scale Data and PCA: Wild Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wml.combined <- ScaleData(wml.combined,verbose=F)
wml.combined <- RunPCA (wml.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Wild Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
wml.combined <- RunUMAP(wml.combined,reduction = "pca", dims = 1:20)
wml.combined <- FindNeighbors(wml.combined,reduction = "pca", dims = 1:20)
wml.combined <- FindClusters(wml.combined,resolution = 0.5)
saveRDS(wml.combined,file="Wild_Medium_Low.rds")
rm(wml.anchors)
rm(wml.combined)

cat("Setting Up Integration Anchors: Mutant Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mln.anchors <- FindIntegrationAnchors(object.list = list(mutant_low,mutant_neg), dims = 1:20)
cat("Integrating Data: Mutant Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mln.combined <- IntegrateData(anchorset = mln.anchors, dims = 1:20)
DefaultAssay(mln.combined) <- "integrated"
cat("Scale Data and PCA: Mutant Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mln.combined <- ScaleData(mln.combined,verbose=F)
mln.combined <- RunPCA (mln.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Mutant Low vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mln.combined <- RunUMAP(mln.combined,reduction = "pca", dims = 1:20)
mln.combined <- FindNeighbors(mln.combined,reduction = "pca", dims = 1:20)
mln.combined <- FindClusters(mln.combined,resolution = 0.5)
saveRDS(mln.combined,file="Mutant_Low_Negative.rds")
rm(mln.anchors)
rm(mln.combined)

cat("Setting Up Integration Anchors: Mutant Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mmn.anchors <- FindIntegrationAnchors(object.list = list(mutant_med,mutant_neg), dims = 1:20)
cat("Integrating Data: Mutant Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mmn.combined <- IntegrateData(anchorset = mmn.anchors, dims = 1:20)
DefaultAssay(mmn.combined) <- "integrated"
cat("Scale Data and PCA: Mutant Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mmn.combined <- ScaleData(mmn.combined,verbose=F)
mmn.combined <- RunPCA (mmn.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Mutant Medium vs Negative\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mmn.combined <- RunUMAP(mmn.combined,reduction = "pca", dims = 1:20)
mmn.combined <- FindNeighbors(mmn.combined,reduction = "pca", dims = 1:20)
mmn.combined <- FindClusters(mmn.combined,resolution = 0.5)
saveRDS(mmn.combined,file="Mutant_Medium_Negative.rds")
rm(mmn.anchors)
rm(mmn.combined)

cat("Setting Up Integration Anchors: Mutant Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mml.anchors <- FindIntegrationAnchors(object.list = list(mutant_med,mutant_low), dims = 1:20)
cat("Integrating Data: Mutant Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mml.combined <- IntegrateData(anchorset = mml.anchors, dims = 1:20)
DefaultAssay(mml.combined) <- "integrated"
cat("Scale Data and PCA: Mutant Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mml.combined <- ScaleData(mml.combined,verbose=F)
mml.combined <- RunPCA (mml.combined,npcs=30,verbose=F)
cat("Run UMAP, find cluster, find neighbors: Mutant Medium vs Low\n",file = "/home/cc8dm/SingleCell/Rlog_Objects.out",append=T)
mml.combined <- RunUMAP(mml.combined,reduction = "pca", dims = 1:20)
mml.combined <- FindNeighbors(mml.combined,reduction = "pca", dims = 1:20)
mml.combined <- FindClusters(mml.combined,resolution = 0.5)
saveRDS(mln.combined,file="Mutant_Medium_Low.rds")
rm(mml.anchors)
rm(mml.combined)


cat("Process complete\n",file = "Rlog_Objects.out",append=T)
