logFile <- file("/home/cc8dm/SingleCell/Rlog.out")
cat("starting Seurat log: Loading Seurat library\n",file=logFile)

library("Seurat")
cat("Loading cowplot library\n",file="/home/cc8dm/SingleCell/Rlog.out",append=T)
library("cowplot")

#Set up Seurat Objects
cat("Setting up v1 Seurat object\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
v1.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_med_mt_1234/outs/filtered_feature_bc_matrix.h5")
v1 <- CreateSeuratObject(counts = v1.data,project = "Conc Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v1[["percent.mt"]] <- PercentageFeatureSet(v1,pattern = "^MT-")
#v1 <- subset(v1,subset = percent.mt < 5)
#insert condition to compare: example: v1$type = "Mutant", v1$concentration = "Low"
#v1$strain = "Wild"
v1$concentration = "Medium"
v1 <- NormalizeData(v1,verbose=F)
v1 <- FindVariableFeatures(v1,selection.method = "vst") #other selection methods exist, help("FindVariableFeatures")

cat("Setting up v2 Seurat object\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
v2.data <- Read10X_h5(filename = "/project/biocomplexity/isentry/djm_analysis/HK353DSXX/S_neg_mt_1234/outs/filtered_feature_bc_matrix.h5")
v2 <- CreateSeuratObject(counts = v2.data,project = "Conc Comparison")
#Filter samples with a high concentration of mitochondrial RNA,, setting threshold at 5%
#v2[["percent.mt"]] <- PercentageFeatureSet(v2,pattern = "^MT-")
#v2 <- subset(v2,subset = percent.mt < 5)
#insert condition to compare: example: v2$type = "Wild", v2$concentration = "High"
#v2$strain = "Mutant"
v2$concentration = "None"
v2 <- NormalizeData(v2,verbose=F)
v2 <- FindVariableFeatures(v2,selection.method = "vst")

#Perform integration
cat("Find integration anchors\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
v.anchors <- FindIntegrationAnchors(object.list = list(v1,v2), dims = 1:20) #Look into dims??
cat("IntegrateData\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
v.combined <- IntegrateData(anchorset = v.anchors, dims = 1:20)

#Perform integrated analysis
cat("Set default assay\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
DefaultAssay(v.combined) <- "integrated" 

#standard workflow for visualisation and clustering
cat("Scale data and run PCA\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
v.combined <- ScaleData(v.combined,verbose=F)
v.combined <- RunPCA (v.combined,npcs=30,verbose=F)
#t-SNE and Clustering
cat("Run UMAP, find cluster, find neighbors\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
v.combined <- RunUMAP(v.combined,reduction = "pca", dims = 1:20)
v.combined <- FindNeighbors(v.combined,reduction = "pca", dims = 1:20)
v.combined <- FindClusters(v.combined,resolution = 0.5)

#Visualization
cat("Visualize first plot\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
jpeg("Mutant_med_low.jpeg")
p1 <- DimPlot(v.combined,reduction = "umap", group.by = "concentration") #include condition from above for group.by
p2 <- DimPlot(v.combined, reduction = "umap", label = T)
plot_grid(p1,p2)

cat("Visualize second plot\n",file = "/home/cc8dm/SingleCell/Rlog.out",append=T)
jpeg("Mutant_med_low_contrast.jpeg")
DimPlot(v.combined,reduction = "umap", split.by = "concentration")
dev.off()
cat("Process complete\n",file = "Rlog.out",append=T)
