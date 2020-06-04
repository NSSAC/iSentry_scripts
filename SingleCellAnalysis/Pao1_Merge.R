library("Seurat")
#library("cowplot")
#library("ggplot2")
library("sctransform")
library("future")

plan("multiprocess",workers=4)
options(future.globals.maxSize = 120 * 1024 ^ 3)

pao1_top_path = "/scratch/cc8dm/SingleCell_PAO1/StarSolo/"
pao1_end_path = "/Solo.out/Gene/filtered/"

mito_genes = c("TRNF","RNR1","TRNV","RNR2","TRNL1","TRNI","TRNQ","TRNM","TRNW","TRNA","TRNN","TRNC","TRNY","TRNS1","TRND","TRNK","TRNG","TRNR","TRNH","TRNS2","TRNL2","TRNE","TRNT","TRNP")

#Pao1 Mock
pao1_mock_dir = paste(pao1_top_path,"Native_A549",pao1_end_path,sep="")
pao1_mock.data <- Read10X(data.dir=mock_dir)
pao1_mock = CreateSeuratObject(counts=mock.data)
pao1_mock[["Pao1_Genes"]] <- PercentageFeatureSet(mock,pattern=".peg.")
pao1_mock[["percent.mt"]] <- PercentageFeatureSet(mock,features=mito_genes)
pao1_mock$experiment = "Mock"

#Pao1 Mu_Hi
pao1_mu_hi_dir = paste(pao1_top_path,"S_high_mt",pao1_end_path,sep="")
pao1_mu_hi.data <- Read10X(data.dir=pao1_mu_hi_dir)
pao1_mu_hi = CreateSeuratObject(counts=pao1_mu_hi.data)
pao1_mu_hi[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_hi,pattern=".peg.")
pao1_mu_hi[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_hi,features=mito_genes)
pao1_mu_hi$experiment = "Mutant High"

#Pao1 WT_Hi
pao1_wt_hi_dir = paste(pao1_top_path,"S_high_wt",pao1_end_path,sep="")
pao1_wt_hi.data <- Read10X(data.dir=pao1_wt_hi_dir)
pao1_wt_hi = CreateSeuratObject(counts=pao1_wt_hi.data)
pao1_wt_hi[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_hi,pattern=".peg.")
pao1_wt_hi[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_hi,features=mito_genes)
pao1_wt_hi$experiment = "Wild High"

#Pao1 Mu_Lo
pao1_mu_lo_dir = paste(pao1_top_path,"S_low_mt",pao1_end_path,sep="")
pao1_mu_lo.data <- Read10X(data.dir=pao1_mu_lo_dir)
pao1_mu_lo = CreateSeuratObject(counts=pao1_mu_lo.data)
pao1_mu_lo[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_lo,pattern=".peg.")
pao1_mu_lo[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_lo,features=mito_genes)
pao1_mu_lo$experiment = "Mutant Low"

#Pao1 WT_Lo
pao1_wt_lo_dir = paste(pao1_top_path,"S_low_wt",pao1_end_path,sep="")
pao1_wt_lo.data <- Read10X(data.dir=pao1_wt_lo_dir)
pao1_wt_lo = CreateSeuratObject(counts=pao1_wt_lo.data)
pao1_wt_lo[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_lo,pattern=".peg.")
pao1_wt_lo[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_lo,features=mito_genes)
pao1_wt_lo$experiment = "Wild Low"

#Pao1 Mu_Neg
pao1_mu_neg_dir = paste(pao1_top_path,"S_neg_mt",pao1_end_path,sep="")
pao1_mu_neg.data <- Read10X(data.dir=pao1_mu_neg_dir)
pao1_mu_neg = CreateSeuratObject(counts=pao1_mu_neg.data)
pao1_mu_neg[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_neg,pattern=".peg.")
pao1_mu_neg[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_neg,features=mito_genes)
pao1_mu_neg$experiment = "Mutant Negative"

#Pao1 WT_Neg
pao1_wt_neg_dir = paste(pao1_top_path,"S_neg_wt",pao1_end_path,sep="")
pao1_wt_neg.data <- Read10X(data.dir=pao1_wt_neg_dir)
pao1_wt_neg = CreateSeuratObject(counts=pao1_wt_neg.data)
pao1_wt_neg[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_neg,pattern=".peg.")
pao1_wt_neg[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_neg,features=mito_genes)
pao1_wt_neg$experiment = "Wild Negative"

#Pao1 Mu_Med
pao1_mu_med_dir = paste(pao1_top_path,"S_med_mt",pao1_end_path,sep="")
pao1_mu_med.data <- Read10X(data.dir=pao1_mu_med_dir)
pao1_mu_med = CreateSeuratObject(counts=pao1_mu_med.data)
pao1_mu_med[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_med,pattern=".peg.")
pao1_mu_med[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_med,features=mito_genes)
pao1_mu_med$experiment =  "Mutant Medium"

#Pao1 WT_Med
pao1_wt_med_dir = paste(pao1_top_path,"S_med_wt",pao1_end_path,sep="")
pao1_wt_med.data <- Read10X(data.dir=pao1_wt_med_dir)
pao1_wt_med = CreateSeuratObject(counts=pao1_wt_med.data)
pao1_wt_med[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_med,pattern=".peg.")
pao1_wt_med[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_med,features=mito_genes)
pao1_wt_med$experiment =  "Wild Medium"

#Merge Objects
merge.object <- merge(pao1_mock,y=c(pao1_mu_hi,pao1_wt_hi,pao1_mu_lo,pao1_wt_lo,pao1_mu_neg,pao1_wt_neg,pao1_mu_med,pao1_wt_med),add.cell.ids=c("M","MH","WH","ML","WL","MN","WN","MM","WM"))

#Split single object into datasets
merge.experiment <- SplitObject(merge.object,split.by="experiment")

#SCTransform on each experiment
for (i in 1:length(merge.experiment)) {
	merge.experiment[[i]] <- SCTransform(merge.experiment[[i]],verbose=FALSE)
}

#Prepare for integration
merge.features <- SelectIntegrationFeatures(object.list = merge.experiment, nfeatures = 3000)
merge.experiment <- PrepSCTIntegration(object.list = merge.experiment,anchor.features=merge.features,verbose=FALSE)

#Find Integration Anchors
merge.anchors <- FindIntegrationAnchors(object.list = merge.experiment,normalization.method="SCT",anchor.features=merge.features,k.filter=50)
merge.integrated <- IntegrateData(anchorset = merge.anchors,normalization.method="SCT")

#switch to integrated assay
DefaultAssay(merge.integrated) <- "integrated"

#Clustering failed, save before that point
saveRDS(merge.integrated,file="Homo_Pao1_Seurat_PreClustering_CR3.rds")

#Standard visualization and clustering workflow
merge.integrated <- RunPCA(merge.integrated,npcs=30,verbose=FALSE)
merge.integrated <- RunUMAP(merge.integrated,reduction = "pca", dims=1:30)
merge.integrated <- FindNeighbors(merge.integrated,reduction="pca",dims=1:20)
merge.integrated <- FindClusters(merge.integrated,resolution=0.5)

#Save data object
saveRDS(merge.integrated,file="Homo_Pao1_Seurat_CR3.rds")
