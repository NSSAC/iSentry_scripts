library("Seurat")
#library("cowplot")
#library("ggplot2")
library("sctransform")
library("future")

plan("multiprocess",workers=2)
options(future.globals.maxSize = 200 * 1024 ^ 3)

staph_top_path = "/scratch/cc8dm/SingleCell_StaphJE2/"
staph_end_path = "/Solo.out/Gene/filtered/"
pao1_top_path = "/scratch/cc8dm/SingleCell_PAO1/StarSolo/"
pao1_end_path = "/Solo.out/Gene/filtered/"

mito_genes = c("TRNF","RNR1","TRNV","RNR2","TRNL1","TRNI","TRNQ","TRNM","TRNW","TRNA","TRNN","TRNC","TRNY","TRNS1","TRND","TRNK","TRNG","TRNR","TRNH","TRNS2","TRNL2","TRNE","TRNT","TRNP")

#Staph Mock
staph_mock_dir = paste(staph_top_path,"Mock",staph_end_path,sep="")
staph_mock.data <- Read10X(data.dir=staph_mock_dir)
staph_mock = CreateSeuratObject(counts=staph_mock.data)
staph_mock[["Staph_Genes"]] <- PercentageFeatureSet(staph_mock,pattern=".peg.")
staph_mock[["percent.mt"]] <- PercentageFeatureSet(staph_mock,features=mito_genes)
staph_mock$experiment = "Mock"
staph_mock$species = "Staph"
staph_mock$condition = "Staph Mock"

#Staph Mu_Hi
staph_mu_hi_dir = paste(staph_top_path,"Mu_Hi",staph_end_path,sep="")
staph_mu_hi.data <- Read10X(data.dir=staph_mu_hi_dir)
staph_mu_hi = CreateSeuratObject(counts=staph_mu_hi.data)
staph_mu_hi[["Staph_Genes"]] <- PercentageFeatureSet(staph_mu_hi,pattern=".peg.")
staph_mu_hi[["percent.mt"]] <- PercentageFeatureSet(staph_mu_hi,features=mito_genes)
staph_mu_hi$experiment = "Mutant High"
staph_mu_hi$species = "Staph"
staph_mu_hi$condition = "Staph Mutant High"

#Staph WT_Hi
staph_wt_hi_dir = paste(staph_top_path,"WT_Hi",staph_end_path,sep="")
staph_wt_hi.data <- Read10X(data.dir=staph_wt_hi_dir)
staph_wt_hi = CreateSeuratObject(counts=staph_wt_hi.data)
staph_wt_hi[["Staph_Genes"]] <- PercentageFeatureSet(staph_wt_hi,pattern=".peg.")
staph_wt_hi[["percent.mt"]] <- PercentageFeatureSet(staph_wt_hi,features=mito_genes)
staph_wt_hi$experiment = "Wild High"
staph_wt_hi$species = "Staph"
staph_wt_hi$condition = "Staph Wild High"

#Staph Mu_Lo
staph_mu_lo_dir = paste(staph_top_path,"Mu_Lo",staph_end_path,sep="")
staph_mu_lo.data <- Read10X(data.dir=staph_mu_lo_dir)
staph_mu_lo = CreateSeuratObject(counts=staph_mu_lo.data)
staph_mu_lo[["Staph_Genes"]] <- PercentageFeatureSet(staph_mu_lo,pattern=".peg.")
staph_mu_lo[["percent.mt"]] <- PercentageFeatureSet(staph_mu_lo,features=mito_genes)
staph_mu_lo$experiment = "Mutant Low"
staph_mu_lo$species = "Staph"
staph_mu_lo$condition = "Staph Mutant Low"

#Staph WT_Lo
staph_wt_lo_dir = paste(staph_top_path,"WT_Lo",staph_end_path,sep="")
staph_wt_lo.data <- Read10X(data.dir=staph_wt_lo_dir)
staph_wt_lo = CreateSeuratObject(counts=staph_wt_lo.data)
staph_wt_lo[["Staph_Genes"]] <- PercentageFeatureSet(staph_wt_lo,pattern=".peg.")
staph_wt_lo[["percent.mt"]] <- PercentageFeatureSet(staph_wt_lo,features=mito_genes)
staph_wt_lo$experiment = "Wild Low"
staph_wt_lo$species = "Staph"
staph_wt_lo$condition = "Staph Wild Low"

#Staph Mu_Neg
staph_mu_neg_dir = paste(staph_top_path,"Mu_Neg",staph_end_path,sep="")
staph_mu_neg.data <- Read10X(data.dir=staph_mu_neg_dir)
staph_mu_neg = CreateSeuratObject(counts=staph_mu_neg.data)
staph_mu_neg[["Staph_Genes"]] <- PercentageFeatureSet(staph_mu_neg,pattern=".peg.")
staph_mu_neg[["percent.mt"]] <- PercentageFeatureSet(staph_mu_neg,features=mito_genes)
staph_mu_neg$experiment = "Mutant Negative"
staph_mu_neg$species = "Staph"
staph_mu_neg$condition = "Staph Mutant Negative"

#Staph WT_Neg
staph_wt_neg_dir = paste(staph_top_path,"WT_Neg",staph_end_path,sep="")
staph_wt_neg.data <- Read10X(data.dir=staph_wt_neg_dir)
staph_wt_neg = CreateSeuratObject(counts=staph_wt_neg.data)
staph_wt_neg[["Staph_Genes"]] <- PercentageFeatureSet(staph_wt_neg,pattern=".peg.")
staph_wt_neg[["percent.mt"]] <- PercentageFeatureSet(staph_wt_neg,features=mito_genes)
staph_wt_neg$experiment = "Wild Negative"
staph_wt_neg$species = "Staph"
staph_wt_neg$condition = "Staph Wild Negative"

#Pao1 Mock
pao1_mock_dir = paste(pao1_top_path,"Native_A549",pao1_end_path,sep="")
pao1_mock.data <- Read10X(data.dir=pao1_mock_dir)
pao1_mock = CreateSeuratObject(counts=pao1_mock.data)
pao1_mock[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mock,pattern=".peg.")
pao1_mock[["percent.mt"]] <- PercentageFeatureSet(pao1_mock,features=mito_genes)
pao1_mock$experiment = "Mock"
pao1_mock$species = "Pao1"
pao1_mock$condition = "Pao1 Mock"

#Pao1 Mu_Hi
pao1_mu_hi_dir = paste(pao1_top_path,"S_high_mt",pao1_end_path,sep="")
pao1_mu_hi.data <- Read10X(data.dir=pao1_mu_hi_dir)
pao1_mu_hi = CreateSeuratObject(counts=pao1_mu_hi.data)
pao1_mu_hi[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_hi,pattern=".peg.")
pao1_mu_hi[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_hi,features=mito_genes)
pao1_mu_hi$experiment = "Mutant High"
pao1_mu_hi$species = "Pao1"
pao1_mu_hi$condition = "Pao1 Mutant High"

#Pao1 WT_Hi
pao1_wt_hi_dir = paste(pao1_top_path,"S_high_wt",pao1_end_path,sep="")
pao1_wt_hi.data <- Read10X(data.dir=pao1_wt_hi_dir)
pao1_wt_hi = CreateSeuratObject(counts=pao1_wt_hi.data)
pao1_wt_hi[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_hi,pattern=".peg.")
pao1_wt_hi[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_hi,features=mito_genes)
pao1_wt_hi$experiment = "Wild High"
pao1_wt_hi$species = "Pao1"
pao1_wt_hi$condition = "Pao1 Wild High"

#Pao1 Mu_Lo
pao1_mu_lo_dir = paste(pao1_top_path,"S_low_mt",pao1_end_path,sep="")
pao1_mu_lo.data <- Read10X(data.dir=pao1_mu_lo_dir)
pao1_mu_lo = CreateSeuratObject(counts=pao1_mu_lo.data)
pao1_mu_lo[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_lo,pattern=".peg.")
pao1_mu_lo[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_lo,features=mito_genes)
pao1_mu_lo$experiment = "Mutant Low"
pao1_mu_lo$species = "Pao1"
pao1_mu_lo$condition = "Pao1 Mutant Low"

#Pao1 WT_Lo
pao1_wt_lo_dir = paste(pao1_top_path,"S_low_wt",pao1_end_path,sep="")
pao1_wt_lo.data <- Read10X(data.dir=pao1_wt_lo_dir)
pao1_wt_lo = CreateSeuratObject(counts=pao1_wt_lo.data)
pao1_wt_lo[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_lo,pattern=".peg.")
pao1_wt_lo[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_lo,features=mito_genes)
pao1_wt_lo$experiment = "Wild Low"
pao1_wt_lo$species = "Pao1"
pao1_wt_lo$condition = "Pao1 Wild Low"

#Pao1 Mu_Neg
pao1_mu_neg_dir = paste(pao1_top_path,"S_neg_mt",pao1_end_path,sep="")
pao1_mu_neg.data <- Read10X(data.dir=pao1_mu_neg_dir)
pao1_mu_neg = CreateSeuratObject(counts=pao1_mu_neg.data)
pao1_mu_neg[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_neg,pattern=".peg.")
pao1_mu_neg[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_neg,features=mito_genes)
pao1_mu_neg$experiment = "Mutant Negative"
pao1_mu_neg$species = "Pao1"
pao1_mu_neg$condition = "Pao1 Mutant Negative"

#Pao1 WT_Neg
pao1_wt_neg_dir = paste(pao1_top_path,"S_neg_wt",pao1_end_path,sep="")
pao1_wt_neg.data <- Read10X(data.dir=pao1_wt_neg_dir)
pao1_wt_neg = CreateSeuratObject(counts=pao1_wt_neg.data)
pao1_wt_neg[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_neg,pattern=".peg.")
pao1_wt_neg[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_neg,features=mito_genes)
pao1_wt_neg$experiment = "Wild Negative"
pao1_wt_neg$species = "Pao1"
pao1_wt_neg$condition = "Pao1 Wild Negative"

#Pao1 Mu_Med
pao1_mu_med_dir = paste(pao1_top_path,"S_med_mt",pao1_end_path,sep="")
pao1_mu_med.data <- Read10X(data.dir=pao1_mu_med_dir)
pao1_mu_med = CreateSeuratObject(counts=pao1_mu_med.data)
pao1_mu_med[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_mu_med,pattern=".peg.")
pao1_mu_med[["percent.mt"]] <- PercentageFeatureSet(pao1_mu_med,features=mito_genes)
pao1_mu_med$experiment =  "Mutant Medium"
pao1_mu_med$species = "Pao1"
pao1_mu_med$condition = "Pao1 Mutant Medium"

#Pao1 WT_Med
pao1_wt_med_dir = paste(pao1_top_path,"S_med_wt",pao1_end_path,sep="")
pao1_wt_med.data <- Read10X(data.dir=pao1_wt_med_dir)
pao1_wt_med = CreateSeuratObject(counts=pao1_wt_med.data)
pao1_wt_med[["Pao1_Genes"]] <- PercentageFeatureSet(pao1_wt_med,pattern=".peg.")
pao1_wt_med[["percent.mt"]] <- PercentageFeatureSet(pao1_wt_med,features=mito_genes)
pao1_wt_med$experiment =  "Wild Medium"
pao1_wt_med$species = "Pao1"
pao1_wt_med$condition = "Pao1 Wild Medium"

#Merge Objects
merge.object <- merge(staph_mock,y=c(staph_mu_hi,staph_wt_hi,staph_mu_lo,staph_wt_lo,staph_mu_neg,staph_wt_neg,pao1_mock,pao1_mu_hi,pao1_wt_hi,pao1_mu_lo,pao1_wt_lo,pao1_mu_neg,pao1_wt_neg,pao1_mu_med,pao1_wt_med),add.cell.ids=c("S_M","S_MH","S_WH","S_ML","S_WL","S_MN","S_WN","P_M","P_MH","P_WH","P_ML","P_WL","P_MN","P_WN","P_MM","P_WM"))

#Split single object into datasets
merge.condition <- SplitObject(merge.object,split.by="condition")

#SCTransform on each condition
for (i in 1:length(merge.condition)) {
	merge.condition[[i]] <- SCTransform(merge.condition[[i]],verbose=FALSE)
}

#Prepare for integration
merge.features <- SelectIntegrationFeatures(object.list = merge.condition, nfeatures = 3000)
merge.condition <- PrepSCTIntegration(object.list = merge.condition,anchor.features=merge.features,verbose=FALSE)

#Find Integration Anchors
merge.anchors <- FindIntegrationAnchors(object.list = merge.condition,normalization.method="SCT",anchor.features=merge.features,k.filter=50)
merge.integrated <- IntegrateData(anchorset = merge.anchors,normalization.method="SCT")

#switch to integrated assay
DefaultAssay(merge.integrated) <- "integrated"

#Clustering failed, save before that point
saveRDS(merge.integrated,file="Merge_Pao1_Staph_PreClustering_Seurat.rds")

#Standard visualization and clustering workflow
merge.integrated <- RunPCA(merge.integrated,npcs=30,verbose=FALSE)
merge.integrated <- RunUMAP(merge.integrated,reduction = "pca", dims=1:30)
merge.integrated <- FindNeighbors(merge.integrated,reduction="pca",dims=1:20)
merge.integrated <- FindClusters(merge.integrated,resolution=0.5)

#Save data object
saveRDS(merge.integrated,file="Merge_Pao1_Staph_Seurat.rds")
