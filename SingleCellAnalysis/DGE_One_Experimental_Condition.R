library(Seurat) 

data <- readRDS("../Merge_Pao1_Staph_Seurat.rds")

DefaultAssay(data) <- "SCT"
condition_param <- commandArgs(trailingOnly=TRUE)
condition = unlist(condition_param[2])
out_file = paste(condition,"CenterCluster",sep="_vs_")
out_file = paste(out_file,"_Staph.txt",sep="")
out_file = gsub(" ","_",out_file)
#Reassign cluster names
data <- subset(data,subset=species=="Staph")
data <- subset(data,subset=experiment==toString(condition))
Idents(data) <- "seurat_clusters"
data <- RenameIdents(data,'0'='A','1'='A','2'='A','3'='A')
#Concentration comparison table
exp.frame <- data.frame(Gene_Name=character(),Gene=character(),Experiment=character(),Cluster=character(),p_val=double(),avg_logFC=double(),pct.1=double(),pct.2=double(),p_val_adj=double())
#Split datasets by experiment
#Differential Expression: Wilcoxon ranked-sum test
for (i in 4:15) {
    base = data.frame(Experiment=condition,Cluster=toString(i))
    result = tryCatch({
        markers <- FindMarkers(data,ident.1=toString(i),ident.2="A",assay="SCT")
        res = cbind(base,markers)
        res = cbind(rownames(markers),res)
        exp.frame = rbind(exp.frame,res)
        remove(res)
        remove(markers)
        remove(Gene)
    }, warning = function(w) {
        print(paste("MY_WARNING: ",w))
    }, error = function(e) {
        print(paste("MY_ERROR: ",e))
    }, finally = {
        remove(base)
    })
} 
write.table(exp.frame,file=out_file,quote=F,sep="\t")
