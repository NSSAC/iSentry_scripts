library(Seurat) 

data <- readRDS("../Merge_Pao1_Staph_Seurat.rds")

DefaultAssay(data) <- "SCT"
#condition_param <- commandArgs(trailingOnly=TRUE)
#condition = unlist(condition_param[2])
condition = "Wild High"
out_file = paste(condition,"_Pao1_vs_Staph.txt",sep="")
#Reassign cluster names
data <- subset(data,subset=experiment==toString(condition))
Idents(data) <- "species"
#Concentration comparison table
exp.frame <- data.frame(Gene_Name=character(),Gene=character(),Experiment=character(),Cluster=character(),p_val=double(),avg_logFC=double(),pct.1=double(),pct.2=double(),p_val_adj=double())
#Differential Expression: Wilcoxon ranked-sum test
for (i in 0:14) {
    base = data.frame(Experiment=condition,Cluster=toString(i))
    curr_data <- subset(data,subset=seurat_clusters==toString(i))
    result = tryCatch({
        markers <- FindMarkers(curr_data,ident.1="Pao1",ident.2="Staph",assay="SCT")
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
