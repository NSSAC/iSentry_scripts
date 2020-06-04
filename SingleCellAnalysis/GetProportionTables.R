data = readRDS("Merge_Pao1_Staph_Seurat.rds")

library(Seurat)

DefaultAssay(data) <- "SCT"

data.staph = subset(data,subset=species=="Pao1")
remove(data)

#conditions = c("Mock","Wild Negative","Mutant Negative","Wild Low","Mutant Low","Wild High","Mutant High")
conditions = c("Mock","Wild Negative","Mutant Negative","Wild Low","Mutant Low","Wild Medium","Mutant Medium","Wild High","Mutant High")

for (cond in conditions) {
    data.cond = subset(data.staph,subset=experiment==cond)
    seurat_clusters = droplevels(unique(data.cond$seurat_clusters))
    for (cluster in seurat_clusters) {
        data.cluster = subset(data.cond,subset=seurat_clusters==cluster)
        DefaultAssay(data.cluster) <- "SCT"
        data.dotplot = DotPlot(data.cluster,features=rownames(data.cluster))
        filename = paste("Pao1",cond,cluster,"DotPlot.txt",sep="_")
        filename = paste("Proportion_Tables",filename,sep="/")
        write.table(data.dotplot$data,filename,sep="\t",quote=F)
        remove(data.dotplot)
        remove(data.cluster)
    }
    remove(data.cond)
}
