#Load Necessary Libraries
library(Seurat)

#assuming the dataset file has been completely merged and is stored as an RDS file
#Final function performed on it was the UMAP function
alldata <- readRDS("All_Datasets.rds")
alldata <- FindNeighbors(alldata,reduction="pca",dims = 1:20)
alldata <- FindClusters(alldata,resolution=0.5)

#Concentration comparison table
conc.frame <- data.frame(Genotype=character(),Cluster=character(),Conc1=character(),Conc2=character(),p_val=double(),avg_logFC=double(),pct.1=double(),pct.2=double(),p_val_adj=double())
#Split datasets by strain
Idents(alldata) <- "strain"
alldata.wild <- subset(alldata,subset = strain == "Wild")
alldata.mutant <- subset(alldata,subset = strain == "Mutant")
lastcluster = length(levels(alldata$seurat_clusters)) - 1
#Mutant type
for (i in 0:lastcluster) {
    #Then, rowbind the dataframe to the growing frame
    currcluster <- subset(alldata.mutant,subset = seurat_clusters == toString(i))
    Idents(currcluster) <- "concentration"
    #low vs negative
    base = data.frame(Genotype="Mutant",Cluster=toString(i),Conc1="Low",Conc2="Negative")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Low",ident.2="Negative")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #medium vs negative
    base = data.frame(Genotype="Mutant",Cluster=toString(i),Conc1="Medium",Conc2="Negative")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Medium",ident.2="Negative")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #medium vs low
    base = data.frame(Genotype="Mutant",Cluster=toString(i),Conc1="Medium",Conc2="Low")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Medium",ident.2="Low")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #high vs negative
    base = data.frame(Genotype="Mutant",Cluster=toString(i),Conc1="High",Conc2="Negative")
    result = tryCatch({
                markers <-  FindMarkers(currcluster,ident.1="High",ident.2="Negative")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #high vs low
    base = data.frame(Genotype="Mutant",Cluster=toString(i),Conc1="High",Conc2="Low")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="High",ident.2="Low")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    remove(base)
    #high vs medium'
    base = data.frame(Genotype="Mutant",Cluster=toString(i),Conc1="High",Conc2="Medium")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="High",ident.2="Medium")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    remove(currcluster)
}

#wild type
for (i in 0:lastcluster) {
    #Then, rowbind the dataframe to the growing frame
    currcluster <- subset(alldata.wild,subset = seurat_clusters == toString(i))
    Idents(currcluster) <- "concentration"
    #low vs negative
    base = data.frame(Genotype="Wild",Cluster=toString(i),Conc1="Low",Conc2="Negative")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Low",ident.2="Negative")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #medium vs negative
    base = data.frame(Genotype="Wild",Cluster=toString(i),Conc1="Medium",Conc2="Negative")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Medium",ident.2="Negative")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #medium vs low
    base = data.frame(Genotype="Wild",Cluster=toString(i),Conc1="Medium",Conc2="Low")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Medium",ident.2="Low")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #high vs negative
    base = data.frame(Genotype="Wild",Cluster=toString(i),Conc1="High",Conc2="Negative")
    result = tryCatch({
                markers <-  FindMarkers(currcluster,ident.1="High",ident.2="Negative")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #high vs low
    base = data.frame(Genotype="Wild",Cluster=toString(i),Conc1="High",Conc2="Low")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="High",ident.2="Low")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    #high vs medium'
    base = data.frame(Genotype="Wild",Cluster=toString(i),Conc1="High",Conc2="Medium")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="High",ident.2="Medium")
                res = cbind(base,markers)
                conc.frame = rbind(conc.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    remove(currcluster)
}

write.table(conc.frame,file ="AllData_Concentration_Comparison.txt",quote=F,sep="\t")

#Strain comparison
strain.frame <- data.frame(Concentration=character(),Cluster=character(),Cond1=character(),Cond2=character(),p_val=double(),avg_logFC=double(),pct.1=double(),pct.2=double(),p_val_adj=double())
Idents(alldata) <- "concentration"
#lastcluster = length(levels(alldata$seurat_clusters)) - 1
alldata.Neg <- subset(alldata,subset = concentration == "Negative")
alldata.Low <- subset(alldata,subset = concentration == "Low")
alldata.Medium <- subset(alldata,subset = concentration == "Medium")
alldata.High <- subset(alldata,subset = concentration == "High")

#Negative
for (i in 0:lastcluster) {
    #Then, rowbind the dataframe to the growing frame
    currcluster <- subset(alldata.Neg,subset = seurat_clusters == toString(i))
    Idents(currcluster) <- "strain"
    #low vs negative
    base = data.frame(Concentration="Negative",Cluster=toString(i),Cond1="Mutant",Cond2="Wild")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Mutant",ident.2="Wild")
                res = cbind(base,markers)
                strain.frame = rbind(strain.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    remove(currcluster)
}

#Low
for (i in 0:lastcluster) {
    #Then, rowbind the dataframe to the growing frame
    currcluster <- subset(alldata.Low,subset = seurat_clusters == toString(i))
    Idents(currcluster) <- "strain"
    #low vs negative
    base = data.frame(Concentration="Low",Cluster=toString(i),Cond1="Mutant",Cond2="Wild")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Mutant",ident.2="Wild")
                res = cbind(base,markers)
                strain.frame = rbind(strain.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    remove(currcluster)
}

#Medium
for (i in 0:lastcluster) {
    #Then, rowbind the dataframe to the growing frame
    currcluster <- subset(alldata.Medium,subset = seurat_clusters == toString(i))
    Idents(currcluster) <- "strain"
    #low vs negative
    base = data.frame(Concentration="Medium",Cluster=toString(i),Cond1="Mutant",Cond2="Wild")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Mutant",ident.2="Wild")
                res = cbind(base,markers)
                strain.frame = rbind(strain.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    remove(currcluster)
}

#High
for (i in 0:lastcluster) {
    #Then, rowbind the dataframe to the growing frame
    currcluster <- subset(alldata.High,subset = seurat_clusters == toString(i))
    Idents(currcluster) <- "strain"
    #low vs negative
    base = data.frame(Concentration="High",Cluster=toString(i),Cond1="Mutant",Cond2="Wild")
    result = tryCatch({
                markers <- FindMarkers(currcluster,ident.1="Mutant",ident.2="Wild")
                res = cbind(base,markers)
                strain.frame = rbind(strain.frame,res)
                remove(res)
                remove(markers)
            }, warning = function(w) {
                print(paste("MY_WARNING:  ",w))
            }, error = function(e) {
                print(paste("MY_ERROR:  ",e))
            }, finally = {
                remove(base)
            })
    remove(currcluster)
}

write.table(strain.frame,file ="AllData_Strain_Comparison.txt",quote=F,sep="\t")

