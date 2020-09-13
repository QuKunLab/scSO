# Title     : TODO
# Objective : TODO
# Created by: jeffery
# Created on: 2019-05-27

library(SingleCellExperiment)
library(Seurat)

# ----------------------------------------------------------------------------------------------
data_sets = c('biase','deng2','goolam', 'kolodziejczyk','pollen',
             'rca','yan','treutlein','usoskin','zeisel','ting2_149','klein2')
data_sets = c('baron')
data_sets = c('Baron_initial')
setwd("E:/artical1/result_cluster_2019_5_29/res_SelfPretreatment")
# ----------------------------------------------------------------------------------------------
n_pca_vector <- c()
for(data_set in data_sets)
{   
    file_dir = paste0('./Seurat/', data_set)
    if (!file.exists(file_dir))
    {
        dir.create(file_dir)
    }
    df = read.csv(paste0("../initial_Data/",data_set, '.csv'),row.names= 1,
    header=TRUE,fileEncoding = "UTF-8",stringsAsFactors=FALSE)
    pbmc <- CreateSeuratObject(counts = df, project = data_set)
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- FindVariableFeatures(pbmc)
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs=20)
    pbmc <- JackStraw(pbmc)
    pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

    p_values_ = pbmc@reductions$pca@jackstraw@overall.p.values[,2]
    count = 0
    for(val in p_values_)
    {
        if(val < 0.05)
        {
            count= count + 1
        }
        else{
            break
        }
    }
    if(count <=1){
        count = 2
    }
    n_pca_vector<-c(n_pca_vector  ,count)
    for(i in 1:50)
    {
      seurat_obj <- FindNeighbors(pbmc, dims = 1:count)
      seurat_obj <- FindClusters(seurat_obj)  # resolution = 0.5
      file_name = paste0('/res_', as.character(i), '.csv')
      write.csv(Idents(seurat_obj), file =paste0(file_dir, file_name))
    }
}
write.csv(n_pca_vector,file = 'seurat_Npca.csv')




