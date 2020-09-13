# Title     : TODO
# Objective : TODO
# Created by: jeffery
# Created on: 2019/5/27
# Title     : Seurat Cluster
# Objective :
# Created by: jeffery
# Created on: 2019/5/23

library(SingleCellExperiment)
library(Seurat)

# -----------------------------------------------------------------------------------------------
setwd("E:/artical1/result_cluster_2019_5_29/res_OurPretreatment")
data_sets = c('biase','deng','goolam','klein', 'kolodziejczyk','pollen',
             'rca','yan','treutlein','usoskin','zeisel','ting_149', 'baron')
# -----------------------------------------------------------------------------------------------
data_sets = 'zeisel'
n_pca_vector = c()
for(data_set in data_sets)
{
    file_dir = paste0('./Seurat/', data_set)
    if (!file.exists(file_dir))
    {
        dir.create(file_dir)
    }
    df = read.csv(paste0("../Our_pretreatment_Data/",data_set, '.csv'),row.names= 1,
    header=TRUE,fileEncoding = "UTF-8",stringsAsFactors=FALSE)
    seurat_obj <- CreateSeuratObject(counts = df, project = data_set, min.cells = 0 , min.features = 0)
    all.genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)

        seurat_obj <- RunPCA(seurat_obj, features = row.names(seurat_obj), npcs=20)
        seurat_obj <- JackStraw(seurat_obj)
        seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
        p_values_ = seurat_obj@reductions$pca@jackstraw@overall.p.values[,2]
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
          seurat_obj <- FindNeighbors(seurat_obj, dims = 1:count)
          seurat_obj <- FindClusters(seurat_obj)  # resolution = 0.5
          file_name = paste0('/res_', as.character(i), '.csv')
          write.csv(Idents(seurat_obj), file =paste0(file_dir, file_name))
        }
}

write.csv(n_pca_vector,file = 'seurat_Npca.csv')

seurat_obj <- RunTSNE(seurat_obj, npcs=20)
DimPlot(seurat_obj, reduction = "tsne")






