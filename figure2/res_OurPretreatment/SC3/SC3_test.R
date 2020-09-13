library(SingleCellExperiment)
library(SC3)
library(stringr)
library(Hmisc)
setwd("/home/math/hyl2016/SC3/")
data_set = c('biase','baron','goolam','kolodziejczyk','klein','pollen',
             'yan','deng','treutlein','rca','usoskin','zeisel','ting_149')
data_set = c('Baron_initial')
  Data_id = 1;
  data = read.csv(paste0('./Data/',data_set[Data_id],".csv"),row.names= 1)
  data = read.csv(paste0('E:/artical1/result_cluster_2019_5_29/Our_pretreatment_Data/Baron_initial.csv'),row.names= 1)
  ## --------------------------------------------------------------------------
  # create a SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(counts = (10^(as.matrix(data))-1),logcounts = as.matrix(data))
  ) 
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  ## --------------------------------------------------------------------------
for(n_times in 2:50)
{
  sce2 <- sc3(sce, ks = 22 ,gene_filter = FALSE,n_cores = 4)
 # sce2 <- sc3_run_svm(sce2, ks = 24)
  resname <- str_c("res_",n_times,".csv")  
  #path2 = paste0("./",data_set[Data_id],"/",resname)
  path2 = paste0("E:/artical1/result_cluster_2019_5_29/res_OurPretreatment/SC3/Baron_initial/",resname)
  write.csv(sce2@colData,file = path2)
}

  