library(SingleCellExperiment)
library(SC3)
library(stringr)
library(Hmisc)
#setwd("/home/math/hyl2016/SC3/")
setwd("E:/artical1/result_cluster_2019_5_29/res_SelfPretreatment/SC3")
data_set = c('biase','goolam','kolodziejczyk','klein2','pollen',
             'yan','deng2','treutlein','rca','usoskin','zeisel','Baron_initial')
Data_id = 12;
#data = read.csv(paste0('./Data/',data_set[Data_id],".csv"),row.names= 1)
data = read.csv(paste0('../../initial_Data/',data_set[Data_id],".csv"),row.names= 1)
## --------------------------------------------------------------------------
# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data),
    logcounts = log2(as.matrix(data) + 1)
  )
) 
file_dir = paste0('./', data_set[Data_id])
if (!file.exists(file_dir))
{
  dir.create(file_dir)
}
# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# define spike-ins
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

## --------------------------------------------------------------------------
for(n_times in 1:50)
{
  sce2 <- sc3(sce, ks = 24 ,n_cores = 12)
  #sce2 <- sc3_run_svm(sce2, ks = 24)
  file_name = paste0('/res_', as.character(n_times), '.csv')
  write.csv(sce2@colData, file =paste0(file_dir, file_name))
}

