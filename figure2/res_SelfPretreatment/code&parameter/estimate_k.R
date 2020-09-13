library(SingleCellExperiment)
library(SC3)
#setwd("/home/math/hyl2016/SC3/")
setwd("E:/artical1/result_cluster_2019_5_29/res_SelfPretreatment/SC3")
data_set = c('biase','goolam','kolodziejczyk','klein2','pollen',
             'yan','deng2','treutlein','rca','usoskin','zeisel','Baron_initial')
Data_id = 12;
data = read.csv(paste0('../../initial_Data/',data_set[Data_id],".csv"),row.names= 1)
## --------------------------------------------------------------------------
# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data),
    logcounts = log2(as.matrix(data) + 1)
  )
) 
# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# define spike-ins
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)

## --------------------------------------------------------------------------
sce <- sc3_estimate_k(sce)
(metadata(sce)$sc3)$k_estimation
