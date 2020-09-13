library(SingleCellExperiment)
library(SC3)
data_set = c('biase','baron','goolam','kolodziejczyk','klein','pollen',
             'yan','deng','treutlein','rca','usoskin','zeisel','ting_149')
setwd("/home/math/hyl2016/SC3/")
Data_id = 13;
data = read.csv(paste0('./Data/',data_set[Data_id],".csv"),row.names= 1)
  ## --------------------------------------------------------------------------
  # create a SingleCellExperiment object
sce <- SingleCellExperiment(
    assays = list(counts = (10^(as.matrix(data))-1),logcounts = as.matrix(data))
  )
  # define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
  ## --------------------------------------------------------------------------
sce <- sc3_estimate_k(sce)
(metadata(sce)$sc3)$k_estimation
