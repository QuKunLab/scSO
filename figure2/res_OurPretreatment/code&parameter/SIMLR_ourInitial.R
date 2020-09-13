# Title     : TODO
# Objective : TODO
# Created by: jeffery
# Created on: 2019-05-26
# Title     : SIMLR Cluster
# Objective :
# Created by: jeffery
# Created on: 2019/5/23
#
library(SIMLR)
setwd("E:/artical1/result_cluster_2019_5_29/res_OurPretreatment")
#'biase','deng','goolam',
data_sets = c('klein','kolodziejczyk','pollen',
              'rca','yan','treutlein','usoskin','zeisel','ting_149','baron')
data_sets = c('biase','yan','treutlein','goolam')
for(data_set in data_sets)
{
    file_dir = paste0('./SIMLR/', data_set)
    if (!file.exists(file_dir))
    {
        dir.create(file_dir)
    }
    df = read.csv(paste0('../Our_pretreatment_Data/', data_set, '.csv'),row.names= 1, header=TRUE,fileEncoding = "UTF-8",stringsAsFactors=FALSE)
    c = 0;
    NUMC = 2:15  # 20 15
    if(dim(df)[2] > 4000)
    {
      id_cell <- sample(1:dim(df)[2],3000)
      res_example = SIMLR_Estimate_Number_of_Clusters(df[,id_cell], NUMC = NUMC)
      c= NUMC[which.min(res_example$K1)]
    }
    else
    {
      res_example = SIMLR_Estimate_Number_of_Clusters(df, NUMC = NUMC)
      c= NUMC[which.min(res_example$K1)]
    }
    

    if(dim(df)[2] < 150)
    {
        # small size data
        for(i in c(1:50))
        {
            res = SIMLR(X = df, c = c)
            file_name = paste0('/res_',as.character(i),'.csv')
            write.csv(res$y$cluster, file =paste0(file_dir, file_name),col.names = F)
        }
    }
    else
    {
        for(i in c(1:50))
        {
            example_large_scale = SIMLR_Large_Scale(X = df, c = c, kk=30)
            file_name = paste0('/res_',as.character(i),'.csv')
            write.csv(example_large_scale$y$cluster, file =paste0(file_dir, file_name),col.names = F)
        }
    }
}

