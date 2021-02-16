import scSO as ss
import numpy as np
import time
from sklearn import metrics
import pandas as pd

test1 = ss.scSO()
data_set = ['biase.csv', "deng.csv", 'goolam.csv', 'yan.csv', 'treutlein.csv', 'kolodziejczyk.csv',
            'usoskin.csv', 'zeisel.csv', 'pollen.csv', 'rca.csv','Baron.csv', 'klein.csv']
ARI_fileName = "./ARI.txt"
with open(ARI_fileName , encoding="utf-8",mode="w") as file:  
    file.write("DataName\tNCluster\tARI\n") 
for data_name in data_set:
    test1.read_data(data_path="../../test/"+data_name)
    print(np.shape(test1.Data))

    test1.filter_genes(bound_low=0.1, bound_up=8.5)
    test1.LogNormalize()
    test1.scSO_SNMF()
    test1.GMM_fit_and_calcu_0Space()
    labels_true =  pd.read_csv(filepath_or_buffer="../../figure2/ref/"+data_name,
                                sep=',',
                                header=None,
                                index_col=None).values
    y = test1.Cluster_func(Lambda=0.15) #default 0.15
    
    
    # np.savetxt(data_name, y, delimiter=',')
    ARI = metrics.adjusted_rand_score(labels_true[:,0], y)
    print(data_name,len(np.unique(y)), ARI)

    with open(ARI_fileName , encoding="utf-8",mode="a") as file:  
        file.write(data_name[:-4]+"\t%d\t%f\n"%(len(np.unique(y)), ARI))   