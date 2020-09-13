#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/5/9 10:55
# @Author  : YinLei Hu
import numpy as np
import pandas as pd
import scanpy as sc
from os.path import join, exists
from os import mkdir

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/pbmc3k.h5ad'  # the file that will store the analysis results
sc.settings.set_figure_params(dpi=80)
npca_set = [2,20,9,19,27,14,12,13,2,22,17,17,6]

np.random.seed(2019)

data_name = ['biase', 'yan', 'kolodziejczyk', 'pollen', 'goolam', 'deng', 'klein', 'zeisel', 'rca',	'usoskin',
             'treutlein',	'ting_149', 'baron']

for data_ID in range(13):
    file_dir= join(r"E:\artical1\result_cluster_2019_5_29\res_OurPretreatment\Scanpy", data_name[data_ID])
    if not exists(file_dir):
        mkdir(file_dir)
    adata_old = pd.read_csv("E:/artical1/result_cluster_2019_5_29/Data/" + data_name[data_ID] + '.csv' , header = 0)
    adata = sc.AnnData(np.transpose(adata_old.values[:,1:]), obs= adata_old.columns[1:], var=adata_old.Row)
    adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
    adata.raw = adata

    for i in range(1, 51):
        # tmp = adata.raw
        sc.tl.pca(adata, svd_solver='arpack',n_comps=20)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
        sc.tl.louvain(adata)
        adata.obs['louvain'].to_csv(join(file_dir, 'res_'+ str(i) + '.csv'))
