# -*- coding:utf-8 _*-  
""" 
@author:jeffery
@time: 2019/05/27 18:01
@contact: jeffery_cpu@163.com
"""
import numpy as np
import pandas as pd
import scanpy as sc
from os.path import join, exists
from os import mkdir

#
sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
data_name = [ 'biase', 'yan', 'kolodziejczyk', 'pollen', 'goolam', 'deng2', 'klein2', 'zeisel', 'rca',	'usoskin',
             'treutlein','ting2_149','baron','Baron_initial']

for data_ID in range(1):
    data_ID = 13
    adata_old = pd.read_csv("E:/artical1/result_cluster_2019_5_29/initial_Data/" + data_name[data_ID]+'.csv', header=0)

    adata = sc.AnnData(np.transpose(adata_old.values[:, 1:]), obs=adata_old.columns[1:], var=adata_old.geneID)
    adata.obs_names = adata_old.columns[1:]
    adata.var_names = list(adata_old.geneID)
    file_dir = join(r'E:\artical1\result_cluster_2019_5_29\res_SelfPretreatment\Scanpy', data_name[data_ID])
    if not exists(file_dir):
        mkdir(file_dir)
    adata.var_names_make_unique()
    if data_ID != 12:
        sc.pp.filter_cells(adata, min_genes=0)
        sc.pp.filter_genes(adata, min_cells=1)
        mito_genes = adata.var_names.str.startswith('MT-')
        adata.obs['percent_mito'] = np.sum(
            adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
        adata.obs['n_counts'] = adata.X.sum(axis=1)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    if data_ID != 12:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var['highly_variable']]
        sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])  # !
    sc.pp.scale(adata, max_value=10)

    for i in range(1,51):
        sc.tl.pca(adata, svd_solver='arpack', n_comps=20)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
        sc.tl.louvain(adata)
        adata.obs['louvain'].to_csv(join(file_dir, 'res_' + str(i) + '.csv'), header=False, index=False)

