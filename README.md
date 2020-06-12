
<!-- TOC -->
- [SMAFS](#smafs)
    - [Installation](#installation)
    - [GUI Tutorial](#gui-tutorial)
    - [Command-Line Tutorial](#command-line-tutorial)
    - [Help](#help)
<!-- TOC -->


SMAFS
===
![](./images/Fig.jpg)

SMAFS is a new algorithm for scRNA-seq data clustering based on the low-rank hypothesis of the gene expression matrix, and the combination of the sparse non-negative matrix factorization and the spectral method.

### Installation

Because SMAFS relies on the [cvx](http://cvxr.com/cvx/download/) optimization toolbox, you need to download the [cvx](http://cvxr.com/cvx/download/) toolbox and then the extraction path is set to "C: \ Software \ cvx". After downloading the project, first set the "icon" folder to the MATLAB search path. Then just double-click "SNMFS_GUI.mlapp" to start running SMAFS

### GUI Tutorial

Note: If it is the first time to cluster scRNA-seq data, it is recommended to select "Coarse". If you want to further subdivide the cell population, we recommend choosing "Fine". If using SMAFS to clustering 10x data, one should select three files matrix.mtx , barcodes.tsv and genes.tsv at the same time, and upload them to SMAFS.

<div align=center>
<img src= './images/SMAFS_app.png'  width="95%" height="50%"  />
</div>

### Command-Line Tutorial

We provided a "demo.m" file to introduce how to use SMAFS. the detail is in the following.
```matlab
    %read data
    path.FileName = 'bigclass_5.csv';path.Path = 'E:/sc-RNA_Clustering/SMAFS_app/sub_data/';
    [geneName, cellName,A_source]= read_data(path);

    % normalization
    app.A_source = A_source;
    app.geneName = geneName;
    app.cellName = cellName;
    app.Upperbound = 8.5;%Upper bounds for filtering genes
    app.Lowerbound = 0.1;%lower bounds for filtering genes
    app = normalize(app);

    %calculate the initial value of W
    app = calculate_initial_W(app);

    %SNMF dimension reduction
    app.Maxiter = 500;
    app.tolerance = 1e-7;
    method = 'Corase';%Used to determine whether to use 'Fine' or 'Corase' method for dimensionality reduction
    MAGICsmoothing = 0;%1 means use MAGIC smooth after dimension reduction, and vice versa
    app = SNMF(app,method, MAGICsmoothing);

    %calculate the simmilarity matrix
    app = caulecte_similarMat(app);

    %Spectral Clustering base on L0
    Lambda = 0.012;%A parameter that controls the number of clusters. If lambda > 1, you will get floor (lambda) clusters. If ilambda < 1, the number of categories will be obtained through the optimization algorithm.
    app.Minimumgroupsize = 6;%Minimum number of cells in a cluster
    app = Clustering(app,Lambda);
    
    %extract cell label
    cell_label = app.celltype_label;
```
### Help
If you have any questions or require assistance using SMAFS, please contact us: hyl2016@mail.ustc.edu.cn .