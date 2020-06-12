clc;
clear;
close all;
%read data
    path.FileName = 'deng.csv';path.Path = '../Data/';
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
    method = 'Fine';%Used to determine whether to use 'Fine' or 'Corase' method for dimensionality reduction
    MAGICsmoothing = 1;%1 means use MAGIC smooth after dimension reduction, and vice versa
    app = SNMF(app,method, MAGICsmoothing);

    %calculate the simmilarity matrix
    app = caulecte_similarMat(app);

    %Spectral Clustering base on L0
    Lambda = 0.012;%A parameter that controls the number of clusters. If lambda > 1, you will get floor (lambda) clusters. If ilambda < 1, the number of categories will be obtained through the optimization algorithm.
    app.Minimumgroupsize = 0;%Minimum number of cells in a cluster. If app.Minimumgroupsize=0£¬ SMAFS will use the default value.

%     app.SpectralNumber = 4;
    app = Clustering(app,Lambda);
    
    %extract cell label
    cell_label = app.celltype_label;
    
    class_(:,1) = cell_label;
    [~,a,class_(:,2)] = unique(cellName);
    sprintf("reference: %d  ,SMAFS: %d",length(a),length(unique(cell_label )))
    test_resoult(class_)
    
    %Store cell label in [output_path 'label'] folder, meanwhile, extract and store different clusters of A_sourse in [output_path 'sub_data']  
    output_path = './';
    Save_file(app,output_path)
    
    
    