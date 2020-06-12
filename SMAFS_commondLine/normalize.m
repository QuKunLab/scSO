function app = normalize(app)
    non_rate = sum(app.A_source>0,2)/size(app.A_source,2);           
    app.A_norm = app.A_source(non_rate>0.00&non_rate<1,:);           
    app.genes_used_for_clustering = app.geneName(non_rate>0.00&non_rate<1,1);

    mean_count = mean(app.A_norm,2);    
    mean_value = mean(mean_count);   
    low_bound = (app.Lowerbound)*mean_value;
    up_bound = app.Upperbound*mean_value;           
    index = mean_count < up_bound & mean_count >low_bound; 
    app.A_norm = app.A_norm(index,:);%extract the feature of suitable mean value  
    app.genes_used_for_clustering = app.genes_used_for_clustering(index,:);   

    sprintf(' After gene filter, the dataset contains %d cells and %d genes.',size(app.A_norm,2),size(app.A_norm,1))

    app.A_norm(isnan(app.A_norm)) = 0;
    app.A_norm = app.A_norm .*  repmat(10000./(sum(app.A_norm) + 0.000000001),size(app.A_norm,1) ,1); 

    app.A_norm = log10(app.A_norm + 1);       
end

