function app = calculate_initial_W(app)
   '**It is calculating initial W and estimating the rank.**'                        
    app.initial_W = [];
    app.A_norm = sparse(app.A_norm);
    [V,D,~] = svds(gpuArray(app.A_norm),100);
    V = gather(V);
    D = gather(D);
    app.initial_W =V(:,1:min(50,size(V,2)));
    latent_new = diag(D);           
    %             app.eigen_vector = latent_new;
    latent_new = latent_new(latent_new > 0);            

    initial_W_P = max(app.initial_W, 0);
    initial_W_N = max(-app.initial_W, 0);
    n_P_big = sum(initial_W_P.* initial_W_P);
    n_N_big = sum(initial_W_N.* initial_W_N);

    if sum(n_P_big >=n_N_big)>0
        app.initial_W(:,n_P_big >=n_N_big) = initial_W_P(:,n_P_big >=n_N_big);
    end
    if sum(n_P_big <n_N_big)>0
        app.initial_W(:,n_P_big <n_N_big) = initial_W_N(:,n_P_big <n_N_big);
    end
    n_svd = min(50,length(latent_new)-1);            
    t = 1:n_svd;
    latent_new = latent_new/latent_new(1);

    latent_new_diff = latent_new(1:end-1)./ latent_new(2:end)-1;

    bound = 0.085;
    
    figure
    hold on
    plot(t ,  latent_new_diff(1:n_svd),'xr');             

    for i =1:length(latent_new_diff) -10  
        app.Rank = i; 
        if (latent_new_diff(i) >= bound) &&all(latent_new_diff(i+1 : i+10)<bound )
            break;
        end
    end 

    plot(t , repmat(bound , 1 ,n_svd),'k');
    hold off
    title('The distribution of the relative change of eigenvalues')
    app.Rank = max(app.Rank , 3);
    app.Gaussiancorenumbre = app.Rank + 1;

sprintf('**It has finished calculating initial W and estimating rank.**\n The estimated rank is %d.',app.Rank)

end

