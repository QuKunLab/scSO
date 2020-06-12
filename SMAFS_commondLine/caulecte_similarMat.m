function app = caulecte_similarMat(app)
            'It is calculating similar matrix'
            beginTime = tic;
            try
                app.similar_matrix = [];   
                app.P_mat = [];
                app.H_new = [];
                app = Gaussion_classify(app); 
                app.similar_matrix = app.P_mat*app.P_mat';                 
                sprintf('Total %f features are selected' , size(app.H_new,1))                        
            catch ErrorInfo %caught error is an MException object
                disp(ErrorInfo);
                'gaussian mixed failed'             
                return;
            end
            sprintf('Calculating similar matrix costs %f s',toc(beginTime));           
            
            sum_row = sum(app.similar_matrix,2);
            Lapace_similar_Matrix = diag(sum_row)-app.similar_matrix;
            La_new = Lapace_similar_Matrix.*(abs(Lapace_similar_Matrix) > 0.0001);
            clear Lapace_similar_Matrix;
            
            'It is svd';
                      
             beginTime = tic;           
             La_2 = max(sum_row)*eye(size(La_new,1)) - La_new;
             [singleValue_U,eigen_SM]  = svds(gpuArray(La_2),45);
             eigen_SM= gather(eigen_SM);
             singleValue_U = gather(singleValue_U); 
             eigenValue_V=  max(sum_row) - diag(eigen_SM);
             app.eigenVector_SM = singleValue_U;
             clear singleValue_U 

             
             app.eigenVector_SM = app.eigenVector_SM(:,~isnan(eigenValue_V));
             eigenValue_V  = eigenValue_V (~isnan(eigenValue_V));
             [eigenValue_V ,index] = sort(abs(eigenValue_V),'ascend');             
             eigenValue_V = eigenValue_V/max(sum_row);
             app.SingleValue = eigenValue_V;
             app.eigenVector_SM = app.eigenVector_SM(:,index);            
        
            value_P = max(app.eigenVector_SM, 1e-16);
            value_N = max(-app.eigenVector_SM, 1e-16);            
            n_P_big = sum(value_P.* value_P);
            n_N_big = sum(value_N.*value_N);
            if sum(n_P_big <n_N_big)>0 
                app.eigenVector_SM(:,n_P_big <n_N_big) = -app.eigenVector_SM(:,n_P_big <n_N_big);
            end
            sprintf('svd costs %f s',toc(beginTime))
            
             app.Fielder_vector = [];
             t = 1:size(app.eigenVector_SM,1);
             figure('Name','Eigenvectors of similar matrix');
             for i = 1 :min(9 , size(app.eigenVector_SM ,2))
                subplot(3,3,i)
                plot(t,sort(app.eigenVector_SM(:,i)),'.k','MarkerSize',10)
                title(['\lambda',num2str(i),' : ',num2str(eigenValue_V(i))]);
             end

            addpath('C:\Software\cvx\builtins');
            addpath('C:\Software\cvx\commands');
            addpath('C:\Software\cvx\functions');
            addpath('C:\Software\cvx\lib');  
            addpath('C:\Software\cvx\structures');
            addpath('C:\Software\cvx\functions\vec_');
            addpath('C:\Software\cvx');
%                 x_ing = 1- eigen_SM_vector(1:49);
            n_eigen =min( 30 , size(app.eigenVector_SM ,2));
            x_ing  = eigenValue_V(1:n_eigen);
            x_old_ = x_ing;
            n_old  = length(x_old_);
            index_i = [1:n_old-1,1:n_old-1];
            index_j = [1:n_old-1,2:n_old];
            value_A = [-ones(1,n_old-1),ones(1,n_old-1)];
            A_diff = sparse(index_i , index_j , value_A);
            clear  index_i index_j value_A;
            error = 1;
            while error > 0.0001       
                y = A_diff * x_old_;
                wight = 1./(y+0.00000001);
                wight = diag(wight);
                A_diff_wight = wight*A_diff;
                x_new = zero_norm_inner(x_old_,A_diff_wight,0.010);%NMF
                error = norm(x_new - x_old_ , 2);
                x_old_ = x_new;                        
            end
            rmpath('C:\Software\cvx\builtins');
            rmpath('C:\Software\cvx\commands');
            rmpath('C:\Software\cvx\functions');
            rmpath('C:\Software\cvx\lib');                
            rmpath('C:\Software\cvx\structures');
            rmpath('C:\Software\cvx\functions\vec_');
            rmpath('C:\Software\cvx');  
            
            figure('Name','Preliminary estimated Fiedler vector')
            hold on
            plot(1:n_eigen,x_ing,'.')
            x_old_ = roundn(x_old_ ,-3);%Rounding off
            plot(1:n_eigen,x_old_)
            hold off
            
            app.singleValue_bound = x_old_(2);
            n_eigenVector = sum(x_ing < (1.5*app.singleValue_bound));%NMF
            F_v = mean(app.eigenVector_SM(:,1:n_eigenVector),2);
            if max(F_v) - min(F_v) < 1e-4
                    app.eigenVector_SM(:,1) = 0;
                    n_eigenVector = n_eigenVector +1;                    
            end

            app.SpectralNumber = min(max(n_eigenVector,2),app.Rank-1);
            app.Fielder_vector = mean(app.eigenVector_SM(:,1:app.SpectralNumber),2); 
            
            [x_old, ~ ] = sort(app.Fielder_vector);            
            y = x_old(2:end) - x_old(1:end-1);
            block_size=max(floor( 0.025*length(x_old)),3);
            for j  =block_size : length(y)-block_size
                if y(j,1) ~= max(y(j-block_size+1:j+block_size,1))
                    y(j) =0;
                end
           end
                 
            title( 'the vector of the smallest no zero eigenvalue' ,'Color', 'm');
            xlabel( 'idx feature');
            ylabel( 'Value');            
            t = 1:length(x_old);
            plot(t ,x_old,'.k');
end

function app = Gaussion_classify(app)
     app.H_new = app.H(var(app.H')>0,:);
     app.GaussiancorenumbreEditField.Value = size(app.H_new,1)+1;
     ClusterNumber = app.GaussiancorenumbreEditField.Value;
     app.P_mat=[];    
     rng(1,'twister')
     Sigma = {'diagonal','full'};
     SharedCovariance = {true,false};
     options = statset('Display','off','MaxIter',1000);
     D = app.H_new - mean(app.H_new,2);
     Maha_D = sum(D.*(diag(var(app.H_new,1,2).^(-1))*D));
     [~,index_cell] = sort(Maha_D,'ascend');
     H_sorted = app.H_new(:,index_cell);
     gmfit = fitgmdist(H_sorted',ClusterNumber,'CovarianceType',Sigma{1},...
                              'SharedCovariance',SharedCovariance{1},'Options',options);
     [~,~,P_posterior] = cluster(gmfit,app.H_new');
     app.P_mat = [app.P_mat,P_posterior];

    for i = 1:49
        rng(i+1,'twister'); % For reproducibility
        gmfit = fitgmdist(H_sorted',ClusterNumber ,'Start','plus','CovarianceType','diagonal',...
                                   'SharedCovariance',true,'Options',options);
        [~,~,P_posterior] = cluster( gmfit,app.H_new');  
        app.P_mat = [app.P_mat,P_posterior];
    end
    app.P_mat = app.P_mat * sqrt(0.02); 
end

        