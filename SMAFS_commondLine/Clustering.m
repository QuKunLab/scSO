function app = Clustering(app,Lambda,n_Spectral)
%Lambda: determint the number of cluster.
%n_Spectral: the number of eigen vectors used to calculate the Fielder vector.

    nn_eigenVector = app.SpectralNumber;
    if  nargin <3            
            app.Fielder_vector = mean(app.eigenVector_SM(:,1:nn_eigenVector),2); 
            [x_old,x_index1] = sort( app.Fielder_vector);            
            num_subClass = floor(Lambda);
            x_new = [];
            state_plus = 0;
            if num_subClass >= 1
                t = 1:length(x_old);
                x_old  = x_old  - mean(x_old );
                x_old = x_old /abs(min(x_old ))*2;
                y = x_old(2:end) - x_old(1:end-1);
                block_size =max(floor( 0.025*length(y)),5);
                for j  = block_size : length(y)-block_size
                    if y(j,1) ~= max(y(j-block_size+1:j+block_size,1))
                        y(j) =0;
                    end
                end
                [~ , index_y]  = sort(y(block_size:end-block_size),'descend');
                index_X_old = index_y(1:num_subClass-1)+ block_size-1;
                subBound = 0.5 * (x_old(index_X_old) + x_old(index_X_old+1));
                subBound = sort(subBound);
                x_new(x_old<subBound(1)) = mean(x_old(x_old<subBound(1)));
                for i = 2: num_subClass -1
                    x_new(x_old>=subBound(i-1)&x_old<subBound(i)) =...
                    mean(x_old(x_old>=subBound(i-1)&x_old<subBound(i)));                    
                end
                x_new(x_old>=subBound(num_subClass -1)) = mean(x_old(x_old>=subBound(num_subClass -1)));
                x_index = x_index1;

            else 
                addpath('C:\Software\cvx\builtins');
                addpath('C:\Software\cvx\commands');
                addpath('C:\Software\cvx\functions');
                addpath('C:\Software\cvx\lib');  
                addpath('C:\Software\cvx\structures');
                addpath('C:\Software\cvx\functions\vec_');
                addpath('C:\Software\cvx');
                x_old  = x_old  - mean(x_old);
                x_old = x_old /abs(min(x_old))*2;
                x_new1 = calcu_FielderVector(app,x_old,Lambda);
                x_new1 = roundn(x_new1 ,-3);%Rounding off
                x_new = x_new1;
                x_index = x_index1;
                %判断nn_eigenVector+1个最小的奇异值向量的平均作为fielder向量是否可以得到更多的类
    %%{
              if (app.SingleValue(app.SpectralNumber + 1)-app.SingleValue(app.SpectralNumber)  )< 0.03...%SNMF
                         &&app.SpectralNumber < max(app.Rank-1,3)                     
    %                  if (app.SingleValue(nn_eigenVector + 1)-app.singleValue_bound )< 0.017...%NMF
    %                         &&nn_eigenVector <= max(app.ReducedDimensionEditField.Value-2,3)
                    x_old = mean(app.eigenVector_SM(:,1:app.SpectralNumber+1),2);
                    x_old  = x_old  - mean(x_old);
                    x_old = x_old /abs(min(x_old ))*2;
                    [x_old,x_index2] = sort(x_old);
                    x_new2 = calcu_FielderVector(app,x_old,Lambda);
                    x_new2= roundn(x_new2 ,-3);%Rounding off     
                    if length(unique(x_new2)) > length(unique(x_new1))              
                        state_plus = 1;
                        x_new = x_new2;  
                        x_index = x_index2;
                        nn_eigenVector = app.SpectralNumber + 1;
                        pause(0.001)
                    end
              end
    %}

    % %{
             x_old = mean(app.eigenVector_SM(:,1:app.SpectralNumber-1),2);

             if (app.SpectralNumber -1)>2&&(max(x_old) - min(x_old) > 1e-5)

                x_old  = x_old  - mean(x_old );
                x_old = x_old /abs(min(x_old ))*2;
                [x_old,x_index2] = sort(x_old);
                x_new2 = calcu_FielderVector(app,x_old,Lambda);
                x_new2= roundn(x_new2 ,-3);%Rounding off 
                 if length(unique(x_new2)) >= length(unique(x_new))                             
                        x_new = x_new2;  
                        x_index = x_index2;
                        nn_eigenVector = app.SpectralNumber -1;
                        pause(0.001)
                 else
                     t = 1 : length(x_new);
                     x_old = mean(app.eigenVector_SM(:,1:nn_eigenVector),2);
                     x_old  = x_old  - mean(x_old );
                     x_old = x_old /abs(min(x_old ))*2;
                     [x_old,x_index] = sort(x_old);      
                 end 
             end
     %}

                rmpath('C:\Software\cvx\builtins');
                rmpath('C:\Software\cvx\commands');
                rmpath('C:\Software\cvx\functions');
                rmpath('C:\Software\cvx\lib');                
                rmpath('C:\Software\cvx\structures');
                rmpath('C:\Software\cvx\functions\vec_');
                rmpath('C:\Software\cvx');  
            end        
             x_new_set = union(x_new ,[]);
             app.celltype_label = zeros(length(app.cellName),1);
             for i=1:length(x_index)
                 app.celltype_label(x_index(i),1) = find(x_new(i)==x_new_set);
             end
    else
            app.Fielder_vector = mean(app.eigenVector_SM(:,1:n_Spectral),2); 
            [x_old,x_index] = sort( app.Fielder_vector);            
            num_subClass = floor(Lambda);
            x_new = [];
            if num_subClass >= 1
                t = 1:length(x_old);
                x_old  = x_old  - mean(x_old );
                x_old = x_old /abs(min(x_old ))*2;
                y = x_old(2:end) - x_old(1:end-1);
                block_size =max(floor( 0.025*length(y)),5);
                for j  = block_size : length(y)-block_size
                    if y(j,1) ~= max(y(j-block_size+1:j+block_size,1))
                        y(j) =0;
                    end
                end
                [~ , index_y]  = sort(y(block_size:end-block_size),'descend');
                index_X_old = index_y(1:num_subClass-1)+ block_size-1;
                subBound = 0.5 * (x_old(index_X_old) + x_old(index_X_old+1));
                subBound = sort(subBound);
                x_new(x_old<subBound(1)) = mean(x_old(x_old<subBound(1)));
                for i = 2: num_subClass -1
                    x_new(x_old>=subBound(i-1)&x_old<subBound(i)) =...
                    mean(x_old(x_old>=subBound(i-1)&x_old<subBound(i)));                    
                end
                x_new(x_old>=subBound(num_subClass -1)) = mean(x_old(x_old>=subBound(num_subClass -1)));

            else 
                addpath('C:\Software\cvx\builtins');
                addpath('C:\Software\cvx\commands');
                addpath('C:\Software\cvx\functions');
                addpath('C:\Software\cvx\lib');  
                addpath('C:\Software\cvx\structures');
                addpath('C:\Software\cvx\functions\vec_');
                addpath('C:\Software\cvx');
                x_old  = x_old  - mean(x_old);
                x_old = x_old /abs(min(x_old))*2;
                x_new = calcu_FielderVector(app,x_old,Lambda);
                x_new = roundn(x_new ,-3);%Rounding off
                
                rmpath('C:\Software\cvx\builtins');
                rmpath('C:\Software\cvx\commands');
                rmpath('C:\Software\cvx\functions');
                rmpath('C:\Software\cvx\lib');                
                rmpath('C:\Software\cvx\structures');
                rmpath('C:\Software\cvx\functions\vec_');
                rmpath('C:\Software\cvx');  
            end        
             x_new_set = union(x_new ,[]);
             app.celltype_label = zeros(length(app.cellName),1);
             for i=1:length(x_index)
                 app.celltype_label(x_index(i),1) = find(x_new(i)==x_new_set);
             end
    end
end

function x_new2 = calcu_FielderVector(app,V_s,lambda)
%This function is used to calculate L0 Approach of V_s
% V_s: is a vector will be approached
% lambda: (scalar)  L0 parameter
% figure_handle: handle of figure 
    n_old = length(V_s);
    index_i = [1:n_old-1,1:n_old-1];
    index_j = [1:n_old-1,2:n_old];
    value_A = [-ones(1,n_old-1),ones(1,n_old-1)];
    A_diff = sparse(index_i , index_j , value_A);
    clear index_i index_j value_A         
     x_old_ = V_s;
    error = 1;
    while error > 0.0001       
        y = A_diff * x_old_;
        block_size = app.Minimumgroupsize;
        if  block_size == 0                
            block_size =max(ceil(0.01*length(y)),6);
        end
        y(1:block_size,1) = 0;
        y(end - block_size+1:end,1) = 0;
        for j  = block_size : length(y)-block_size
            if y(j,1) ~= max(y(j-block_size+1:j+block_size,1))
                y(j) =0.0;
            end
        end
        wight = 1./(y.^1.4+0.0000000001);
        wight = diag(wight);
        A_diff_wight = wight*A_diff;      
        x_new2 = zero_norm_inner(x_old_,A_diff_wight,lambda);
        error = norm(x_new2 - x_old_ , 2);
        x_old_ = x_new2;
    end
    if app.Corase == 1
        x_new2 = dealWrongCluster(x_new2, V_s);
    end
end

function y = dealWrongCluster(x,x_old)
            x_old = roundn(x_old ,-4);
            x= roundn(x ,-3);%Rounding off    
            y = x;
            x_unique = unique(x,'stable');
                %%{
            for i = 1:length(x_unique)
                    if length(unique(x_old(x== x_unique(i),1))) <= 3 
                             continue
                    end      
                    x_ing = x_old(x== x_unique(i),1);
                    length_x = sum(x== x_unique(i));
                    index_1 = max(floor(0.25*length_x),2);
                    index_2 = max(floor(0.75*length_x),4);
                    x_ing = x_ing(index_1:index_2) - x_ing(index_1);
                    sig = 2*mean(x_ing)/length(x_ing);
                    if abs(sig) > 0.012
                        if i == 1
                            y(x==  x_unique(i),1) = x_unique(i+1);
                            x_unique(i) = x_unique(i+1);
                            continue;
                        end
                        if i==length(x_unique)
                            y(x==  x_unique(i),1) = x_unique(i-1);
                            x_unique(i) = x_unique(i-1);
                            continue;
                        end
                        if x_unique(i) - x_unique(i-1) >x_unique(i+1) - x_unique(i)
                            y(x==  x_unique(i),1) = x_unique(i+1);
                            x_unique(i) = x_unique(i+1);
                        else
                            y(x==  x_unique(i),1) = x_unique(i-1);
                            x_unique(i) = x_unique(i-1);
                        end
            
                    end
                            
            end 
                    %}   
                    
                    
                    
end 
