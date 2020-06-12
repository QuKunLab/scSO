function app = SNMF(app,method, MAGICsmoothing,alpha,beta)
    %method:'Corase':
    %               'Fine':
    %MAGICsmoothing= 1,

    if strcmpi(method,'Corase')
        app.Fine =0;
        app.Corase = 1;
    end
    
     if strcmpi(method,'Fine')
        app.Fine =1;
        app.Corase = 0;
     end
    if  nargin <4
         if strcmpi(method,'Corase')
                app.Alpha = 0.05;
                app.Beta = 0.085;
         end
    
         if strcmpi(method,'Fine')
                app.Alpha = 0.01;
                app.Beta = 0.017;
         end
    else
            app.Alpha = alpha;
            app.Beta = beta;
    end
    app.W = [];
    app.H = []; 
    [app.W , app.H ,error ,state] = SNMF_saveMemery_cpu(app);            

    if MAGICsmoothing == 1
        app.H = MAGIC_PCA(app)';
    end

  sprintf('The relative error of the SNMF is : %f ' ,error);         
  '**It is Finished SNMF.**'
    
    var_H = std(app.H , 0 , 2);        
    bound = 0.4*mean(var_H);
    app.Gaussiancorenumbre = sum(var_H >= bound) + 1;
end



function [W , H ,error_new,state] = SNMF_saveMemery_cpu(app)
            %This program is used to implement sparse matrix non-negative decomposition algorithm
            [m , n] = size(app.A_norm);
           
            if app.Fine == 1
                W = app.initial_W(:,2:app.Rank + 1);%+ 0.000001 * rand(m ,k);
            end
            
            if app.Corase == 1
                W = app.initial_W(:,1:app.Rank);                
            end
            H =rand(app.Rank ,n);
            A_F = sum(sum(app.A_norm.^2));
            error_new = sum(sum(H.*((W'*W)*H))) - 2*sum(sum((H.*(W'*app.A_norm))));
            error_old = 1 + error_new/A_F;
            error_min = error_old;
            error_new = error_old;
            
            alph = app.Alpha;
            beta = app.Beta;
            
            gama = 1 - alph - beta ;
            
            iter = 1;
            count = 0;
            W_old = W;
            H_old = H;
            count_rise = 0;%calculate the number of raising steps;
            state = 0;
            bound_error =app.tolerance;
            ones_1_k = ones(1,app.Rank);
            A_new_H =[gama*app.A_norm;zeros(1,n)];  
            A_new_W = [gama* (app.A_norm') ; zeros(1,m)];
            sprintf('Iter£¬ error£¬ difference, count, count_rise')
            while count < 5 && iter<=app.Maxiter
                    W_new = [gama*W ; beta * ones_1_k];        
                    C_mat = (W_new')*W_new;
                    B_mat = (W_new')*A_new_H;
            
                    tol = 10*norm(W_new,1)*length(W_new)*eps;
                    H = SNMF_inner_mex(C_mat,B_mat,tol);
             
                    H_new  = [gama*(H');alph * ones_1_k];      
                    C_mat = (H_new')*H_new;
                    B_mat = (H_new')*A_new_W;
            
                    tol = 10*norm(H_new,1)*length(H_new)*eps;
                    W = SNMF_inner_mex(C_mat,B_mat,tol)';    
                    error_new = sum(sum(H.*((W'*W)*H))) - 2*sum(sum((H.*(W'*app.A_norm))));
                    error_new = 1 + error_new/A_F;
                    error = abs(error_new - error_old);
                    if error_new > (error_old + bound_error)
                        count_rise = count_rise + 1;
                    else
                        count_rise = 0;
                    end
                
                    if error_new < error_min
                         error_min = error_new;
                         W_old = W;
                         H_old = H;
                    end
                    if count_rise>3
                        'The rank selection of the matrix is larger and trying a smaller rank.' 
                        error = error_new;
                        W = W_old;
                        H = H_old;
                        return;
                    end
                     
                    if (error<= bound_error)
                        count = count +1;                        
                    else
                        count =0;
                    end
                    
                    sprintf('Iter %d£¬ %f£¬ %0.5g, %d, %d', iter, error_new, abs(error_new - error_old),count,count_rise)
                    error_old = error_new;
            
                    iter = iter+1;
                     pause(0.00001);
                     if iter == app.Maxiter
                        return;
                     end
            end
            W = W_old;
            H = H_old;
            state = 1;
end