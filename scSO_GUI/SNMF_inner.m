%#coder.gpu.kernelfun
function H = SNMF_inner(C_mat , B_mat ,tol)
        n = size(B_mat,2);
        H = zeros(size(B_mat));
        parfor i = 1 : double(n)                            
                H(:,i) = lsqnonneg_our( C_mat,B_mat(:,i) ,tol);
        end
end