function x_new = zero_norm_inner(x_old,A_diff,Lambda)
        %This  function is  used to L0 optimization
        n_old = length(x_old);
        abs_x_old = abs(x_old);
        n_old_1 = n_old -1;
        cvx_begin quiet
            variable x_new(n_old)
            variable y(n_old_1)
            minimize((1 - Lambda)*norm(x_old - x_new , 2)+Lambda* norm(y , 1))
            subject to
            y == A_diff*x_new;
        cvx_end
end
        