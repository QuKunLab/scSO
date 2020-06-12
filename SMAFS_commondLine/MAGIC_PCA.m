 function pc_imputed = MAGIC_PCA(app, varargin)
            % run_magic  Run MAGIC for imputing and denoising of single-cell data
            %   [pc_imputed, U, pc] = run_magic(data, varargin) runs MAGIC on data (rows:
            %   cells, columns: genes) with default parameter settings and returns the
            %   imputed data in a compressed format.
            %
            %   The compressed format consists of loadings (U) and imputed principal
            %   components (pc_imputed). To obtain gene values form this compressedf format
            %   either run project_genes.m or manually project (pc_imputed * U') either all
            %   genes or a subset (pc_imputed * U(idx,:)'). Also returned are the original
            %   principal components (pc);
            %
            %   Since pc_imputed and U are both narrow matrices the imputed data can be
            %   stored in a memory efficient way, without having to store the dense
            %   matrix.
            %
            %   Supplied data can be a sparse matrix, in which case MAGIC will be more
            %   memory efficient.
            %
            %   [...] = phate(data, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to
            %   specify optional parameter name/value pairs that control further details
            %   of PHATE.  Parameters are:
            %
            %   'npca' - number of PCA components to do MAGIC on. Defaults to 100.
            %
            %   'k' - number of nearest neighbors for bandwidth of adaptive alpha
            %   decaying kernel or, when a=[], number of nearest neighbors of the knn
            %   graph. For the unweighted kernel we recommend k to be a bit larger,
            %   e.g. 10 or 15. Defaults to 10.
            %
            %   'a' - alpha of alpha decaying kernel. when a=[] knn (unweighted) kernel
            %   is used. Defaults to 15.
            %
            %   't' - number of diffusion steps. Defaults to [] wich autmatically picks
            %   the optimal t.
            %
            %   'distfun' - Distance function used to compute kernel. Defaults to
            %   'euclidean'.
            %
            %   'make_plot_opt_t' - Boolean flag for plotting the optimal t analysis.
            %   Defaults to true.
            k = 15;
            a = 15;
            t = [];
            distfun = 'euclidean';
            make_plot_opt_t = true;
            
            % get input parameters
            for i=1:length(varargin)
                % k for knn adaptive sigma
                if(strcmp(varargin{i},'k'))
                   k = lower(varargin{i+1});
                end
                % a (alpha) for alpha decaying kernel
                if(strcmp(varargin{i},'a'))
                   a = lower(varargin{i+1});
                end
                % diffusion time
                if(strcmp(varargin{i},'t'))
                   t = lower(varargin{i+1});
                end
                % npca
                if(strcmp(varargin{i},'npca'))
                   npca = lower(varargin{i+1});
                end
                % make plot optimal t
                if(strcmp(varargin{i},'make_plot_opt_t'))
                   make_plot_opt_t = lower(varargin{i+1});
                end
            end
            
            % compute kernel
            disp 'computing kernel'
            K = compute_kernel(app, app.H', 'k', k, 'a', a, 'distfun', distfun);
            
            % row stochastic
            P = bsxfun(@rdivide, K, sum(K,2));
            
            % optimal t
            if isempty(t)
                disp 'imputing using optimal t'
                pc_imputed = compute_optimal_t(app, app.H', P, 'make_plot', make_plot_opt_t);
            else
                disp 'imputing using provided t'
                pc_imputed = app.H';
                for I=1:t
                    disp(['t = ' num2str(I)]);
                    pc_imputed = P * pc_imputed;
                end
            end
            disp 'done.'
end


function K = compute_kernel(app,data, varargin)
    % K = computer_alpha_kernel_sparse(data, varargin)
    % Computes sparse alpha-decay kernel
    % varargin: 
    %   'npca' (default = [], no PCA)
    %       Perform fast random PCA before computing distances
    %   'k' (default = 5)
    %       k for the knn distances for the locally adaptive bandwidth
    %   'a' (default = 10)
    %       The alpha exponent in the alpha-decaying kernel
    %   'distfun' (default = 'euclidean')
    %       Input distance function
    k = 5;
    a = 10;
    npca = [];
    distfun = 'euclidean';
    % get the input parameters
    if ~isempty(varargin)
        for j = 1:length(varargin)
            % k nearest neighbora
            if strcmp(varargin{j}, 'k')
                k = varargin{j+1};
            end
            % alpha
            if strcmp(varargin{j}, 'a')
                a = varargin{j+1};
            end
            % npca to project data
            if strcmp(varargin{j}, 'npca')
                npca = varargin{j+1};
            end
            % distfun
            if strcmp(varargin{j}, 'distfun')
                distfun = varargin{j+1};
            end
        end
    end

    th = 1e-4;

    k_knn = k * 20;

    bth=(-log(th))^(1/a);

    disp 'Computing alpha decay kernel:'

    N = size(data, 1); % number of cells

    if ~isempty(npca)
        disp '   PCA'
        data_centr = bsxfun(@minus, data, mean(data,1));
        [U,~,~] = randPCA(data_centr', npca); % fast random svd
        %[U,~,~] = svds(data', npca);
        data_pc = data_centr * U; % PCA project
    else
        data_pc = data;
    end

    disp(['Number of samples = ' num2str(N)])

    % Initial knn search and set the radius
    disp(['First iteration: k = ' num2str(k_knn)])
    [idx, kdist]=knnsearch(data_pc,data_pc,'k',k_knn,'Distance',distfun);
    epsilon=kdist(:,k+1);

    % Find the points that have large enough distance to be below the kernel
    % threshold
    below_thresh=kdist(:,end)>=bth*epsilon;

    idx_thresh=find(below_thresh);

    if ~isempty(idx_thresh) 
        K=exp(-(kdist(idx_thresh,:)./epsilon(idx_thresh)).^a);
        K(K<=th)=0;
        K=K(:);
        i = repmat(idx_thresh',1,size(idx,2));
        i = i(:);
        idx_temp=idx(idx_thresh,:);
        j = idx_temp(:);
    end

    disp(['Number of samples below the threshold from 1st iter: ' num2str(length(idx_thresh))])

    % Loop increasing k by factor of 20 until we cover 90% of the data
    while length(idx_thresh)<.9*N
        k_knn=min(20*k_knn,N);
        data_pc2=data_pc(~below_thresh,:);
        epsilon2=epsilon(~below_thresh);
        disp(['Next iteration: k= ' num2str(k_knn)])
        [idx2, kdist2]=knnsearch(data_pc,data_pc2,'k',k_knn,'Distance',distfun);

    %     Find the points that have large enough distance
        below_thresh2=kdist2(:,end)>=bth*epsilon2;
        idx_thresh2=find(below_thresh2);

        if ~isempty(idx_thresh2)
            K2=exp(-(kdist2(idx_thresh2,:)./epsilon2(idx_thresh2)).^a);
            K2(K2<=th)=0;
            idx_notthresh=find(~below_thresh);
            i2=repmat(idx_notthresh(idx_thresh2)',1,size(idx2,2));
            i2=i2(:);
            idx_temp=idx2(idx_thresh2,:);
            j2=idx_temp(:);

            i=[i; i2];
            j=[j; j2];
            K=[K; K2(:)];
    %         Add the newly thresholded points to the old ones
            below_thresh(idx_notthresh(idx_thresh2))=1;
            idx_thresh=find(below_thresh);
        end
        disp(['Number of samples below the threshold from the next iter: ' num2str(length(idx_thresh))])
    end

    % Radius search for the rest
    if length(idx_thresh)<N
        disp(['Using radius based search for the rest'])
        data_pc2=data_pc(~below_thresh,:);
        epsilon2=epsilon(~below_thresh);
        [idx2, kdist2]=rangesearch(data_pc,data_pc2,bth*max(epsilon2),'Distance',distfun);
        idx_notthresh=find(~below_thresh);
        for m=1:length(idx2)
            i=[i; idx_notthresh(m)*ones(length(idx2{m}),1)];
            j=[j; idx2{m}'];
            K2=exp(-(kdist2{m}./epsilon2(m)).^a);
            K2(K2<=th)=0;
            K=[K; K2(:)];
        end

    end

    % Build the kernel
    K = sparse(i, j, K);

    disp '   Symmetrize affinities'
    K = K + K';
    disp '   Done computing kernel'
end
function W = compute_operator(app,data, varargin)
    % W = compute_operator(data, varargin)
    %   computes diffusion operator W
    % varargin:
    %   'npca' (default = 20)
    %       perform fast random PCA before computing distances
    %   'ka' (default = 10)
    %       k of adaptive kernel
    %       0 for non-adaptive (standard gaussian) kernel with bandwidth
    %       epsilon
    %   'k' (default = 30)
    %       k of kNN graph
    %   'epsilon' (default = 1)
    %       kernel bandwith, if epsilon = 0 kernel will be uniform, i.e.
    %       unweighted kNN graph (ka will be ignored)

    % set up default parameters
    k = 12;
    ka = 4;
    npca = 100;
    epsilon = 1;
    lib_size_norm = true;
    log_transform = false;

    % get the input parameters
    if ~isempty(varargin)
        for j = 1:length(varargin)
            % k nearest neighbor
            if strcmp(varargin{j}, 'ka')
                ka = varargin{j+1};
            end
            % for knn-autotune
            if strcmp(varargin{j}, 'k')
                k = varargin{j+1};
            end
            % epsilon
            if strcmp(varargin{j}, 'epsilon')
                epsilon = varargin{j+1};
            end
            % npca to project data
            if strcmp(varargin{j}, 'npca')
                npca = varargin{j+1};
            end
            % library size normalization
            if strcmp(varargin{j}, 'lib_size_norm')
                lib_size_norm = varargin{j+1};
            end
            % log transform
            if strcmp(varargin{j}, 'log_transform')
                log_transform = varargin{j+1};
            end
        end
    end

    % library size normalization
    if lib_size_norm
        disp 'Library size normalization'
        libsize  = sum(data,2);
        data = bsxfun(@rdivide, data, libsize) * median(libsize);
    end

    % log transform
    if log_transform
        disp(['Log transform, with pseudo count ' num2str(pseudo_count)]);
        data = log(data + pseudo_count);
    end

    N = size(data, 1); % number of cells

    if ~isempty(npca)
        disp 'PCA'
        data_centr = bsxfun(@minus, data, mean(data,1));
        [U,~,~] = randPCA(data_centr', npca); % fast random svd
        %[U,~,~] = svds(data', npca);
        data_pc = data_centr * U; % PCA project
    else
        data_pc = data;
    end

    disp 'Computing distances'
    [idx, dist] = knnsearch(data_pc, data_pc, 'k', k);

    disp 'Adapting sigma'
    dist = bsxfun(@rdivide, dist, dist(:,ka));

    i = repmat((1:N)',1,size(idx,2));
    i = i(:);
    j = idx(:);
    s = dist(:);
    if epsilon > 0
        W = sparse(i, j, s);
    else
        W = sparse(i, j, ones(size(s))); % unweighted kNN graph
    end

    disp 'Symmetrize distances'
    W = W + W';

    if epsilon > 0
        disp 'Computing kernel'
        [i,j,s] = find(W);
        i = [i; (1:N)'];
        j = [j; (1:N)'];
        s = [s./(epsilon^2); zeros(N,1)];
        s = exp(-s);
        W = sparse(i,j,s);
    end

    disp 'Markov normalization'
    W = bsxfun(@rdivide, W, sum(W,2)); % Markov normalization

    disp 'done'
end
function [data_opt_t, t_opt]  = compute_optimal_t(app,data, DiffOp, varargin)

    t_max = 32;
    make_plot = true;
    th = 1e-3;
    data_opt_t = [];

    if ~isempty(varargin)
        for j = 1:length(varargin)
            if strcmp(varargin{j}, 't_max')
                t_max = varargin{j+1};
            end
            if strcmp(varargin{j}, 'make_plot')
                make_plot = varargin{j+1};
            end
            if strcmp(varargin{j}, 'th')
                th = varargin{j+1};
            end
        end
    end

    data_prev = data;
    if make_plot
        error_vec = nan(t_max,1);
        for I=1:t_max
            disp(['t = ' num2str(I)]);
            data_curr = DiffOp * data_prev;
            error_vec(I) = procrustes(data_prev, data_curr);
            if error_vec(I) < th && isempty(data_opt_t)
                data_opt_t = data_curr;
            end
            data_prev = data_curr;
        end
        t_opt = find(error_vec < th, 1, 'first');

        figure;
        hold all;
        plot(1:t_max, error_vec, '*-');
        plot(t_opt, error_vec(t_opt), 'or', 'markersize', 10);
        xlabel 't'
        ylabel 'error'
        axis tight
        ylim([0 ceil(max(error_vec)*10)/10]);
        plot(xlim, [th th], '--k');
        legend({'y' 'optimal t' ['y=' num2str(th)]});
        set(gca,'xtick',1:t_max);
        set(gca,'ytick',0:0.1:1);
    else
        for I=1:t_max
            disp(['t = ' num2str(I)]);
            data_curr = DiffOp * data_prev;
            error = procrustes(data_prev, data_curr);
            if error < th
                t_opt = I;
                data_opt_t = data_curr;
                break
            end
            data_prev = data_curr;
        end
    end

    disp(['optimal t = ' num2str(t_opt)]);
end