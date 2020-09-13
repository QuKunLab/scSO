function [value_ARI,value_NMI,value_P,value_E ] = test_resoult(class_)
%This function is used to calculate the various indicators for evaluating the classification results.
[n , ~] = size(class_);%number of cell
idx_class_1 = unique(class_(:,1),'stable');
idx_class_2 = unique(class_(:,2),'stable');

num_cell_type_1 = length(idx_class_1);
num_cell_type_2 = length(idx_class_2);
matrix_n = zeros(num_cell_type_1 , num_cell_type_2);
for i = 1 : num_cell_type_1
    A_set = find(class_(:,1)==idx_class_1(i));
    for j = 1 : num_cell_type_2
        B_set = find(class_(:,2)==idx_class_2(j));
        matrix_n(i,j) = length(intersect(A_set,B_set));
    end 
end
value_P = 0.0;
value_E = 0.0;
for i = 1 : num_cell_type_1
    for j = 1 : num_cell_type_2
        if matrix_n(i , j)~=0
            value_E = value_E + matrix_n(i , j)* log2(matrix_n(i , j) / sum(matrix_n(: , j)));
        end
    end 
    value_P = value_P + max(matrix_n(i , :));
end
value_P = value_P/n;
value_E = value_E/(n*log2(num_cell_type_1));

part_1 = 0.0;
part_1_1 = 0.0;
part_1_2 = 0.0;
for i = 1 : num_cell_type_1
    if sum(matrix_n(i , :))>=2
        part_1_1 = part_1_1 + nchoosek(sum(matrix_n(i , :)) , 2);
    end
    for j = 1 : num_cell_type_2
        if i == 1
            if sum(matrix_n(:, j))>=2
                part_1_2 = part_1_2 + nchoosek(sum(matrix_n(: ,j)) , 2);
            end
        end
        if matrix_n(i , j)>=2
            part_1 = part_1 + nchoosek(matrix_n(i , j),2);
        end
    end
end
part_1 = part_1 - part_1_1 *part_1_2 / (nchoosek(n , 2));
part_2 = 0.5 * (part_1_1 + part_1_2) - part_1_1 *part_1_2 / (nchoosek(n , 2));
value_ARI =part_1 / part_2;
value_NMI = Cal_NMI(class_(:,1), class_(:,2)) ;


end

function score = Cal_NMI(true_labels, cluster_labels)

true_labels = double(true_labels);
cluster_labels = double(cluster_labels);
true_labels = true_labels-min(true_labels(:))+1;
cluster_labels = cluster_labels-min(cluster_labels(:))+1;


if size(true_labels,2)>size(true_labels,1)
    true_labels = true_labels';
end

if size(cluster_labels,2)>size(cluster_labels,1)
    cluster_labels = cluster_labels';
end

n = length(true_labels);
cat = spconvert([(1:n)' true_labels ones(n,1)]);
cls = spconvert([(1:n)' cluster_labels ones(n,1)]);
cls = cls';
cmat = full(cls * cat);

n_i = sum(cmat, 1); % Total number of data for each true label (CAT), n_i
n_j = sum(cmat, 2); % Total number of data for each cluster label (CLS), n_j

% Calculate n*n_ij / n_i*n_j
[row, col] = size(cmat);
product = repmat(n_i, [row, 1]) .* repmat(n_j, [1, col]);
index = find(product > 0);
n = sum(cmat(:));
product(index) = (n*cmat(index)) ./ product(index);
% Sum up n_ij*log()
index = find(product > 0);
product(index) = log(product(index));
product = cmat .* product;
score = sum(product(:));
% Divide by sqrt( sum(n_i*log(n_i/n)) * sum(n_j*log(n_j/n)) )
index = find(n_i > 0);
n_i(index) = n_i(index) .* log(n_i(index)/n);
index = find(n_j > 0);
n_j(index) = n_j(index) .* log(n_j(index)/n);
denominator = sqrt(sum(n_i) * sum(n_j));

% Check if the denominator is zero
if denominator == 0
  score = 0;
else
  score = score / denominator;
end
end

