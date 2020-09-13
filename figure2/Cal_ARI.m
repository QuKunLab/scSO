function value_ARI = Cal_ARI(class_)
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
up_part = part_1 - part_1_1 *part_1_2 / (nchoosek(n , 2));
low_part = 0.5 * (part_1_1 + part_1_2) - part_1_1 *part_1_2 / (nchoosek(n , 2));
value_ARI =up_part / low_part;
end