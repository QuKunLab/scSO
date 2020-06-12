function [geneName, cellName,A_source]= read_data(path)
            %This function is used to read the gene expression matrix
            %Note: this Function only used to read .csv, .tsv and 10x
            %for .csv and .tsv file, you need to set path.FileName and path.Path. For example, path.FileName= 'a.csv';path.Path = './';
            %for 10x format, you only need to set path.Path that include 'matrix.mtx', 'barcodes.tsv' and 'genes.tsv'. For example, path.Path = './';
            if ~isempty(path.FileName)
                %% read csv file        
                if strcmpi(path.FileName(end-2:end),'csv')
                    Data_no_change = importdata([path.Path,path.FileName],',',1);
                    A_source_head = Data_no_change.textdata;
                    A_source =Data_no_change.data;
                    
                    % Summing the expression of genes of the same name
                    [geneName,index_ing,b] = unique(A_source_head(2:end,1),'stable');
                    unique_b= unique(b);
                    A_source = A_source(index_ing,:);
                    for i = 1 : length(unique_b)
                        if sum(b == unique_b(i))>1
                            A_source(i,:) = sum(Data_no_change.data(b == unique_b(i),:));
                        end
                    end
                    cellName = A_source_head(1,2:end);
                    return;
                end
                
                %% read tsv file
                if  strcmpi(path.FileName(end-2:end),'tsv')
                    Data_no_change = importdata([path.Path,path.FileName],'\t',1);
                    A_source_head = Data_no_change.textdata;
                    A_source =Data_no_change.data;
                
                    % Summing the expression of genes of the same name
                    [geneName,index_ing,b] = unique(A_source_head(2:end,1),'stable');
                    unique_b= unique(b);
                    A_source = A_source(index_ing,:);
                    for i = 1 : length(unique_b)
                        if sum(b == unique_b(i))>1
                            A_source(i,:) = sum(Data_no_change.data(b == unique_b(i),:));
                        end
                    end
                    cellName = A_source_head(1,2:end);
                    return
                end
            else    
            
                %% read mtx file
                fid = fopen([path.Path,'barcodes.tsv']);
                cellName = textscan(fid, '%s','delimiter', '\t');
                fclose(fid);
                cellName = cellName{1};
                
                fileID = fopen([path.Path,'matrix.mtx']);
                C = textscan(fileID,'%n %n %n','CommentStyle','%');
                fclose(fileID);
                A_source = sparse(C{1}(2:end),C{2}(2:end),C{3}(2:end),C{1}(1),C{2}(1));
                clear C;
                
                fid = fopen([path.Path,'genes.tsv']);
                geneName = textscan(fid, '%s %s','delimiter', '\t');
                geneName = geneName{1,1};
                fclose(fid);
                
                % Sum the expression of genes of the same name
                [geneName,index_ing,b] = unique(geneName(2:end,1),'stable');
                unique_b= unique(b);
                A_source = A_source(index_ing,:);
                for i = 1 : length(unique_b)
                    if sum(b == unique_b(i))>1
                        A_source(i,:) = sum(A_source(b == unique_b(i),:));
                    end
                end
             end
        end