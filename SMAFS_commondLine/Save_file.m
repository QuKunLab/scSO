function Save_file(app,output_path)
          
            output_path = [output_path,'/'];
            [CelltypeName, ~,ref_Label] = unique(app.cellName,'stable');
            
            'It is saving results'
            mkdir([output_path,'label'])
            mkdir([output_path,'sub_data'])
            
            csvwrite([output_path,'label/ref.csv'],ref_Label);
            csvwrite([output_path,'label/0_label.csv'],app.celltype_label);
            fileID = fopen([output_path, '/sub_data/gene_name.csv'],'w');
            formatSpec1 = '%s\n'; 
            fprintf(fileID,formatSpec1, app.geneName{:});
            fclose(fileID);
%           gene_name = app.A_source_head(2:end,1);
            A_save = [(1:1:size(app.A_source,1))',app.A_source];
            for i = 1:max(app.celltype_label)
                    fileID = fopen([output_path, '/sub_data/bigclass_',num2str(i),'.csv'],'w');
                    index_ing = find(app.celltype_label==i);
                    formatSpec1 = [repmat('%s,' ,1,length(index_ing)),'%s\n'];                 
                    formatSpec2 = ['gene%d,',repmat('%f,' ,1,length(index_ing)-1),'%f\n'];   
                    colum_name = ['geneName',app.cellName(1,index_ing)];
                    fprintf(fileID,formatSpec1,colum_name{1,1:length(index_ing)+1});
                    fprintf(fileID,formatSpec2,A_save(:, [1 ; 1+index_ing])');
                    fclose(fileID);                    
            end
            sprintf('The results have saved in %s.',output_path)
        end