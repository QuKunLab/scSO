data_set = {'biase','deng','goolam','klein','kolodziejczyk','pollen','rca','treutlein','usoskin','yan','zeisel','Baron'};
path = './res_OurPretreatment/';
j = 1 %change this number to chose dataset
ref = csvread(['./ref/', data_set{j}, '.csv']);
label = csvread([path,'Our/' ,data_set{j}, '.csv']);
our  = Cal_ARI([ref,label]);
A= data_set{j};
for i = 1:50
   label = readtable([path,'SC3/' ,[upper(A(1)),A(2:end)],'/res_',num2str(i), '.csv']); 
   label = label(:,2).Variables;
    if iscell(label) 
        label = string(label);
        label =str2double(label);
    end
   SC3(i)  = Cal_ARI([ref,label]);    
    label = csvread([path,'Scanpy/'  ,[upper(A(1)),A(2:end)],'/res_',num2str(i),'.csv']);
    Scanpy(i)  = Cal_ARI([ref,label(:,2)]);
    
    label = readtable([path,'Seurat/' ,A,'/res_',num2str(i),'.csv']);
    label = label(:,2).Variables;
    label = str2double(label);
    Seurat(i)  = Cal_ARI([ref,label]);
    
    label = readtable([path,'SIMLR/' ,[upper(A(1)),A(2:end)],'/res_',num2str(i), '.csv']);
    SIMLR(i)  = Cal_ARI([ref,label.x]);
end

A = [repmat(our,50,1),SC3',Scanpy',Seurat', SIMLR'];
figure
X = categorical({'Our','SC3','Scanpy','Seurat', 'SIMLR'});
X = reordercats(X,{'Our','SC3','Scanpy','Seurat', 'SIMLR'});  
plot_bar_plus_dot(gca,[1,2,3,4,5] ,A,0.2 , 0)
title(data_set{j})
