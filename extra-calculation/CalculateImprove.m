% filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/May17/output';
% downpad = [20 26 19    23    18    20    19    24];
downpad = [0 0 0 0 0 0 0 0];
filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/FinalResults/3T/slab/';
subjectid_list = ['115825'; '208024'; '208226'; '211316'; '385046'; '412528'; '555651'; '902242'];
% subjectid_list = ['100610'; '102311'; '102816'; '104416'; '105923'; '108323'; '109123'; '111312'];

for i=1:8
    % load the B0 map and B0shim map
    for slab_num = 1:2
        
       load([fullfile(filepath, label), '/B0map/', num2str(slab_num), '_',subjectid_list(i,:), '_B0.mat']);
        load([fullfile(filepath, label), '/B0map/', num2str(slab_num), '_',subjectid_list(i,:), 'B0shim.mat']);   
        B00(:,:,(slab_num-1)*5+1:slab_num*5) = B0;
        B00shim(:,:,(slab_num-1)*5+1:slab_num*5) = B0shim;       
    end
    clear B0 B0shim
    B0 = B00;
    B0shim = B00shim;
    mask = ~isnan(B0);
    B0whole = B0(mask);
    whole_std(i) = std(B0whole);
    
    B0slab = B0(:,:,:);
    B0shimslab = B0shim(:,:,:);
    B0slabf = B0slab(mask(:,:,:) ~= 0);
    B0slab_STD = std(B0slabf);
    B0shimslabf = B0shimslab(mask(:,:,:) ~= 0);
    B0slabshim_STD = std(B0shimslabf);
    improve_STD = B0slabshim_STD/B0slab_STD;
    result(i,1) = improve_STD;
    B0slab_STD_temporal_list(i,1) = B0slab_STD;
    B0slabshim_STD_temporal_list(i,1) = B0slabshim_STD;
    
    clear B0 B0shim B00 B00shim
    for slab_num = 4:5
       load([fullfile(filepath, label), '/B0map/', num2str(slab_num), '_',subjectid_list(i,:), '_B0.mat']);
        load([fullfile(filepath, label), '/B0map/', num2str(slab_num), '_',subjectid_list(i,:), 'B0shim.mat']);   
        B00(:,:,(slab_num-4)*5+1:(slab_num-3)*5) = B0;
        B00shim(:,:,(slab_num-4)*5+1:(slab_num-3)*5) = B0shim;         
    end
    clear B0 B0shim
    B0 = B00;
    B0shim = B00shim;
    mask = ~isnan(B0);
    B0whole = B0(mask);
    whole_std(i) = std(B0whole);
    
    B0slab = B0(:,:,:);
    B0shimslab = B0shim(:,:,:);
    B0slabf = B0slab(mask(:,:,:) ~= 0);
    B0slab_STD = std(B0slabf);
    B0shimslabf = B0shimslab(mask(:,:,:) ~= 0);
    B0slabshim_STD = std(B0shimslabf);
    improve_STD = B0slabshim_STD/B0slab_STD;
    result(i,2) = improve_STD;
    B0slab_STD_frontal_list(i,1) = B0slab_STD;
    B0slabshim_STD_frontal_list(i,1) = B0slabshim_STD;

    result(i,3) = std(B0whole);
end

result(9,:) = mean(result, 1);
result(10,:) = std(result, 1);

writematrix(result,fullfile(fullfile(filepath, label),'B0temporal_frontal.xlsx'));

writematrix(B0slab_STD_temporal_list,fullfile(fullfile(filepath, label),'B0slab_STD_temporal_list.xlsx'));
writematrix(B0slabshim_STD_temporal_list,fullfile(fullfile(filepath, label),'B0slabshim_STD_temporal_list.xlsx'));

writematrix(B0slab_STD_frontal_list,fullfile(fullfile(filepath, label),'B0slab_STD_frontal_list.xlsx'));
writematrix(B0slabshim_STD_frontal_list,fullfile(fullfile(filepath, label),'B0slabshim_STD_frontal_list.xlsx'));