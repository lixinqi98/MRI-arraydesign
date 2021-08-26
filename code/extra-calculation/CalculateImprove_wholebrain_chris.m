filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/May17/output';
% downpad = [20 26 19    23    18    20    19    24];
% downpad = [0 0 0 0 0 0 0 0];
% filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/FinalResults/3T/wholebrain';
subjectid_list = ['115825'];
slab(:,1) = 1:5;
slab(:,2) = 6:10;
slab(:,3) = 11:15;
slab(:,4) = 16:20;
slab(:,5) = 21:25;
slab(:,6) = 26:30;
slab(:,7) = 31:35;
slab(:,8) = 36:40;
slab(:,9) = 41:45;
slab(:,10) = 46:50;
slab(:,11) = 51:55;
slab(:,12) = 56:60;
slab(:,13) = 61:65;
slab(:,14) = 66:70;
for i=1
    % load the B0 map and B0shim map
    for slab_num = 1:14
        
       load([fullfile(filepath, label), '/B0map/', num2str(slab_num), '_',subjectid_list(i,:), '_B0.mat']);
        load([fullfile(filepath, label), '/B0map/', num2str(slab_num), '_',subjectid_list(i,:), '_B0shim.mat']);   
        B00(:,:,(slab_num-1)*5+1:slab_num*5) = B0;
        B00shim(:,:,(slab_num-1)*5+1:slab_num*5) = B0shim;       
    end
    mask = ~isnan(B00);
    B0whole = B00(mask);
    B0shimwhole = B00shim(mask);
    whole_std(i) = std(B0whole);
    result(i) = std(B0shimwhole);
end

writematrix(result,fullfile(fullfile(filepath, label),'B0shimstd.xlsx'));