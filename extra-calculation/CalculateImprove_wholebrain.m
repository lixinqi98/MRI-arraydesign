% filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/May17/output';
% downpad = [20 26 19    23    18    20    19    24];
downpad = [0 0 0 0 0 0 0 0];
filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/FinalResults/7T/wholebrain';
% subjectid_list = ['115825'; '208024'; '208226'; '211316'; '385046'; '412528'; '555651'; '902242'];
subjectid_list = ['132118'; '134627'; '134829'; '135124'; '137128'; '140117'; '144226'; '145834'];
% subjectid_list = ['100610'; '102311'; '102816'; '104416'; '105923'; '108323'; '109123'; '111312'];
slab(:,1) = 52:56;
slab(:,2) = 57:61;
slab(:,3) = 62:66;
slab(:,4) = 67:71;
slab(:,5) = 72:76;
slab(:,6) = 77:81;
temporal = 52:61;
frontal = 67:71;
% slab(:,1) = 16:20;
% slab(:,2) = 21:25;
% slab(:,3) = 26:30;
% slab(:,4) = 31:35;
% slab(:,5) = 36:40;
% slab(:,6) = 41:45;
% temporal = 16:25;
% frontal = 31:40;
for i=1:8
    % load the B0 map and B0shim map
    load([fullfile(filepath, label), '/ImageReport/', subjectid_list(i,:), '_B0.mat']);
    load([fullfile(filepath, label), '/ImageReport/', subjectid_list(i,:), '_B0shim.mat']);

    mask = ~isnan(B0);
    B0whole = B0(mask);
    B0shimwhole = B0shim(mask);
    whole_std(i) = std(B0whole);
    
    B0slab = B0(:,:,temporal+downpad(i));
    B0shimslab = B0shim(:,:,temporal+downpad(i));
    B0slabf = B0slab(mask(:,:,temporal+downpad(i)) ~= 0);
    B0slab_STD = std(B0slabf);
    B0shimslabf = B0shimslab(mask(:,:,temporal+downpad(i)) ~= 0);
    B0slabshim_STD = std(B0shimslabf);
    improve_STD = B0slabshim_STD/B0slab_STD;
    result(i,1) = improve_STD;
    B0slab_STD_temporal_list(i,1) = B0slab_STD;
    B0slabshim_STD_temporal_list(i,1) = B0slabshim_STD;
    
    clear B0slab B0shimslab B0slabf B0shimslabf
    B0slab = B0(:,:,frontal+downpad(i));
    B0shimslab = B0shim(:,:,frontal+downpad(i));
    B0slabf = B0slab(mask(:,:,frontal+downpad(i)) ~= 0);
    B0slab_STD = std(B0slabf);
    B0shimslabf = B0shimslab(mask(:,:,frontal+downpad(i)) ~= 0);
    B0slabshim_STD = std(B0shimslabf);
    improve_STD = B0slabshim_STD/B0slab_STD;
    result(i,2) = improve_STD;
    B0slab_STD_frontal_list(i,2) = B0slab_STD;
    B0slabshim_STD_frontal_list(i,2) = B0slabshim_STD;

    result(i,3) = std(B0shimwhole);
end

result(9,:) = mean(result, 1);
result(10,:) = std(result, 1);

writematrix(result,fullfile(fullfile(filepath, label),'B0temporal_frontal.xlsx'));
writematrix(B0slab_STD_temporal_list,fullfile(fullfile(filepath, label),'B0slab_STD_temporal_list.xlsx'));
writematrix(B0slabshim_STD_temporal_list,fullfile(fullfile(filepath, label),'B0slabshim_STD_temporal_list.xlsx'));

writematrix(B0slab_STD_frontal_list,fullfile(fullfile(filepath, label),'B0slab_STD_frontal_list.xlsx'));
writematrix(B0slabshim_STD_frontal_list,fullfile(fullfile(filepath, label),'B0slabshim_STD_frontal_list.xlsx'));