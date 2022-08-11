% SH downpad 20
filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/May17/output';
% filepath = '/Users/mona/codeWorkSpace/Cedars-Sinai/FinalResults/3T/wholebrain';
% downpad = [20 26 19    23    18    20    19    24];
downpad = [0 0 0 0 0 0 0 0];
% filepath = '/Users/mona/OneDrive/Cedars-sinai/2021/ShimmingResults/Jul2';
% subjectid_list = ['115825'; '208024'; '208226'; '211316'; '385046'; '412528'; '555651'; '902242'];
subjectid_list = ['132118'; '134627'; '134829'; '135124'; '137128'; '140117'; '144226'; '145834'];
slab(:,1) = 52:56;
slab(:,2) = 57:61;
slab(:,3) = 62:66;
slab(:,4) = 67:71;
slab(:,5) = 72:76;
slab(:,6) = 77:81;
% slab(:,1) = 16:20;
%     slab(:,2) = 21:25;
%     slab(:,3) = 26:30;
%     slab(:,4) = 31:35;
%     slab(:,5) = 36:40;
%     slab(:,6) = 41:45;
for i=1:8
    % load the B0 map and B0shim map
    load([fullfile(filepath, label), '/ImageReport/', subjectid_list(i,:), '_B0.mat']);
    load([fullfile(filepath, label), '/ImageReport/', subjectid_list(i,:), '_B0shim.mat']);
    mask = ~isnan(B0);
    B0whole = B0(mask);
    whole_std(i) = std(B0whole);
    for j = 1:6
        B0slab = B0(:,:,slab(:,j)+downpad(i));
        B0shimslab = B0shim(:,:,slab(:,j)+downpad(i));
        B0slabf = B0slab(mask(:,:,slab(:,j)+downpad(i)) ~= 0);
        B0slab_STD = std(B0slabf);
        B0shimslabf = B0shimslab(mask(:,:,slab(:,j)+downpad(i)) ~= 0);
        B0slabshim_STD = std(B0shimslabf);
        improve_STD = B0slabshim_STD/B0slab_STD;
        result(i,j) = improve_STD;
        B0slab_STD_list(i,j) = B0slab_STD;
        B0slabshim_STD_list(i,j) = B0slabshim_STD;
    end
    result(i,7) = std(B0whole);
end

mean_slab = mean(result, 1);
std_slab = std(result, 1);

mean_slabShimming = mean(B0slabshim_STD_list, 1);
std_slabShimming = std(B0slabshim_STD_list, 1);
result_1 = cat(1, mean_slabShimming, std_slabShimming);
result = cat(1, mean_slab, std_slab);

writematrix(result,fullfile(fullfile(filepath, label),'std_slab_summary.xlsx'));
writematrix(result_1,fullfile(fullfile(filepath, label),'std_summary.xlsx'));

writematrix(B0slab_STD_list,fullfile(fullfile(filepath, label),'B0slab_STD_list.xlsx'));
writematrix(B0slabshim_STD_list,fullfile(fullfile(filepath, label),'B0slabshim_STD_list.xlsx'));