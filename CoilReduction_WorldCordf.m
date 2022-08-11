function [] = CoilReduction_WorldCordf(method,B0_allbrain,Bz_mapped,DC_limit,resultprefix,label)
% find the redundant coils
% Created. 20/10/2019 Mona
% 
% use function: ChooseCoil_WorldCordf.m, PlotReductionCurve.m
% Params: 
%       @method       : choose optimize methods
%       @B0_allbrain  : B0 field maps of all brains
%       @Bz_mapped    : Generated Bz map
%       @DC_limit     : DC current limit (A)
%       @resultprefix : specify path to save results
%       @label        : used to distinguish different coil settings

    fid = fopen(fullfile(resultprefix,[label,'_iteration_resultslog.txt']),'a+');
    fprintf(fid,"\r\n\r\n\r\n=================================%s===================================\r\n",datestr(now));
    [nx, ny, nz, nc] = size(Bz_mapped);
    subjectnum = size(B0_allbrain,4);
    mask_neg = (B0_allbrain <= -40*2*pi);
    B0_allbrain(mask_neg) = nan;
    
%     do not use this if shimming on 7T brains
%     B0_allbrain = B0_allbrain./(2*pi);    
    
    mask_allbrain = ~isnan(B0_allbrain); 
    B0_allbrain = B0_allbrain.*mask_allbrain;
    union_mask = false([nx ny nz]);
    for i = 1:subjectnum
        mask = mask_allbrain(:,:,:,i);
        union_mask = union_mask | mask;
    end
    for i=1:size(Bz_mapped,4)
        Bz_ROI = Bz_mapped(:,:,:,i).*union_mask;
        Bzf_unionmask(:,i) = Bz_ROI(union_mask);
        Bzf_iter(:,i) = Bz_ROI(:);
    end
    N_coils = size(Bz_mapped,4);
    % set up required parameters
    threshold = subjectnum / 2;
    iteration = 1;
    flag = 1;
    indexlist = 1:nc;
    deletelist = 1:nc;

    % do the rank for Bzf
    % [~,reserved_index] = rref(Bzf_unionmask,0.1);
    reserved_index = 1:nc;
    totaldeletenum = nc - length(reserved_index);
    deletelist(reserved_index) = [];

    Bzf_iter(:,deletelist) = [];
    indexlist(deletelist) = [];
    tic
    while (flag == 1)
        clear delete_index improve_STD improve_mean delete_index
        currentnum = nc-totaldeletenum;
        if(currentnum == 0)
            fprintf(fid,'Coil number is 0. Can not Reduce');
            fprintf('Coil number is 0. Can not Reduce');
            break
        end
        [delete_index,improve_STD,improve_mean,flag] = ChooseCoil_WorldCordf(method,Bzf_iter,B0_allbrain,mask_allbrain,DC_limit,iteration,resultprefix,threshold,fid);

        if(~flag)
            fprintf('dimension can not be reduced,compelete the coil reduction');
            fprintf(fid,'dimension can not be reduced,compelete the coil reduction');
            break
        end
        improve_STD_iter(iteration,:) = improve_STD;
        improve_mean_iter(iteration,:) = improve_mean;
        iteration = iteration + 1;    
        Bzf_iter(:,delete_index) = [];
        deletelist(totaldeletenum+1:totaldeletenum+length(delete_index)) = indexlist(delete_index);
        indexlist(delete_index) = [];
        totaldeletenum = totaldeletenum + length(delete_index);
    end
    toc
    STD_curve = mean(improve_STD_iter,2);
    mean_curve = mean(improve_mean_iter,2);
    save(fullfile(resultprefix,'STD_curve.mat'),'STD_curve');
    save(fullfile(resultprefix,'mean_curve.mat'),'mean_curve');
    save(fullfile(resultprefix,'deletelist.mat'),'deletelist');
    save(fullfile(resultprefix,'improve_STD_iter.mat'),'improve_STD_iter');
    save(fullfile(resultprefix,'improve_mean_iter.mat'),'improve_mean_iter');
    PlotReductionCurve(mean_curve,STD_curve,resultprefix,N_coils);
end

