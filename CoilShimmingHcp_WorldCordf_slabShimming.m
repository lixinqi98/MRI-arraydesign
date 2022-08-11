function [improve_STD,improve_mean,DC,B0shim_STD,B0shim_mean,B0_STD,B0_mean] = CoilShimmingHcp_WorldCordf_slabShimming(B0,Bzf_temp,DC_limit,resultprefix,subjectid,slab_num)
% CoilShimmingHcp_WorldCordf.m
% calculate the shimming result for each subject and plot the shimming
% results
% Created. 20/06/2021 Mona
% 
% Params: 
%       @B0              : resized B0 field map
%       @Bzf_tmp         : reshaped Bz (flattened)
%       @DC_limit        : DC current limit (A)
%       @resultprefix    : path prefix to save results
%       @subjectid       : 
%       @downpad         : added padding for the brain along z-axis(only
%                          used to plot figures)


    % Zero the negative B0 that is <-40*2pi Hz
    % *todo : find out why scale down 2pi for hcp data, maybe the problem of unwrap*
    mask_neg = (B0 <= -40*2*pi);
    B0(mask_neg) = nan;
    B0 = B0./(2*pi); % (do not use this when shimming 7T brain)   
    mask = ~isnan(B0);  
    B0 = B0.*mask;
    
    B0f = double(B0(mask)); 
    Bzf = Bzf_temp(mask(:),:);
    [~,nc] = size(Bzf_temp);
    lb0 = -ones(nc,1)*DC_limit;  %lower bound  dc limit
    ub0 = ones(nc,1)*DC_limit;   % upper bound dc limit 
    X0 = zeros(nc,1);                 % initial value, 0

    options7 = optimset('Algorithm','trust-region-reflective','MaxFunEvals',1e14, 'MaxIter',...
        1e12,'TolFun',1e-20,'TolX',1e-16,'Disp','off');%,'Largescale','off',);
    [X,resnorm,residual] = lsqlin(Bzf,B0f,[],[],[],[],lb0,ub0,X0,options7);
%     X = mldivide(Bzf,B0f);
    DC = reshape(X,[nc 1]); % in (A)

    % find the shimming B field of coils (whole brain)
    shimming = sum(Bzf*DC,2);
    B0shim = B0f - shimming;
    shimming= reshape(sum(Bzf_temp*DC,2), size(B0));  % in (Hz)

    B0shim_temp = (B0 - shimming).* mask;
    
    % improve_std: the std of the brain after shimming / the std of the initial brain
    % improve_mean: the abs mean of the brain after shimming / the abs mean of the initial brain
    B0_STD = std(B0f);
    B0shim_STD = std(B0shim);
    improve_STD = B0shim_STD/B0_STD;
    B0_mean = mean(abs(B0f));
    B0shim_mean = mean(abs(B0shim));
    improve_mean = B0shim_mean/B0_mean;
    
%     plot the shimming results
%     PlotShimmingResult(B0,B0shim_temp,mask,resultprefix,subjectid,downpad);
    PlotShimmingResult_slabShimming_grey(B0,B0shim_temp,mask,resultprefix,subjectid,slab_num);
    PlotShimmingResult_slabShimming(B0,B0shim_temp,mask,resultprefix,subjectid,slab_num);
    clear B0shim
    B0shim = B0shim_temp;
    imgpath = fullfile(resultprefix,'B0map');
%     mkdir(imgpath)
    save(fullfile(imgpath,[num2str(slab_num), '_', subjectid,'_B0shim.mat']), 'B0shim');
    save(fullfile(imgpath,[num2str(slab_num), '_', subjectid,'_B0.mat']), 'B0');
    close all
end

