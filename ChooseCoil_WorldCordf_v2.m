function [delete_index,mean_STD,reduce_flag] = ChooseCoil_WorldCordf_v2(method,Bz_iter,B0_allbrain,mask_allbrain,DC_limit,iteration,resultprefix,threshold,fid)
% using PCA to find the redundant coil index in this iteration. Test the
% pca results on each subject, if the number of converged subjects is
% larger than threshould, the loop should stop.
% Created. 20/10/2019 Mona
% 
% Params: 
%       @method         : if 1, use lsqlin; others, use mldivide
%       @Bzf_iter       : flattened Bz map, updated in each iteration
%       @B0_allbrain    : B0 field maps of all brains after masking
%       @mask_allbrain  : masks
%       @DC_limit       : DC current limit (A)
%       @iteration      :
%       @resultprefix   :
%       @threshold      : converge requirement
%       @fid            : specify the file to save iteration results


% can not use batch because of function lsqlin
    tic
    [nx,ny,nz,subjectnum] = size(B0_allbrain);
    [~,~,~,nc] = size(Bz_iter);
    
    
    slab(:,1) = 16:20;
    slab(:,2) = 21:25;
    slab(:,3) = 26:30;
    slab(:,4) = 31:35;
    slab(:,5) = 36:40;
    slab(:,6) = 41:45;
%     mean_STD(7*1) record the 'average improve %' on subjectnum(16) subjects of 6 slab and whole brain shimming
%     slab_improve_STD(subjectnum*7) record the STD improvement of each
%     subject in each slab
%     coil_weight(16*nc*7) record the each coil's weight on every subject
%     in each slab
    % slab shimming
    for i = 1:6
        [slab_coil_weight, improve_STD, improve_mean, flag] = PCA_getWeight(B0_allbrain(:,:,slab(:,i),:), Bz_iter(:,:,slab(:,i),:), mask_allbrain(:,:,slab(:,i),:), 15, subjectnum, nc, method);
         
        mean_STD(i) = mean(improve_STD); 
        slab_improve_STD(:,i) = improve_STD;
        slab_improve_mean(:,i) = improve_mean;
        coil_weight(:,:,i) = slab_coil_weight;
    end
    % whole brain shimming
    [slab_coil_weight, improve_STD, improve_mean, flag] = PCA_getWeight(B0_allbrain, Bz_iter, mask_allbrain, 10, subjectnum, nc, method);
    mean_STD(7) = mean(improve_STD);
    slab_improve_STD(:,7) = improve_STD;
    slab_improve_mean(:,7) = improve_mean;
    coil_weight(:,:,7) = slab_coil_weight;
    
    save(fullfile(resultprefix,['iteration_',num2str(iteration),'_coilWeights_',num2str(nc),'.mat']),'coil_weight');
    
    save(fullfile(resultprefix,['iteration_',num2str(iteration),'_slabImproveSTD_',num2str(nc),'.mat']),'slab_improve_STD');
    save(fullfile(resultprefix,['iteration_',num2str(iteration),'_slabImproveMean_',num2str(nc),'.mat']),'slab_improve_mean');

    ensemble_coil_weight = sum(coil_weight(:,:,1),1)./12 + sum(coil_weight(:,:,2),1)./12+ sum(coil_weight(:,:,3),1)./12 + sum(coil_weight(:,:,4),1)./12 + sum(coil_weight(:,:,5),1)./12 + sum(coil_weight(:,:,6),1)./12 + sum(coil_weight(:,:,7),1)./2;        
    
    elapsedTime = toc;
    fprintf('%12s   iteration: %3d|improve_STD slab1:%6.4f|slab2:%6.4f|slab:%6.4f|slab4:%6.4f|slab5:%6.4f|slab6:%6.4f|whole brain:%6.4f|  curr_num:%3d|   time: %5.2f\r\n',datestr(now),iteration,mean_STD(1),mean_STD(2),mean_STD(3),mean_STD(4),mean_STD(5),mean_STD(6),mean_STD(7),nc,elapsedTime);
    fprintf(fid,'%12s   iteration: %3d|improve_STD slab1:%6.4f  slab2:%6.4f  slab3:%6.4f  slab4:%6.4f  slab5:%6.4f  slab6:%6.4f  whole brain:%6.4f|  curr_num:%5d|   time: %5.2f\r\n',datestr(now),iteration,mean_STD(1),mean_STD(2),mean_STD(3),mean_STD(4),mean_STD(5),mean_STD(6),mean_STD(7),nc,elapsedTime); 
    
    if (sum(flag) >= threshold)
        reduce_flag = 0;
        delete_index = [];
        return;
    end
    reduce_flag = 1;
    sumweight = sum(ensemble_coil_weight,1);
    [~,sortindex] = sort(sumweight);
    reducedper = round(nc/100)+1;
    delete_index = sortindex(1:reducedper);
end

function [coil_weight, improve_STD, improve_mean, flag] = PCA_getWeight(B0, Bz, mask, DC_limit, subjectnum, nc, method)
    flag = zeros([subjectnum 1]);
    improve_STD = zeros([subjectnum 1]);
    improve_mean = zeros([subjectnum 1]);
    for i = 1:nc
        Bz_ROI = Bz(:,:,:,i);
        Bzf(:,i) = Bz_ROI(:);
    end
    for i = 1:subjectnum
        tempB0 = B0(:,:,:,i);
        tempMask = mask(:,:,:,i);
        B0f_mask = double(tempB0(tempMask));
        Bzf_mask = double(Bzf(tempMask,:));
        lb0 = -ones(nc,1)*DC_limit;         %lower bound  dc limit
        ub0 = ones(nc,1)*DC_limit;          % upper bound dc limit 
        X0 = zeros(nc,1);                   % initial value, 0
        if(method == 1)
            options7 = optimset('Algorithm','trust-region-reflective','MaxFunEvals',1e14, 'MaxIter',...
                1e12,'TolFun',1e-20,'TolX',1e-16,'Disp','off');%,'Largescale','off',);
            X = lsqlin(Bzf_mask,B0f_mask,[],[],[],[],lb0,ub0,X0,options7);
        else
            X = mldivide(Bzf,B0f);   
        end
        [improve_STD(i),improve_mean(i),~] = OutputResultforEachBrain_WorldCordf(X,Bzf_mask,B0f_mask);
        DC(:,i) = reshape(X,[nc 1]); % in (A)
    end

%     Using PCA
    for i = 1:subjectnum
        Bzf_ebrain = Bzf(mask(:,:,:,i),:);
        Bz_applycurrent = Bzf_ebrain.*repmat(DC(:,i)',size(Bzf_ebrain,1),1);
        [coeff,score,latent,tsquared,explained,mu] = pca(Bz_applycurrent);
        percent = cumsum(latent)./sum(latent);
        index = find(percent>0.9999,1);
        % check converge
        if(index == nc)
            flag(i) = 1;
        end
        tempweight = sum(abs(coeff(:,1:index)'));
        coil_weight(i,:) = tempweight/norm(tempweight);
    end

end

