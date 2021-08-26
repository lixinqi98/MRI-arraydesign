function [delete_index,improve_STD,improve_mean,reduce_flag] = ChooseCoil_WorldCordf(method,Bzf_iter,B0_allbrain,mask_allbrain,DC_limit,iteration,resultprefix,threshold,fid)
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
    [~,nc] = size(Bzf_iter);
    for i = 1:subjectnum
        clear B0f Bzf lb0 ub0 X0
        temp_B0 = B0_allbrain(:,:,:,i);
        temp_mask = mask_allbrain(:,:,:,i);
        B0f = double(temp_B0(temp_mask));
        Bzf= double(Bzf_iter(temp_mask,:));        
        lb0 = -ones(nc,1)*DC_limit;         %lower bound  dc limit
        ub0 = ones(nc,1)*DC_limit;          % upper bound dc limit 
        X0 = zeros(nc,1);                   % initial value, 0
        if(method == 1)
            options7 = optimset('Algorithm','trust-region-reflective','MaxFunEvals',1e14, 'MaxIter',...
                1e12,'TolFun',1e-20,'TolX',1e-16,'Disp','off');%,'Largescale','off',);
            X = lsqlin(Bzf,B0f,[],[],[],[],lb0,ub0,X0,options7);
        else
            X = mldivide(Bzf,B0f);   
        end
        [improve_STD(i),improve_mean(i),~] = OutputResultforEachBrain_WorldCordf(X,Bzf,B0f);
        DC(:,i) = reshape(X,[nc 1]); % in (A)
    end
    if (size(DC,1) == nc && size(DC,2) == subjectnum)
        save(fullfile(resultprefix,['iteration_',num2str(iteration),'_coils_',num2str(nc),'DC.mat']),'DC');
        shimming_img= reshape(Bzf_iter*DC, [nx ny nz subjectnum]);
        B0shim_img = (B0_allbrain - shimming_img).* mask_allbrain;    
        f1 = figure;
        row = round(sqrt(subjectnum))+1;
        for i = 1:subjectnum
            subplot(row,row,i)
            B0shim_sli = B0shim_img(:,:,20,i);
            B0_sli = B0_allbrain(:,:,20,i);
            imshowpair(B0_sli,B0shim_sli,'montage')
        end
        savefig(f1,fullfile(resultprefix,['iteration',num2str(iteration),'_shimcompare.fig']));
        saveas(f1,fullfile(resultprefix,['iteration',num2str(iteration),'_shimcompare.png']));
    end
    close all
    improve_STD_allmean = mean(improve_STD);
    improve_mean_allmean = mean(improve_mean);
    elapsedTime = toc;
    fprintf('%12s   iteration: %5d|    improve_STD:%5.2f|  improve_mean:%5.2f|  curr_num:%5d|   time: %5.2f\r\n',datestr(now),iteration,improve_STD_allmean,improve_mean_allmean,nc,elapsedTime);
    fprintf(fid,'%12s   iteration: %5d|    improve_STD:%5.2f|  improve_mean:%5.2f|  curr_num:%5d|  time: %5.2f\r\n',datestr(now),iteration,improve_STD_allmean,improve_mean_allmean,nc,elapsedTime);
    
    flag = zeros([subjectnum 1]);
    for i = 1:subjectnum
        Bzf_iter_ebrain = Bzf_iter(mask_allbrain(:,:,:,i),:);
        Bz_applycurrent = Bzf_iter_ebrain.*repmat(DC(:,i)',size(Bzf_iter_ebrain,1),1);
        [coeff,score,latent,tsquared,explained,mu] = pca(Bz_applycurrent);
        percent = cumsum(latent)./sum(latent);
        index = find(percent>0.9999,1);
        % check converge
        if(index == nc)
            flag(i) = 1;
        end
        tempweight = sum(abs(coeff(:,1:index)'));
        weight_coil(i,:) = tempweight/norm(tempweight);
    end
    if (sum(flag) >= threshold)
        reduce_flag = 0;
        delete_index = [];
        return;
    end
    reduce_flag = 1;
    sumweight = sum(weight_coil,1);
    [~,sortindex] = sort(sumweight);
    reducedper = round(nc/100)+1;
    delete_index = sortindex(1:reducedper);
end

