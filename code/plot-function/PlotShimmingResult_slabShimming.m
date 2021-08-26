function [] = PlotShimmingResult_slabShimming(B0,B0shim,mask,path,subjectid,slab_num)
    for i = 1:size(B0,3)
        B0_slice = B0(:,:,i);           % get B0 of each slice in ROI
        B0_slicemask = B0_slice(mask(:,:,i) ~= 0);   % remove the non-zero term
        B0_slice_STD(i) = std(B0_slicemask);  % scale by 2*pi  
    end

    % Calculate the stdv of shimmed B0 of each slice (Hz)
    for i = 1:size(B0shim,3)
        B0shim_slice = B0shim(:,:,i);   % get B0shim of each slice in ROI
        B0shim_slicemask = B0shim_slice(mask(:,:,i) ~= 0);  % remove the non-zero term
        B0shim_slice_STD(i) = std(B0shim_slicemask);        
    end
    %% Plot the figure
    
    imgpath = fullfile(path,'ImageReport');
    m = 100;   % upper & lower limits for the B0 maps (Hz)
    
    figure1 = figure('Name',['Slab Shimming slab - ', num2str(slab_num)]);
    set(figure1,'position',[200,100,700,300])
    % ------------------------ Plot B0 before shimming ---------------------%
    temp = B0(:,:,:);
    minB0 = min(temp(:)); maxB0 = max(temp(:));
    for i=1:5
        subplot(2,5,i), axis off,imshow(B0(:,:,i),[minB0 maxB0]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i)));
    end
    colorbar('Position',[0.95 0.6 0.01 0.3]);
    % ------------------------ Plot B0 after shimming ---------------------%
    temp = B0shim(:,:,:);
    minB0 = min(temp(:)); maxB0 = max(temp(:));
    for i=1:5
        subplot(2,5,i+5), imshow(B0shim(:,:,i),[minB0 maxB0]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i)));
    end
    colorbar('Position',[0.95 0.15 0.01 0.3]);

    savefig(figure1,fullfile(imgpath,[num2str(slab_num), '_', subjectid, '_slabShimming.fig']));
    saveas(figure1,fullfile(imgpath,[num2str(slab_num),'_', subjectid,'_slabShimming.png']));
   
%     save(fullfile(imgpath,[num2str(slab_num), '_', subjectid,'B0shim.mat']), 'B0shim');
%     save(fullfile(imgpath,[num2str(slab_num), '_', subjectid,'_B0.mat']), 'B0');
end
