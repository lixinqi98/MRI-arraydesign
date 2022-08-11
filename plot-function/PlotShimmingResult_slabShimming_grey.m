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
    
    for i=1:5
        temp = B0(:,:,i);
        maxB0 = max(temp(:));
        temp = temp./maxB0;
        temp(isnan(temp)) = 0.5;
        subplot(2,5,i), axis off,imshow(temp, [-1 1])
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i)));
    end
    colorbar('Position',[0.95 0.6 0.01 0.3]);
    % ------------------------ Plot B0 after shimming ---------------------%
    for i=1:5
        temp = B0shim(:,:,i);
        maxB0 = max(temp(:));
        temp = temp./maxB0;
        temp(isnan(temp)) = 0.5;
        subplot(2,5,i+5), imshow(temp, [-1 1])
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i)));
    end
    colorbar('Position',[0.95 0.15 0.01 0.3]);

    savefig(figure1,fullfile(imgpath,[num2str(slab_num), '_', subjectid, '_slabShimming_grey.fig']));
    saveas(figure1,fullfile(imgpath,[num2str(slab_num),'_', subjectid,'_slabShimming_grey.png']));
end
