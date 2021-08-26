function [] = PlotShimmingResult(B0,B0shim,mask,path,subjectid,downpad)
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
    % ------------------------ Plot B0 before shimming ---------------------%
    imgpath = fullfile(path,'ImageReport');
    m = 100;   % upper & lower limits for the B0 maps (Hz)
    
    figure1 = figure('Name','Initial B0 Map 1 to 36');
    set(figure1,'position',[200,100,700,400])
    for i=1:36
        subplot(5,8,i), axis off,imshow(B0(:,:,i+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i+downpad)));
    end
    savefig(figure1,fullfile(imgpath,[subjectid,'_initialB0_sli1to36.fig']));
    saveas(figure1,fullfile(imgpath,[subjectid,'_initialB0_sli1to36.png']));

    figure2 = figure('Name','Initial B0 Map 37 to 72');
    set(figure2,'position',[200,100,700,400])
    for i=1:36
        subplot(5,8,i), imshow(B0(:,:,i+36+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i+36+downpad)));
    end
    savefig(figure2,fullfile(imgpath,[subjectid,'_initialB0_sli37to72.fig']));
    saveas(figure2,fullfile(imgpath,[subjectid,'_initialB0_sli37to72.png']));

    figure3 = figure('Name','B0 aftershimming Map 1 to 36');
    set(figure3,'position',[200,100,700,400])
    for i=1:36
        subplot(5,8,i), imshow(B0shim(:,:,i+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i+downpad)));
    end
    savefig(figure3,fullfile(imgpath,[subjectid,'_B0shimming_sli1to36.fig']));
    saveas(figure3,fullfile(imgpath,[subjectid,'_B0shimming_sli1to36.png']));

    figure4 = figure('Name','B0 aftershimming Map 37 to 72');
    set(figure4,'position',[200,100,700,400])
    for i=1:36
        subplot(5,8,i), imshow(B0shim(:,:,i+36+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i+36+downpad)));
    end
    savefig(figure4,fullfile(imgpath,[subjectid,'_B0shimming_sli37to72.fig']));
    saveas(figure4,fullfile(imgpath,[subjectid,'_B0shimming_sli37to72.png']));
    
    save(fullfile(imgpath,[subjectid,'_B0shim.mat']), 'B0shim');
    save(fullfile(imgpath,[subjectid,'_B0.mat']), 'B0');
    
end
