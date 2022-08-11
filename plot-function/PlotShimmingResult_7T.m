function [] = PlotShimmingResult_7T(B0,B0shim,mask,path,subjectid,downpad)
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
    
    figure1 = figure('Name','Initial B0 Map 1 to 33');
    set(figure1,'position',[200,100,700,400])
    for i=1:33
        subplot(5,7,i), axis off,imshow(B0(:,:,i+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i+downpad)));
    end
    savefig(figure1,fullfile(imgpath,[subjectid,'_initialB0_sli1to33.fig']));
    saveas(figure1,fullfile(imgpath,[subjectid,'_initialB0_sli1to33.png']));

    figure2 = figure('Name','Initial B0 Map 34 to 66');
    set(figure2,'position',[200,100,700,400])
    for i=1:33
        subplot(5,7,i), imshow(B0(:,:,i+33+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i+33+downpad)));
    end
    savefig(figure2,fullfile(imgpath,[subjectid,'_initialB0_sli34to66.fig']));
    saveas(figure2,fullfile(imgpath,[subjectid,'_initialB0_sli34to66.png']));
    
    figure3 = figure('Name','Initial B0 Map 67 to 99');
    set(figure3,'position',[200,100,700,400])
    for i=1:33
        subplot(5,7,i), imshow(B0(:,:,i+66+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i+66+downpad)));
    end
    savefig(figure3,fullfile(imgpath,[subjectid,'_initialB0_sli67to99.fig']));
    saveas(figure3,fullfile(imgpath,[subjectid,'_initialB0_sli67to99.png']));
    
    figure4 = figure('Name','Initial B0 Map 100 to 131');
    set(figure2,'position',[200,100,700,400])
    for i=1:32
        subplot(5,7,i), imshow(B0(:,:,i+99+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0_slice_STD(i+99+downpad)));
    end
    savefig(figure4,fullfile(imgpath,[subjectid,'_initialB0_sli100to131.fig']));
    saveas(figure4,fullfile(imgpath,[subjectid,'_initialB0_sli100to131.png']));

%     after shimming   
    figure5 = figure('Name','B0 aftershimming Map 1 to 36');
    set(figure5,'position',[200,100,700,400])
    for i=1:33
        subplot(5,7,i), imshow(B0shim(:,:,i+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i+downpad)));
    end
    savefig(figure5,fullfile(imgpath,[subjectid,'_B0shimming_sli1to33.fig']));
    saveas(figure5,fullfile(imgpath,[subjectid,'_B0shimming_sli1to33.png']));

    figure6 = figure('Name','B0 aftershimming Map 34 to 66');
    set(figure6,'position',[200,100,700,400])
    for i=1:33
        subplot(5,7,i), imshow(B0shim(:,:,i+33+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i+33+downpad)));
    end
    savefig(figure6,fullfile(imgpath,[subjectid,'_B0shimming_sli34to66.fig']));
    saveas(figure6,fullfile(imgpath,[subjectid,'_B0shimming_sli34to66.png']));
    
    figure7 = figure('Name','B0 aftershimming Map 67 to 99');
    set(figure7,'position',[200,100,700,400])
    for i=1:33
        subplot(5,7,i), imshow(B0shim(:,:,i+66+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i+66+downpad)));
    end
    savefig(figure7,fullfile(imgpath,[subjectid,'_B0shimming_sli67to99.fig']));
    saveas(figure7,fullfile(imgpath,[subjectid,'_B0shimming_sli67to99.png']));

    figure8 = figure('Name','B0 aftershimming Map 100 to 131');
    set(figure8,'position',[200,100,700,400])
    for i=1:32
        subplot(5,7,i), imshow(B0shim(:,:,i+99+downpad),[-m m]),  colormap jet 
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%.1f',B0shim_slice_STD(i+99+downpad)));
    end
    savefig(figure8,fullfile(imgpath,[subjectid,'_B0shimming_sli100to131.fig']));
    saveas(figure8,fullfile(imgpath,[subjectid,'_B0shimming_sli100to131.png']));
    
    save(fullfile(imgpath,[subjectid,'_B0shim.mat']), 'B0shim');
    save(fullfile(imgpath,[subjectid,'_B0.mat']), 'B0');
end
