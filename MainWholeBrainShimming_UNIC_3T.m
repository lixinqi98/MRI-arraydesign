% Reproduce Carissa's results(based on Carissa's brain position)
% Brief conclusion:
%                 STEP 1: prepare the data path
%                 STEP 2: load data and generate figures
%                 STEP 3: generate aligned Bz
%                 STEP 4: shim 
%                 STEP 5: coil reduction

% Created. 24/05/2021 Mona

% use function B0_multi_Edgecut.m, BiotSavart100Points.m, CoilShimmingHcp_WorldCordf.m, CoilReduction_WorldCordf.m, PlotReservedReducedCoil_sphereHelmet.m

%% STEP 1 : choose the required data and data path 
% for fast implement and on different pc
addpath('tool');

fask_impl = 0;
if (fask_impl)
    datapath = uigetdir('','Select the data Folder');
    resultpath = uigetdir('','Select the result Folder');
    [hcpfile,hcppath] = uigetfile('*Phase.nii','Select the HCP data');
    if (isequal(hcpfile,0)||isequal(datapath,0)||isequal(resultpath,0))
        error('User selected Cancel');
    end
else
    datapath = '/Users/mona/codeWorkSpace/Cedars-Sinai/7TShimming/CoilData';
    resultpath = '/Users/mona/codeWorkSpace/Cedars-Sinai/May17/output';
    hcppath = '/Users/mona/codeWorkSpace/Cedars-Sinai/Data/HCPdataset_test';    
end
resultprefix = fullfile(resultpath,[label,'Results']);
mkdir(resultprefix);
imgpath = fullfile(resultprefix,'ImageReport');
mkdir(imgpath);
filelist_ph = dir(fullfile(hcppath,'*Phase.nii'));
filelist_mag = dir(fullfile(hcppath,'*Magnitude.nii'));
subjectnum = length(filelist_ph);


%% STEP 2: load coil position and generate the figure
% load in the position of coil and helmet
% using the *Original HCP data*

% change the input mat file for different design
load(fullfile(datapath,coilposition));
rr = 0.24/2; % radius of semi-sphere
height = 2*rr;
offsety = -0.01; % in m
offsetz = 0.02; % in m
% set resolution and thickness in unit(m)
ratio = 1000; % convert mm to m

% change the helmet and coil position to fit the world coordinate of the
% brain
coilx = coilcordf(1,:,:);
coily = coilcordf(3,:,:);
coilz = coilcordf(2,:,:);
mapped_coilcordf = cat(1,coilx,coily,coilz);
N_coils = size(coilcordf,3);
f1 = figure('Name','check_brain'); 
set(gcf,'position',[100,0,1000,800]);
t = title('Check Brain Position(Z)');
t.FontSize = 15;
t.FontWeight = 'bold';
for i = 1:subjectnum
    subplot(2,4,i)
    example_brain = load_nii(fullfile(filelist_ph(i).folder,filelist_ph(i).name));
%     allbrain_img(:,:,:,i) = mask_HCP(fullfile(filelist_mag(i).folder,filelist_mag(i).name),fullfile(filelist_ph(i).folder,filelist_ph(i).name));
    allbrain_img(:,:,:,i) = B0_multi_Edgecut(fullfile(filelist_mag(i).folder,filelist_mag(i).name),fullfile(filelist_ph(i).folder,filelist_ph(i).name));
    subjectname = filelist_ph(i).name;
    subjectid(i,:) = subjectname(1:6);

    [nx0,ny0,nz0] = size(allbrain_img(:,:,:,i));

    res = example_brain.hdr.dime.pixdim(2:4)./ratio;   % resolution of the slice (mm -> m),in HCP dataset the resolution is about 2mm
    affine_matrix = cat(1,example_brain.hdr.hist.srow_x,example_brain.hdr.hist.srow_y,example_brain.hdr.hist.srow_z,[0 0 0 1]);
    
%     since we don't need the brain's motion info, change the affine matrix
    affine_matrix(1:3,1:3) = eye(3).*repmat(res',1,3).*ratio;
    orig = [-nx0/2*res(1) -ny0/2*res(2)+offsety -nz0/2*res(3)+offsetz]';
    affine_matrix(1:3,4) = orig.*ratio;
    
    [x11,y11,z11] = ind2sub(size(allbrain_img(:,:,:,i)),find(~isnan(allbrain_img(:,:,:,i))));

    index = cat(1,x11',y11',z11',ones(size(x11')));
    world_cordf = affine_matrix*index;
    x11f = world_cordf(1,:)/ratio;
    y11f = world_cordf(2,:)/ratio;
    z11f = world_cordf(3,:)/ratio;
%     
%     zoffset(i) = abs(max(helmetz(:)))-abs(max(z11f))-0.0065; % make sure the top of head 6.5mm far away from the helmet
%     z11f = z11f + zoffset(i);
    for j = 1:N_coils
        % attention the transfer of the coordinate
        x1f = mapped_coilcordf(1,:,j);
        y1f = mapped_coilcordf(2,:,j);
        z1f = mapped_coilcordf(3,:,j);
        plot3(x1f,y1f,z1f,'linewidth',1.5), hold on
        xlabel x
        ylabel y
        zlabel z
    end
    
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])
    zlim([-0.15 0.15])
    plot3(x11f,y11f,z11f,'b.','MarkerSize',6),hold on;
    R = rr;   
    [X,Y,Z] = cylinder(R,100);
    Z = (Z - 1)*height/2+0.01;
    test1 = surf(X,Y,Z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    set(test1,'FaceColor',[0 0 0],'FaceAlpha',0.3)
    hold on
    R1=rr;
    [X1,Y1,Z1] = sphere;
    Z1(Z1<-0.001) = NaN;
    X1 = X1*R1; Y1 = Y1*R1; Z1 = Z1*R1;
    test11 = surf(X1,Y1,Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    set(test11 ,'FaceColor',[0 0 0],'FaceAlpha',0.3)
    hold on
end
for view_index = 1:4
    for i = 1:subjectnum
        subplot(2,4,i),
        set(gcf,'Color',[1 1 1])
        h = title(sprintf('%s',subjectid(i,:)));
        if(view_index == 1)
            view(180,0);
        elseif(view_index == 2)
            view(90,0);
        elseif(view_index == 3)
            view(0,90);
        elseif(view_index == 4)
            view(0,0);
        end
    end
    
    xlim([-0.15 0.15])
    ylim([-0.15 0.15])
    zlim([-0.15 0.15])
    savefig(f1,fullfile(resultprefix,[label,'_viewindex',num2str(view_index),'_check_position.fig']));
end

% generate a figure for last brain for detail
f2 = figure('Name','detail look');
set(gcf,'position',[0,100,2000,300])
for view_index = 1:4
    subplot(1,4,view_index)
    for j = 1:N_coils
        % attention the transfer of the coordinate
        x1f = mapped_coilcordf(1,:,j);
        y1f = mapped_coilcordf(2,:,j);
        z1f = mapped_coilcordf(3,:,j);
        xc1 = mean(x1f);
        yc1 = mean(y1f);
        zc1 = mean(z1f);
        text(double(xc1),double(yc1),double(zc1),num2str(j),'FontSize',10,'color','w')
        plot3(x1f,y1f,z1f,'linewidth',1.5), hold on
        xlabel x
        ylabel y
        zlabel z
    end
    if(view_index == 1)
        view(180,0);
    elseif(view_index == 2)
        view(90,0);
    elseif(view_index == 3)
        view(0,90);
    elseif(view_index == 4)
        view(0,0);
    end
    plot3(x11f,y11f,z11f,'b.','MarkerSize',6),hold on;
    R = rr;   
    [X,Y,Z] = cylinder(R,100);
    Z = (Z - 1)*height/2+0.01;
    test1 = surf(X,Y,Z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    set(test1,'FaceColor',[0 0 0],'FaceAlpha',0.3)
    hold on
    R1=rr;
    [X1,Y1,Z1] = sphere;
    Z1(Z1<-0.001) = NaN;
    X1 = X1*R1; Y1 = Y1*R1; Z1 = Z1*R1;
    test11 = surf(X1,Y1,Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    set(test11 ,'FaceColor',[0 0 0],'FaceAlpha',0.3)
    hold on
end
savefig(f2,fullfile(resultprefix,[label,'helmet_brain_coil.fig']));

%% STEP 3: generate the Bz field map
pad_allbrain_img = allbrain_img;
[nx0, ny0, nz0, ~] = size(pad_allbrain_img);
brain_point(:,1) = affine_matrix*[1 1 1 1]'./ratio;
brain_point(:,2) = affine_matrix*[nx0 ny0 nz0 1]'./ratio; % get two diagonal endpoints to use meshgrid
[y,x,z] = meshgrid(brain_point(2,1):res(1):brain_point(2,2),brain_point(1,1):res(2):brain_point(1,2),brain_point(3,1):res(3):brain_point(3,2));
bz_point(:,1) = [x(1,1,1),y(1,1,1),z(1,1,1)]';
[n1,n2,n3] = size(x);
bz_point(:,2) = [x(n1,n2,n3),y(n1,n2,n3),z(n1,n2,n3)]';
current = 1;
N_coils = size(coilcordf,3);

Bz_mapped = BiotSavart100Points(coilcordf,N_coils,x,y,z,current).*42.6*1e6;
save(fullfile(resultprefix,[label,'_Bz_mapped.mat']),'Bz_mapped');
% load(Bz_path);

[x1,y1,z1,~] = size(Bz_mapped(:,:,:,:));
resize_allbrain_img = pad_allbrain_img(1:x1,1:y1,1:z1,:); 

%!!!!!!!!!!!!!! check whether the brain and the bz mapped !!!!!!!!!!!!!!!%
% - 1st: check the shimming region
disp(brain_point)
disp(bz_point)
%%
% - 2nd: check the relative position of coil and the brain, especillay the
% ear coils and the forehead position
% slicer(resize_allbrain_img,allbrain_img);
f4 = figure('Name','Check the Bz and B0');
row = round(sqrt(N_coils))+1;
for i = 1:N_coils
    subplot(row,row,i)
    bz_test = Bz_mapped(:,:,40,i);
    b0_test = resize_allbrain_img(:,:,40,1);
    imshowpair(bz_test,b0_test,'montage')
    h = title(['CoilNo.',num2str(i)]);
end



%% STEP 4: Shimming for each brain
close all
fid = fopen(fullfile(resultprefix,'Shimmingresultlog.txt'),'a+');
fprintf(fid,"\r\n\r\n\r\n\r\n=============================%s===================================\r\n",datestr(now));

% DC_limit = 10;   % DC current limit (A), from paper "Modeling Real Shim Fields 
                        % for Very High Degree (and Order) B0 Shimming of the Human Brain at 9.4 T"
[nx, ny, nz, nc] = size(Bz_mapped); 
for i=1:size(Bz_mapped,4)
    Bz_ROI = Bz_mapped(:,:,:,i);
    Bzf_temp (:,i) = Bz_ROI(:);
end
for i = 1:subjectnum   
    subjectname = filelist_ph(i).name;
    subjectid = subjectname(1:6);
    % calculate the shimming result
    [improve_STD,improve_mean,DC,B0shim_STD,B0shim_mean,B0_STD,B0_mean] = ...
        CoilShimmingHcp_WorldCordf(resize_allbrain_img(:,:,:,i),double(Bzf_temp),DC_limit,resultprefix,subjectid,0);% output the result and save the result
    fprintf("%12s|  %8s|     %5.4f| %5.4f|  %5.4f\r\n",datestr(now),subjectid,improve_STD,improve_mean,max(DC(:)));
    fprintf(fid,"%12s    %8s    %5.4f   %5.4f   %5.4f\r\n",datestr(now),subjectid,improve_STD,improve_mean,max(DC(:)));
    coilcurrent(i,:) = DC(:)';
    improve_STD_record(i) = improve_STD;
    improve_mean_record(i) = improve_mean;
    B0shim_STD_list(i) = B0shim_STD;
    B0_STD_list(i) = B0_STD;
    B0shim_mean_list(i) = B0shim_mean;
    B0_mean_list(i) = B0_mean;
end
save(fullfile(resultprefix,'coilcurrent_allsubject.mat'),'coilcurrent');
save(fullfile(resultprefix,'improve_STD_allsubject.mat'),'improve_STD_record');
save(fullfile(resultprefix,'improve_mean_allsubject.mat'),'improve_mean_record');


tran_improve_STD = improve_STD_record';
tran_improve_mean = improve_mean_record';
max_improve_STD = max(tran_improve_STD);
min_improve_STD = min(tran_improve_STD);
mean_improve_STD = mean(tran_improve_STD);
std_improve_STD = std(tran_improve_STD);

max_current = max(abs(coilcurrent),[],2);
mean_current = mean(abs(coilcurrent),2);
result = cat(2,tran_improve_STD,max_current,tran_improve_mean,mean_current);
writematrix(result,fullfile(resultprefix,'result.xlsx'));
std_result = cat(2,B0shim_STD_list',B0_STD_list',B0shim_mean_list',B0_mean_list');
writematrix(std_result,fullfile(resultprefix,'std_result.xlsx'));
fprintf(fid,'\r\n\r\n\r\nmax_improve_STD: %3.4f min_improve_STD: %3.4f  mean_improve_STD: %3.4f std_improve_STD: %3.4f\r\n\r\n\r\n',...
    max_improve_STD,min_improve_STD,mean_improve_STD,std_improve_STD);


% % %% STEP 5 : Coil Reduction
% method = 1;
% resultprefix = fullfile(resultprefix,'coilreduction');
% mkdir(resultprefix);
% 
% CoilReduction_WorldCordf(method,resize_allbrain_img,Bz_mapped,DC_limit,resultprefix,label)
% 
% PlotReservedReducedCoil_sphereHelmet(resultprefix,label,mapped_coilcordf,x11f,y11f,z11f,rr)