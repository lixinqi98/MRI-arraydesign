% MainShimming_WorldCordf_SHorder.m 
% Shimming the brain in the world coordinate and ouput the shimming results
% Use SH_order to generate the shimming field

% Created. 17/10/2019 Mona

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

% change the input excel file for different design
load(fullfile(datapath, 'helmetposition.mat'));

% set resolution and thickness in unit(m)
ratio = 1000; % convert mm to m

% change the helmet and coil position to fit the world coordinate of the
% brain
helmetx = helmetcordf(1,:);
helmety = helmetcordf(3,:);
helmetz = helmetcordf(2,:);
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
    [nx0,ny0,nz0] = size(allbrain_img(:,:,:,i));
    res = example_brain.hdr.dime.pixdim(2:4)./ratio;   % resolution of the slice (mm -> m),in HCP dataset the resolution is 2mm
%     affine_matrix = cat(1,example_brain.hdr.hist.srow_x,example_brain.hdr.hist.srow_y,example_brain.hdr.hist.srow_z,[0 0 0 1]);
    affine_matrix(1:3,1:3) = eye(3).*repmat(res',1,3).*ratio;
    orig = [-nx0/2*res(1) -ny0/2*res(2) -nz0/2*res(3)]';
    affine_matrix(1:3,4) = orig.*ratio;
    [x11,y11,z11] = ind2sub(size(allbrain_img(:,:,:,i)),find(~isnan(allbrain_img(:,:,:,i))));
    index = cat(1,x11',y11',z11',ones(size(x11')));
    world_cordf = affine_matrix*index;
    x11f = world_cordf(1,:)/ratio;
    y11f = world_cordf(2,:)/ratio;
    z11f = world_cordf(3,:)/ratio;
    zoffset(i) = abs(max(helmetz(:)))-abs(max(z11f))-0.0065;
    z11f = z11f + zoffset(i);
% 
%     plot3(x11f,y11f,z11f,'b.','MarkerSize',6),hold on;
%     scatter3(helmetx,helmety,helmetz),hold on;
%     hold on
end
%% STEP 3: generate the Bz field map
[nx0, ny0, nz0, ~] = size(allbrain_img);
brain_point(:,1) = affine_matrix*[1 1 1 1]'./ratio;
brain_point(:,2) = affine_matrix*[nx0 ny0 nz0 1]'./ratio;
[y,x,z] = meshgrid(brain_point(2,1):res(1):brain_point(2,2),brain_point(1,1):res(2):brain_point(1,2),brain_point(3,1):res(3):brain_point(3,2));
bz_point(:,1) = [x(1,1,1),y(1,1,1),z(1,1,1)]';
[n1,n2,n3] = size(x);
bz_point(:,2) = [x(n1,n2,n3),y(n1,n2,n3),z(n1,n2,n3)]';

[B_collect, S_collect, title_collect] = SH_order_20210527(order,x.*100,y.*100,z.*100); 
N_coils = length(B_collect);
% Create 4-D array of Bz field including the (x,y,z,# of coil)
Bz_mapped = zeros([size(x,1) size(x,2) size(x,3) N_coils]); % preallocation

for jr = 1:N_coils
    temp = cell2mat(B_collect(jr));
    Bz_mapped(:,:,:,jr) = temp .* S_collect(jr);
end
save(fullfile(resultprefix,[label,'_Bz_mapped.mat']),'Bz_mapped');

[x1,y1,z1,~] = size(Bz_mapped(:,:,:,:));
resize_allbrain_img = allbrain_img(1:x1,1:y1,1:z1,:); 



%% STEP 4: Shimming for each brain
close all
fid = fopen(fullfile(resultprefix,'Shimmingresultlog.txt'),'a+');
fprintf(fid,"\r\n\r\n\r\n\r\n=============================%s===================================\r\n",datestr(now));
% DC_limit = 10;   % DC current limit (A), from paper "Modeling Real Shim Fields 
%                         % for Very High Degree (and Order) B0 Shimming of the Human Brain at 9.4 T"
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
        CoilShimmingHcp_WorldCordf(resize_allbrain_img(:,:,:,i),Bzf_temp,DC_limit,resultprefix,subjectid,0);
    % output the result and save the result
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
std_result = cat(2,B0_STD_list',B0_mean_list',B0shim_STD_list',B0shim_mean_list');
writematrix(std_result,fullfile(resultprefix,'std_result.xlsx'));
fprintf(fid,'\r\n\r\n\r\nmax_improve_STD: %3.4f min_improve_STD: %3.4f  mean_improve_STD: %3.4f std_improve_STD: %3.4f\r\n\r\n\r\n',...
    max_improve_STD,min_improve_STD,mean_improve_STD,std_improve_STD);
