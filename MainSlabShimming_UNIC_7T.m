% Try to reproduce Carissa's results
% Brief conclusion:
%                 STEP 1: prepare the data path
%                 STEP 2: load data and generate the figure(
%                 STEP 3: generate Bz
%                 STEP 4: shim 
%                 STEP 5: reduce the redundant coil


% Created. 24/05/2020 Mona


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
    hcppath = '/Users/mona/codeWorkSpace/Cedars-Sinai/Data/Processed_Data_test';    
end
resultprefix = fullfile(resultpath,[label,'Results']);
mkdir(resultprefix);
imgpath = fullfile(resultprefix,'ImageReport');
mkdir(imgpath);
filelist_ph = dir(fullfile(hcppath,'*.nii'));
subjectnum = length(filelist_ph);


%% STEP 2: load coil position and generate the figure
% load in the position of coil and helmet
% using the *Original HCP data*

% change the input mat file for different design
load(fullfile(datapath,coilposition));
rr = 0.24/2; % radius of semi-sphere
height = 2*rr;
% set resolution and thickness in unit(m)
ratio = 1000; % convert mm to m

for i = 1:subjectnum
    example_brain = load_nii(fullfile(filelist_ph(i).folder,filelist_ph(i).name));
    brain_size(:,i) = size(example_brain.img);
end

max_x = max(brain_size(1,:));
max_y = max(brain_size(2,:));
max_z = max(brain_size(3,:));
for i = 1:subjectnum
    example_brain = load_nii(fullfile(filelist_ph(i).folder,filelist_ph(i).name));
    x_pad = max_x - brain_size(1,i);
    y_pad = max_y - brain_size(2,i);
    z_pad = max_z - brain_size(3,i);
    temp_img = padarray(example_brain.img, [round(x_pad/2),round(y_pad/2),round(z_pad/2)], nan, 'post');
    temp_img = padarray(temp_img, [x_pad-round(x_pad/2),y_pad-round(y_pad/2),z_pad-round(z_pad/2)], nan, 'pre');
    allbrain_img(:,:,:,i) = temp_img;
end
% use the same affine matrix for all brains, because we ignore the motion
% of brain here.
offsetx = 0;
offsety = -0.020;
offsetz = 0.03;
[nx0,ny0,nz0] = size(allbrain_img(:,:,:,i));
res = example_brain.hdr.dime.pixdim(2:4)./ratio;   % resolution of the slice (mm -> m),in HCP dataset the resolution is about 2mm
affine_matrix = cat(1,example_brain.hdr.hist.srow_x,example_brain.hdr.hist.srow_y,example_brain.hdr.hist.srow_z,[0 0 0 1]);
affine_matrix(1:3,1:3) = eye(3).*repmat(res',1,3).*ratio;
orig = [-nx0/2*res(1)+offsetx -ny0/2*res(2)+offsety -nz0/2*res(3)+offsetz]';
affine_matrix(1:3,4) = orig.*ratio;

allbrain_img(:,:,:,1:22) = [];
subjectnum = 8;

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
    subjectname = filelist_ph(i).name;
    subjectid(i,:) = subjectname(1:6);

    [nx0,ny0,nz0] = size(allbrain_img(:,:,:,i));
    
    [x11,y11,z11] = ind2sub(size(allbrain_img(:,:,:,i)),find(~isnan(allbrain_img(:,:,:,i))));

    index = cat(1,x11',y11',z11',ones(size(x11')));
    world_cordf = affine_matrix*index;
    x11f = world_cordf(1,:)/ratio;
    y11f = world_cordf(2,:)/ratio;
    z11f = world_cordf(3,:)/ratio;

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
    saveas(f1,fullfile(resultprefix,[label,'_viewindex',num2str(view_index),'_check_position.png']));
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
saveas(f2,fullfile(resultprefix,[label,'helmet_brain_coil.png']));

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
slab(:,1) = 52:56;
slab(:,2) = 57:61;
slab(:,3) = 62:66;
slab(:,4) = 67:71;
slab(:,5) = 72:76;
slab(:,6) = 77:81;
for j = 1:6
    for i=1:size(Bz_mapped,4)
        Bz_ROI = Bz_mapped(:,:,slab(:,j),i);
        Bzf_temp (:,i) = Bz_ROI(:);
    end
    for i = 1:subjectnum   
        subjectname = filelist_ph(i).name;
        subjectid = subjectname(1:6);
        % calculate the shimming result
        [improve_STD,improve_mean,DC,B0shim_STD,B0shim_mean,B0_STD,B0_mean] = ...
            CoilShimmingHcp_WorldCordf_slabShimming(resize_allbrain_img(:,:,slab(:,j),i),double(Bzf_temp),DC_limit,resultprefix,subjectid,j);% output the result and save the result
        fprintf("%12s| slab %d |  %8s|     %5.4f| %5.4f|  %5.4f\r\n",datestr(now),j,subjectid,improve_STD,improve_mean,max(DC(:)));
        fprintf(fid,"%12s slab %d    %8s    %5.4f   %5.4f   %5.4f\r\n",datestr(now),j,subjectid,improve_STD,improve_mean,max(DC(:)));
        coilcurrent(i,:,j) = DC(:)';
        improve_STD_record(i,j) = improve_STD;
        improve_mean_record(i,j) = improve_mean;
        B0shim_STD_list(i,j) = B0shim_STD;
        B0_STD_list(i,j) = B0_STD;
        B0shim_mean_list(i,j) = B0shim_mean;
        B0_mean_list(i,j) = B0_mean;
    end
end
save(fullfile(resultprefix,'coilcurrent_allsubject.mat'),'coilcurrent');
save(fullfile(resultprefix,'improve_STD_allsubject.mat'),'improve_STD_record');
save(fullfile(resultprefix,'improve_mean_allsubject.mat'),'improve_mean_record');


tran_improve_STD = improve_STD_record';
tran_improve_mean = improve_mean_record';

max_current = max(abs(coilcurrent),[],2);
mean_current = mean(abs(coilcurrent),2);
for i = 1:6
     result(subjectnum*(i-1)+1:subjectnum*i,:) = cat(2,tran_improve_STD(i,:)',max_current(:,1,i),tran_improve_mean(i,:)',mean_current(:,1,i));    
end
writematrix(result,fullfile(resultprefix,'result.xlsx'));
std_result = cat(2,B0shim_STD_list,B0_STD_list,B0shim_mean_list,B0_mean_list);
writematrix(std_result,fullfile(resultprefix,'std_result.xlsx'));

mean_STD = mean(improve_STD_record, 1);
fprintf(fid,'%12s  improve_STD slab1:%6.4f  slab2:%6.4f  slab3:%6.4f  slab4:%6.4f  slab5:%6.4f  slab6:%6.4f\r\n',datestr(now),mean_STD(1),mean_STD(2),mean_STD(3),mean_STD(4),mean_STD(5),mean_STD(6)); 