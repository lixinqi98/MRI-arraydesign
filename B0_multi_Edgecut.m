function B0_improvement = B0_multi_Edgecut(B0_M,B0_P)
% 2019/06/10 Chris Huang
% This function is used for B0 improvement after alignment
% The alignment_read_in part need only run once
%%
%[B0, mask] = mask_HCP_Downsample(B0_M.name,B0_P.name);%mask
[B0, mask] = mask_HCP(B0_M,B0_P);%mask
%% Morphological Processing to change mask
[row,col,z] = size(B0);
for slice = 1:z
mask1 = mask(:,:,slice);
%figure;
%imagesc(mask1_I);
mask_I = Morphological_Process(mask1,row,col);
B0_O = B0(:,:,slice);
%% B0 improvement
B0_I = B0_O.*mask_I;
%subplot(222);
%imagesc(B0_I); colormap jet;caxis([-100,100]);title(["Edge Cutting(without Thresholding),slice:",num2str(slice)]);%colorbar;
%% Patial Edge Reduction
ite = 2;
B0_I = Patial_Edge_Reduction(B0_I,ite);
%% Thresholding
t_value = -40*2*pi; % Threshold value,I use 50Hz here instead of 40Hz, because 50Hz gives better outcome
Mask_T(:) = B0_I(:) > t_value;
Mask_T = reshape(Mask_T,[row,col]);
B0_I_T = B0_I.*Mask_T;
%subplot(223);
%imagesc(B0_I_T); colormap jet;caxis([-100,100]);title(["Edge Cutting(with Threshold),slice",num2str(slice)]);%colorbar;
%% Further process for scatter points
ite2 = 20;
B0_I_T_S = scatter_point(B0_I_T,slice,ite2,row,col);
%subplot(122);
%imagesc(rot90(B0_I_T_S,-1)); colormap jet;caxis([-100,100]);title(["Futher Edge Cutting(with Threshold),slice",num2str(slice)]);%colorbar;
B0_improvement(:,:,slice) = B0_I_T_S;
end
% for j=1:size(B0_improvement,3)
%     B0_outcome(:,:,j) = rot90(B0_improvement(:,:,j),3); 
% end
B0_improvement(find(B0_improvement==0)) = nan;
end

%% Morphological Processing function
function mask1 = Morphological_Process(mask1,row,col)
if sum((mask1(:)==1))>0.3*row*col
    num = 1;
else 
    num = 4;
end
for i = 1:num
mask1_I = bwmorph(mask1,'remove',1);
mask1 = mask1-mask1_I;
end
%figure;
%imagesc(mask1)
end
%% Patial Edge Reduction function
function B0_I = Patial_Edge_Reduction(B0_I,ite)
for i = 1:ite
mask2 = edge(B0_I,'Sobel',200);
B0_I = B0_I.*(~mask2);
end
end
%% Further process for scatter points function
function B0_I_T = scatter_point(B0_I_T,slice,ite2,row,col)
if slice>=21&&slice<=30 % This range is based on the most of brain pattern 
B0_I_T_C = B0_I_T;
for k = 1:ite2
for i = 2:row-1 % along readout direction
    for j = 2: col-1
        if (B0_I_T_C(i,j-1)==0)&&(B0_I_T_C(1,j+1)==0)&&(B0_I_T(i,j)>20*2*pi)
            B0_I_T(i,j) = 0;
        end
    end
end

for i = 2:row-1 % along phase encoding direction
    for j = 2: col-1
        if (B0_I_T_C(i-1,j)==0)&&(B0_I_T_C(i+1,j)==0)&&(B0_I_T(i,j)>20*2*pi)
            B0_I_T(i,j) = 0;
        end
    end
end

B0_I_T_C = B0_I_T;
end
end
end