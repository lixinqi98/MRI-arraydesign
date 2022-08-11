function [rawNII_ph_masked, Mask] = mask_HCP(filename_mag,filename_ph)
% ======================================================================= %
% mask_HCP.m
% ======================================================================= %
% This function is going to mask the phase of raw data file of HCP, and
% return the .mat file of B0 phase.
%
%
% Input: 
% =======
%  filename_mag     : filename of magnitude nii file, an example could be: 
%                     '/Users/apple/Documents/MATLAB/Human Connectome Project/
%                      Raw data/100408_strc_FieldMap_Magnitude.nii'; 
%  filename_ph      : filename of phases nii file, an example could be: 
%                     '/Users/apple/Documents/MATLAB/Human Connectome Project/
%                      Raw data/100408_strc_FieldMap_Phase.nii'; 
% Output:
% =======
%  rawNII_ph_masked : aligned phase of brain map in rad/s
% 
%


% -------------------------------- Load data ---------------------------- %
rawNII_mag_struct = load_nii(filename_mag);
rawNII_mag = double(rawNII_mag_struct.img(:,:,:,1));
rawNII_ph_struct = load_nii(filename_ph);
rawNII_ph = double(rawNII_ph_struct.img);
rawNII_ph = rawNII_ph/(2*pi); % convert the unit from rad/s to Hz


% ------------------------------ Create mask --------------------------- %
% scale the magnitude by dividing its maximum value
iMag=rawNII_mag/max(rawNII_mag(:));  

n=1;
pixdim = rawNII_ph_struct.hdr.dime.pixdim;       % pixel dimension
voxel_size = [pixdim(2), pixdim(3), pixdim(4)];  % in mm, (was [2,2,2]);
matrix_size = size(iMag);

AllPhasemap(n).Mask=BET(iMag,matrix_size,voxel_size);
%one more time to trim the FOV
%AllPhasemap(n).Mask(:,:,20:end)=0;%!!!!hard code for 1 time use 
[Mi,Mj,Mk]=ind2sub(size(AllPhasemap(n).Mask),find(AllPhasemap(n).Mask));
Maskind=[Mi,Mj,Mk];
tempFOVstart=[min(Mi) min(Mj) min(Mk)];
tempFOVend=[max(Mi) max(Mj) max(Mk)];
tempiMag=iMag(tempFOVstart(1):tempFOVend(1), tempFOVstart(2):tempFOVend(2), tempFOVstart(3):tempFOVend(3));
tempMask=BET(tempiMag,size(tempiMag),voxel_size);
%Mask=zeros(size(iMag));
Mask=NaN(size(iMag));
Mask(tempFOVstart(1):tempFOVend(1), tempFOVstart(2):tempFOVend(2), tempFOVstart(3):tempFOVend(3))=tempMask;
Mask(Mask==0)=NaN;
AllPhasemap(n).Mask=Mask;

% ------------------------- Find the masked B0 phase -------------------- %
rawNII_ph_masked = rawNII_ph .* Mask;
