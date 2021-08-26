% run three coil settings

% !!! Coil setup, label
% !!! 7T/3T : scale down by 2pi!!
% !!! ? change to plot 7T? 
% !!! ? cut 2 pixels
% !!! current
% 
% 
% 
% clear all
% label = 'iPRES_32_5A_7T_testset_slabShimming_';
% coilposition = 'iPRES_32_coilxyz_360_m.mat';
% DC_limit = 5;
% SlabShimmingShimming_reproduce_Carissa_7T;
% 
% clear all
% label = 'iPRES_32_5A_7T_testset_wholebrain_';
% coilposition = 'iPRES_32_coilxyz_360_m.mat';
% DC_limit = 5;
% Shimming_reproduce_Carissa_7T;
% 
% % 
% clear all
% label = '68to32_30A_7T_testset_wholebrain_';
% coilposition = '68to32_v2_7T_coilxyz_360_m.mat';
% DC_limit = 30;
% Shimming_reproduce_Carissa_7T;
% 
% clear all
% label = '99to53_30A_7T_testset_wholebrain_';
% coilposition = '99to53_v5_7T_coilxyz_360_m.mat';
% Bz_path = '/Users/mona/codeWorkSpace/Cedars-Sinai/May17/output/99to53_30A_7T_trainset_wholebrain_Results/99to53_30A_7T_trainset_wholebrain__Bz_mapped.mat'
% DC_limit = 30;
% Shimming_reproduce_Carissa_7T;
% % 
% label = 'SH2_10A_7T_testset_wholebrain_';
% DC_limit = 10;
% order = 2;
% MainShimming_WorldCordf_SHorder_7T;
% 
% clear all
% label = 'SH3_10A_7T_testset_wholebrain_';
% DC_limit = 10;
% order = 3;
% MainShimming_WorldCordf_SHorder_7T;
% clear all
% label = 'SH4_10A_7T_testset_wholebrain_';
% DC_limit = 10;
% order = 4;
% MainShimming_WorldCordf_SHorder_7T;
% 
% clear all
% label = 'SH5_10A_7T_testset_wholebrain_';
% DC_limit = 10;
% order = 5;
% MainShimming_WorldCordf_SHorder_7T;
% 
% clear all
% label = 'SH6_10A_7T_testset_wholebrain_';
% DC_limit = 10;
% order = 6;
% MainShimming_WorldCordf_SHorder_7T;
% 
% clear all
% label = 'SH7_10A_7T_testset_wholebrain_';
% DC_limit = 10;
% order = 7; 
% MainShimming_WorldCordf_SHorder_7T;
% /Users/mona/codeWorkSpace/Cedars-Sinai/May17/output/SH7_10A_3T_testset_wholebrain_Results
% label = '99to53_3T_15A_v4_slabShimming_testset_Results';
% CalculateWholeBrainSlabImprove;
% 
% label = '99to32_3T_15A_v4_slabShimming_testset_Results';
% CalculateWholeBrainSlabImprove;

% clear all
% label = '99to32_10A_3T_testset_wholebrain_Results';
% CalculateImprove;
% % 
% 
clear all
label = '99to53_3T_15A_v4_slabShimming_testset_Results';
CalculateImprove;
% % 
clear all
label = '99to32_3T_15A_v4_slabShimming_testset_Results';
CalculateImprove;

clear all
label = 'iPRES_32_5A_3T_testset_slabShimming_Results';
CalculateImprove;

clear all 
label = 'SH3_10A_3T_testset_slabShimming_Results';
CalculateImprove;
clear all
label = 'SH4_10A_3T_testset_slabShimming_Results';
CalculateImprove;
clear all
label = 'SH5_10A_3T_testset_slabShimming_Results';
CalculateImprove;
clear all
label = 'SH6_10A_3T_testset_slabShimming_Results';
CalculateImprove;
clear all
label = 'SH7_10A_3T_testset_slabShimming_Results';
CalculateImprove;
