close all
clear all
%drg_cell_check
  
% this_imageFullFileName='E:\SFTP\Nikon\20211028-fcm7-MoM\FCM7-MoM_ncorr.mat';
% this_outFullFileName='E:\SFTP\Nikon\20211028-fcm7-MoM\FCM7-MoM_ncorr_ext.mat';

% this_imageFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_H_H_ncorr.mat';
% this_outFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_H_H_ncorr_ext.mat';

% this_imageFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_HR_CHR2_ML_ncorr.mat';
% this_outFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_HR_CHR2_ML_ncorr_ext.mat';
% 
% this_imageFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_HR_OFL_ncorr.mat';
% this_outFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_HR_OFL_ncorr_ext.mat';
% 
% this_imageFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_HR_PCDH2_ML_ncorr.mat';
% this_outFullFileName='E:\SFTP\Nikon\20211105-fcm07-home-odor-task\20211105_FCM7_HR_PCDH2_ML_ncorr_ext.mat';

% this_imageFullFileName='E:\SFTP\ncorr_ext\20210617\Grin1_fsds_home_otherPcdh2_XY1623961175_Z0_T00000_C0_ncorr.mat';
% this_outFullFileName='E:\SFTP\ncorr_ext\20210617\Grin1_fsds_home_otherPcdh2_XY1623961175_Z0_T00000_C0_ncorr_ext.mat';



[this_outFileName,thisoutPathName] = uigetfile({'*.mat'},'Select the EXTRACT output file');

this_outFullFileName=[thisoutPathName this_outFileName];

load(this_outFullFileName)



[this_ImageFileName,thisImagePathName] = uigetfile({'*.mat'},'Select the NorRMCorre file');

this_imageFullFileName=[thisImagePathName this_ImageFileName];


%.mat from NoRMCorre
M=load(this_imageFullFileName);

try
    M_in=M.M2; %drg_lp_low_RAM
catch
    M_in=M.imDataMc; %This is for the EXTRACT version of NoRMCorre drg_batch_NoRMCorreExt.m
end

cell_check(output,M_in)

