%drg_batch_extract.m
close all 
clear all

%These are config settings for Fabio
%Initialize config

config =[];
config = get_defaults(config);

config.trace_output_option='baseline_adjusted'; %'baseline_adjusted', 'raw', 'noneg'
 
% config.avg_cell_radius=5;
% config.num_partitions_x=2;
% config.num_partitions_y=2;
% config.spatial_highpass_cutoff = 5;
% config.downsample_time_by = 6;
% config.thresholds.low_ST_index_thresh = -1;
% config.thresholds.spatial_corrupt_thresh = 3;
% config.cellfind_max_steps = 800;
% config.second_df = [];
% config.max_iter = 10;
% config.trace_quantile = 0.25; % This is because you wanted to see fluctuations, I would keep this ~ 0.25 to not sample too much noise. In the ss below, it is 0.25.
% config.use_gpu=0; %0 for my Mac, 1 for sphgpu
% 
% % change these as needed
% config.cellfind_min_snr = 5; %4
% config.thresholds.T_min_snr=10; %7

config.avg_cell_radius=7;
config.num_partitions_x=1;
config.num_partitions_y=1;
config.spatial_highpass_cutoff = 5;
config.downsample_time_by = 6;
config.thresholds.low_ST_index_thresh = -1;
config.thresholds.spatial_corrupt_thresh = 3;
config.cellfind_max_steps = 800;
config.second_df = [];
config.max_iter = 10;
config.trace_quantile = 0.1; % This is because you wanted to see fluctuations, I would keep this ~ 0.25 to not sample too much noise. In the ss below, it is 0.25.
config.use_gpu=0; %0 for my Mac, 1 for sphgpu

% change these as needed
config.cellfind_min_snr = 5;
config.thresholds.T_min_snr=10;
 

% imageFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210302/20210302_Grin3_homeodor_otherGRIN1_1_XY1614724413_Z0_T0000_C0_motioncorrection.tiff';
% outFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210302/20210302_Grin3_homeodor_otherGRIN1_1_XY1614724413_Z0_T0000_C0_extract2.mat';

% imageFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210302/20210302_Grin1_homeodor_otherGRIN3_XY1614719876_Z0_T00000_C0_motioncorrection.tiff';
% outFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210302/20210302_Grin1_homeodor_otherGRIN3_XY1614719876_Z0_T00000_C0_extract2.mat';

% imageFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210308/20210308_Grin4_fsds_amyl_acetophenone_XY1615243443_Z0_T0000_C0_clean_motioncorrection.tiff';
% outFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210308/20210308_Grin4_fsds_amyl_acetophenone_XY1615243443_Z0_T0000_C0_extract2.mat';

% imageFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210228/20210228_Grin4_fsds1_amyl_acetophen_XY1614548423_Z0_T0000_C0-1_fixed_motioncorrected.tif';
% outFullFileName='/Users/restrepd/Documents/Projects/SFTP/CalmAn_20210228/20210228_Grin4_fsds1_amyl_acetophen_XY1614548423_Z0_T0000_C0-1_extract2.mat';
% 


% %File 1
imageFullFileName{1}='/Users/restrepd/Documents/Projects/Basal_Forebrain/290424/0_fiji_ncorre.mat';
outFullFileName{1}='/Users/restrepd/Documents/Projects/Basal_Forebrain/290424/0_fiji_ncorre_extract.mat';

%Simulation 1
% imageFullFileName{1}='/Users/restrepd/Documents/Projects/SFTP/Simulations/Simulation1.mat';
% outFullFileName{1}='/Users/restrepd/Documents/Projects/SFTP/Simulations/Simulation1ex_bla_hpc15.mat';


% 
% % File 2
% imageFullFileName{2}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20210617\20210617\CalmAn\Grin1_fsds_home_otherPcdh2_XY1623961175_Z0_T00000_C0_ncorr.mat';
% outFullFileName{2}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20210617\20210617\CalmAn\Grin1_fsds_home_otherPcdh2_XY1623961175_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 3
% imageFullFileName{3}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20210803\20210803_FCM2_Home_OtherGrin1_2_XY1628020657_Z0_T00000_C0_ncorr.mat';
% outFullFileName{3}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20210803\20210803_FCM2_Home_OtherGrin1_2_XY1628020657_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 4
% imageFullFileName{4}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210422\20210422_Grin4_home_home_control_XY1619120309_Z0_T00000_C0_ncorr.mat';
% outFullFileName{4}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210422\20210422_Grin4_home_home_control_XY1619120309_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 5
% imageFullFileName{5}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210420\20210420_Grin3_home_home_control_XY1618950608_Z0_T00000_C0_ncorr.mat';
% outFullFileName{5}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210420\20210420_Grin3_home_home_control_XY1618950608_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 6
% imageFullFileName{6}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210415\20210415_Grin1_home_home_control_XY1618521884_Z0_T00000_C0_ncorr.mat';
% outFullFileName{6}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210415\20210415_Grin1_home_home_control_XY1618521884_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 7
% imageFullFileName{7}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210409\20210409_Grin4_home_male_XY1617998540_Z0_T00000_C0_ncorr.mat';
% outFullFileName{7}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210409\20210409_Grin4_home_male_XY1617998540_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 8
% imageFullFileName{8}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210330\20210330_Grin3_fsds_Home_otherMale_XY1617138791_Z0_T00000_C0_ncorr.mat';
% outFullFileName{8}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210330\20210330_Grin3_fsds_Home_otherMale_XY1617138791_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 9
% imageFullFileName{9}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210408\20210408_Grin3_Home_othermale_XY1617911391_Z0_T00000_C0_ncorr.mat';
% outFullFileName{9}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210408\20210408_Grin3_Home_othermale_XY1617911391_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 10
% imageFullFileName{10}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210415\20210415_Grin1_home_othermale_XY1618518579_Z0_T00000_C0_ncorr.mat';
% outFullFileName{10}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210415\20210415_Grin1_home_othermale_XY1618518579_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 11
% imageFullFileName{11}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210402\20210402_Grin1_fsds_home_othermale_XY1617394264_Z0_T00000_C0_ncorr.mat';
% outFullFileName{11}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210402\20210402_Grin1_fsds_home_othermale_XY1617394264_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 12
% imageFullFileName{12}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210409\20210409_Grin4_home_otherGRIN3_KX_XY1618002992_Z0_T00000_C0_ncorr.mat';
% outFullFileName{12}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210409\20210409_Grin4_home_otherGRIN3_KX_XY1618002992_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 13
% imageFullFileName{13}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210408\20210408_Grin3_Home_otheGRIN4_KX_XY1617916254_Z0_T00000_C0_ncorr.mat';
% outFullFileName{13}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210408\20210408_Grin3_Home_otheGRIN4_KX_XY1617916254_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 14
% imageFullFileName{14}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210317\20210317_Grin1_fsds_Home_OtherGrin3_kxanestesia_XY1616018201_Z0_T00000_C0_ncorr.mat';
% outFullFileName{14}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210317\20210317_Grin1_fsds_Home_OtherGrin3_kxanestesia_XY1616018201_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 15
% imageFullFileName{15}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210308\CalmAn_20210308\20210308_Grin4_fsds_homeodor_otherGRIN3_XY1615240058_Z0_T00000_C0_ncorr.mat';
% outFullFileName{15}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210308\20210308_Grin4_fsds_homeodor_otherGRIN3_XY1615240058_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 16
% imageFullFileName{16}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210302\CalmAn_20210302\20210302_Grin1_homeodor_otherGRIN3_XY1614719876_Z0_T00000_C0_ncorr.mat';
% outFullFileName{16}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210302\CalmAn_20210302\20210302_Grin1_homeodor_otherGRIN3_XY1614719876_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 17
% imageFullFileName{17}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210302\CalmAn_20210302\20210302_Grin3_homeodor_otherGRIN1_2_XY1614725527_Z0_T00000_C0_ncorr.mat';
% outFullFileName{17}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210302\CalmAn_20210302\20210302_Grin3_homeodor_otherGRIN1_2_XY1614725527_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 18
% imageFullFileName{18}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210323\20210323_Grin3_fsds_Amyl_Acethop_XY1616532966_Z0_T00000_C0_ncorr.mat';
% outFullFileName{18}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210323\20210323_Grin3_fsds_Amyl_Acethop_XY1616532966_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 19
% imageFullFileName{19}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210317\20210317_Grin1_fsds_AmylSp_AcetopheSm_XY0_Z0_T00000_C0_ncorr.mat';
% outFullFileName{19}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210317\20210317_Grin1_fsds_AmylSp_AcetopheSm_XY0_Z0_T00000_C0_nn_ncorr_ext.mat';
% 
% % File 20
% imageFullFileName{20}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210308\CalmAn_20210308\20210308_Grin4_fsds_amyl_acetophenone_XY1615243443_Z0_T0000_C0_ncorr.mat';
% outFullFileName{20}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\CalmAn_20210308\CalmAn_20210308\20210308_Grin4_fsds_amyl_acetophenone_XY1615243443_Z0_T0000_C0_nn_ncorr_ext.mat';

%% Are the files here?
all_files_here=1;
for ii_file=1:length(imageFullFileName)
    if exist(imageFullFileName{ii_file})~=0
        fprintf(1, ['\nRead file number %d \n'],ii_file);
    else
        all_files_here=0;
        fprintf(1, ['\nFile number %d could not be read\n'],ii_file);
    end
end

%% now run extractor
if all_files_here==1
    for ii_file=1:length(imageFullFileName)
        this_imageFullFileName=imageFullFileName{ii_file};
        M=load(this_imageFullFileName);
        
        tic
   
        %Perform the extraction
        output=[];
        try
            M_in=M.M2; %drg_lp_low_RAM
        catch
            M_in=M.imDataMc; %This is for the EXTRACT version of NoRMCorre drg_batch_NoRMCorreExt.m
        end
        output=extractor(M_in,config);
        % output=extractor('example.h5:/data',config); % If movie is large, do not pre-load. Use this!
        
        toc
        
        %Save the data
        this_outFullFileName=outFullFileName{ii_file};
        save(this_outFullFileName,'output','config','imageFullFileName','outFullFileName','-v7.3')
    end
end

%% Save all figures
% FolderName = tempdir;   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = [FolderName 'Figure' num2str(iFig) '.fig'];
%   savefig(FigHandle, FigName);
% end

% Perform post-processing such as cell checking and further data analysis.
% run drg_cell_check
% Check example_tutorial.m for more in depth tutorial!