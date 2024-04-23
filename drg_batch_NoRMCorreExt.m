%drg_batch_NoRMCorreExt
clear

warning('off')
 
%Figures are saved in this folder
% FolderName='/home/restrepd/minian/20220824_FCM22/';

% name = '/Users/restrepd/Documents/Projects/Andrew/Area2_hizoom_Z9482_500frames_nostim_XY0_Z0_T000_C0.tiff';

% name='/Users/restrepd/Documents/Projects/Ming/20190116_mmPD21_6f01_OB/2-newodorspm-150um-26pwer_XY1547680461_Z0_T0000_C0.tif'

% name='/Users/restrepd/Documents/Projects/SFTP/20210617/CalmAn/Grin1_fsds_home_otherPcdh2_XY1623961175_Z0_T00000_C0.tiff';

%File 1
% name='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20210621\20210621\CalmAn\20210621_Grin4_fsds_home_pcdh2_XY1624311429_Z0_T00000_C0.tiff';
% out_name='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20210621\20210621\CalmAn\20210621_Grin4_fsds_home_pcdh2_XY1624311429_Z0_T00000_C0_ncorr.mat';

% name='E:\SFTP\20210621\20210621\CalmAn\20210621_Grin4_fsds_home_pcdh2_XY1624311429_Z0_T00000_C0.tiff';
% out_name='E:\SFTP\20210621\20210621\CalmAn\20210621_Grin4_fsds_home_pcdh2_XY1624311429_Z0_T00000_C0_ncorr.mat';
% 
% name{1}='/home/restrepd/minian/20220824_FCM22/20220824_FCM22_012.avi';
% out_name{1}='/home/restrepd/minian/20220824_FCM22/20220824_FCM22_012_ncorre.mat';

name{1}='/Users/restrepd/Documents/Grants/BRAIN Technology R01 June 30th 2023/Figures/20230518_PVG8f001_40Hz.sld - 20230518_PVG8f_PMT60_PWR2_NM920_40HZ_E_30_Iter_3501_output.tif';
out_name{1}='/Users/restrepd/Documents/Grants/BRAIN Technology R01 June 30th 2023/Figures/20230518_PVG8f001_40Hz.sld - 20230518_PVG8f_PMT60_PWR2_NM920_40HZ_E_30_Iter_3501_output_ncorr.mat';


% % File 21
% name{1}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201115\20201115_Grin1_Ca1_spm1_secondday_XY1605476190_Z0_T0000_C0.tiff';
% out_name{1}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201115\20201115_Grin1_Ca1_spm1_secondday_XY1605476190_Z0_T0000_C0_ncorr.mat';
% 
% % File 22
% name{2}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201115\20201115_Grin1_Ca1_spm2_secondday_XY1605477484_Z0_T0000_C0.tiff';
% out_name{2}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201115\20201115_Grin1_Ca1_spm2_secondday_XY1605477484_Z0_T0000_C0_ncorr.mat';

% % File 23
% name{3}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201116\20201116_Grin1_Ca1_spm1_power30_5fms_XY1605565652_Z0_T0000_C0.tiff';
% out_name{3}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201116\20201116_Grin1_Ca1_spm1_power30_5fms_XY1605565652_Z0_T0000_C0_ncorr.mat';
% 
% % File 24
% name{4}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201116\20201116_Grin1_Ca1_spm2_power30_5fms_XY1605566913_Z0_T0000_C0.tiff';
% out_name{4}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201116\20201116_Grin1_Ca1_spm2_power30_5fms_XY1605566913_Z0_T0000_C0_ncorr.mat';
% 
% % File 25
% name{5}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201117\20201117_Grin1_Ca1_spm1_XY1605649383_Z0_T000_C0.tiff';
% out_name{5}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201117\20201117_Grin1_Ca1_spm1_XY1605649383_Z0_T000_C0_ncorr.mat';
% 
% % File 26
% name{6}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201117\20201117_Grin2_Ca1_spm2_XY1605657246_Z0_T0000_C0.tiff';
% out_name{6}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201117\20201117_Grin2_Ca1_spm2_XY1605657246_Z0_T0000_C0_ncorr.mat';
% 
% % File 27
% name{7}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201117\20201117_Grin2_Ca1_spm3_XY1605659478_Z0_T0000_C0.tiff';
% out_name{7}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201117\20201117_Grin2_Ca1_spm3_XY1605659478_Z0_T0000_C0_ncorr.mat';
% 
% % File 28
% name{8}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin1_Ca1_spm2_XY1605807596_Z0_T0000_C0.tiff';
% out_name{8}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin1_Ca1_spm2_XY1605807596_Z0_T0000_C0_ncorr.mat';
% 
% % File 29
% name{9}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin1_Ca1_spm3_XY1605820670_Z0_T0000_C0.tiff';
% out_name{9}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin1_Ca1_spm3_XY1605820670_Z0_T0000_C0_ncorr.mat';
% 
% % File 30
% name{10}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin2_Ca1_spm1_XY1605812694_Z0_T0000_C0.tiff';
% out_name{10}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin2_Ca1_spm1_XY1605812694_Z0_T0000_C0_ncorr.mat';
% 
% % File 31
% name{11}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin2_Ca1_spm3_XY1605816795_Z0_T0000_C0.tiff';
% out_name{11}='S:\CDB\Restrepo Lab\SFTP\Fabio_Good_Data\20201119\20201119_Grin2_Ca1_spm3_XY1605816795_Z0_T0000_C0_ncorr.mat';

all_files_here=1;
for ii_file=1:length(name)
    try
        Y = read_file_nmcorr(name{ii_file},1,1);
        fprintf(1, ['\nRead file number %d \n'],ii_file);
    catch
       all_files_here=0;
       fprintf(1, ['\nFile number %d could not be read\n'],ii_file);
    end
end

config=[];
config = get_defaults_mc(config);
config.use_gpu = 1;
config.file_type='tif';

if all_files_here==1
    figNo=0;
    
    for ii_file=1:length(name)
        input=name{ii_file};
        output=out_name{ii_file};
        
        %Run normcorre pipeline
        tic; drg_run_normcorre_pipeline(input,output,config); toc;        
    
    end
end

