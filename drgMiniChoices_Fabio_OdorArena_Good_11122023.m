function handles=drgMiniChoices_Fabio_OdorArena_Good_11122023


handles.this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
handles.outputFile='test.mat';

handles.first_file=1;

handles.training_fraction=0.9;


%The user can define what time period to use spikes from (with respect to the output).
handles.bins_before=10;
handles.bins_current=1;
handles.bins_after=0;


%Note: The data brought into the Kording lab jupyter notebbok seems to be
%binned in 200 msec bins
handles.dt=0.2;
handles.dt_miniscope=1/30;
handles.n_shuffle=10;
handles.which_training_algorithm=[2 3];

%Now list the files
handles.dFF_file{1}='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
handles.arena_file{1}='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';
handles.group(1)=1;


handles.dFF_file{2}='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
handles.arena_file{2}='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';
handles.group(2)=1;


handles.dFF_file{3}='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
handles.arena_file{3}='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';
handles.group(3)=2;