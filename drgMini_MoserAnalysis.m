function handles_out=drgMini_MoserAnalysis(handles_choices)
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020

warning('off')



if exist('handles_choices')==0
    clear all
    close all

    handles_choices.resume_processing=1; %Make this 1 if you want the program to keep all fits already calculated

    handles_choice.is_sphgpu=0;
    is_sphgpu=handles_choice.is_sphgpu;

    % if gpuDeviceCount>0
    %     handles_choice.is_gpu=1;
    %     gpuDevice(1)
    % else
    %     handles_choice.is_gpu=0;
    % end
    % 
    % is_gpu=handles_choice.is_gpu;

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 100 msec bins
    dt=0.1;
    handles_choices.dt=dt;
    dt_miniscope=1/30;
    handles_choices.dt_miniscope=dt_miniscope;
    %Note that n_shuffle is changed to a maximum of ii_n_training
    % n_shuffle=5;
    n_shuffle_SI=50; %number of shuffles for SD normalization of information content
    % handles_choices.n_shuffle=n_shuffle;
    handles_choices.n_shuffle_SI=n_shuffle_SI;

    % handles_choices.distances_in_mm=1; %If the distance is already in mm there will not be morphing of space
    handles_choices.z=1; %Convert dFF to z
    handles_choices.save_tag='deczdFFopt2'; %This will be used in the name of the save file
    %Files for dFF data

    %First troubleshooting file
    if is_sphgpu==1
        this_path='/data/SFTP/PreProcessedDR/20220804_FCM22/';
    else
        this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    end
    dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';

    % if handles_choices.distances_in_mm==0
    %     arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';
    % else
    arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';
    % end
    handles_choices.group=1;


    %Lane 4 only
    % this_path='/data/SFTP/PreProcessedDR/20220526_FCM6_withodor_lane4/';
    % 
    % dFF_file='20220526_FCM6_withodor_miniscope_sync_L4_ncorre_ext_nonneg.mat';
    % 
    % arena_file='20220526_FCM6withodor_odorarena_L4_sync_mm.mat';
    % 
    % handles_choices.group=2;

    %Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4
    %Group 2 is rewarded, with odor lane 4, no odor in lane 1
    %Group 3 is rewarded, with odor lane 1, no odor in lane 4
    %Group 4 is rewarded, with no odor in lane 1 and lane 4


    % handles_choices.process_these_ROIs=[5 8 29 55 98 107 137 138 204];

    % %This was run with handles_choices.save_tag='deczdFFopt'
    % handles_choices.dt_decoding_op=2; %All bins from -dt_decoding to 0 will be included in the decoder
    % %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    % %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    % handles_choices.dt_decoding_xy=0.6;

    %This was run with handles_choices.save_tag='deczdFFopt2'
    handles_choices.dt_decoding_op=1; %All bins from -dt_decoding to 0 will be included in the decoder
    %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    handles_choices.dt_decoding_xy=0.3;

    %Second troubleshooting files
    % this_path='/data/SFTP/PreProcessedDR/20220713_FCM6/';
    % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % 
    % if handles_choices.distances_in_mm==0
    %     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';
    % else
    %     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat'; %Distances converted to mm
    % end
    % 
    % handles_choices.process_these_ROIs=[36 43 53 73 86 93 95 97 106 111 113 117 125 131 146 158 160];
    % 
    % handles_choices.dt_decoding_op=2; %All bins from -dt_decoding to 0 will be included in the decoder
    % %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    % %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    % handles_choices.dt_decoding_xy=0.6;

    %Third troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220727_FCM19/';
    % dFF_file='20220727_FCM19_withodor_miniscope_sync_L1andL4_ncorre_ext_nonneg.mat';
    % arena_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm.mat'; %Distances converted to mm
   
    % handles_choices.dt_decoding_op=2; %All bins from -dt_decoding to 0 will be included in the decoder
    % %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    % %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    % handles_choices.dt_decoding_xy=0.6;


    %No odor troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    % dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';
    % 
    % handles_choices.dt_decoding_op=2; %All bins from -dt_decoding to 0 will be included in the decoder
    % %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    % %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    % handles_choices.dt_decoding_xy=0.6;


    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

    % handles_choices.which_odor_plume=2;   %1=use odor plume from the publication, 2=use odor plume simulated by Aaron
    handles_choices.cm_from_floor=2;

    handles_choices.weber_fechner=1;
    handles_choices.alpha=1;
    handles_choices.multiplier=1;
    %0 is Stevens Law, R proportional to C^alpha
    %1 is Weber-Flechner law R proportional to log(C)
    %See Copelli et al DOI: 10.1103/PhysRevE.65.060901

    handles_choices.algo=7;
    %1 ann
    %2 glm
    %3 tree
    %4 svm with radial basis
    %5 Gaussian process regression
    %6 GPR with 'KernelFunction', 'squaredexponential', 'Standardize', true
    %7 ann optimize hyper parameters, this is VERY slow!




   


    % handles_choices.repeats=1;
    handles_choices.sh_repeats=1; %Number of repeats for shuffled decoding
    % handles_choices.sh_repeats=10;
    % handles_choices.best_repeat=1; %If multiple runs are performed choose the best result

    handles_choices.max_overlap=3; %Maximum overlap of the shuffled segments
    handles_choices.displayFigures=1;




   
    no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
    fprintf(1,['Number of time bins for op decoding ' num2str(no_dec_time_bins_op) '\n\n'])


    
    no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
    fprintf(1,['Number of time bins for xy decoding ' num2str(no_dec_time_bins_xy) '\n\n'])


    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    % handles_choices.n_shuffle=n_shuffle;



else

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=handles_choices.dt;
    dt_miniscope=handles_choices.dt_miniscope;
    % n_shuffle=handles_choices.n_shuffle;

    % distances_in_mm=handles_choices.distances_in_mm;

    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;



    % which_odor_plume=handles_choices.which_odor_plume;   %1=use odor plume from the publication, 2=use odor plume simulated by Aaron
    cm_from_floor=handles_choices.cm_from_floor;

    weber_fechner=handles_choices.weber_fechner;
    alpha=handles_choices.alpha;
    multiplier=handles_choices.multiplier;
    %0 is Stevens Law, R proportional to C^alpha
    %1 is Weber-Flechner law R proportional to log(C)
    %See Copelli et al DOI: 10.1103/PhysRevE.65.060901

    algo=handles_choices.algo;
    %1 ann
    %2 glm
    %3 tree
    %4 svm with radial basis
    %5 Gaussian process regression
    %6 GPR with 'KernelFunction', 'squaredexponential', 'Standardize', true
    %7 ann optimize hyper parameters, this is VERY slow!


    group=handles_choices.group; %1=odor plume, 2=spatial


    save_tag=handles_choices.save_tag; %This will be used in the name of the save file


    % repeats=handles_choices.repeats;

    max_overlap=handles_choices.max_overlap; %Maximum overlap of the shuffled segments
    displayFigures=handles_choices.displayFigures;




    dt_decoding_op=handles_choices.dt_decoding_op; %All bins from -dt_decoding to 0 will be included in the decoder
    %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
    fprintf(1,['Number of time bins for op decoding ' num2str(no_dec_time_bins_op) '\n\n'])


    dt_decoding_xy=handles_choices.dt_decoding_xy;
    no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
    fprintf(1,['Number of time bins for xy decoding ' num2str(no_dec_time_bins_xy) '\n\n'])
    is_sphgpu=handles_choices.is_sphgpu;
    n_shuffle_SI=handles_choices.n_shuffle_SI;

end

if is_sphgpu==1
    addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
    addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
    addpath('/home/restrepd/Documents/MATLAB/drgMaster')
end

max_overlap=handles_choices.max_overlap;

figNo=0;

%Load odor plume and calculate the odor plume for lane 1 and lane 4

%Use Aaron's simulation
%Plot odor plume
if is_sphgpu==1
    load('/data/SFTP/PreProcessedDR/Odor Arena Plumes/odorArenaPlumesDR.mat')
else
    load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumesDR.mat')
end

%Now make mean_plume matrices for lanes 1 and 4
%Determine which lanes have odorants
%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4
%Group 2 is rewarded, with odor lane 4, no odor in lane 1
%Group 3 is rewarded, with odor lane 1, no odor in lane 4
%Group 4 is rewarded, with no odor in lane 1 and lane 4

switch handles_choices.group
    case 1
        %Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4
        handles_choices.lane1_odor_on=1;
        handles_choices.lane4_odor_on=1;
    case 2
        %Group 2 is rewarded, with odor lane 4, no odor in lane 1
        handles_choices.lane1_odor_on=0;
        handles_choices.lane4_odor_on=1;
    case 3
        %Group 3 is rewarded, with odor lane 1, no odor in lane 4
        handles_choices.lane1_odor_on=1;
        handles_choices.lane4_odor_on=0;
    case 4
        %Group 4 is rewarded, with no odor in lane 1 and lane 4
        handles_choices.lane1_odor_on=0;
        handles_choices.lane4_odor_on=0;
end

%Calculate mean for lane 4
this_lane=4; %Note that I had mistakenly made lane 1 close to y=0, this is =1, not 4 even though I am calcualting lane 4
for ii_source=1:length(odor_plumes.source)
    if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==handles_choices.cm_from_floor)
        this_source=ii_source;
    end
end
x_for_plume=10*odor_plumes.source(this_source).x;
y_for_plume=10*(odor_plumes.source(this_source).y-min(odor_plumes.source(this_source).y));
mean_plume_l4=odor_plumes.source(this_source).mean_plume;

%Calculate mean lane 1
this_lane=1; %Note that I had mistakenly made lane 4 close to y=480, this is =4, not 1 even though I am calcualting lane 4
for ii_source=1:length(odor_plumes.source)
    if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==handles_choices.cm_from_floor)
        this_source=ii_source;
    end
end
mean_plume_l1=odor_plumes.source(this_source).mean_plume;
 
mean_plume_l4=mean_plume_l4-min(mean_plume_l4(:));
mean_plume_l1=mean_plume_l1-min(mean_plume_l1(:));
 
min_nonzero=min([min(mean_plume_l4(mean_plume_l4~=0)) min(mean_plume_l1(mean_plume_l1~=0))]);
 
if handles_choices.lane4_odor_on==1
    mean_plume_l4(mean_plume_l4==0)=min_nonzero;
else
    mean_plume_l4(:,:)=min_nonzero;
end

%Shift the plume to 7 cm (20 mm)
if handles_choices.weber_fechner==0
    mean_plume_l4=handles_choices.multiplier*mean_plume_l4.^handles_choices.alpha;
else
    %Weber-Frechner
    mean_plume_l4=handles_choices.multiplier*log10(mean_plume_l4);
end

if handles_choices.lane1_odor_on==1
    mean_plume_l1(mean_plume_l1==0)=min_nonzero;
else
    mean_plume_l1(:,:)=min_nonzero;
end

if handles_choices.weber_fechner==0
    mean_plume_l1=handles_choices.multiplier*mean_plume_l1.^handles_choices.alpha;
else
    %Weber-Frechner
    mean_plume_l1=handles_choices.multiplier*log10(mean_plume_l1);
end

minC=min([min(mean_plume_l4(:)) min(mean_plume_l1(:))]);
maxC=max([max(mean_plume_l4(:)) max(mean_plume_l1(:))]);

if handles_choices.displayFigures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l4')
    colormap fire
    shading interp
    caxis([minC maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')

    title('Mean odor plume lane 4 before shift')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l1')
    colormap fire
    shading interp
    caxis([minC maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')

    title('Mean odor plume lane 1 before shift')
end
 
%Now shift to the actual dimensions of the odorant arena
new_mean_plume_l1=mean_plume_l1(y_for_plume>=20,:);
mean_plume_l1=[];
mean_plume_l1=new_mean_plume_l1;
new_mean_plume_l4=mean_plume_l4(y_for_plume<=480,:);
mean_plume_l4=new_mean_plume_l4;

new_y_for_plume=y_for_plume(:,y_for_plume<=480);
y_for_plume=[];
y_for_plume=new_y_for_plume;

%Now shift the center of lane 4 to 70 mm
new_mean_plume_l4(y_for_plume>=20,:)=mean_plume_l4(y_for_plume<=460,:);
from_y_ii=find(y_for_plume<70,1,'first');
new_mean_plume_l4(length(y_for_plume):-1:from_y_ii,:)=mean_plume_l4((y_for_plume>50)&(y_for_plume<=120),:);
mean_plume_l4=new_mean_plume_l4;

% dy=abs(y_for_plume(2)-y_for_plume(1));
% ii_dy_shift=20/dy;
% mean_plume_l4=zeros(size(new_mean_plume_l4,1),size(new_mean_plume_l4,2));
% mean_plume_l4(ii_dy_shift+1:end,:)=new_mean_plume_l4(1:size(new_mean_plume_l4,1)-ii_dy_shift,:);
% ii_y_center=find(y_for_plume==410);
% mean_plume_l4(1:ii_y_center-1,:)=mean_plume_l4(ii_y_center+ii_y_center-1:-1:ii_y_center+1,:);
% 

if handles_choices.displayFigures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

 
    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l4')
    colormap fire
    shading interp
    caxis([minC maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')

    title('Mean odor plume shifted to lane 4')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l1')
    colormap fire
    shading interp
    caxis([minC maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')

    title('Mean odor plume shifted to lane 1')
end

pffft=1;
% end

%Restart random seeds
rng('shuffle');


dFF=[];
if strcmp(dFF_file(end-3:end),'.mat')
    %This reads the extract file
    load([this_path dFF_file])
    dFF=zeros(size(output.temporal_weights,1),size(output.temporal_weights,2));
    for traceNo=1:size(output.temporal_weights,2)
        dFF(:,traceNo)=output.temporal_weights(:,traceNo);
    end
else
    if strcmp(dFF_file(end-4:end),'.hdf5')
        %This is an hdf5 generated by CaImAn
        dFF=h5read([this_path dFF_file],'/estimates/F_dff' );
    else
        %This is a csv file created from ImageJ
        dFF=readmatrix([this_path dFF_file]);
    end
end

%Transform to z if desired
if handles_choices.z==1
    for trNo=1:size(dFF,2)
        this_dFF=zeros(size(dFF,1),1);
        this_dFF(:,1)=dFF(:,trNo);
        this_dFF=this_dFF-min(this_dFF);
        this_dFF=this_dFF/std(this_dFF);
        dFF(:,trNo)=this_dFF;
    end
end
% dFF=readmatrix([this_path dFF_file]); %Timepoints x ROIs

load([this_path arena_file])
 
%Note that Ryan divided by 4
% if handles_choices.distances_in_mm~=1
%     arena.xsync=4*arena.xsync;
%     arena.ysync=4*arena.ysync;
% end

%Extract trials
trials=[];


%Extract odor on using the camera sync
at_end=0;
ii=0;
jj=0;
jj_l1=0;
jj_l4=0;
while at_end==0
    next_ii=find(arena.odorsync(ii+1:end)==1,1,'first');
    if ~isempty(next_ii)
        jj=jj+1;
        trials.ii_odor(jj)=ii+next_ii;
        trials.x_odor(jj)=arena.xsync(ii+next_ii);
        trials.y_odor(jj)=arena.ysync(ii+next_ii);

        ii=ii+next_ii;
        ii_mini=arena.index_flirsynctominiscope(ii);

        if sum(arena.laneodor1(ii_mini-5:ii_mini+5)==1)>0
            %Note: laneodor4 is 1 only for one time point
            jj_l1=jj_l1+1;
            trials.ii_laneodor1(jj_l1)=ii;
            trials.x_laneodor1(jj_l1)=trials.x_odor(jj);
            trials.y_laneodor1(jj_l1)=trials.y_odor(jj);
        end

        if sum(arena.laneodor4(ii_mini-5:ii_mini+5)==1)>0
            %Note: laneodor4 is 1 only for one time point
            jj_l4=jj_l4+1;
            trials.ii_laneodor4(jj_l4)=ii;
            trials.x_laneodor4(jj_l4)=trials.x_odor(jj);
            trials.y_laneodor4(jj_l4)=trials.y_odor(jj);
        end

        next_ii=find(arena.odorsync(ii+1:end)==0,1,'first');
        if ~isempty(next_ii)
            ii=ii+next_ii;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end

trials.jj_l1=jj_l1;
trials.jj_l4=jj_l4;

%Extract lanewater1
at_end=0;
ii_flir=0;
jj=0;
while at_end==0
    next_ii_flir=find(arena.lanewater1(ii_flir+1:end)==1,1,'first');
    if ~isempty(next_ii_flir)
        next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
        jj=jj+1;
        trials.ii_lanewater1(jj)=next_ii;
        trials.x_lanewater1(jj)=arena.xsync(next_ii);
        trials.y_lanewater1(jj)=arena.ysync(next_ii);
        next_ii_flir2=find(arena.lanewater1(ii_flir+next_ii_flir+1:end)==0,1,'first');
        %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
        if ~isempty(next_ii_flir2)
            ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end


%Extract lanewater4
at_end=0;
ii_flir=0;
jj=0;
while at_end==0
    next_ii_flir=find(arena.lanewater4(ii_flir+1:end)==1,1,'first');
    if ~isempty(next_ii_flir)
        next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
        jj=jj+1;
        trials.ii_lanewater4(jj)=next_ii;
        trials.x_lanewater4(jj)=arena.xsync(next_ii);
        trials.y_lanewater4(jj)=arena.ysync(next_ii);
        next_ii_flir2=find(arena.lanewater4(ii_flir+next_ii_flir+1:end)==0,1,'first');
        %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
        if ~isempty(next_ii_flir2)
            ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end


if handles_choices.displayFigures==1
    %Display the location of trial start and water delivery
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on

    if trials.jj_l1>0
        plot(trials.x_laneodor1,trials.y_laneodor1,'or')
        plot(trials.x_lanewater1,trials.y_lanewater1,'xr')
    end

    if trials.jj_l4>0
        plot(trials.x_laneodor4,trials.y_laneodor4,'ob')
        plot(trials.x_lanewater4,trials.y_lanewater4,'xb')
    end

    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Location of trial start (o) and water delivery (x)')


    %Plot mouse pathway
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on
    % these_no_points=150;
    % plot(arena.xsync(end),arena.ysync(end),'ok')
    % plot(arena.xsync(end-these_no_points:end),arena.ysync(end-these_no_points:end),'-b')
    plot(arena.xsync,arena.ysync,'-b')
    set(gca, 'YDir', 'reverse');
    % title(['End of trajectory for ' arena_file])
    title(['Trajectory for ' arena_file])
    ylabel('y (mm)')
    xlabel('x (mm)')
    % legend('end of trajectory','trajectory')

    %Display the lane 1 and lane 4 trials
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on
    seconds_to_plot=5;
    ii_seconds_to_plot=ceil(seconds_to_plot/dt_miniscope);
    %Lane 1
    if trials.jj_l1>0
        for iil1=1:length(trials.ii_laneodor1)
            plot(arena.xsync(trials.ii_laneodor1(iil1):trials.ii_laneodor1(iil1)+ii_seconds_to_plot),arena.ysync(trials.ii_laneodor1(iil1):trials.ii_laneodor1(iil1)+ii_seconds_to_plot),'-b')
            plot(arena.xsync(trials.ii_laneodor1(iil1)),arena.ysync(trials.ii_laneodor1(iil1)),'ob')
            plot(arena.xsync(trials.ii_laneodor1(iil1)+ii_seconds_to_plot),arena.ysync(trials.ii_laneodor1(iil1)),'xb')
            pffft=1;
        end
    end

    %Lane 4
    if trials.jj_l4>0
        for iil4=1:length(trials.ii_laneodor4)
            plot(arena.xsync(trials.ii_laneodor4(iil4):trials.ii_laneodor4(iil4)+ii_seconds_to_plot),arena.ysync(trials.ii_laneodor4(iil4):trials.ii_laneodor4(iil4)+ii_seconds_to_plot),'-r')
            plot(arena.xsync(trials.ii_laneodor4(iil4)),arena.ysync(trials.ii_laneodor4(iil4)),'or')
            plot(arena.xsync(trials.ii_laneodor4(iil4)+ii_seconds_to_plot),arena.ysync(trials.ii_laneodor4(iil4)),'xr')
            pffft=1;
        end
    end
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Lane 1 and 4 trajectories lane 1 (blue) and lane 4 (red)')

end

% if handles_choices.distances_in_mm==1
    %The distances are already in mm

    %Bin positions into dt time bins
    pos=[];
    pos(:,1)=arena.xsync;
    pos(:,2)=arena.ysync;
    no_time_points=size(pos,1);

% else
%     %The physical odor arena is 50 cm along the air flow and 48 cm wide,
%     %with lanes 1 and 4, 5 cm from the wall
%     %We need to morph the video camera locations to the physical location
%     %For now I will do this isometrically
%     y_mean_video_lane1=mean(trials.y_lanewater1);
%     y_mean_video_lane4=mean(trials.y_lanewater4);
%     y_mean_video_lanes14=mean([y_mean_video_lane1 y_mean_video_lane4]);
%     delta_y_physical_lane14=380;
%     y_half_physical=480/2;
% 
% 
%     %Bin positions into dt time bins
%     pos=[];
%     pos(:,1)=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(arena.xsync-95);
%     pos(:,2)=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(arena.ysync-y_mean_video_lanes14)+y_half_physical;
%     no_time_points=size(pos,1);
% 
%     %Plot morphed end locations
%     if handles_choices.displayFigures==1
%         %Display the location of trial start and water delivery
%         figNo=figNo+1;
%         try
%             close(figNo)
%         catch
%         end
% 
%         hFig = figure(figNo);
% 
%         set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% 
%         hold on
% 
%         morphed_x1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_laneodor1-95);
%         morphed_y1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_laneodor1-y_mean_video_lanes14)+y_half_physical;
%         plot(morphed_x1,morphed_y1,'or')
% 
%         morphed_x4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_laneodor4-95);
%         morphed_y4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_laneodor4-y_mean_video_lanes14)+y_half_physical;
%         plot(morphed_x4,morphed_y4,'ob')
% 
%         % plot([95 95],[200 250],'-r')
%         % plot([95 95],[0 50],'-b')
% 
%         morphed_x1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_lanewater1-95);
%         morphed_y1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_lanewater1-y_mean_video_lanes14)+y_half_physical;
%         plot(morphed_x1,morphed_y1,'xr')
% 
%         morphed_x4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_lanewater4-95);
%         morphed_y4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_lanewater4-y_mean_video_lanes14)+y_half_physical;
%         plot(morphed_x4,morphed_y4,'xb')
%         set(gca, 'YDir', 'reverse');
%         xlabel('x (mm)')
%         ylabel('y (mm)')
%         title('Location of trial start (o) and water delivery (x), x and y morphed')
%         %Note: This is clearly not correct, the mice move to 25 mm and 225 mm
%         %This is half as large as the physical odor arena!!
% 
% 
%     end
% 
% end

dFF_times=[1:no_time_points]*dt_miniscope;

no_neurons=size(dFF,2);
if ~isfield(handles_choices,'process_these_ROIs')
    handles_choices.process_these_ROIs=[1:no_neurons];
end
process_these_ROIs=handles_choices.process_these_ROIs;
if ~isfield(handles_choices,'resume_processing')
    resume_processing=0;
else
    resume_processing=handles_choices.resume_processing;
end

 
if resume_processing==1
    load([this_path arena_file(1:end-4) '_' handles_choices.save_tag '_' num2str(handles_choices.algo) ...
    num2str(handles_choices.weber_fechner) num2str(handles_choices.alpha)...
    num2str(no_dec_time_bins_op) num2str(no_dec_time_bins_xy) '.mat'])
end

handles_choices.resume_processing=resume_processing;
handles_choices.process_these_ROIs=process_these_ROIs;

no_time_bins=ceil(dFF_times(end)/dt);
time_binned=[1:no_time_bins]*dt-dt/2;
neural_data=zeros(no_time_bins,no_neurons);
pos_binned=zeros(no_time_bins,2);

for ii_time_bin=1:no_time_bins
    time_from=time_binned(ii_time_bin)-dt/2;
    time_to=time_binned(ii_time_bin)+dt/2;
    pos_binned(ii_time_bin,:)=mean(pos((dFF_times>=time_from)&(dFF_times<time_to),:),1);
    for ii_neuron=1:no_neurons
        neural_data(ii_time_bin,ii_neuron)=mean(dFF((dFF_times>=time_from)&(dFF_times<time_to),ii_neuron),1);
    end
end

%Now calculate the behavioral performance
trials.hit1=0;
trials.miss1=0;
trials.lane1=0;
trials.lane4=0;
trials.hit4=0;
trials.miss4=0;
trials.odor_trNo=0;

trim_factor=no_time_bins/no_time_points;

dii_bin_trial=[];
dii_trial=[];
for trNo=1:length(trials.ii_odor)
    if trials.jj_l1>0
        this_ii_laneodor1=find(abs(trials.ii_odor(trNo)-trials.ii_laneodor1)<3);
        if ~isempty(this_ii_laneodor1)
            if trNo==length(trials.ii_odor)
                ii_next=length(arena.xsync);
            else
                ii_next=trials.ii_odor(trNo+1);
            end
            this_water=find((trials.ii_lanewater1>trials.ii_odor(trNo))&(trials.ii_lanewater1<ii_next));
            if ~isempty(this_water)
                trials.hit1=trials.hit1+1;
                trials.hit1_ii_bin_start(trials.hit1)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.hit1_ii_bin_end(trials.hit1)=ceil(trim_factor*trials.ii_lanewater1(this_water));
                trials.hit1_ii_start(trials.hit1)=trials.ii_odor(trNo);
                trials.hit1_ii_end(trials.hit1)=trials.ii_lanewater1(this_water);
                trials.lane1=trials.lane1+1;
                trials.lane1_ii_bin_start(trials.lane1)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.lane1_ii_bin_end(trials.lane1)=ceil(trim_factor*trials.ii_lanewater1(this_water));
                trials.lane1_ii_start(trials.lane1)=trials.ii_odor(trNo);
                trials.lane1_ii_end(trials.lane1)=trials.ii_lanewater1(this_water);
                dii_bin_trial=[dii_bin_trial trials.hit1_ii_bin_end(trials.hit1)-trials.hit1_ii_bin_start(trials.hit1)];
                dii_trial=[dii_trial trials.hit1_ii_end(trials.hit1)-trials.hit1_ii_start(trials.hit1)];
                trials.odor_trNo=trials.odor_trNo+1;
                trials.odor_ii_bin_start(trials.odor_trNo)=trials.hit1_ii_bin_start(trials.hit1);
                trials.odor_ii_bin_end(trials.odor_trNo)=trials.hit1_ii_bin_end(trials.hit1);
                trials.odor_ii_start(trials.odor_trNo)=trials.hit1_ii_start(trials.hit1);
                trials.odor_ii_bin_end(trials.odor_trNo)=trials.hit1_ii_bin_end(trials.hit1);
                trials.odor_trial_type(trials.odor_trNo)=1;
                trials.odor_lane(trials.odor_trNo)=1;
                pfft=1;
            else
                trials.miss1=trials.miss1+1;
                trials.miss1_ii_bin_start(trials.miss1)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.miss1_ii_start(trials.miss1)=trials.ii_odor(trNo);
                trials.lane1=trials.lane1+1;
                trials.lane1_ii_bin_start(trials.lane1)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.lane1_ii_start(trials.lane1)=trials.ii_odor(trNo);
                trials.lane1_ii_bin_end(trials.lane1)=-1;
                trials.odor_trNo=trials.odor_trNo+1;
                trials.odor_ii_bin_start(trials.odor_trNo)=trials.miss1_ii_bin_start(trials.miss1);
                trials.odor_ii_start(trials.odor_trNo)=trials.miss1_ii_start(trials.miss1);
                trials.odor_trial_type(trials.odor_trNo)=2;
                trials.odor_lane(trials.odor_trNo)=1;
                pfft=1;
            end
        end
    end

    if trials.jj_l4>0
        this_ii_laneodor4=find(abs(trials.ii_odor(trNo)-trials.ii_laneodor4)<3);
        if ~isempty(this_ii_laneodor4)
            if trNo==length(trials.ii_odor)
                ii_next=length(arena.xsync);
            else
                ii_next=trials.ii_odor(trNo+1);
            end
            this_water=find((trials.ii_lanewater4>trials.ii_odor(trNo))&(trials.ii_lanewater4<ii_next));
            if ~isempty(this_water)
                trials.hit4=trials.hit4+1;
                trials.hit4_ii_bin_start(trials.hit4)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.hit4_ii_bin_end(trials.hit4)=ceil(trim_factor*trials.ii_lanewater4(this_water));
                trials.hit4_ii_start(trials.hit4)=trials.ii_odor(trNo);
                trials.hit4_ii_end(trials.hit4)=trials.ii_lanewater4(this_water);
                trials.lane4=trials.lane4+1;
                trials.lane4_ii_bin_start(trials.lane4)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.lane4_ii_bin_end(trials.lane4)=ceil(trim_factor*trials.ii_lanewater4(this_water));
                trials.lane4_ii_start(trials.lane4)=trials.ii_odor(trNo);
                trials.lane4_ii_end(trials.lane4)=trials.ii_lanewater4(this_water);
                dii_bin_trial=[dii_bin_trial trials.hit4_ii_bin_end(trials.hit4)-trials.hit4_ii_bin_start(trials.hit4)];
                dii_trial=[dii_trial trials.hit4_ii_end(trials.hit4)-trials.hit4_ii_start(trials.hit4)];
                trials.odor_trNo=trials.odor_trNo+1;
                trials.odor_ii_bin_start(trials.odor_trNo)=trials.hit4_ii_bin_start(trials.hit4);
                trials.odor_ii_bin_end(trials.odor_trNo)=trials.hit4_ii_bin_end(trials.hit4);
                trials.odor_ii_start(trials.odor_trNo)=trials.hit4_ii_start(trials.hit4);
                trials.odor_ii_end(trials.odor_trNo)=trials.hit4_ii_end(trials.hit4);
                trials.odor_trial_type(trials.odor_trNo)=3;
                trials.odor_lane(trials.odor_trNo)=4;
                pfft=1;
            else
                trials.miss4=trials.miss4+1;
                trials.miss4_ii_bin_start(trials.miss4)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.miss4_ii_start(trials.miss4)=trials.ii_odor(trNo);
                trials.lane4=trials.lane4+1;
                trials.lane4_ii_bin_start(trials.lane4)=ceil(trim_factor*trials.ii_odor(trNo));
                trials.lane4_ii_start(trials.lane4)=trials.ii_odor(trNo);
                trials.lane4_ii_bin_end(trials.lane4)=-1;
                trials.lane4_ii_end(trials.lane4)=-1;
                trials.odor_trNo=trials.odor_trNo+1;
                trials.odor_ii_bin_start(trials.odor_trNo)=trials.miss4_ii_bin_start(trials.miss4);
                trials.odor_ii_start(trials.odor_trNo)=trials.miss4_ii_start(trials.miss4);
                trials.odor_trial_type(trials.odor_trNo)=4;
                trials.odor_lane(trials.odor_trNo)=4;
                pffft=1;
            end
        end
    end
end

for ii_miss=1:trials.miss1
    trials.miss1_ii_bin_end(ii_miss)=trials.miss1_ii_bin_start(ii_miss)+ceil(mean(dii_bin_trial));
    trials.miss1_ii_end(ii_miss)=trials.miss1_ii_start(ii_miss)+ceil(mean(dii_trial));
end

for ii_miss=1:trials.miss4
    trials.miss4_ii_bin_end(ii_miss)=trials.miss4_ii_bin_start(ii_miss)+ceil(mean(dii_bin_trial));
    trials.miss4_ii_end(ii_miss)=trials.miss4_ii_start(ii_miss)+ceil(mean(dii_trial));
end

for ii_odor=1:trials.odor_trNo
    if trials.odor_trial_type(ii_odor)==4
        trials.odor_ii_bin_end(ii_odor)=trials.odor_ii_bin_start(ii_odor)+ceil(mean(dii_bin_trial));
        trials.odor_ii_end(ii_odor)=trials.odor_ii_start(ii_odor)+ceil(mean(dii_trial));
    end
    if trials.odor_trial_type(ii_odor)==2
        trials.odor_ii_bin_end(ii_odor)=trials.odor_ii_bin_start(ii_odor)+ceil(mean(dii_bin_trial));
        trials.odor_ii_end(ii_odor)=trials.odor_ii_start(ii_odor)+ceil(mean(dii_trial));
    end
end

if trials.jj_l1>0
    for iil1=1:length(trials.lane1_ii_bin_end)
        if trials.lane1_ii_bin_end(iil1)==-1
            trials.lane1_ii_bin_end(iil1)=trials.lane1_ii_bin_start(iil1)+ceil(mean(dii_bin_trial));
            trials.lane1_ii_end(iil1)=trials.lane1_ii_start(iil1)+ceil(mean(dii_trial));
        end
    end
end

if trials.jj_l4>0
    for iil1=1:length(trials.lane4_ii_bin_end)
        if trials.lane4_ii_bin_end(iil1)==-1
            trials.lane4_ii_bin_end(iil1)=trials.lane4_ii_bin_start(iil1)+ceil(mean(dii_bin_trial));
            trials.lane4_ii_end(iil1)=trials.lane4_ii_start(iil1)+ceil(mean(dii_trial));
        end
    end
end

percent_correct=100*(trials.hit4+trials.hit1)/(trials.hit4+trials.hit1+trials.miss4+trials.miss1);
percent_correct1=100*(trials.hit1)/(trials.hit1+trials.miss1);
percent_correct4=100*(trials.hit4)/(trials.hit4+trials.miss4);
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct lane 1 ' num2str(percent_correct1) ' percent correct lane 4 ' num2str(percent_correct4) '\n\n'])

lane1_trials=trials.hit1+trials.miss1;
lane4_trials=trials.hit4+trials.miss4;

[p_hit_miss, Q]= chi2test([trials.hit1, trials.miss1; trials.hit4, trials.miss4]);
fprintf(1,['Chi squared testing difference in hit vs miss between lane 1 and 4 ' num2str(p_hit_miss) '\n\n'])

[p_hit_miss_lane1, Q]= chi2test([trials.hit1, trials.miss1; (trials.hit1+trials.miss1)/2, (trials.hit1+trials.miss1)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for lane 1 ' num2str(p_hit_miss_lane1) '\n\n'])

[p_hit_miss_lane4, Q]= chi2test([trials.hit4, trials.miss4; (trials.hit4+trials.miss4)/2, (trials.hit4+trials.miss4)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for lane 4 ' num2str(p_hit_miss_lane4) '\n\n'])

[p_hit_miss_both_lanes, Q]= chi2test([trials.hit1+trials.hit4, trials.miss1+trials.miss4; (trials.hit4+trials.hit1+trials.miss4+trials.miss1)/2, (trials.hit4+trials.hit1+trials.miss4+trials.miss1)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for both lanes ' num2str(p_hit_miss_both_lanes) '\n\n'])


if handles_choices.displayFigures==1

    %Display the lane 1 trials
    if trials.jj_l1>0
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

        hold on

        %These are hits
        for iil1=1:trials.hit1
            plot(pos_binned(trials.hit1_ii_bin_start(iil1):trials.hit1_ii_bin_end(iil1),1),pos_binned(trials.hit1_ii_bin_start(iil1):trials.hit1_ii_bin_end(iil1),2),'-r')
            plot(pos_binned(trials.hit1_ii_bin_start(iil1),1),pos_binned(trials.hit1_ii_bin_start(iil1),2),'or')
            plot(pos_binned(trials.hit1_ii_bin_end(iil1),1),pos_binned(trials.hit1_ii_bin_end(iil1),2),'xr')
        end

        %These are miss
        for iil1=1:trials.miss1
            plot(pos_binned(trials.miss1_ii_bin_start(iil1):trials.miss1_ii_bin_end(iil1),1),pos_binned(trials.miss1_ii_bin_start(iil1):trials.miss1_ii_bin_end(iil1),2),'-b')
            plot(pos_binned(trials.miss1_ii_bin_start(iil1),1),pos_binned(trials.miss1_ii_bin_start(iil1),2),'ob')
            plot(pos_binned(trials.miss1_ii_bin_end(iil1),1),pos_binned(trials.miss1_ii_bin_end(iil1),2),'xb')
        end

        set(gca, 'YDir', 'reverse');
        xlabel('x (mm)')
        ylabel('y (mm)')
        title('Lane 1 trajectories hit (red) and miss (blue)')
    end

    %Display the lane 4 trials
    if trials.jj_l4>0
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

        hold on

        %These are hits
        for iil4=1:trials.hit4
            plot(pos_binned(trials.hit4_ii_bin_start(iil4):trials.hit4_ii_bin_end(iil4),1),pos_binned(trials.hit4_ii_bin_start(iil4):trials.hit4_ii_bin_end(iil4),2),'-r')
            plot(pos_binned(trials.hit4_ii_bin_start(iil4),1),pos_binned(trials.hit4_ii_bin_start(iil4),2),'or')
            plot(pos_binned(trials.hit4_ii_bin_end(iil4),1),pos_binned(trials.hit4_ii_bin_end(iil4),2),'xr')
        end

        %These are miss
        for iil4=1:trials.miss4
            plot(pos_binned(trials.miss4_ii_bin_start(iil4):trials.miss4_ii_bin_end(iil4),1),pos_binned(trials.miss4_ii_bin_start(iil4):trials.miss4_ii_bin_end(iil4),2),'-b')
            plot(pos_binned(trials.miss4_ii_bin_start(iil4),1),pos_binned(trials.miss4_ii_bin_start(iil4),2),'ob')
            plot(pos_binned(trials.miss4_ii_bin_end(iil4),1),pos_binned(trials.miss4_ii_bin_end(iil4),2),'xb')
        end

        set(gca, 'YDir', 'reverse');
        xlabel('x (mm)')
        ylabel('y (mm)')
        title('Lane 4 trajectories hit (red) and miss (blue)')
    end
end

% nan_mask=logical(ones(1,no_time_bins));
% for ii_neuron=1:no_neurons
%     this_nan_mask=ones(1,no_time_bins);
%     this_nan_mask(1,:)=~isnan(neural_data(:,ii_neuron));
%     nan_mask=nan_mask&this_nan_mask;
% end
% neural_data_trimmed=neural_data(nan_mask,:);
% pos_binned_trimmed=pos_binned(nan_mask,:);
% no_time_bins=sum(nan_mask);
no_time_bins=size(neural_data,1);

% %Do z scores
% mean_neural_data_col=mean(neural_data,1);
% mean_neural_data=repmat(mean_neural_data_col,no_time_bins,1);
% 
% std_neural_data_col=std(neural_data,1);
% std_neural_data=repmat(std_neural_data_col,no_time_bins,1);
% 
% neural_data=(neural_data-mean_neural_data)./std_neural_data;

pffft=1;

% % Format for Wiener Filter, Wiener Cascade, XGBoost, and Dense Neural Network
% % Put in "flat" format, so each "neuron / time" is a single feature
% % i.e. each time point in the before and after window becomes a different
% % "neuron"
% all_bins_per_window=bins_before+bins_after+bins_current;
% no_X_dFF_neurons=no_neurons*all_bins_per_window;
% X_dFF=zeros(no_time_bins,no_X_dFF_neurons);
%
% for ii_t=bins_before+1:no_time_bins-bins_after
%     ii_n=0;
%     for no_win=1:all_bins_per_window
%         ii_this_t=ii_t-bins_before+no_win-1;
%         X_dFF(ii_t,ii_n+1:ii_n+no_neurons)=neural_data(ii_this_t,:);
%         ii_n=ii_n+no_neurons;
%     end
% end
% 
% %Now do the neural networks
% tic
% 
% 
% %Train with per trial data using a leave one out approach
% if resume_processing==0
%     per_ROI=[];
%     for this_ROI=handles_choices.process_these_ROIs
%         per_ROI(this_ROI).results=[];
%     end
% else
%     per_ROI=handles_out.per_ROI;
% 
% end
% 
% if resume_processing==0
%     per_ROI_sh=[];
%     for this_ROI=handles_choices.process_these_ROIs
%         per_ROI_sh(this_ROI).repeats=[];
%     end
% else
%     per_ROI_sh=handles_out.per_ROI_sh;
% end
% 
for this_ROI=handles_choices.process_these_ROIs
%     %We have this while try catch loop here because we were getting this error:
%     %         This parallel pool has been shut down.
%     %
%     % Caused by:
%     %     The parallel pool shut down because the client lost connection to worker 6. Check the network
%     %     connection or restart the parallel pool with 'parpool'.
%     % if resume_processing==0
%     %     per_ROI(this_ROI).results=[];
%     % end
%     if isempty(per_ROI(this_ROI).results)
%         thisROI_done=0;
%         while thisROI_done==0
%             try
% 
%                 %training_range_template has all the trials
%                 training_range_template=zeros(1,no_time_bins);
%                 lane_template=zeros(1,no_time_bins);
%                 odor_plume_template=zeros(1,no_time_bins);
%                 for trNo=1:trials.odor_trNo
%                     dFF_pred(trNo).dFFpred_xy=[];
%                     % dFF_pred(trNo).dFFpred_xyl=[];
%                     % dFF_pred(trNo).dFFpred_xyop=[];
%                     dFF_pred(trNo).dFFpred_op=[];
%                     dFF_pred(trNo).dFF=[];
%                     dFF_pred(trNo).dFFpred_xyop=[];
% 
%                     % x_pred(trNo).data=[];
% 
%                     x_predictedstart=trials.odor_ii_bin_start(trNo);
%                     x_predictedend=trials.odor_ii_bin_end(trNo);
%                     training_range_template(x_predictedstart:x_predictedend)=1;
%                     if trials.odor_lane(trNo)==1
%                         lane_template(x_predictedstart:x_predictedend)=1;
%                         for x_ii=x_predictedstart:x_predictedend
%                             this_x=pos_binned(x_ii,1);
%                             this_y=pos_binned(x_ii,2);
%                             [minabx,ii_minax]=min(abs(x_for_plume-this_x));
%                             [minaby,ii_minay]=min(abs(y_for_plume-this_y));
%                             this_ca=mean_plume_l1(ii_minay,ii_minax);
%                             odor_plume_template(x_ii)=this_ca;
%                         end
%                     else
%                         lane_template(x_predictedstart:x_predictedend)=4;
%                         for x_ii=x_predictedstart:x_predictedend
%                             this_x=pos_binned(x_ii,1);
%                             this_y=pos_binned(x_ii,2);
%                             [minabx,ii_minax]=min(abs(x_for_plume-this_x));
%                             [minaby,ii_minay]=min(abs(y_for_plume-this_y));
%                             this_ca=mean_plume_l4(ii_minay,ii_minax);
%                             odor_plume_template(x_ii)=this_ca;
%                         end
%                     end
%                 end
% 
%                 start_toc=toc;
%                 % Create the parallel pool once before the main loop
%                 % if isempty(gcp('nocreate'))
%                 % parpool;  % This will use the default number of workers
%                 % end
%                 % parfor trNo=1:trials.odor_trNo
                for trNo=1:trials.odor_trNo
% 
                    this_test_range=zeros(1,no_time_bins);
%                     % if trNo==1
%                     %     ii_test_range_start=1;
%                     % else
%                     %     ii_test_range_start=trials.odor_ii_bin_end(trNo-1)+15;
%                     % end
% 
                    ii_test_range_start=trials.odor_ii_bin_start(trNo);
                    ii_test_range_end=trials.odor_ii_bin_end(trNo);
% 
%                     % if trNo==trials.odor_trNo
%                     %     ii_test_range_end=no_time_bins;
%                     % else
%                     %     ii_test_range_end=trials.odor_ii_bin_end(trNo)+15;
%                     % end
% 
                    this_test_range(ii_test_range_start:ii_test_range_end)=1;
%                     this_training_range=logical(training_range_template)&(~logical(this_test_range));
% 
%                     % ii_valid_range=ceil(valid_range*no_time_bins);
%                     %     ii_test_range=ceil(test_range*no_time_bins);
%                     this_dFFtrain=zeros(sum(this_training_range),1);
%                     this_dFFtrain=neural_data(this_training_range,this_ROI);
%                     % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
                    this_dFFtest=zeros(sum(this_test_range),1);
                    this_dFFtest(:,1)=neural_data(logical(this_test_range),this_ROI);
                    dFF_pred(trNo).dFF=this_dFFtest;
% 
%                     no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
%                     XYtrain=zeros(sum(this_training_range),size(pos_binned,2)*no_dec_time_bins_xy);
%                     XYtest=zeros(sum(logical(this_test_range)),size(pos_binned,2)*no_dec_time_bins_xy);
% 
% 
%                     for ii_shift=0:no_dec_time_bins_xy-1
%                         shifted_this_training_range=zeros(1,length(this_training_range));
%                         shifted_this_training_range(1+ii_shift:end)=this_training_range(1:end-ii_shift);
% 
%                         shifted_this_test_range=zeros(1,length(this_test_range));
%                         shifted_this_test_range(1+ii_shift:end)=this_test_range(1:end-ii_shift);
% 
%                         XYtrain(:,size(pos_binned,2)*ii_shift+1:size(pos_binned,2)*(ii_shift+1))=pos_binned(logical(shifted_this_training_range),:);
%                         XYtest(:,size(pos_binned,2)*ii_shift+1:size(pos_binned,2)*(ii_shift+1))=pos_binned(logical(shifted_this_test_range),:);
%                     end
% 
                    dFF_pred(trNo).this_XY=pos_binned(logical(this_test_range),:);
% 
%                     no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
%                     odor_plume_train=zeros(sum(this_training_range),no_dec_time_bins_op);
%                     odor_plume_test=zeros(sum(logical(this_test_range)),no_dec_time_bins_op);
% 
%                     for ii_shift=0:no_dec_time_bins_op-1
%                         shifted_this_training_range=zeros(1,length(this_training_range));
%                         shifted_this_training_range(1+ii_shift:end)=this_training_range(1:end-ii_shift);
% 
%                         shifted_this_test_range=zeros(1,length(this_test_range));
%                         shifted_this_test_range(1+ii_shift:end)=this_test_range(1:end-ii_shift);
% 
%                         odor_plume_train(:,ii_shift+1)=odor_plume_template(logical(shifted_this_training_range));
%                         odor_plume_test(:,ii_shift+1)=odor_plume_template(logical(shifted_this_test_range));
%                     end
% 
%                     %Decode using neural network
% 
%                     %Decode using x and y only
%                     % opts = struct('UseParallel', true);
%                     switch handles_choices.algo
%                         case 1
%                             Mdl_dFFxy = fitrnet(XYtrain,this_dFFtrain,'Standardize',true);
%                         case 2
%                             Mdl_dFFxy = fitglm(XYtrain,this_dFFtrain);
%                         case 3
%                             Mdl_dFFxy = fitrtree(XYtrain,this_dFFtrain);
%                         case 4
%                             Mdl_dFFxy = fitrsvm(XYtrain,this_dFFtrain, 'KernelFunction','rbf','Standardize', true);
%                         case 5
%                             Mdl_dFFxy = fitrgp(XYtrain,this_dFFtrain);
%                         case 6
%                             Mdl_dFFxy = fitrgp(XYtrain,this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
%                         case 7
%                             opts = struct('ShowPlots', false, ...
%                                 'Verbose', 0, ...
%                                 'MaxObjectiveEvaluations', 15,...
%                                 'UseParallel',true);
%                             Mdl_dFFxy = fitrnet(XYtrain,this_dFFtrain,'OptimizeHyperparameters', 'auto',...
%                                 'HyperparameterOptimizationOptions', opts);
%                     end
%                     dFF_pred(trNo).dFFpred_xy=predict(Mdl_dFFxy,XYtest);
%                     dFF_pred(trNo).Mdl_dFFxy=Mdl_dFFxy;
% 
%                     %Decode using odor plume only
%                     switch handles_choices.algo
%                         case 1
%                             Mdl_dFFop = fitrnet(odor_plume_train,this_dFFtrain,'Standardize',true);
%                         case 2
%                             Mdl_dFFop = fitglm(odor_plume_train,this_dFFtrain);
%                         case 3
%                             Mdl_dFFop = fitrtree(odor_plume_train,this_dFFtrain);
%                         case 4
%                             Mdl_dFFop = fitrsvm(odor_plume_train,this_dFFtrain, 'KernelFunction','rbf','Standardize', true);
%                         case 5
%                             Mdl_dFFop = fitrgp(odor_plume_train,this_dFFtrain);
%                         case 6
%                             Mdl_dFFop = fitrgp(odor_plume_train,this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
%                         case 7
%                             opts = struct('ShowPlots', false, ...
%                                 'Verbose', 0, ...
%                                 'MaxObjectiveEvaluations', 15,...
%                                 'UseParallel',true);
%                             Mdl_dFFop = fitrnet(odor_plume_train,this_dFFtrain,'OptimizeHyperparameters', 'auto',...
%                                 'HyperparameterOptimizationOptions', opts);
%                     end
%                     dFF_pred(trNo).dFFpred_op=predict(Mdl_dFFop,odor_plume_test);
%                     dFF_pred(trNo).Mdl_dFFop=Mdl_dFFop;
% 
%                     %Decode using both xy and odor plume only
%                     switch handles_choices.algo
%                         case 1
%                             Mdl_dFFxyop = fitrnet([XYtrain odor_plume_train],this_dFFtrain,'Standardize',true);
%                         case 2
%                             Mdl_dFFxyop = fitglm([XYtrain odor_plume_train],this_dFFtrain);
%                         case 3
%                             Mdl_dFFxyop = fitrtree([XYtrain odor_plume_train],this_dFFtrain);
%                         case 4
%                             Mdl_dFFxyop = fitrsvm([XYtrain odor_plume_train],this_dFFtrain, 'KernelFunction','rbf','Standardize', true);
%                         case 5
%                             Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train],this_dFFtrain);
%                         case 6
%                             Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train],this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
%                         case 7
%                             opts = struct('ShowPlots', false, ...
%                                 'Verbose', 0, ...
%                                 'MaxObjectiveEvaluations', 15,...
%                                 'UseParallel',true);
%                             Mdl_dFFxyop = fitrnet([XYtrain odor_plume_train],this_dFFtrain,'OptimizeHyperparameters', 'auto',...
%                                 'HyperparameterOptimizationOptions', opts);
%                     end
%                     dFF_pred(trNo).dFFpred_xyop=predict(Mdl_dFFxyop,[XYtest odor_plume_test]);
%                     dFF_pred(trNo).Mdl_dFFxyop=Mdl_dFFxyop;
% 
%                     pffft=1;
%                     % y_pred(trNo).data=predict(MdlY2,XdFFtest);
                end
%                 % delete(gcp('nocreate'));
% 
% 
%                 %Parse out the parfor loop output
%                 all_dFFpred_xy=[];
%                 all_dFFpred_xyop=[];
%                 all_dFFpred_op=[];
                all_dFF=[];
                all_lanes=[];
                all_XY=[];
% 
%                 all_dFFpredl1_xy=[];
%                 all_dFFpredl1_xyop=[];
%                 all_dFFpredl1_op=[];
                all_dFFl1=[];
% 
%                 all_dFFpredl4_xy=[];
%                 all_dFFpredl4_xyop=[];
%                 all_dFFpredl4_op=[];
                all_dFFl4=[];
% 
                all_XYl1=[];
                all_XYl4=[];
% 
                all_tr=[];
                all_trl1=[];
                all_trl4=[];
% 
                ii_trl1=0;
                ii_trl4=0;
% 
                for trNo=1:trials.odor_trNo
%                     per_ROI(this_ROI).results.trial(trNo).dFFpred_xy=dFF_pred(trNo).dFFpred_xy;
%                     per_ROI(this_ROI).results.trial(trNo).dFFpred_xyop=dFF_pred(trNo).dFFpred_xyop;
%                     per_ROI(this_ROI).results.trial(trNo).dFFpred_op=dFF_pred(trNo).dFFpred_op;
%                     per_ROI(this_ROI).results.trial(trNo).dFF=dFF_pred(trNo).dFF;
% 
%                     % per_ROI(this_ROI).results.trial(trNo).Mdl_dFFxy=dFF_pred(trNo).Mdl_dFFxy;
%                     % per_ROI(this_ROI).results.trial(trNo).Mdl_dFFxyop=dFF_pred(trNo).Mdl_dFFxyop;
%                     % per_ROI(this_ROI).results.trial(trNo).Mdl_dFFop=dFF_pred(trNo).Mdl_dFFop;
% 
%                     all_dFFpred_xy=[all_dFFpred_xy dFF_pred(trNo).dFFpred_xy'];
%                     all_dFFpred_xyop=[all_dFFpred_xyop dFF_pred(trNo).dFFpred_xyop'];
%                     all_dFFpred_op=[all_dFFpred_op dFF_pred(trNo).dFFpred_op'];
                    all_dFF=[all_dFF dFF_pred(trNo).dFF'];
                    all_tr(trNo)=size(all_XY,1)+1;
                    all_XY=[all_XY; dFF_pred(trNo).this_XY];
% 
% 
                    if trials.odor_lane(trNo)==1
                        ii_trl1=ii_trl1+1;
                        all_trl1(ii_trl1)=size(all_XYl1,1)+1;
                        all_XYl1=[all_XYl1; dFF_pred(trNo).this_XY];
                        % all_dFFpredl1_xy=[all_dFFpredl1_xy dFF_pred(trNo).dFFpred_xy'];
                        % all_dFFpredl1_xyop=[all_dFFpredl1_xyop dFF_pred(trNo).dFFpred_xyop'];
                        % all_dFFpredl1_op=[all_dFFpredl1_op dFF_pred(trNo).dFFpred_op'];
                        all_dFFl1=[all_dFFl1 dFF_pred(trNo).dFF'];
                        all_lanes=[all_lanes ones(1,length(dFF_pred(trNo).dFF))];
                    else
                        ii_trl4=ii_trl4+1;
                        all_trl4(ii_trl4)=size(all_XYl4,1)+1;
                        all_XYl4=[all_XYl4; dFF_pred(trNo).this_XY];
                        % all_dFFpredl4_xy=[all_dFFpredl4_xy dFF_pred(trNo).dFFpred_xy'];
                        % all_dFFpredl4_xyop=[all_dFFpredl4_xyop dFF_pred(trNo).dFFpred_xyop'];
                        % all_dFFpredl4_op=[all_dFFpredl4_op dFF_pred(trNo).dFFpred_op'];
                        all_dFFl4=[all_dFFl4 dFF_pred(trNo).dFF'];
                        all_lanes=[all_lanes 4*ones(1,length(dFF_pred(trNo).dFF))];
                    end
                end
% 
                per_ROI(this_ROI).results.all_lanes=all_lanes;
% 
                per_ROI(this_ROI).results.all_XY=all_XY;
                per_ROI(this_ROI).results.all_XYl1=all_XYl1;
                per_ROI(this_ROI).results.all_XYl4=all_XYl4;

                per_ROI(this_ROI).results.all_tr=all_tr;
                per_ROI(this_ROI).results.all_trl1=all_trl1;
                per_ROI(this_ROI).results.all_trl4=all_trl4;
% 
                per_ROI(this_ROI).results.all_dFF=all_dFF;
%                 per_ROI(this_ROI).results.all_dFFpred_xy=all_dFFpred_xy;
%                 per_ROI(this_ROI).results.all_dFFpred_xyop=all_dFFpred_xyop;
%                 per_ROI(this_ROI).results.all_dFFpred_op=all_dFFpred_op;
% 
                per_ROI(this_ROI).results.all_dFFl1=all_dFFl1;
%                 per_ROI(this_ROI).results.all_dFFpredl1_xy=all_dFFpredl1_xy;
%                 per_ROI(this_ROI).results.all_dFFpredl1_xyop=all_dFFpredl1_xyop;
%                 per_ROI(this_ROI).results.all_dFFpredl1_op=all_dFFpredl1_op;
% 
                per_ROI(this_ROI).results.all_dFFl4=all_dFFl4;
%                 per_ROI(this_ROI).results.all_dFFpredl4_xy=all_dFFpredl4_xy;
%                 per_ROI(this_ROI).results.all_dFFpredl4_xyop=all_dFFpredl4_xyop;
%                 per_ROI(this_ROI).results.all_dFFpredl4_op=all_dFFpredl4_op;
% 
%                 thisRxy=corrcoef([all_dFF' all_dFFpred_xy']);
%                 per_ROI(this_ROI).results.Rxy=thisRxy(1,2);
% 
%                 if (~isempty(all_dFFl1))&(~isempty(all_dFFpredl1_xy))
%                     thisRxyl1=corrcoef([all_dFFl1' all_dFFpredl1_xy']);
%                     per_ROI(this_ROI).results.Rxyl1=thisRxyl1(1,2);
%                 else
%                     per_ROI(this_ROI).results.Rxyl1=NaN;
%                 end
% 
%                 if (~isempty(all_dFFl4))&(~isempty(all_dFFpredl4_xy))
%                     thisRxyl4=corrcoef([all_dFFl4' all_dFFpredl4_xy']);
%                     per_ROI(this_ROI).results.Rxyl4=thisRxyl4(1,2);
%                 else
%                     per_ROI(this_ROI).results.Rxyl4=NaN;
%                 end
% 
%                 thisRxyop=corrcoef([all_dFF' all_dFFpred_xyop']);
%                 per_ROI(this_ROI).results.Rxyop=thisRxyop(1,2);
% 
%                 if (~isempty(all_dFFl1))&(~isempty(all_dFFpredl1_xyop))
%                     thisRxyopl1=corrcoef([all_dFFl1' all_dFFpredl1_xyop']);
%                     per_ROI(this_ROI).results.Rxyopl1=thisRxyopl1(1,2);
%                 else
%                     per_ROI(this_ROI).results.Rxyopl1=NaN;
%                 end
% 
%                 if (~isempty(all_dFFl4))&(~isempty(all_dFFpredl4_xyop))
%                     thisRxyopl4=corrcoef([all_dFFl4' all_dFFpredl4_xyop']);
%                     per_ROI(this_ROI).results.Rxyopl4=thisRxyopl4(1,2);
%                 else
%                     per_ROI(this_ROI).results.Rxyopl4=NaN;
%                 end
% 
% 
%                 thisRop=corrcoef([all_dFF' all_dFFpred_op']);
%                 per_ROI(this_ROI).results.Rop=thisRop(1,2);
% 
%                 if (~isempty(all_dFFl1))&(~isempty(all_dFFpredl1_op))
%                     thisRopl1=corrcoef([all_dFFl1' all_dFFpredl1_op']);
%                     per_ROI(this_ROI).results.Ropl1=thisRopl1(1,2);
%                 else
%                     per_ROI(this_ROI).results.Ropl1=NaN;
%                 end
% 
%                 if (~isempty(all_dFFl4))&(~isempty(all_dFFpredl4_op))
%                     thisRopl4=corrcoef([all_dFFl4' all_dFFpredl4_op']);
%                     per_ROI(this_ROI).results.Ropl4=thisRop(1,2);
%                 else
%                     per_ROI(this_ROI).results.Ropl4=NaN;
%                 end
% 
%                 fprintf(1, 'ROI No %d, rho for prediction dFFxy, dFFop dFFxyop %d %d %d\n'...
%                     ,this_ROI, per_ROI(this_ROI).results.Rxy,per_ROI(this_ROI).results.Rop...
%                     ,per_ROI(this_ROI).results.Rxyop);
% 
%                 thisROI_done=1;
%                 fprintf(1,['Elapsed time ' num2str((toc-start_toc)/(60)) ' mins \n\n'])
%             catch
%                 delete(gcp('nocreate'));
%             end
%         end
%         handles_out.per_ROI=per_ROI;
%         % save([this_path arena_file(1:end-4) '_' handles_choices.save_tag '_' num2str(handles_choices.algo) ...
%         %     num2str(handles_choices.weber_fechner) num2str(handles_choices.alpha) num2str(handles_choices.which_odor_plume)...
%         %     num2str(no_dec_time_bins_op) num2str(no_dec_time_bins_xy) '.mat'],'handles_out','handles_choices','trials','-v7.3')
% 
% 
% 
% 
%         %Train with per trial with shuffled data using a leave one out approach
% 
%         % per_ROI_sh=[];
%         for ii_repeat=1:handles_choices.sh_repeats
% 
% 
%             % if skip_this_ii==0
%             %training_range_template has all the trials
%             training_range_template=zeros(1,no_time_bins);
%             lane_template=zeros(1,no_time_bins);
%             odor_plume_template=zeros(1,no_time_bins);
% 
% 
%             for trNo=1:trials.odor_trNo
%                 dFF_pred(trNo).dFFpred_xy=[];
%                 % dFF_pred(trNo).dFFpred_xyl=[];
%                 dFF_pred(trNo).dFF=[];
%                 % x_pred(trNo).data=[];
% 
%                 Mdl_dFFxy_pars(trNo).pars=[];
%                 Mdl_dFFop_pars(trNo).pars=[];
%                 Mdl_dFFxyop_pars(trNo).pars=[];
% 
%                 x_predictedstart=trials.odor_ii_bin_start(trNo);
%                 x_predictedend=trials.odor_ii_bin_end(trNo);
%                 training_range_template(x_predictedstart:x_predictedend)=1;
%                 if trials.odor_lane(trNo)==1
%                     lane_template(x_predictedstart:x_predictedend)=1;
%                     for x_ii=x_predictedstart:x_predictedend
%                         this_x=pos_binned(x_ii,1);
%                         this_y=pos_binned(x_ii,2);
%                         [minabx,ii_minax]=min(abs(x_for_plume-this_x));
%                         [minaby,ii_minay]=min(abs(y_for_plume-this_y));
%                         this_ca=mean_plume_l1(ii_minay,ii_minax);
%                         odor_plume_template(x_ii)=this_ca;
%                     end
%                 else
%                     lane_template(x_predictedstart:x_predictedend)=4;
%                     for x_ii=x_predictedstart:x_predictedend
%                         this_x=pos_binned(x_ii,1);
%                         this_y=pos_binned(x_ii,2);
%                         [minabx,ii_minax]=min(abs(x_for_plume-this_x));
%                         [minaby,ii_minay]=min(abs(y_for_plume-this_y));
%                         this_ca=mean_plume_l4(ii_minay,ii_minax);
%                         odor_plume_template(x_ii)=this_ca;
%                     end
%                 end
%             end
% 
% 
%             start_toc=toc;
% 
%             % Create the parallel pool once before the main loop
%             parfor_is_done=0;
%             % if isempty(gcp('nocreate'))
%             %     parpool;  % This will use the default number of workers
%             % end
%             while parfor_is_done==0
%                 try
%                     if isempty(gcp('nocreate'))
%                         parpool;  % This will use the default number of workers
%                     end
%                     parfor trNo=1:trials.odor_trNo
%                         % for trNo=1:trials.odor_trNo
% 
%                         this_test_range=zeros(1,no_time_bins);
%                         % if trNo==1
%                         %     ii_test_range_start=1;
%                         % else
%                         %     ii_test_range_start=trials.odor_ii_bin_end(trNo-1)+15;
%                         % end
%                         %
%                         % if trNo==trials.odor_trNo
%                         %     ii_test_range_end=no_time_bins;
%                         % else
%                         %     ii_test_range_end=trials.odor_ii_bin_end(trNo)+15;
%                         % end
% 
%                         ii_test_range_start=trials.odor_ii_bin_start(trNo);
%                         ii_test_range_end=trials.odor_ii_bin_end(trNo);
% 
%                         this_test_range(ii_test_range_start:ii_test_range_end)=1;
%                         this_training_range=logical(training_range_template)&(~logical(this_test_range));
% 
%                         % ii_valid_range=ceil(valid_range*no_time_bins);
%                         %     ii_test_range=ceil(test_range*no_time_bins);
%                         this_dFFtrain=zeros(sum(this_training_range),1);
%                         this_dFFtrain=neural_data(this_training_range,this_ROI);
%                         % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
%                         this_dFFtest=zeros(sum(this_test_range),1);
%                         this_dFFtest(:,1)=neural_data(logical(this_test_range),this_ROI);
%                         dFF_pred(trNo).dFF=this_dFFtest;
% 
%                         %Setup the pre variables
%                         no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
%                         XYtrain_pre=zeros(sum(this_training_range),size(pos_binned,2)*no_dec_time_bins_xy);
%                         XYtest=zeros(sum(logical(this_test_range)),size(pos_binned,2)*no_dec_time_bins_xy);
% 
%                         for ii_shift=0:no_dec_time_bins_xy-1
%                             shifted_this_training_range=zeros(1,length(this_training_range));
%                             shifted_this_training_range(1+ii_shift:end)=this_training_range(1:end-ii_shift);
% 
%                             shifted_this_test_range=zeros(1,length(this_test_range));
%                             shifted_this_test_range(1+ii_shift:end)=this_test_range(1:end-ii_shift);
% 
%                             XYtrain_pre(:,size(pos_binned,2)*ii_shift+1:size(pos_binned,2)*(ii_shift+1))=pos_binned(logical(shifted_this_training_range),:);
%                             XYtest(:,size(pos_binned,2)*ii_shift+1:size(pos_binned,2)*(ii_shift+1))=pos_binned(logical(shifted_this_test_range),:);
% 
%                         end
% 
%                         no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
%                         odor_plume_train_pre=zeros(sum(this_training_range),no_dec_time_bins_op);
%                         odor_plume_test=zeros(sum(logical(this_test_range)),no_dec_time_bins_op);
% 
%                         for ii_shift=0:no_dec_time_bins_op-1
%                             shifted_this_training_range=zeros(1,length(this_training_range));
%                             shifted_this_training_range(1+ii_shift:end)=this_training_range(1:end-ii_shift);
% 
%                             shifted_this_test_range=zeros(1,length(this_test_range));
%                             shifted_this_test_range(1+ii_shift:end)=this_test_range(1:end-ii_shift);
% 
%                             odor_plume_train_pre(:,ii_shift+1)=odor_plume_template(logical(shifted_this_training_range));
%                             odor_plume_test(:,ii_shift+1)=odor_plume_template(logical(shifted_this_test_range));
%                         end
% 
%                         %Now shuffle
%                         % XYtrain_pre=pos_binned(this_training_range,:);
%                         shuffle_length=floor(size(XYtrain_pre,1)/(trials.odor_trNo*2));
% 
%                         this_overlap=10;
%                         while this_overlap>max_overlap
%                             perm_indices=randperm(trials.odor_trNo*2);
%                             this_overlap=0;
%                             for ii_p=1:trials.odor_trNo*2
%                                 if perm_indices(ii_p)==ii_p
%                                     this_overlap=this_overlap+1;
%                                 end
%                             end
%                         end
% 
%                         % % Yvalid=pos_binned(ii_valid_range(1):ii_valid_range(2),:);
%                         % XYtest=pos_binned(logical(this_test_range),:);
%                         %
%                         % lane_train_pre=lane_template(this_training_range);
%                         % lane_train_pre=lane_train_pre';
%                         % lane_test=lane_template(logical(this_test_range));
%                         %
%                         % odor_plume_train_pre=odor_plume_template(this_training_range);
%                         % odor_plume_train_pre=odor_plume_train_pre';
%                         % odor_plume_test=odor_plume_template(logical(this_test_range));
% 
% 
% 
% 
%                         %Now shuffle
%                         shuffle_length=floor(size(XYtrain_pre,1)/(trials.odor_trNo*2));
% 
%                         this_overlap=10;
%                         while this_overlap>max_overlap
%                             perm_indices=randperm(trials.odor_trNo*2);
%                             this_overlap=0;
%                             for ii_p=1:trials.odor_trNo*2
%                                 if perm_indices(ii_p)==ii_p
%                                     this_overlap=this_overlap+1;
%                                 end
%                             end
%                         end
% 
%                         XYtrain=zeros(size(XYtrain_pre,1),size(XYtrain_pre,2));
%                         % lane_train=zeros(size(XYtrain_pre,1),1);
%                         odor_plume_train=zeros(size(odor_plume_train_pre,1),size(odor_plume_train_pre,2));
%                         for sh_ii=1:trials.odor_trNo*2
%                             XYtrain((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
%                                 XYtrain_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
%                             % lane_train((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
%                             %     lane_train_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
%                             odor_plume_train((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
%                                 odor_plume_train_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
%                         end
% 
%                         %Make the last few equal to the first few
%                         if trials.odor_trNo*2*shuffle_length<size(XYtrain,1)
%                             delta_ii=size(XYtrain,1)-trials.odor_trNo*2*shuffle_length;
%                             XYtrain(trials.odor_trNo*2*shuffle_length+1:end,:)=...
%                                 XYtrain_pre(1:delta_ii,:);
%                             % lane_train(trials.odor_trNo*2*shuffle_length+1:end,:)=...
%                             %     lane_train_pre(1:delta_ii,:);
%                             odor_plume_train(trials.odor_trNo*2*shuffle_length+1:end,:)=...
%                                 odor_plume_train_pre(1:delta_ii,:);
%                         end
% 
% 
%                         %Decode using neural network
% 
%                         % opts = struct('UseParallel', true);
% 
%                         %Decode using x and y only
%                         switch handles_choices.algo
%                             case 1
%                                 Mdl_dFFxy = fitrnet(XYtrain,this_dFFtrain','Standardize',true);
%                             case 2
%                                 Mdl_dFFxy = fitglm(XYtrain,this_dFFtrain');
%                             case 3
%                                 Mdl_dFFxy = fitrtree(XYtrain,this_dFFtrain');
%                             case 4
%                                 Mdl_dFFxy = fitrsvm(XYtrain,this_dFFtrain', 'KernelFunction','rbf','Standardize', true);
%                             case 5
%                                 Mdl_dFFxy = fitrgp(XYtrain,this_dFFtrain');
%                             case 6
%                                 Mdl_dFFxy = fitrgp(XYtrain,this_dFFtrain','KernelFunction', 'squaredexponential', 'Standardize', true);
%                             case 7
%                                 % opts = struct('ShowPlots', false, ...
%                                 %     'Verbose', 0, ...
%                                 %     'MaxObjectiveEvaluations', 15);
%                                 % Mdl_dFFxy = fitrnet(XYtrain,this_dFFtrain','OptimizeHyperparameters', 'auto',...
%                                 %     'HyperparameterOptimizationOptions', opts);
% 
%                                 bestHyperparameters = dFF_pred(trNo).Mdl_dFFxy.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
%                                 activationsCell = cellstr(bestHyperparameters.Activations);
%                                 standardizeCell = cellstr(bestHyperparameters.Standardize);
%                                 layer_sizes=bestHyperparameters.Layer_1_Size;
%                                 if ~isnan(bestHyperparameters.Layer_2_Size)
%                                     layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
%                                 end
%                                 if ~isnan(bestHyperparameters.Layer_3_Size)
%                                     layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
%                                 end
%                                 Mdl_dFFxy = fitrnet(XYtrain,this_dFFtrain,'LayerSizes', layer_sizes, ...
%                                     'Activations', activationsCell{1}, ...
%                                     'Lambda', bestHyperparameters.Lambda, ...
%                                     'Standardize', strcmpi(standardizeCell{1},'true'));
% 
%                                 Mdl_dFFxy_pars(trNo).pars.activations=activationsCell{1};
%                                 Mdl_dFFxy_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
%                                 Mdl_dFFxy_pars(trNo).pars.Standardize=standardizeCell{1};
% 
%                         end
%                         dFF_pred(trNo).dFFpred_xy=predict(Mdl_dFFxy,XYtest);
% 
% 
% 
%                         %Decode using odor plume
%                         switch handles_choices.algo
%                             case 1
%                                 Mdl_dFFop = fitrnet([odor_plume_train],this_dFFtrain,'Standardize',true);
%                             case 2
%                                 Mdl_dFFop = fitglm([odor_plume_train],this_dFFtrain);
%                             case 3
%                                 Mdl_dFFop = fitrtree([odor_plume_train],this_dFFtrain);
%                             case 4
%                                 Mdl_dFFop = fitrsvm([odor_plume_train],this_dFFtrain, 'KernelFunction','rbf','Standardize', true);
%                             case 5
%                                 Mdl_dFFop = fitrgp([odor_plume_train],this_dFFtrain);
%                             case 6
%                                 Mdl_dFFop = fitrgp([odor_plume_train],this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
%                             case 7
%                                 bestHyperparameters = dFF_pred(trNo).Mdl_dFFop.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
%                                 activationsCell = cellstr(bestHyperparameters.Activations);
%                                 standardizeCell = cellstr(bestHyperparameters.Standardize);
%                                 layer_sizes=bestHyperparameters.Layer_1_Size;
%                                 if ~isnan(bestHyperparameters.Layer_2_Size)
%                                     layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
%                                 end
%                                 if ~isnan(bestHyperparameters.Layer_3_Size)
%                                     layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
%                                 end
%                                 Mdl_dFFop = fitrnet([odor_plume_train],this_dFFtrain,'LayerSizes', layer_sizes, ...
%                                     'Activations', activationsCell{1}, ...
%                                     'Lambda', bestHyperparameters.Lambda, ...
%                                     'Standardize', strcmpi(standardizeCell{1},'true'));
% 
%                                 Mdl_dFFop_pars(trNo).pars.activations=activationsCell{1};
%                                 Mdl_dFFop_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
%                                 Mdl_dFFop_pars(trNo).pars.Standardize=standardizeCell{1};
%                         end
%                         dFF_pred(trNo).dFFpred_op=predict(Mdl_dFFop,odor_plume_test);
% 
%                         %Decode using xy and odor plume
%                         switch handles_choices.algo
%                             case 1
%                                 Mdl_dFFxyop = fitrnet([XYtrain odor_plume_train],this_dFFtrain,'Standardize',true);
%                             case 2
%                                 Mdl_dFFxyop = fitglm([XYtrain odor_plume_train],this_dFFtrain);
%                             case 3
%                                 Mdl_dFFxyop = fitrtree([XYtrain odor_plume_train],this_dFFtrain);
%                             case 4
%                                 Mdl_dFFxyop = fitrsvm([XYtrain odor_plume_train],this_dFFtrain, 'KernelFunction','rbf','Standardize', true);
%                             case 5
%                                 Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train],this_dFFtrain);
%                             case 6
%                                 Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train],this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
%                             case 7
%                                 % opts = struct('ShowPlots', false, ...
%                                 %     'Verbose', 0, ...
%                                 %     'MaxObjectiveEvaluations', 15);
%                                 % Mdl_dFFxyop = fitrnet([XYtrain odor_plume_train],this_dFFtrain,'OptimizeHyperparameters', 'auto',...
%                                 %     'HyperparameterOptimizationOptions', opts);
% 
%                                 bestHyperparameters = dFF_pred(trNo).Mdl_dFFxy.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
%                                 activationsCell = cellstr(bestHyperparameters.Activations);
%                                 standardizeCell = cellstr(bestHyperparameters.Standardize);
%                                 layer_sizes=bestHyperparameters.Layer_1_Size;
%                                 if ~isnan(bestHyperparameters.Layer_2_Size)
%                                     layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
%                                 end
%                                 if ~isnan(bestHyperparameters.Layer_3_Size)
%                                     layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
%                                 end
%                                 Mdl_dFFxyop = fitrnet([XYtrain odor_plume_train],this_dFFtrain,'LayerSizes', layer_sizes, ...
%                                     'Activations', activationsCell{1}, ...
%                                     'Lambda', bestHyperparameters.Lambda, ...
%                                     'Standardize', strcmpi(standardizeCell{1},'true'));
% 
%                                 Mdl_dFFxyop_pars(trNo).pars.activations=activationsCell{1};
%                                 Mdl_dFFxyop_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
%                                 Mdl_dFFxyop_pars(trNo).pars.Standardize=standardizeCell{1};
%                         end
%                         dFF_pred(trNo).dFFpred_xyop=predict(Mdl_dFFxyop,[XYtest odor_plume_test]);
% 
% 
% 
%                         pffft=1;
%                         % y_pred(trNo).data=predict(MdlY2,XdFFtest);
%                     end
%                     parfor_is_done=1;
%                 catch
%                 end
%             end
%             delete(gcp('nocreate'));
% 
% 
% 
%             %Parse out the parfor loop output
%             all_dFFpred_xy=[];
%             all_dFFpred_xyl=[];
%             all_dFFpred_xyop=[];
%             all_dFFpred_op=[];
%             all_dFF=[];
%             for trNo=1:trials.odor_trNo
%                 per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).dFFpred_xy=dFF_pred(trNo).dFFpred_xy;
%                 all_dFFpred_xy=[all_dFFpred_xy dFF_pred(trNo).dFFpred_xy'];
%                 % per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).dFFpred_xyl=dFF_pred(trNo).dFFpred_xyl;
%                 % all_dFFpred_xyl=[all_dFFpred_xyl dFF_pred(trNo).dFFpred_xyl'];
%                 all_dFFpred_xyop=[all_dFFpred_xyop dFF_pred(trNo).dFFpred_xyop'];
%                 all_dFFpred_op=[all_dFFpred_op dFF_pred(trNo).dFFpred_op'];
%                 per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).dFF=dFF_pred(trNo).dFF;
%                 all_dFF=[all_dFF dFF_pred(trNo).dFF'];
%             end
% 
%             per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFF=all_dFF;
%             per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_xy=all_dFFpred_xy;
%             % per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_xyl=all_dFFpred_xyl;
%             per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_xyop=all_dFFpred_xyop;
%             per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_op=all_dFFpred_op;
% 
%             thisRxy=corrcoef([all_dFF' all_dFFpred_xy']);
%             per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy=thisRxy(1,2);
%             % thisRxyl=corrcoef([all_dFF' all_dFFpred_xyl']);
%             % per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyl=thisRxyl(1,2);
%             thisRxyop=corrcoef([all_dFF' all_dFFpred_xyop']);
%             per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop=thisRxyop(1,2);
%             thisRop=corrcoef([all_dFF' all_dFFpred_op']);
%             per_ROI_sh(this_ROI).repeats(ii_repeat).Rop=thisRop(1,2);
% 
%             per_ROI(this_ROI).results.Mdl_dFFxy_pars=Mdl_dFFxy_pars;
%             per_ROI(this_ROI).results.Mdl_dFFop_pars=Mdl_dFFop_pars;
%             per_ROI(this_ROI).results.Mdl_dFFxyop_pars=Mdl_dFFxyop_pars;
% 
% 
%             % fprintf(1, 'Repeat %d ROI No %d, rho for shuffled prediction dFFxy, dFFxyl dFFxyop, dFFop %d %d %d %d\n\n',ii_repeat,...
%             %     this_ROI,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyl,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop,per_ROI_sh(this_ROI).repeats(ii_repeat).Rop);
%             fprintf(1, 'Repeat %d ROI No %d, rho for shuffled prediction dFFxy, dFFop, dFFxyop %d %d %d\n',ii_repeat,...
%                 this_ROI,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rop...
%                 ,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop);
%             fprintf(1,['Elapsed time ' num2str((toc-start_toc)/(60)) ' mins shuffled\n\n'])
%             %         thisROI_done=1;
%             %     catch
%             %         delete(gcp('nocreate'));
%             %     end
%             % end
% 
%         end
%         handles_out.per_ROI_sh=per_ROI_sh;
%         save([this_path arena_file(1:end-4) '_' handles_choices.save_tag '_' num2str(handles_choices.algo) ...
%             num2str(handles_choices.weber_fechner) num2str(handles_choices.alpha)...
%             num2str(no_dec_time_bins_op) num2str(no_dec_time_bins_xy) '.mat'],'handles_out','handles_choices','trials','-v7.3')
%     end
end
% 
% fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])
% 
% handles_out.per_ROI=per_ROI;
% handles_out.per_ROI_sh=per_ROI_sh;
handles_out.trials=trials;


% figNo=figNo+1;
% try
%     close(figNo)
% catch
% en
%Calculate 95th percentile
% figNo=0;
% noROIs=length(handles_out.per_ROI);
% if isfield(handles_out.per_ROI(1),'sh_repeats')
%     noshRepeats=length(handles_out.per_ROI(1).sh_repeats);
% else
%     noshRepeats=length(handles_out.per_ROI(1).repeats);
% end
% 
% noshRepeats=handles_choices.sh_repeats;
% 
% shRops=[];
% shRxys=[];
% shRxyops=[];
% for ii_ROI=handles_choices.process_these_ROIs
%     for ii_sh_repeats=1:noshRepeats
%         shRops=[shRops handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rop];
%         shRxys=[shRxys handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rxy];
%         shRxyops=[shRxyops handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rxyop];
%     end
% end
% 
% % Bootstrap resampling op
% try
%     bootstrap_samples_op = bootstrp(1000, @prctile, shRops, 95);
% catch
%     pffft=1;
% end
% 
% % Compute the 95th percentile of the bootstrap samples
% Ropsh_bootstrap_95th_percentile = prctile(bootstrap_samples_op, 95);
% 
% % Bootstrap resampling xy
% bootstrap_samples_xy = bootstrp(1000, @prctile, shRxys, 95);
% 
% % Compute the 95th percentile of the bootstrap samples
% Rxysh_bootstrap_95th_percentile = prctile(bootstrap_samples_xy, 95);
% 
% % Bootstrap resampling xyop
% bootstrap_samples_xyop = bootstrp(1000, @prctile, shRxyops, 95);
% 
% % Compute the 95th percentile of the bootstrap samples
% Rxyopsh_bootstrap_95th_percentile = prctile(bootstrap_samples_xyop, 95);
% 
% 
% %Find the Rops that are above 95th percentile
% 
% 
% Rops=[];
% Rxys=[];
% Rxyops=[];
% iiROI_Rops=[];
% iiROI_Rxys=[];
% iiROI_Rxyops=[];
% Rops_above_95=[];
% Rxys_above_95=[];
% Rxyops_above_95=[];
% Rops_for_Rxys_above_95=[];
% Roxyl1s_for_Rxys_above_95=[];
% Roxyl4s_for_Rxys_above_95=[];
% Rxys_for_Rops_above_95=[];
% Ropl1s_for_Rops_above_95=[];
% Ropl4s_for_Rops_above_95=[];
% Rxys_for_Rxyops_above_95=[];
% Rxyopl1s_for_Rxyops_above_95=[];
% Rxyopl4s_for_Rxyops_above_95=[];
% Rxyops_for_Rxys_above_95=[];
% 
% for ii_ROI=handles_choices.process_these_ROIs
% 
%     Rops=[Rops handles_out.per_ROI(ii_ROI).results.Rop];
%     Rxys=[Rxys handles_out.per_ROI(ii_ROI).results.Rxy];
%     Rxyops=[Rxyops handles_out.per_ROI(ii_ROI).results.Rxyop];
% 
%     if handles_out.per_ROI(ii_ROI).results.Rop>Ropsh_bootstrap_95th_percentile
%         iiROI_Rops=[iiROI_Rops ii_ROI];
%         % fileNo_Rops=[fileNo_Rops fileNo];
%         Rops_above_95=[Rops_above_95 handles_out.per_ROI(ii_ROI).results.Rop];
%         Rxys_for_Rops_above_95=[Rxys_for_Rops_above_95 handles_out.per_ROI(ii_ROI).results.Rxy];
%         Ropl1s_for_Rops_above_95=[Ropl1s_for_Rops_above_95 handles_out.per_ROI(ii_ROI).results.Ropl1];
%         Ropl4s_for_Rops_above_95=[Ropl4s_for_Rops_above_95 handles_out.per_ROI(ii_ROI).results.Ropl4];
%     end
%     if handles_out.per_ROI(ii_ROI).results.Rxy>Rxysh_bootstrap_95th_percentile
%         iiROI_Rxys=[iiROI_Rxys ii_ROI];
%         % fileNo_Rxys=[fileNo_Rxys fileNo];
%         Rxys_above_95=[Rxys_above_95 handles_out.per_ROI(ii_ROI).results.Rxy];
%         Rops_for_Rxys_above_95=[Rops_for_Rxys_above_95 handles_out.per_ROI(ii_ROI).results.Rop];
%         Rxyops_for_Rxys_above_95=[Rxyops_for_Rxys_above_95 handles_out.per_ROI(ii_ROI).results.Rxyop];
%         Roxyl1s_for_Rxys_above_95=[Roxyl1s_for_Rxys_above_95 handles_out.per_ROI(ii_ROI).results.Rxyl1];
%         Roxyl4s_for_Rxys_above_95=[Roxyl4s_for_Rxys_above_95 handles_out.per_ROI(ii_ROI).results.Rxyl4];
%     end
%     if handles_out.per_ROI(ii_ROI).results.Rxyop>Rxyopsh_bootstrap_95th_percentile
%         iiROI_Rxyops=[iiROI_Rxyops ii_ROI];
%         % fileNo_Rxys=[fileNo_Rxys fileNo];
%         Rxyops_above_95=[Rxyops_above_95 handles_out.per_ROI(ii_ROI).results.Rxyop];
%         Rxys_for_Rxyops_above_95=[Rxys_for_Rxyops_above_95 handles_out.per_ROI(ii_ROI).results.Rxy];
%         Rxyopl1s_for_Rxyops_above_95=[Rxyopl1s_for_Rxyops_above_95 handles_out.per_ROI(ii_ROI).results.Rxyopl1];
%         Rxyopl4s_for_Rxyops_above_95=[Rxyopl4s_for_Rxyops_above_95 handles_out.per_ROI(ii_ROI).results.Rxyopl4];
%     end
% 
% end
% if handles_choices.displayFigures==1
%     %Display Rho op vs Rho xy figure for this ii_wf and fileNo
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
% 
%     hFig = figure(figNo);
% 
%     set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% 
% 
%     hold on
% 
%     %Plot Rops above 95
%     for ii_rho1=1:length(Rops_above_95)
%         match_found=0;
%         for ii_rho2=1:length(Rxys_above_95)
%             if (iiROI_Rxys(ii_rho2)==iiROI_Rops(ii_rho1))
%                 match_found=1;
%             end
%         end
%         if match_found==1
%             plot(Rxys_for_Rops_above_95(ii_rho1),Rops_above_95(ii_rho1),'ok')
%         else
%             plot(Rxys_for_Rops_above_95(ii_rho1),Rops_above_95(ii_rho1),'ob')
%         end
%     end
% 
%     %Plot Rxys above 95
%     for ii_rho1=1:length(Rxys_above_95)
%         match_found=0;
%         for ii_rho2=1:length(Rops_above_95)
%             if (iiROI_Rops(ii_rho2)==iiROI_Rxys(ii_rho1))
%                 match_found=1;
%             end
%         end
%         if match_found==1
%             plot(Rxys_above_95(ii_rho1),Rops_for_Rxys_above_95(ii_rho1),'ok')
%         else
%             plot(Rxys_above_95(ii_rho1),Rops_for_Rxys_above_95(ii_rho1),'or')
%         end
%     end
% 
%     this_xlim=xlim;
%     this_ylim=ylim;
%     text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.15*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rop>95%','Color',[0 0 1])
%     text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.10*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy>95%','Color',[1 0 0])
%     text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.05*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy and Rop>95%','Color',[0 0 0])
% 
% 
%     plot([this_xlim(1) this_xlim(2)],[this_xlim(1) this_xlim(2)],'-k')
% 
%     xlabel('Rho xy')
%     ylabel('Rho op')
%     title(['Rho op vs. Rho xy log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha)])
% 
%     %Display Rho xyop vs Rho xy figure for this ii_wf and fileNo
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
% 
%     hFig = figure(figNo);
% 
%     set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% 
% 
%     hold on
% 
%     %Plot Rxyops above 95
%     for ii_rho1=1:length(Rxyops_above_95)
%         match_found=0;
%         for ii_rho2=1:length(Rxys_above_95)
%             if (iiROI_Rxys(ii_rho2)==iiROI_Rxyops(ii_rho1))
%                 match_found=1;
%             end
%         end
%         if match_found==1
%             plot(Rxys_for_Rxyops_above_95(ii_rho1),Rxyops_above_95(ii_rho1),'ok')
%         else
%             plot(Rxys_for_Rxyops_above_95(ii_rho1),Rxyops_above_95(ii_rho1),'ob')
%         end
%     end
% 
%     %Plot Rxys above 95
%     for ii_rho1=1:length(Rxys_above_95)
%         match_found=0;
%         for ii_rho2=1:length(Rxyops_above_95)
%             if (iiROI_Rxyops(ii_rho2)==iiROI_Rxys(ii_rho1))
%                 match_found=1;
%             end
%         end
%         if match_found==1
%             plot(Rxys_above_95(ii_rho1),Rxyops_for_Rxys_above_95(ii_rho1),'ok')
%         else
%             plot(Rxys_above_95(ii_rho1),Rxyops_for_Rxys_above_95(ii_rho1),'or')
%         end
%     end
% 
%     this_xlim=xlim;
%     this_ylim=ylim;
%     text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.15*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxyop>95%','Color',[0 0 1])
%     text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.10*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy>95%','Color',[1 0 0])
%     text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.05*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy and Rxyop>95%','Color',[0 0 0])
% 
%     plot([this_xlim(1) this_xlim(2)],[this_xlim(1) this_xlim(2)],'-k')
% 
%     xlabel('Rho xy')
%     ylabel('Rho xyop')
%     title(['Rho xyop vs. Rho xy log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha) ])
% 
% 
% end
% 
% 
% 
% %Plot the traces and predicted traces
% %Rxys above 95 percentile
% iiROI_Rxys_to_plot=[];
% if ~isempty(iiROI_Rxys)
%     to_sort=[Rxys_above_95' iiROI_Rxys'];
%     sorted_rows=sortrows(to_sort,1);
%     iiROI_Rxys_to_plot=zeros(1,length(iiROI_Rxys));
%     iiROI_Rxys_to_plot(1,:)=sorted_rows(:,2);
% 
%     these_alldFF=[];
% 
%     for ii=1:length(iiROI_Rxys)
%         these_alldFF=[these_alldFF per_ROI(iiROI_Rxys(ii)).results.all_dFF];
%     end
% 
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
% 
%     hFig = figure(figNo);
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
%     hold on
% 
%     % Determine the y spacing of the traces
%     y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));
% 
% 
%     for iiROI=1:length(iiROI_Rxys)
%         this_ROI=iiROI_Rxys_to_plot(iiROI);
%         no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
%         no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);
% 
%         plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
%         plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)
% 
%         plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xy+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
%         plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xy+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
%     end
% 
%     %Show the last few
%     ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])
% 
%     xlabel('time(sec)')
%     title(['All dFF timecourses for ROIs with xy predictions above 95 percentile'])
% end
% 
% 
% %Show the fits for xyrops
% iiROI_Rxyops_to_plot=[];
% if ~isempty(iiROI_Rxyops)
%     to_sort=[Rxyops_above_95' iiROI_Rxyops'];
%     sorted_rows=sortrows(to_sort,1);
%     iiROI_Rxyops_to_plot=zeros(1,length(iiROI_Rxyops));
%     iiROI_Rxyops_to_plot(1,:)=sorted_rows(:,2);
% 
% 
% 
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
% 
%     hFig = figure(figNo);
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
%     hold on
% 
%     % Determine the y spacing of the traces
%     y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));
% 
% 
%     for iiROI=1:length(iiROI_Rxyops_to_plot)
%         this_ROI=iiROI_Rxyops_to_plot(iiROI);
%         no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
%         no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);
% 
%         plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
%         plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)
% 
%         plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xyop+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
%         plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xyop+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
%     end
% 
%     %Show the last few
%     ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])
% 
%     xlabel('time(sec)')
%     title(['All dFF timecourses for ROIs with xyop predictions above 95 percentile'])
% end
% 
% 
% 
% %Show the fits for ops
% iiROI_Rops_to_plot=[];
% if ~isempty(iiROI_Rops)
%     to_sort=[Rops_above_95' iiROI_Rops'];
%     sorted_rows=sortrows(to_sort,1);
%     iiROI_Rops_to_plot=zeros(1,length(iiROI_Rops));
%     iiROI_Rops_to_plot(1,:)=sorted_rows(:,2);
% 
% 
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
% 
%     hFig = figure(figNo);
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
%     hold on
% 
%     % Determine the y spacing of the traces
%     y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));
% 
% 
%     for iiROI=1:length(iiROI_Rops_to_plot)
%         this_ROI=iiROI_Rops_to_plot(iiROI);
%         no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
%         no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);
% 
%         plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
%         plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)
% 
%         plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_op+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
%         plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_op+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
%     end
% 
%     %Show the last few
%     ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])
% 
%     xlabel('time(sec)')
%     title(['All dFF timecourses for ROIs with op predictions above 95 percentile'])
% end
% 

%Now do space activity maps and calculate information content
%Now show the pseudocolor activity plots
x=25:50:475;
y=24:48:456;

X_square=zeros(10,10);
Y_square=zeros(10,10);
for ii_x=1:10
    for ii_y=1:10
        X_square(ii_x,ii_y)=x(ii_x);
        Y_square(ii_x,ii_y)=y(ii_y);
    end
end

% %Keep track of all iiROIs
% all_iiROIs_above_95=unique([iiROI_Rxyops_to_plot iiROI_Rxys_to_plot iiROI_Rops_to_plot]);
% all_Rs_above_95=zeros(1,length(all_iiROIs_above_95));
% ii_all_iiROIs=0;
% groups_all_iiROIs_above_95=zeros(length(all_iiROIs_above_95),3); %First column op, second xy, third xyop



%Calculate activity maps and information content
figNo=figNo+1;
information_content=zeros(no_neurons,1);
sparsity=zeros(no_neurons,1);
information_contentl1=zeros(no_neurons,1);
sparsityl1=zeros(no_neurons,1);
information_contentl4=zeros(no_neurons,1);
sparsityl4=zeros(no_neurons,1);
information_contentl1l4=zeros(no_neurons,1);
sparsityl4=zeros(no_neurons,1);

spatial_rhol1l4=zeros(no_neurons,1);
delta_center_of_mass=zeros(no_neurons,1);

sh_information_content=zeros(no_neurons,n_shuffle_SI);
sh_sparsity=zeros(no_neurons,n_shuffle_SI);
sh_information_contentl1=zeros(no_neurons,n_shuffle_SI);
sh_sparsityl1=zeros(no_neurons,n_shuffle_SI);
sh_information_contentl4=zeros(no_neurons,n_shuffle_SI);
sh_sparsityl4=zeros(no_neurons,n_shuffle_SI);
sh_information_contentl1l4=zeros(no_neurons,n_shuffle_SI);
sh_sparsityl1l4=zeros(no_neurons,n_shuffle_SI);

sh_spatial_rhol1l4=zeros(no_neurons,n_shuffle_SI);
sh_delta_center_of_mass=zeros(no_neurons,n_shuffle_SI);

glm_pvalues=[];
for this_ROI=handles_choices.process_these_ROIs

    %Initialize variables
    this_dFF_activity=zeros(10,10);
    this_dFF_activity_n=zeros(10,10);
    sum_dFF_activity=0;

    this_dFFl1_activity=zeros(10,10);
    this_dFFl1_activity_n=zeros(10,10);
    sum_dFFl1_activity=0;

    this_dFFl4_activity=zeros(10,10);
    this_dFFl4_activity_n=zeros(10,10);
    sum_dFFl4_activity=0;

    ii_repeat=1;
    these_all_dFF=(per_ROI(this_ROI).results.all_dFF)-min(per_ROI(this_ROI).results.all_dFF);
    these_all_lanes=per_ROI(this_ROI).results.all_lanes;

    %GLM
    glm_dFF=[];
    glm_dFF_ii=0;

    glm_dFFl1=[];
    glm_dFFl1_ii=0;

    glm_dFFl4=[];
    glm_dFFl4_ii=0;

    glm_dFFl14=[];
    glm_dFFl14_ii=0;

    %Now calculate activity for the odor arena activity map
    for ii_t=1:length(per_ROI(this_ROI).results.all_dFF)

        this_x_ii=ceil(per_ROI(this_ROI).results.all_XY(ii_t,1)/50);
        if this_x_ii==11
            this_x_ii=10;
        end

        this_y_ii=ceil(per_ROI(this_ROI).results.all_XY(ii_t,2)/48);
        if this_y_ii==11
            this_y_ii=10;
        end

        this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
        this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
        sum_dFF_activity=sum_dFF_activity+these_all_dFF(ii_t);

        glm_dFF_ii=glm_dFF_ii+1;
        glm_dFF.data(glm_dFF_ii)=these_all_dFF(ii_t);
        glm_dFF.x(glm_dFF_ii)=this_x_ii;
        glm_dFF.y(glm_dFF_ii)=this_y_ii;
        glm_dFF.xy(glm_dFF_ii)=(this_y_ii-1)*10+this_x_ii;

        if this_y_ii<=3
            glm_dFF.lane(glm_dFF_ii)=0; %This is the lane 4 area in the arena
        else
            if this_y_ii>=7
                glm_dFF.lane(glm_dFF_ii)=1; %This is lane 1 area in the areana
            else
                glm_dFF.lane(glm_dFF_ii)=2; %This is the center
            end
        end

        if these_all_lanes(ii_t)==4
            glm_dFF.lane_trial(glm_dFF_ii)=0; %Lane 4 trial
            glm_dFFl14_ii=glm_dFFl14_ii+1;
            glm_dFFl14.data(glm_dFFl14_ii)=these_all_dFF(ii_t);
            glm_dFFl14.lane_trial(glm_dFFl14_ii)=0; %Lane 4 trial
            glm_dFFl4_ii=glm_dFFl4_ii+1;
            glm_dFFl4.data(glm_dFFl4_ii)=these_all_dFF(ii_t);
            glm_dFFl4.xy(glm_dFFl4_ii)=(this_y_ii-1)*10+this_x_ii;
        else
            glm_dFF.lane_trial(glm_dFF_ii)=1; %Lane 1 trial
            glm_dFFl14_ii=glm_dFFl14_ii+1;
            glm_dFFl14.data(glm_dFFl14_ii)=these_all_dFF(ii_t);
            glm_dFFl14.lane_trial(glm_dFFl14_ii)=1; %Lane 4 trial
            glm_dFFl1_ii=glm_dFFl1_ii+1;
            glm_dFFl1.data(glm_dFFl1_ii)=these_all_dFF(ii_t);
            glm_dFFl1.xy(glm_dFFl1_ii)=(this_y_ii-1)*10+this_x_ii;
        end

        if these_all_lanes(ii_t)==1
            this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
            this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFFl1_activity=sum_dFFl1_activity+these_all_dFF(ii_t);
        else
            this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
            this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFFl4_activity=sum_dFFl4_activity+these_all_dFF(ii_t);
        end
    end

    % %Perform the glm
    % fprintf(1, ['glm for lane 1 vs lane 4 trials with x and y\n'])
    %
    % fprintf(1, ['\n\nglm for dFF\n'])
    % tbl = table(glm_dFF.data',glm_dFF.x',glm_dFF.y',glm_dFF.lane_trial',...
    %     'VariableNames',{'dFF','x','y','trial'});
    % mdl = fitglm(tbl,'dFF~x+y+trial'...
    %     ,'CategoricalVars',[4])



    %Perform the glm
    fprintf(1, ['glm for lane 1 vs lane 4 trials with xy\n'])

    fprintf(1, ['\n\nglm for dFF\n'])
    tbl = table(glm_dFF.data',glm_dFF.xy',glm_dFF.lane_trial',...
        'VariableNames',{'dFF','xy','trial'});
    mdl = fitglm(tbl,'dFF~xy+trial'...
        ,'CategoricalVars',[3])

    glm_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;

     %Perform the glm
    fprintf(1, ['glm for lane 1 vs lane 4 trials without xy\n'])

    fprintf(1, ['\n\nglm for dFF lane 1 vs 4, no xy\n'])
    tbl = table(glm_dFFl14.data',glm_dFFl14.lane_trial',...
        'VariableNames',{'dFF','trial'});
    mdl = fitglm(tbl,'dFF~trial'...
        ,'CategoricalVars',[2])

    glm_l14_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;

    %Perform the lane 1 glm
    fprintf(1, ['glm for lane 1 trials with xy\n'])

    fprintf(1, ['\n\nglm for dFF lane 1\n'])
    tbl = table(glm_dFFl1.data',glm_dFFl1.xy',...
        'VariableNames',{'dFF','xy'});
    mdl = fitglm(tbl,'dFF~xy')

    glm_l1_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;

    %Perform the lane 4 glm
    fprintf(1, ['glm for lane 4 trials with xy\n'])

    fprintf(1, ['\n\nglm for dFF lane 4\n'])
    tbl = table(glm_dFFl4.data',glm_dFFl4.xy',...
        'VariableNames',{'dFF','xy'});
    mdl = fitglm(tbl,'dFF~xy')

    glm_l4_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;

    % %Do the ranksum/t-test
    % fprintf(1, ['\n\nRanksum or t-test p values for average MI for each electrode calculated per mouse for PAC theta' freq_names{pacii+1} '\n'])
    % [output_data] = drgMutiRanksumorTtest(input_data);

    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)~=0
                this_dFF_activity(ii_x,ii_y)=this_dFF_activity(ii_x,ii_y)/this_dFF_activity_n(ii_x,ii_y);
                if this_dFFl1_activity_n(ii_x,ii_y)>0
                    this_dFFl1_activity(ii_x,ii_y)=this_dFFl1_activity(ii_x,ii_y)/this_dFFl1_activity_n(ii_x,ii_y);
                end
                if this_dFFl4_activity_n(ii_x,ii_y)>0
                    this_dFFl4_activity(ii_x,ii_y)=this_dFFl4_activity(ii_x,ii_y)/this_dFFl4_activity_n(ii_x,ii_y);
                end
            end
        end
    end


    %
    % this_dFF_activity=this_dFF_activity/sum_dFF_activity;
    % this_dFFl1_activity=((no_trials/2)/trials.lane1)*this_dFFl1_activity/(sum_dFFl1_activity+sum_dFFl4_activity);
    % this_dFFl4_activity=((no_trials/2)/trials.lane4)*this_dFFl4_activity/(sum_dFFl1_activity+sum_dFFl4_activity);


    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .8 .6])

    %Plot the fits
    no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
    no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);

    y_shift=max(per_ROI(this_ROI).results.all_dFF);
    if y_shift==0
        y_shift=1;
    end

    %Plot dFF timecourses
    subplot(2,3,1:3)
    hold on
    ii_plot=0;
    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*ii_plot,'-k','LineWidth',1.5)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*ii_plot,'-k','LineWidth',1.5)

    % plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_op+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    % plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_op+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)

    % ii_plot=1;
    % plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*ii_plot,'-k','LineWidth',1.5)
    % plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*ii_plot,'-k','LineWidth',1.5)
    % 
    % % plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xy+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    % % plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xy+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    % 
    % ii_plot=2;
    % plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*ii_plot,'-k','LineWidth',1.5)
    % plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*ii_plot,'-k','LineWidth',1.5)

    % plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xyop+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    % plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xyop+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)

    %Lane 1
    for ii_tr=1:length(per_ROI(this_ROI).results.all_trl1)
        this_x=per_ROI(this_ROI).results.all_trl1(ii_tr);
        plot([this_x this_x], [-0.2*y_shift y_shift*(ii_plot+1)],'-','Color',[0.6 0.6 0.6])
    end

    %Lane 4
    for ii_tr=1:length(per_ROI(this_ROI).results.all_trl4)
        this_x=no_time_bins1+100+per_ROI(this_ROI).results.all_trl4(ii_tr);
        plot([this_x this_x], [-0.2*y_shift y_shift*(ii_plot+1)],'-','Color',[0.6 0.6 0.6])
    end


    ylim([-0.2*y_shift 1.2*y_shift*(ii_plot+1)])


    % 
    % %Show Rhos
    % 
    % %Rop
    % if sum(this_ROI==iiROI_Rops_to_plot)==1
    %     text(no_time_bins1+10,0,['op ' num2str(per_ROI(this_ROI).results.Rop)], 'Color','k')
    % else
    %     text(no_time_bins1+10,0,['op ' num2str(per_ROI(this_ROI).results.Rop)], 'Color',[0.6 0.6 0.6])
    % end
    % 
    % %Rxy
    % if sum(this_ROI==iiROI_Rxys_to_plot)==1
    %     text(no_time_bins1+10,y_shift,['xy ' num2str(per_ROI(this_ROI).results.Rxy)], 'Color','k')
    % else
    %     text(no_time_bins1+10,y_shift,['xy ' num2str(per_ROI(this_ROI).results.Rxy)], 'Color',[0.6 0.6 0.6])
    % end
    % 
    % %Rxyop
    % if sum(this_ROI==iiROI_Rxyops_to_plot)==1
    %     text(no_time_bins1+10,2*y_shift,['xyop ' num2str(per_ROI(this_ROI).results.Rxyop)], 'Color','k')
    % else
    %     text(no_time_bins1+10,2*y_shift,['xyop ' num2str(per_ROI(this_ROI).results.Rxyop)], 'Color',[0.6 0.6 0.6])
    % end

    %Entire arena
    subplot(2,3,4)
    max_this_dFF_activity=max(this_dFF_activity(:));
    max_this_dFFl4_activity=max(this_dFFl4_activity(:));
    max_this_dFFl1_activity=max(this_dFFl1_activity(:));
    max_activity=max([max_this_dFF_activity max_this_dFFl1_activity max_this_dFFl4_activity]);
    if max_activity==0
        max_activity=1;
    end
    % min_this_dFF_activity=min(this_dFF_activity(:));
    delta_ac=max_activity/255;
    this_masked_dFF_activity=this_dFF_activity;
    this_masked_dFF_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
    drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFF_activity)
    colormap fire
    this_cmap=colormap;
    this_cmap(1,:)=[0.2 0.2 0.2];
    colormap(this_cmap)
    clim([-1.5*delta_ac max_activity])

    shading interp
    set(gca, 'YDir', 'reverse');

    % caxis([minC maxC]);
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    yticks(0:48:480)
    xticks(0:50:500)
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('All trials')

    %Lane 1
    subplot(2,3,5)
    % max_this_dFFl1_activity=max(this_dFFl1_activity(:));
    % min_this_dFFl1_activity=min(this_dFFl1_activity(:));
    % delta_ac=(max_this_dFFl1_activity-min_this_dFFl1_activity)/255;
    this_masked_dFFl1_activity=this_dFFl1_activity;
    this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
    % drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity')
    drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
    colormap(this_cmap)
    clim([-1.5*delta_ac max_activity])
    shading interp
    set(gca, 'YDir', 'reverse');

    % caxis([minC maxC]);
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    yticks(0:48:480)
    xticks(0:50:500)
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Lane 1')

    %Lane 4
    subplot(2,3,6)

    % min_this_dFFl4_activity=min(this_dFFl4_activity(:));
    % delta_ac=(max_this_dFFl4_activity-min_this_dFFl4_activity)/255;
    this_masked_dFFl4_activity=this_dFFl4_activity;
    this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
    drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
    colormap(this_cmap)
    clim([-1.5*delta_ac max_activity])
    shading interp
    set(gca, 'YDir', 'reverse');

    % caxis([minC maxC]);
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    yticks(0:48:480)
    xticks(0:50:500)
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Lane 4')


    sgtitle(['dFF activity map for ROI No ' num2str(this_ROI)])

    pffft=1;

    sh_maps=[];
    for ii_sh=1:n_shuffle_SI
        %Now calculate the activity map after shuffling
        shuffle_length=floor(length(per_ROI(this_ROI).results.all_dFF)/(trials.odor_trNo*2));

        this_overlap=10;
        while this_overlap>max_overlap
            perm_indices=randperm(trials.odor_trNo*2);
            this_overlap=0;
            for ii_p=1:trials.odor_trNo*2
                if perm_indices(ii_p)==ii_p
                    this_overlap=this_overlap+1;
                end
            end
        end

        XYsh=zeros(length(per_ROI(this_ROI).results.all_dFF),2);
        % lane_train=zeros(size(XYtrain_pre,1),1);
        thisXY=per_ROI(this_ROI).results.all_XY;
        these_all_dFF=per_ROI(this_ROI).results.all_dFF;
        max_ii=0;
        for sh_ii=1:trials.odor_trNo*2
            XYsh((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
                thisXY((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
        end

        if sh_ii*shuffle_length+1<=size(thisXY,1)
            XYsh(sh_ii*shuffle_length+1:end,:)=thisXY(1:size(thisXY,1)-sh_ii*shuffle_length,:);
        end

        %Initialize variables
        sh_maps(ii_sh).this_dFF_activity=zeros(10,10);
        sh_maps(ii_sh).this_dFF_activity_n=zeros(10,10);
        sh_maps(ii_sh).sum_dFF_activity=0;

        sh_maps(ii_sh).this_dFFl1_activity=zeros(10,10);
        sh_maps(ii_sh).this_dFFl1_activity_n=zeros(10,10);
        sh_maps(ii_sh).sum_dFFl1_activity=0;

        sh_maps(ii_sh).this_dFFl4_activity=zeros(10,10);
        sh_maps(ii_sh).this_dFFl4_activity_n=zeros(10,10);
        sh_maps(ii_sh).sum_dFFl4_activity=0;

        these_all_dFF=(per_ROI(this_ROI).results.all_dFF)-min(per_ROI(this_ROI).results.all_dFF);
        these_all_lanes=per_ROI(this_ROI).results.all_lanes;

        %Now calculate shuffled activity for the odor arena activity map
        for ii_t=1:length(per_ROI(this_ROI).results.all_dFF)

            this_x_ii=ceil(XYsh(ii_t,1)/50);
            if this_x_ii==11
                this_x_ii=10;
            end

            this_y_ii=ceil(XYsh(ii_t,2)/48);
            if this_y_ii==11
                this_y_ii=10;
            end

            sh_maps(ii_sh).this_dFF_activity(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFF_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
            sh_maps(ii_sh).this_dFF_activity_n(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFF_activity_n(this_x_ii,this_y_ii)+1;
            sh_maps(ii_sh).sum_dFF_activity=sum_dFF_activity+these_all_dFF(ii_t);

            if these_all_lanes(ii_t)==1
                sh_maps(ii_sh).this_dFFl1_activity(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFFl1_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
                sh_maps(ii_sh).this_dFFl1_activity_n(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
                sh_maps(ii_sh).sum_dFFl1_activity=sum_dFFl1_activity+these_all_dFF(ii_t);
            else
                sh_maps(ii_sh).this_dFFl4_activity(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFFl4_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
                sh_maps(ii_sh).this_dFFl4_activity_n(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
                sh_maps(ii_sh).sum_dFFl4_activity=sum_dFFl4_activity+these_all_dFF(ii_t);
            end
        end

        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)~=0
                    sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)=sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)/sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y);
                    if sh_maps(ii_sh).this_dFFl1_activity_n(ii_x,ii_y)>0
                        sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)=sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)/sh_maps(ii_sh).this_dFFl1_activity_n(ii_x,ii_y);
                    end
                    if sh_maps(ii_sh).this_dFFl4_activity_n(ii_x,ii_y)>0
                        sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)=sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)/sh_maps(ii_sh).this_dFFl4_activity_n(ii_x,ii_y);
                    end
                end
            end
        end


    end

    %Calculate information content
    %I will use information theory Markus et al 1994 https://doi.org/10.1002/hipo.450040404,


    %information content = sum(Pi(Ri/R)log2(Ri/R))
    %sparsity = sum((Pi*Ri^2)/R^2)



    %All trials
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    overall_mean_activity=mean(this_dFF_activity(:));
    this_information_content=0;
    this_sparsity=0;
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFF_activity_n(ii_x,ii_y)/sum(this_dFF_activity_n(:));
                if this_dFF_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFF_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFF_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFF_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end
    information_content(this_ROI)=this_information_content;
    sparsity(this_ROI)=this_sparsity;

    %Lane 1
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    this_information_content=0;
    this_sparsity=0;
    overall_mean_activity=mean(this_dFFl1_activity(:));
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFFl1_activity_n(ii_x,ii_y)/sum(this_dFFl1_activity_n(:));
                if this_dFFl1_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFFl1_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end
    information_contentl1(this_ROI)=this_information_content;
    sparsityl1(this_ROI)=this_sparsity;

    %Lane 4
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    this_information_content=0;
    this_sparsity=0;
    overall_mean_activity=mean(this_dFFl4_activity(:));
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFFl4_activity_n(ii_x,ii_y)/sum(this_dFFl4_activity_n(:));
                if this_dFFl4_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFFl4_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end
    information_contentl4(this_ROI)=this_information_content;
    sparsityl4(this_ROI)=this_sparsity;

    %Information content lane 1 lane 4
    %Lane 1
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    this_information_content=0;
    this_sparsity=0;
    overall_mean_activity=mean([this_dFFl1_activity(:); this_dFFl4_activity(:)]);
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFFl1_activity_n(ii_x,ii_y)/sum(this_dFF_activity_n(:));
                if this_dFFl1_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFFl1_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end


    %Lane 4
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFFl4_activity_n(ii_x,ii_y)/sum(this_dFF_activity_n(:));
                if this_dFFl4_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFFl4_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end
    information_contentl1l4(this_ROI)=this_information_content;
    sparsityl1l4(this_ROI)=this_sparsity;


    %Spatial correlation between lanes 1 and 4
    this_rho=corrcoef([this_dFFl1_activity(this_dFF_activity_n>0) this_dFFl4_activity(this_dFF_activity_n>0)]);
    spatial_rhol1l4(this_ROI)=this_rho(1,2);

    %Center of mass
    cmx1=sum(this_dFFl1_activity(this_dFF_activity_n>0).*X_square(this_dFF_activity_n>0))/sum(this_dFFl1_activity(this_dFF_activity_n>0));
    cmx4=sum(this_dFFl4_activity(this_dFF_activity_n>0).*X_square(this_dFF_activity_n>0))/sum(this_dFFl4_activity(this_dFF_activity_n>0));

    cmy1=sum(this_dFFl1_activity(this_dFF_activity_n>0).*Y_square(this_dFF_activity_n>0))/sum(this_dFFl1_activity(this_dFF_activity_n>0));
    cmy4=sum(this_dFFl4_activity(this_dFF_activity_n>0).*Y_square(this_dFF_activity_n>0))/sum(this_dFFl4_activity(this_dFF_activity_n>0));

    delta_center_of_mass(this_ROI)=sqrt((cmx1-cmx4)^2 +(cmy1-cmy4)^2);

    %Now calculate the shuffled values
    these_sh_information_content=zeros(n_shuffle_SI,1);
    these_sh_sparsity=zeros(n_shuffle_SI,1);
    these_sh_information_contentl1=zeros(n_shuffle_SI,1);
    these_sh_sparsityl1=zeros(n_shuffle_SI,1);
    these_sh_information_contentl4=zeros(n_shuffle_SI,1);
    these_sh_sparsityl4=zeros(n_shuffle_SI,1);
    these_sh_spatial_rhol1l4=zeros(n_shuffle_SI,1);
    these_sh_delta_center_of_mass=zeros(n_shuffle_SI,1);

    for ii_sh=1:n_shuffle_SI



        %All trials
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory
        overall_mean_activity=mean(sh_maps(ii_sh).this_dFF_activity(:));
        this_information_content=0;
        this_sparsity=0;
        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end
        these_sh_information_content(ii_sh,1)=this_information_content;
        these_sh_sparsity(ii_sh,1)=this_sparsity;

        %Lane 1
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory
        this_information_content=0;
        this_sparsity=0;
        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFFl1_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end
        these_sh_information_contentl1(ii_sh,1)=this_information_content;
        these_sh_sparsityl1(ii_sh,1)=this_sparsity;

        %Lane 4
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory
        this_information_content=0;
        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFFl4_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end
        these_sh_information_contentl4(ii_sh,1)=this_information_content;
        these_sh_sparsityl4(ii_sh,1)=this_sparsity;

        %Lane 1 and lane 4
        %Lane 1
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory
        this_information_content=0;
        this_sparsity=0;
        overall_mean_activity=mean([sh_maps(ii_sh).this_dFFl1_activity(:); sh_maps(ii_sh).this_dFFl4_activity(:)]);
        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFFl1_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFFl1_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end


        %Lane 4
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory

        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFFl4_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFFl4_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end
        these_sh_information_contentl4(ii_sh,1)=this_information_content;
        these_sh_sparsityl4(ii_sh,1)=this_sparsity;


        %Spatial correlation between lanes 1 and 4
        this_rho=corrcoef([sh_maps(ii_sh).this_dFFl1_activity(sh_maps(ii_sh).this_dFF_activity_n>0) sh_maps(ii_sh).this_dFFl4_activity(sh_maps(ii_sh).this_dFF_activity_n>0)]);
        these_sh_spatial_rhol1l4(ii_sh)=this_rho(1,2);

        %Delta center of mass
        cmx1=sum(sh_maps(ii_sh).this_dFFl1_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*X_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFFl1_activity(sh_maps(ii_sh).this_dFF_activity_n>0));
        cmx4=sum(sh_maps(ii_sh).this_dFFl4_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*X_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFFl4_activity(this_dFF_activity_n>0));

        cmy1=sum(sh_maps(ii_sh).this_dFFl1_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*Y_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFFl1_activity(sh_maps(ii_sh).this_dFF_activity_n>0));
        cmy4=sum(sh_maps(ii_sh).this_dFFl4_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*Y_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFFl4_activity(sh_maps(ii_sh).this_dFF_activity_n>0));

        these_sh_delta_center_of_mass(ii_sh)=sqrt((cmx1-cmx4)^2 +(cmy1-cmy4)^2);
    end

    sh_information_content(this_ROI,:)=these_sh_information_content;
    sh_sparsity(this_ROI,:)=these_sh_sparsity;
    sh_information_contentl1(this_ROI,:)=these_sh_information_contentl1;
    sh_sparsityl1(this_ROI,:)=these_sh_sparsityl1;
    sh_information_contentl4(this_ROI,:)=these_sh_information_contentl4;
    sh_sparsityl4(this_ROI,:)=these_sh_sparsityl4;
    sh_spatial_rhol1l4(this_ROI,:)=these_sh_spatial_rhol1l4;
    sh_delta_center_of_mass(this_ROI,:)=these_sh_delta_center_of_mass;

    pffft=1;
end

% handles_out.Rops=Rops;
% handles_out.Rxys=Rxys;
% handles_out.Rxyops=Rxyops;
% handles_out.iiROI_Rops=iiROI_Rops;
% handles_out.iiROI_Rxys=iiROI_Rxys;
% handles_out.iiROI_Rxyops=iiROI_Rxyops;
% handles_out.Rops_above_95=Rops_above_95;
% handles_out.Rxys_above_95=Rxys_above_95;
% handles_out.Rxyops_above_95=Rxyops_above_95;
% handles_out.Rops_for_Rxys_above_95=Rops_for_Rxys_above_95;
% handles_out.Roxyl1s_for_Rxys_above_95=Roxyl1s_for_Rxys_above_95;
% handles_out.Roxyl4s_for_Rxys_above_95=Roxyl4s_for_Rxys_above_95;
% handles_out.Rxys_for_Rops_above_95=Rxys_for_Rops_above_95;
% handles_out.Ropl1s_for_Rops_above_95=Ropl1s_for_Rops_above_95;
% handles_out.Ropl4s_for_Rops_above_95=Ropl4s_for_Rops_above_95;
% handles_out.Rxys_for_Rxyops_above_95=Rxys_for_Rxyops_above_95;
% handles_out.Rxyopl1s_for_Rxyops_above_95=Rxyopl1s_for_Rxyops_above_95;
% handles_out.Rxyopl4s_for_Rxyops_above_95=Rxyopl4s_for_Rxyops_above_95;
% handles_out.Rxyops_for_Rxys_above_95=Rxyops_for_Rxys_above_95;


handles_out.information_content=information_content;
handles_out.sparsity=sparsity;
handles_out.information_contentl1=information_contentl1;
handles_out.sparsityl1=sparsityl1;
handles_out.information_contentl4=information_contentl4;
handles_out.sparsityl4=sparsityl4;
handles_out.information_contentl1l4=information_contentl4;
handles_out.sparsityl4=sparsityl1l4;
handles_out.spatial_rhol1l4=spatial_rhol1l4;
handles_out.delta_center_of_mass=delta_center_of_mass;
handles_out.sh_information_content=sh_information_content;
handles_out.sh_sparsity=sh_sparsity;
handles_out.sh_information_contentl1=sh_information_contentl1;
handles_out.sh_sparsityl1=sh_sparsityl1;
handles_out.sh_information_contentl4=sh_information_contentl4;
handles_out.sh_sparsityl4=sh_sparsityl4;
handles_out.sh_information_contentl1l4=sh_information_contentl1l4;
handles_out.sh_sparsityl1l4=sh_sparsityl1l4;
handles_out.sh_spatial_rhol1l4=sh_spatial_rhol1l4;
handles_out.sh_delta_center_of_mass=sh_delta_center_of_mass;

handles_out.glm_pvalues=glm_pvalues;
handles_out.glm_l14_pvalues=glm_l14_pvalues;
handles_out.glm_l1_pvalues=glm_l1_pvalues;
handles_out.glm_l4_pvalues=glm_l4_pvalues;



save([this_path arena_file(1:end-4) '_' handles_choices.save_tag '_' num2str(handles_choices.algo) ...
    num2str(handles_choices.weber_fechner) num2str(handles_choices.alpha)...
    num2str(no_dec_time_bins_op) num2str(no_dec_time_bins_xy) '.mat'],'handles_out','handles_choices','trials','-v7.3')


pffft=1;