function handles_out=drgMini_Inspect_conc_xy_traces(handles_choices)
%Does decoding of the odor concentration for a mouse undergoing odor plume navigation 
%following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020


close all

if exist('handles_choices')==0
    clear all

    handles_choices.is_sphgpu=0;
    is_sphgpu=handles_choices.is_sphgpu;

    %Troubleshooting in sphgpu with file 20
    % this_path='/data2/SFTP/PreProcessed/20221117_FCM22_lanes_1_4/';
    % dFF_file='20221117_FCM22_withodor_nearfloor_miniscope_sync_L1andL4_ncorre_fix_ext.mat';
    % arena_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm.mat';
    % 
    % handles_choices.save_path='/data2/SFTP/PreProcessed/test/';

    % %Troubleshooting in Mac with file 3
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220727_FCM19/';
    % dFF_file='20220727_FCM19_withodor_miniscope_sync_L1andL4_ncorre_ext_nonneg.mat';
    % arena_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm.mat';

    %Troubleshooting in Mac with file 20
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20221117_FCM22_lanes_1_4/';
    dFF_file='20221117_FCM22_withodor_nearfloor_miniscope_sync_L1andL4_ncorre_fix_ext.mat';
    arena_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm.mat'
    handles_choices.save_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/test/';

    %Troubleshooting Fabio's files May 14th
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    % dFF_file='20220729_FCM22_withodor_miniscope_sync_L4_ncorre_ext_nonneg.mat';
    % arena_file='20220729_FCM22withodor_odorarena_L4_sync.mat';


    % %First troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    % dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';
    % 
    % handles_choices.save_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/test/';

    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';

    %     %Second troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';

    %First troubleshooting file odor in lanes 1 and 4
    % this_path='/data/SFTP/PreProcessedDR/20220713_FCM6/';
    % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';

    % % Second troubleshooting files odor in lanes 1 and 4
    % this_path='/data/SFTP/PreProcessedDR/20220804_FCM22/';
    % dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';

    %Third troubleshooting odor in lane 4 only ISO1
    % this_path='/data/SFTP/PreProcessedDR/20220526_FCM6_withodor_lane4/'
    % dFF_file='20220526_FCM6_withodor_miniscope_sync_L4_ncorre_ext_nonneg.mat';
    % arena_file='20220526_FCM6withodor_odorarena_L4_sync_mm.mat';


    handles_choices.group=1;

    %Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4
    %Group 2 is rewarded, with odor lane 4, no odor in lane 1
    %Group 3 is rewarded, with odor lane 1, no odor in lane 4
    %Group 4 is rewarded, with no odor in lane 1 and lane 4
    %No odor troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';

    

    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

   

    %     isKording=0;

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    %     dt=0.2;

    %Define the different ranges (training, valid and testing)
    % training_fraction=0.9;
    % handles_choices.training_fraction=training_fraction;

    %     training_range=[0, 0.5];
    %     valid_range=[0.5,0.65];
    %     test_range=[0.5, 1];

    %The user can define what time period to use spikes from (with respect to the output).
    bins_before=0; %How many bins of neural data prior to the output are used for decoding, 10
    bins_current=1; %Whether to use concurrent time bin of neural data, 1
    bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    handles_choices.bins_before=bins_before;
    handles_choices.bins_current=bins_current;
    handles_choices.bins_after=bins_after;

    %Speed for background air flow in mm/sec
    air_flow_speed=50; %5 cm/sec = 50 mm/sec
    handles.air_flow_speed=air_flow_speed;

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.1; %Time bins for decoding
    dt_miniscope=1/30;
    n_shuffle=5; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    handles_choices.n_shuffle=n_shuffle;

    which_training_algorithm=3;
    handles_choices.which_training_algorithm=which_training_algorithm;
    %1=fitrnet
    %2=fitrgp
    %3=fitrtree
    %4=fitglm
    %5=fitrsvm
    %6=fitrtree gpu

    handles_choices.cm_from_floor=2;

    % handles_choices.weber_fechner=1;
    % handles_choices.alpha=1;
    % handles_choices.multiplier=1;
    % handles_choices.lowest_conc=-20;
    % %0 is Stevens Law, R proportional to C^alpha
    % %1 is Weber-Flechner law R proportional to multiplier*log(C)
    % %See Copelli et al DOI: 10.1103/PhysRevE.65.060901


   %Weber-Frechner or Stevens
    handles_choices.weber_fechner=1; 
    %0 is Stevens Law, R proportional to C^alpha
    %1 is Weber-Flechner law R proportional to multiplier*log(C)
    %See Copelli et al DOI: 10.1103/PhysRevE.65.060901

    handles_choices.alpha=1;
    handles_choices.multiplier=1;
    handles_choices.lowest_conc=-200;

    %Hill transform
    handles_choices.hill=0; %0=no Hill transform, 1=Hill transform
    handles_choices.k_half=10^-8; %Hill equation K1/2
    % handles_choices.actual_maxC=(handles_choices.k_half^handles_choices.n_hill); %Measured from the simulated data 0.0043
    
    handles_choices.maxC=10^-6.5;%maxC=0.0043 maximum of simulated odor plume
    handles_choices.n_hill=2; %Hill coefficient


    handles_choices.displayFigures=1;

    handles_choices.trial_start_offset=-15; %This was -10
    handles_choices.trial_end_offset=15;

    handles_choices.train_with_hits=1; %0=train with all trials, 1=train with hit trials, 2=train with miss trials

    handles_choices.showTheseROIs=[1:10];

    % handles_choices.save_tag='OdorConc';


else

    is_sphgpu=handles_choices.is_sphgpu;
    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;
    % training_fraction=handles_choices.training_fraction;
    bins_before=handles_choices.bins_before;
    bins_current=handles_choices.bins_current;
    bins_after=handles_choices.bins_after;
    dt=handles_choices.dt;
    dt_miniscope=handles_choices.dt_miniscope;
    n_shuffle=handles_choices.n_shuffle;
    which_training_algorithm=handles_choices.which_training_algorithm;
    air_flow_speed=handles_choices.air_flow_speed;

end

if is_sphgpu==1
    addpath('/data2/DRMatlab/drgMiniscope')
    addpath('/data2/DRMatlab/m new/Chi Squared')
    addpath('/data2/DRMatlab/drgMaster')
    addpath(genpath('/data2/DRMatlab/m new/kakearney-boundedline-pkg-32f2a1f'))
end

if ~exist(handles_choices.save_path(1:end-1),'dir')
    mkdir(handles_choices.save_path(1:end-1))
end

try
    delete(gcp('nocreate'));
catch
end

setenv('MW_PCT_TRANSPORT_HEARTBEAT_INTERVAL', '700')

switch which_training_algorithm
    case 1
        fprintf(1,['\nArtificial neural network\n\n'])
    case 2
        fprintf(1,['\nGaussian process regression (inherently Bayesian)\n\n'])
    case 3
        fprintf(1,['\nTree\n\n'])
end

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
    case {1,5}
        %Groups 1 and 5 are rewarded, odor ISO1 in both lane 1 and lane 4
        %They differe in cm from floor
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
        %Note that we are testing whether we can "predict" the odor when
        %there is no odor applied
        handles_choices.lane1_odor_on=1;
        handles_choices.lane4_odor_on=1;
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

%Hill transform?
if handles_choices.hill==1
    mean_plume_l4=handles_choices.maxC*(mean_plume_l4.^handles_choices.n_hill)./...
        ((mean_plume_l4.^handles_choices.n_hill)+(handles_choices.k_half.^handles_choices.n_hill));
     mean_plume_l1=handles_choices.maxC*(mean_plume_l1.^handles_choices.n_hill)./...
        ((mean_plume_l1.^handles_choices.n_hill)+(handles_choices.k_half.^handles_choices.n_hill));
end

%Shift the plume to 7 cm (70 mm)
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

% minC=min([min(mean_plume_l4(:)) min(mean_plume_l1(:))]);
minC=prctile([mean_plume_l4(:); mean_plume_l1(:)],0.5);
if minC<handles_choices.lowest_conc
    minC=handles_choices.lowest_conc;
end
maxC=max([max(mean_plume_l4(:)) max(mean_plume_l1(:))]);
if minC==maxC
    maxC=minC+0.1;
end

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

%Old code
% new_mean_plume_l1=mean_plume_l1(y_for_plume>=20,:);
% mean_plume_l1=[];
% mean_plume_l1=new_mean_plume_l1;
% new_mean_plume_l4=mean_plume_l4(y_for_plume<=480,:);
% mean_plume_l4=new_mean_plume_l4;
% 
% new_y_for_plume=y_for_plume(:,y_for_plume<=480);
% y_for_plume=[];
% y_for_plume=new_y_for_plume;
% 
%Now shift the center of lane 4 to 70 mm
% new_mean_plume_l4(y_for_plume>=20,:)=mean_plume_l4(y_for_plume<=460,:);
% from_y_ii=find(y_for_plume<70,1,'first');
% new_mean_plume_l4(length(y_for_plume):-1:from_y_ii,:)=mean_plume_l4((y_for_plume>50)&(y_for_plume<=120),:);
% mean_plume_l4=new_mean_plume_l4;

%There are slight differences between lane 1 and 4 in the low odorant areas
%Here I will merge the data for the simulated lane 1 and lane 4 into a
%lane14 and then I will use those data for the shifted lanes

%The dimensions of the chamber are 50 cm for x and 48 cm for y
%Cut the simulated plumes to 480 mm in the y dimension
new_mean_plume_l1=mean_plume_l1(y_for_plume>=20,:);
mean_plume_l1=[];
mean_plume_l1=new_mean_plume_l1;
new_mean_plume_l4=mean_plume_l4(y_for_plume<=480,:);
mean_plume_l4=new_mean_plume_l4;

new_y_for_plume=y_for_plume(:,y_for_plume<=480);
y_for_plume=[];
y_for_plume=new_y_for_plume;

%Now make a merged lane 1 and 4

%Reflect plume l1
y_reflected_mean_plume_l1=zeros(size(mean_plume_l1,1),size(mean_plume_l1,2));
for ii_y=1:size(mean_plume_l1,1)
    y_reflected_mean_plume_l1(ii_y,:)=mean_plume_l1(size(mean_plume_l1,1)-ii_y+1,:);
end

%Mean14 plume l4
mean14_plume_l4=(y_reflected_mean_plume_l1+mean_plume_l4)/2;

%Reflect plume l4
mean14_plume_l1=zeros(size(mean_plume_l1,1),size(mean_plume_l1,2));
for ii_y=1:size(mean_plume_l1,1)
    mean14_plume_l1(ii_y,:)=mean14_plume_l4(size(mean_plume_l1,1)-ii_y+1,:);
end

%Now shift the center of lane 4 to 70 mm
new_mean14_plume_l4(y_for_plume>=20,:)=mean14_plume_l4(y_for_plume<=460,:);
from_y_ii=find(y_for_plume<70,1,'first');
new_mean14_plume_l4(length(y_for_plume):-1:from_y_ii,:)=mean14_plume_l4((y_for_plume>50)&(y_for_plume<=120),:);
mean14_plume_l4=new_mean14_plume_l4;

%Replace with mean
mean_plume_l4=mean14_plume_l4;
mean_plume_l1=mean14_plume_l1;

%Note that I am replacing minC with prctile because there are a small
%number of points with very low min
minC=prctile([mean14_plume_l4(:); mean14_plume_l1(:)],0.5);
maxC=min([max(mean14_plume_l4(:)) max(mean14_plume_l1(:))]);
if maxC==minC
    maxC=minC+0.1;
end

mean_plume_l4(mean_plume_l4<minC)=minC;
mean_plume_l1(mean_plume_l1<minC)=minC;

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

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.83 .5 .05 .3])

    prain=[minC:(maxC-minC)/99:maxC];
    drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')
end


%Restart random seeds
rng('shuffle');

% gcp;

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


load([this_path arena_file])


%Extract trials
trials=[];

no_time_points=length(arena.xsync);

%Extract lanes using FLIR data
at_end=0;
ii=0;
trNo=0;
trNo_l1=0;
trNo_l4=0;
while at_end==0
    next_ii=find(arena.odor(ii+1:end)==1,1,'first');
    if ~isempty(next_ii)
        trNo=trNo+1;
        % trials.odor_ii(trNo)=ii+next_ii;
        % trials.x_odor(trNo)=arena.xsync(ii+next_ii);
        % trials.y_odor(trNo)=arena.ysync(ii+next_ii);

        ii=ii+next_ii;
        % ii_mini=arena.index_flirsynctominiscope(ii);

        if sum(arena.laneodor1(ii-3:ii+3)==1)>0
            %Note: laneodor1 is 1 only for one time point
            trials.odor_lane(trNo)=1;
        end

        if sum(arena.laneodor4(ii-3:ii+3)==1)>0
            %Note: laneodor4 is 1 only for one time point
            trials.odor_lane(trNo)=4;
        end

        next_ii=find(arena.odor(ii+1:end)==0,1,'first');
        if ~isempty(next_ii)
            ii=ii+next_ii;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end


%Extract odor on using the miniscope sync data
at_end=0;
ii=0;
trNo=0;
trNo_l1=0;
trNo_l4=0;
while at_end==0
    next_ii=find(arena.odorsync(ii+1:end)==1,1,'first');
    if ~isempty(next_ii)
        trNo=trNo+1;
        trials.odor_ii(trNo)=ii+next_ii;
        trials.x_odor(trNo)=arena.xsync(ii+next_ii);
        trials.y_odor(trNo)=arena.ysync(ii+next_ii);

        ii=ii+next_ii;
       
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

trials.odor_trNo=trNo;

%Extract lanewater using miniscope sync
for trNo=1:trials.odor_trNo
    if trNo<trials.odor_trNo
        is_water=find(arena.watersync(trials.odor_ii(trNo):trials.odor_ii(trNo+1)-1)==1,1,'first');
        if isempty(is_water)
            trials.water_per_trial(trNo)=0;
        else
            trials.water_per_trial(trNo)=1;
            trials.water_per_trial_ii(trNo)=trials.odor_ii(trNo)+is_water;
            trials.water_per_trial_x(trNo)=arena.xsync(trials.odor_ii(trNo)+is_water);
            trials.water_per_trial_y(trNo)=arena.ysync(trials.odor_ii(trNo)+is_water);
        end
    else
        is_water=find(arena.watersync(trials.odor_ii(trNo):end)==1,1,'first');
        if isempty(is_water)
            trials.water_per_trial(trNo)=0;
        else
            trials.water_per_trial(trNo)=1;
            trials.water_per_trial_ii(trNo)=trials.odor_ii(trNo)+is_water;
            trials.water_per_trial_x(trNo)=arena.xsync(trials.odor_ii(trNo)+is_water);
            trials.water_per_trial_y(trNo)=arena.ysync(trials.odor_ii(trNo)+is_water);
        end
    end

end



%Display the location of trial start and water delivery
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on

trials.hit1=zeros(1,trials.odor_trNo);
trials.miss1=zeros(1,trials.odor_trNo);
trials.hit4=zeros(1,trials.odor_trNo);
trials.miss4=zeros(1,trials.odor_trNo);
delta_ii_water=[];

for trNo=1:trials.odor_trNo
    if trials.odor_lane(trNo)==1
        plot(trials.x_odor(trNo),trials.y_odor(trNo),'or')
        if trials.water_per_trial(trNo)==1
            plot(trials.water_per_trial_x(trNo),trials.water_per_trial_y(trNo),'xr')
            trials.hit1(trNo)=1;
            trials.odor_trial_type(trNo)=1;
            delta_ii_water=[delta_ii_water trials.water_per_trial_ii(trNo)-trials.odor_ii(trNo)];
            trials.water_ii(trNo)=trials.water_per_trial_ii(trNo);
        else
            trials.miss1(trNo)=1;
            trials.odor_trial_type(trNo)=2;
        end
    else
        plot(trials.x_odor(trNo),trials.y_odor(trNo),'ob')
        if trials.water_per_trial(trNo)==1
            plot(trials.water_per_trial_x(trNo),trials.water_per_trial_y(trNo),'xb')
            trials.hit4(trNo)=1;
            trials.water_ii(trNo)=trials.water_per_trial_ii(trNo);
            trials.odor_trial_type(trNo)=3;
            delta_ii_water=[delta_ii_water trials.water_per_trial_ii(trNo)-trials.odor_ii(trNo)];
        else
            trials.miss4(trNo)=1;
            trials.odor_trial_type(trNo)=4;
        end
    end
end



plot([10 10],[385 435],'-r')
plot([10 10],[25 75],'-b')


xlabel('x')
ylabel('y')
set(gca, 'YDir', 'reverse');
title('Trial start (o) and water delivery (x), red lane 1, blue lane 4')

for trNo=1:trials.odor_trNo
    if trials.odor_lane(trNo)==1
        plot(trials.x_odor(trNo),trials.y_odor(trNo),'or')
        if trials.water_per_trial(trNo)==0
            trials.water_ii(trNo)=trials.odor_ii(trNo)+round(mean(delta_ii_water));
        end
    else
        plot(trials.x_odor(trNo),trials.y_odor(trNo),'ob')
        if trials.water_per_trial(trNo)==0
            trials.water_ii(trNo)=trials.odor_ii(trNo)+round(mean(delta_ii_water));
        end
    end
    if trials.water_ii(trNo)>no_time_points
        trials.water_ii(trNo)=no_time_points;
    end
end


%Bin positions into dt time bins
pos=[];
pos(:,1)=arena.xsync;
pos(:,2)=arena.ysync;



dFF_times=[1:no_time_points]*dt_miniscope;

no_neurons=size(dFF,2)-1;
no_time_bins=round(dFF_times(end)/dt);
time_binned=[1:no_time_bins]*dt-dt/2;
neural_data=zeros(no_time_bins,no_neurons);
pos_binned=zeros(no_time_bins,2);

for ii_time_bin=1:no_time_bins
    time_from=time_binned(ii_time_bin)-dt/2;
    time_to=time_binned(ii_time_bin)+dt/2;
    pos_binned(ii_time_bin,:)=mean(pos((dFF_times>=time_from)&(dFF_times<time_to),:),1);
    for ii_neuron=1:no_neurons
        neural_data(ii_time_bin,ii_neuron)=mean(dFF((dFF_times>=time_from)&(dFF_times<time_to),ii_neuron+1),1);
    end
end

trim_factor=no_time_bins/no_time_points;

for trNo=1:trials.odor_trNo
    trials.odor_ii_start(trNo)=round(trim_factor*trials.odor_ii(trNo));
    trials.odor_ii_end(trNo)=round(trim_factor*trials.water_ii(trNo));
end



percent_correct=100*(sum(trials.hit4)+sum(trials.hit1))/trials.odor_trNo;
percent_correct1=100*sum(trials.hit1)/(sum(trials.hit1)+sum(trials.miss1));
percent_correct4=100*sum(trials.hit4)/(sum(trials.hit4)+sum(trials.miss4));
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct for lane 1 ' num2str(percent_correct1) ' percent correct for lane 4 ' num2str(percent_correct4) '\n\n'])



%Do z scores
mean_neural_data_col=mean(neural_data,1);
mean_neural_data=repmat(mean_neural_data_col,no_time_bins,1);

std_neural_data_col=std(neural_data,1);
std_neural_data=repmat(std_neural_data_col,no_time_bins,1);

neural_data=(neural_data-mean_neural_data)./std_neural_data;

max_these_neural_data=-2000;
min_these_neural_data=20000;
for ii_ROI=handles_choices.showTheseROIs
    max_these_neural_data=max([max_these_neural_data max(neural_data(ii_ROI,:))]);
    min_these_neural_data=min([min_these_neural_data min(neural_data(ii_ROI,:))]);
end

y_shift=2*(max_these_neural_data-min_these_neural_data);
no_these_traces=length(handles_choices.showTheseROIs);



figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on

%Plot trial start and reward
for trialNo=1:length(trials.hit1)

    if (trials.hit1(trialNo)==1)||(trials.hit4(trialNo)==1)
        % This is a hit
        plot([time_binned(trials.odor_ii_start(trialNo)) time_binned(trials.odor_ii_start(trialNo))], [0 (no_these_traces+3)*y_shift],...
            'Color',[86/255 180/255 233/255],'LineWidth',1.5)
        plot([time_binned(trials.odor_ii_end(trialNo)) time_binned(trials.odor_ii_end(trialNo))], [0 (no_these_traces+3)*y_shift],...
            'Color',[86/255 180/255 233/255],'LineWidth',1.5)
    else
        %This is a miss
        plot([time_binned(trials.odor_ii_start(trialNo)) time_binned(trials.odor_ii_start(trialNo))], [0 (no_these_traces+3)*y_shift],...
            'Color',[0/255 158/255 115/255],'LineWidth',1.5)
    end

end


ROInum=0;
for ii_ROI=handles_choices.showTheseROIs
    ROInum=ROInum+1;
    plot(decimate(time_binned,5),decimate(neural_data(:,ii_ROI),5)+y_shift*ROInum,'-k','LineWidth',1)
end

ylim([-y_shift*0.2 (no_these_traces+2)*y_shift])
xlabel('time(sec)')
title(['dFF timecourses  ' no_these_traces ' ROIs'])

%Calculate the crosscorrelations
croscorr_traces=abs(corrcoef(neural_data)); %please note that I am using the absolute value

%Set autocorrelations to zero
for ii=1:size(croscorr_traces,1)
    croscorr_traces(ii,ii)=0;
end
Z = linkage(croscorr_traces,'complete','correlation');
figNo=figNo+1;
[H,T,outperm]=dendrogram(Z,0,'Orientation','left');
set(H,'LineWidth',2)
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.05 .1 .1 .8])

%re-sort the matrix
for ii=1:size(croscorr_traces,1)
    for jj=1:size(croscorr_traces,1)
        perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
    end
end


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .8 .8])
hold on
pcolor(perm_croscorr_traces)
colormap hot

caxis([0    0.6])
title(['Cross correlations for ' dFF_file])

%Plot rainbow
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


prain=[0:0.6/99:0.6];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
%             colormap jet
colormap hot
shading interp
ax=gca;
set(ax,'XTickLabel','')

%Plot the hieararchically sorted traces
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on

max_neural_data=-2000;
min_neural_data=20000;
for ii_ROI=1:no_neurons
    max_neural_data=max([max_neural_data max(neural_data(ii_ROI,:))]);
    min_neural_data=min([min_neural_data min(neural_data(ii_ROI,:))]);
end

y_shift=0.8*(max_neural_data-min_neural_data);

%Plot trial start and reward
for trialNo=1:length(trials.hit1)

    if (trials.hit1(trialNo)==1)||(trials.hit4(trialNo)==1)
        % This is a hit
        plot([time_binned(trials.odor_ii_start(trialNo)) time_binned(trials.odor_ii_start(trialNo))], [0 (no_neurons+3)*y_shift],...
            'Color',[86/255 180/255 233/255],'LineWidth',1.5)
        plot([time_binned(trials.odor_ii_end(trialNo)) time_binned(trials.odor_ii_end(trialNo))], [0 (no_neurons+3)*y_shift],...
            'Color',[86/255 180/255 233/255],'LineWidth',1.5)
    else
        %This is a miss
        plot([time_binned(trials.odor_ii_start(trialNo)) time_binned(trials.odor_ii_start(trialNo))], [0 (no_neurons+3)*y_shift],...
            'Color',[0/255 158/255 115/255],'LineWidth',1.5)
    end

end


for trNo=1:no_neurons
    % for trNo=1:20
    plot(decimate(time_binned,5),decimate(neural_data(:,outperm(trNo)),5)+y_shift*trNo,'-k','LineWidth',1)
end

% ylim([-y_shift*0.2 (no_neurons+2)*y_shift])
ylim([0 592])
xlim([0 380])
xlabel('time(sec)')
title(['All dFF timecourses hierarchical order '])

%Plot hierarchically sorted traces in pcolor
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on

hier_neural_data=zeros(size(neural_data,1),size(neural_data,2));

for trNo=1:no_neurons
    hier_neural_data(:,trNo)=neural_data(:,outperm(trNo));
end

drg_pcolor(repmat(time_binned',1,no_neurons),repmat([1:no_neurons],length(time_binned),1),hier_neural_data)



colormap hot
shading interp
cmax=prctile(neural_data(:),99);
cmin=prctile(neural_data(:),1);
caxis([cmin cmax]);

hold on

%Plot trial start and reward
for trialNo=1:length(trials.hit1)

    if (trials.hit1(trialNo)==1)||(trials.hit4(trialNo)==1)
        % This is a hit
        plot([time_binned(trials.odor_ii_start(trialNo)) time_binned(trials.odor_ii_start(trialNo))], [0 no_neurons+3],...
            'Color',[86/255 180/255 233/255],'LineWidth',1.5)
        plot([time_binned(trials.odor_ii_end(trialNo)) time_binned(trials.odor_ii_end(trialNo))], [0 no_neurons+3],...
            'Color',[86/255 180/255 233/255],'LineWidth',1.5)
    else
        %This is a miss
        plot([time_binned(trials.odor_ii_start(trialNo)) time_binned(trials.odor_ii_start(trialNo))], [0 no_neurons+3],...
            'Color',[0/255 158/255 115/255],'LineWidth',1.5)
    end

end


xlim([0 380])
xlabel('Time (sec)')
ylabel('ROI number');
ylim([1 no_neurons])
title(['All dFF timecourses hierarchical order'])


pffft=1;
