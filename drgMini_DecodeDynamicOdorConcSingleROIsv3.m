function handles_out=drgMini_DecodeDynamicOdorConcSingleROIsv3(handles_choices)
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
    arena_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm.mat';
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
    bins_before=[1 0 0 0 0 0 0 0 0 0]; %Which bins before the current bin are used?
    %[1 0 0] means that only one bin before the current bin is used
    %[0 0 1] means the third bin before is used
    bins_current=1; %Whether to use concurrent time bin of neural data, 1
    % bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    handles_choices.bins_before=bins_before;
    handles_choices.bins_current=bins_current;
    % handles_choices.bins_after=bins_after;

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
    handles_choices.save_tag='odorconctr1';
    % handles_choices.lowest_conc=-9.6332;
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

    handles_choices.train_with_hits=0; %0=train with all trials, 1=train with hit trials, 2=train with miss trials

    handles_choices.no_sub_ROIs=[1 5 10 20 500];
    handles_choices.no_runs__per_sub_ROI=[1 1 1 1 1];

    % handles_outic.ii_ROI((handles_outic.fileNo==20)&(handles_outic.ssi_all_op_info>=3))


else

    is_sphgpu=handles_choices.is_sphgpu;
    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;
    % training_fraction=handles_choices.training_fraction;
    bins_before=handles_choices.bins_before;
    bins_current=handles_choices.bins_current;
    % bins_after=handles_choices.bins_after;
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

minC=min([min(mean_plume_l4(:)) min(mean_plume_l1(:))]);
% minC=prctile([mean_plume_l4(:); mean_plume_l1(:)],0.5);
% minC=-9.6332;
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

% %Note that I am replacing minC with prctile because there are a small
% %number of points with very low min
% minC=prctile([mean14_plume_l4(:); mean14_plume_l1(:)],0.5);
% maxC=min([max(mean14_plume_l4(:)) max(mean14_plume_l1(:))]);
% if maxC==minC
%     maxC=minC+0.1;
% end

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

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on
    histogram ([mean_plume_l1(:) mean_plume_l4(:)])

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

pffft=1;

% Format for Wiener Filter, Wiener Cascade, XGBoost, and Dense Neural Network
% Put in "flat" format, so each "neuron / time" is a single feature
% i.e. each time point in the before and after window becomes a different
% "neuron"
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


no_X_dFF_neurons=no_neurons*(sum(bins_before)+1);
X_dFF=zeros(no_time_bins,no_X_dFF_neurons);

for ii_t=1:no_time_bins
    %Enter current bin
    ii_n=0;
    X_dFF(ii_t,ii_n+1:ii_n+no_neurons)=neural_data(ii_t,:);
    ii_n=ii_n+no_neurons;

    for no_before=1:length(bins_before)
        if bins_before(no_before)==1
            ii_this_t=ii_t-no_before;
            if ii_this_t<=0
                ii_this_t=1;
            end
            X_dFF(ii_t,ii_n+1:ii_n+no_neurons)=neural_data(ii_this_t,:);
            ii_n=ii_n+no_neurons;
        end
    end
end


handles_out.no_neurons=no_neurons;

%Now do the neural networks
tic


%Train with trials using a leave one out approach

%training_range_template has all the trials
training_range_template=zeros(1,no_time_bins);
hit_and_or_miss=zeros(1,no_time_bins);
lane_template=zeros(1,no_time_bins);
odor_plume_template=minC*ones(1,no_time_bins);
for trNo=1:trials.odor_trNo
    % y_pred(trNo).data=[];
    op_pred(trNo).data=[];

    op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
    if op_predictedend>size(X_dFF,1)
        op_predictedend=size(X_dFF,1);
    end
    trials.odor_encounter_ii(trNo)=op_predictedend;
    if op_predictedstart<=0
        op_predictedstart=1;
    end

    training_range_template(op_predictedstart:op_predictedend)=1;

    %trials.odor_trial_type(trNo)
    %1 Hit Lane 1
    %2 Miss Lane 1
    %3 Hit Lane 4
    %4 Miss Lane 1
    switch handles_choices.train_with_hits
        case 0
            %Train with all trials
            hit_and_or_miss(op_predictedstart:op_predictedend)=1;
        case 1
            %Train with hit trials
            if (trials.odor_trial_type(trNo)==1)||(trials.odor_trial_type(trNo)==3)
                hit_and_or_miss(op_predictedstart:op_predictedend)=1;
            end
        case 2
            %Train with miss trials
            if (trials.odor_trial_type(trNo)==2)||(trials.odor_trial_type(trNo)==4)
                hit_and_or_miss(op_predictedstart:op_predictedend)=1;
            end
    end
    %Now calculate the odor plume concentration
    these_x=[];
    these_y=[];
    these_ca=[];
    found_odor=0;
    ii=0;
    if trials.odor_lane(trNo)==1
        %Lane 1
        lane_template(op_predictedstart:op_predictedend)=1;
        for x_ii=op_predictedstart:op_predictedend
            if x_ii<trials.odor_ii_start(trNo)
                odor_plume_template(x_ii)=minC;
            else
                if x_ii>trials.odor_ii_end(trNo)
                    odor_plume_template(x_ii)=minC;
                else
                    this_x=pos_binned(x_ii,1);
                    this_y=pos_binned(x_ii,2);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_l1(ii_minay,ii_minax);
                    %If the odor plume has not reached the animal then
                    %this_ca=minC
                    x_on=(x_ii-trials.odor_ii_start(trNo))*dt*air_flow_speed;
                    if (this_x<=x_on)&(found_odor==0)
                        found_odor=1;
                        trials.odor_encounter_ii(trNo)=x_ii;
                    end
                    if this_x>x_on
                        this_ca=minC;
                    end
                    if this_ca<handles_choices.lowest_conc
                        this_ca=handles_choices.lowest_conc;
                    end

                    odor_plume_template(x_ii)=this_ca;
                    ii=ii+1;
                    these_x(ii)=this_x;
                    these_y(ii)=this_y;
                    these_ca(ii)=this_ca;
                end
            end
        end
        pffft=1;
    else
        %Lane 4
        lane_template(op_predictedstart:op_predictedend)=4;
        for x_ii=op_predictedstart:op_predictedend

            if x_ii<trials.odor_ii_start(trNo)
                odor_plume_template(x_ii)=minC;
            else
                if x_ii>trials.odor_ii_end(trNo)
                    odor_plume_template(x_ii)=minC;
                else
                    this_x=pos_binned(x_ii,1);
                    this_y=pos_binned(x_ii,2);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_l4(ii_minay,ii_minax);
                    %If the odor plume has not reached the animal then
                    %this_ca=minC
                    x_on=(x_ii-trials.odor_ii_start(trNo))*dt*air_flow_speed;
                    if (this_x<=x_on)&(found_odor==0)
                        found_odor=1;
                        trials.odor_encounter_ii(trNo)=x_ii;
                    end
                    if this_x>x_on
                        this_ca=minC;
                    end
                    if this_ca<handles_choices.lowest_conc
                        this_ca=handles_choices.lowest_conc;
                    end
                    odor_plume_template(x_ii)=this_ca;
                    ii=ii+1;
                    these_x(ii)=this_x;
                    these_y(ii)=this_y;
                    these_ca(ii)=this_ca;
                end
            end
        end
        pffft=1;
    end
end


%Perform decoding for each subset of ROIs
% handles_choices.no_sub_ROIs=[1 5 10 20 500];
%  handles_choices.no_runs__per_sub_ROI=[1 1 1 1 1];
ROI_indices=1:no_neurons;


start_toc=toc;

for ii_sub=1:no_neurons

    X_dFF_trimmed=zeros(size(X_dFF,1),sum(bins_before)+1);
    for ii_h=1:sum(bins_before)+1
        X_dFF_trimmed(:,ii_h)=X_dFF(:,ii_sub+(ii_h-1)*no_neurons);
    end

    for trNo=1:trials.odor_trNo

        this_test_range=zeros(1,no_time_bins);
        if trNo==1
            ii_test_range_start=1;
        else
            ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
        end

        if trNo==trials.odor_trNo
            ii_test_range_end=no_time_bins;
        else
            ii_test_range_end=trials.odor_ii_end(trNo)+15;
        end

        this_test_range(ii_test_range_start:ii_test_range_end)=1;
        this_training_range=logical(training_range_template)&(~logical(this_test_range)&(logical(hit_and_or_miss)));

        % ii_valid_range=ceil(valid_range*no_time_bins);
        %     ii_test_range=ceil(test_range*no_time_bins);

        XdFFtrain=X_dFF_trimmed(this_training_range,:);
        % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
        XdFFtest=X_dFF_trimmed(logical(this_test_range),:);

        % XYtrain=pos_binned(this_training_range,:);

        OPtrain=odor_plume_template(:,this_training_range)';
        % Yvalid=pos_binned(ii_valid_range(1):ii_valid_range(2),:);
        %     XYtest=pos_binned(logical(this_test_range),:);
        switch which_training_algorithm
            case {1,2,3,4,5}
                switch which_training_algorithm
                    case 1
                        %Decode using neural network
                        opts = struct('ShowPlots', false, ...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);

                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        parpool;

                        MdlY1 = fitrnet(XdFFtrain,OPtrain,'OptimizeHyperparameters','auto',...
                            'HyperparameterOptimizationOptions', opts);
                        close all
                    case 2
                        %Gaussian process regression, inherently Bayesian
                        opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);
                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        fig = figure('Visible','off');
                        MdlY1 = fitrgp(XdFFtrain,OPtrain,'KernelFunction','squaredexponential',...
                            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                            opts);
                        close all
                    case 3
                        %Binary decision tree
                        opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);
                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        fig = figure('Visible','off');
                        MdlY1 = fitrtree(XdFFtrain,OPtrain,...
                            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                            opts);
                        close(fig)
                        imp=predictorImportance(MdlY1);

                        handles_out.imp.trial(trNo).imp=imp;
                    case 4
                        %GLM
                        MdlY1 = fitglm(XdFFtrain,OPtrain);
                    case 5
                        %SVM
                        opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);
                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        fig = figure('Visible','off');
                        MdlY1 = fitrsvm(XdFFtrain,OPtrain,...
                            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                            opts);

                end

            case 6
                %fitrtree gpu

                %Use gpu
                XdFFtrain_gpu=gpuArray(XdFFtrain);
                OPtrain_gpu=gpuArray(OPtrain);

                MdlY1_gpu = fitrtree(XdFFtrain_gpu,OPtrain_gpu);
                imp_gpu=predictorImportance(MdlY1_gpu);

                %Gather back imp from the gpu
                imp=gather(imp_gpu);
                MdlY1=gather(MdlY1_gpu);

                handles_out.imp.trial(trNo).imp=imp;

        end
        %
        % try
        %     delete(gcp('nocreate'));
        % catch
        % end
        % parpool;
        % MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2),'OptimizeHyperparameters','auto',...
        %     'HyperparameterOptimizationOptions', opts);
        %
        % %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);
        %
        % %     op_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
        % %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

        op_pred(trNo).data=predict(MdlY1,XdFFtest);
        % y_pred(trNo).data=predict(MdlY2,XdFFtest);

        op_pred(trNo).MdlY1=MdlY1;
        % y_pred(trNo).MdlY2=MdlY2;

        % fprintf(1,['Elapsed time ' num2str(toc-start_toc) ' for trial number ' num2str(trNo) ' \n\n'])
    end


    %Parse out the parfor loop output
    op_predicted=zeros(no_time_bins,1);
    % y_predicted=zeros(no_time_bins,1);
    for trNo=1:trials.odor_trNo
        this_test_range=zeros(1,no_time_bins);

        if trNo==1
            ii_test_range_start=1;
        else
            ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
        end

        if trNo==trials.odor_trNo
            ii_test_range_end=no_time_bins;
        else
            ii_test_range_end=trials.odor_ii_end(trNo)+15;
        end

        this_test_range(ii_test_range_start:ii_test_range_end)=1;

        op_predicted(logical(this_test_range),1)=op_pred(trNo).data;


    end

    handles_out.ii_sub(ii_sub).op_predicted=op_predicted;

    no_conv_points=11;
    % conv_win=ones(1,no_conv_points)/no_conv_points;
    conv_win_gauss = gausswin(no_conv_points);
    conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

    op_predicted_conv=conv(op_predicted,conv_win_gauss,'same');

    %Now limit the x and y to max and min
    minop=min(odor_plume_template);
    op_predicted_conv(op_predicted_conv<minop)=6;
    maxop=max(odor_plume_template);
    op_predicted_conv(op_predicted_conv>maxop)=maxop;

    handles_out.ii_sub(ii_sub).op_predicted_conv=op_predicted_conv;

    fprintf(1,['Done with ROI number ' num2str(ii_sub) '\n'])
end

fprintf(1,['Elapsed time for all ROIs ' num2str((toc-start_toc)/(60*60)) ' hrs\n\n'])


%
% %Now do predictions for reversed/permuted training periods

start_toc=toc;

for ii_sub=1:no_neurons

    %We will do a reversal and a circular permutation
    odor_plume_reversed=zeros(size(odor_plume_template,1),size(odor_plume_template,2));
    % pos_binned_reversed=zeros(size(pos_binned,1),size(pos_binned,2));
    offset_ii=(ii_sub-1)*floor(size(pos_binned,1)/no_neurons);
    for ii_trl=1:length(odor_plume_template)
        this_ii_trl=ii_trl+offset_ii;
        if this_ii_trl>length(odor_plume_template)
            offset_ii=-ii_trl+1;
            this_ii_trl=ii_trl+offset_ii;
        end
        odor_plume_reversed(1,length(odor_plume_template)-ii_trl+1,:)=odor_plume_template(1,this_ii_trl);
    end


    X_dFF_trimmed=zeros(size(X_dFF,1),sum(bins_before)+1);
    for ii_h=1:sum(bins_before)+1
        X_dFF_trimmed(:,ii_h)=X_dFF(:,ii_sub+(ii_h-1)*no_neurons);
    end

    for trNo=1:trials.odor_trNo


        this_test_range=zeros(1,no_time_bins);
        if trNo==1
            ii_test_range_start=1;
        else
            ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
        end

        if trNo==trials.odor_trNo
            ii_test_range_end=no_time_bins;
        else
            ii_test_range_end=trials.odor_ii_end(trNo)+15;
        end

        this_test_range(ii_test_range_start:ii_test_range_end)=1;
        this_training_range=logical(training_range_template)&(~logical(this_test_range)&(logical(hit_and_or_miss)));

        % ii_valid_range=ceil(valid_range*no_time_bins);
        %     ii_test_range=ceil(test_range*no_time_bins);

        XdFFtrain=X_dFF_trimmed(this_training_range,:);
        % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
        XdFFtest=X_dFF_trimmed(logical(this_test_range),:);

        % XYtrain=pos_binned(this_training_range,:);

        OPtrain=odor_plume_reversed(:,this_training_range)';
        % Yvalid=pos_binned(ii_valid_range(1):ii_valid_range(2),:);
        %     XYtest=pos_binned(logical(this_test_range),:);
        switch which_training_algorithm
            case {1,2,3,4,5}
                switch which_training_algorithm
                    case 1
                        %Decode using neural network
                        opts = struct('ShowPlots', false, ...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);

                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        parpool;

                        MdlY1 = fitrnet(XdFFtrain,OPtrain,'OptimizeHyperparameters','auto',...
                            'HyperparameterOptimizationOptions', opts);
                        close all
                    case 2
                        %Gaussian process regression, inherently Bayesian
                        opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);
                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        fig = figure('Visible','off');
                        MdlY1 = fitrgp(XdFFtrain,OPtrain,'KernelFunction','squaredexponential',...
                            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                            opts);
                        close all
                    case 3
                        %Binary decision tree
                        opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);
                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        fig = figure('Visible','off');
                        MdlY1 = fitrtree(XdFFtrain,OPtrain,...
                            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                            opts);
                        close(fig)
                        imp=predictorImportance(MdlY1);

                        handles_out.imp.trial(trNo).imp=imp;
                    case 4
                        %GLM
                        MdlY1 = fitglm(XdFFtrain,OPtrain);
                    case 5
                        %SVM
                        opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                            'Verbose', 0, ...
                            'MaxObjectiveEvaluations', 15,...
                            'UseParallel',true);
                        try
                            delete(gcp('nocreate'));
                        catch
                        end
                        fig = figure('Visible','off');
                        MdlY1 = fitrsvm(XdFFtrain,OPtrain,...
                            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                            opts);

                end

            case 6
                %fitrtree gpu

                %Use gpu
                XdFFtrain_gpu=gpuArray(XdFFtrain);
                OPtrain_gpu=gpuArray(OPtrain);

                MdlY1_gpu = fitrtree(XdFFtrain_gpu,OPtrain_gpu);
                imp_gpu=predictorImportance(MdlY1_gpu);

                %Gather back imp from the gpu
                imp=gather(imp_gpu);
                MdlY1=gather(MdlY1_gpu);

                handles_out.imp.trial(trNo).imp=imp;

        end
        %
        % try
        %     delete(gcp('nocreate'));
        % catch
        % end
        % parpool;
        % MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2),'OptimizeHyperparameters','auto',...
        %     'HyperparameterOptimizationOptions', opts);
        %
        % %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);
        %
        % %     op_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
        % %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

        op_pred(trNo).data=predict(MdlY1,XdFFtest);
        % y_pred(trNo).data=predict(MdlY2,XdFFtest);

        op_pred(trNo).MdlY1=MdlY1;
        % y_pred(trNo).MdlY2=MdlY2;

        % fprintf(1,['Elapsed time ' num2str(toc-start_toc) ' for trial number ' num2str(trNo) ' \n\n'])
    end


    %Parse out the parfor loop output
    op_predicted=zeros(no_time_bins,1);
    % y_predicted=zeros(no_time_bins,1);
    for trNo=1:trials.odor_trNo
        this_test_range=zeros(1,no_time_bins);

        if trNo==1
            ii_test_range_start=1;
        else
            ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
        end

        if trNo==trials.odor_trNo
            ii_test_range_end=no_time_bins;
        else
            ii_test_range_end=trials.odor_ii_end(trNo)+15;
        end

        this_test_range(ii_test_range_start:ii_test_range_end)=1;

        op_predicted(logical(this_test_range),1)=op_pred(trNo).data;


    end

    handles_out.ii_sub(ii_sub).op_predicted_sh=op_predicted;
    % handles_out.ROIsubset(ROIsubset).ii_sub(ii_sub).op_reversed=odor_plume_reversed;

    no_conv_points=11;
    % conv_win=ones(1,no_conv_points)/no_conv_points;
    conv_win_gauss = gausswin(no_conv_points);
    conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

    op_predicted_conv=conv(op_predicted,conv_win_gauss,'same');

    %Now limit the x and y to max and min
    minop=min(odor_plume_reversed);
    op_predicted_conv(op_predicted_conv<minop)=6;
    maxop=max(odor_plume_reversed);
    op_predicted_conv(op_predicted_conv>maxop)=maxop;

    handles_out.ii_sub(ii_sub).op_predicted_conv_sh=op_predicted_conv;

    fprintf(1,['Done with shuffled ROI number  ' num2str(ii_sub) '\n'])
end

fprintf(1,['Elapsed time for shuffled ' num2str((toc-start_toc)/(60*60)) ' hrs\n\n'])



handles_out.odor_plume_template=odor_plume_template;
handles_out.trials=trials;

save([handles_choices.save_path arena_file(1:end-4) handles_choices.save_tag '_singleROI.mat'],'handles_out','handles_choices','-v7.3')

pffft=1;
