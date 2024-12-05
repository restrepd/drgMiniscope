function handles_out=drgMini_DecodeOdorConc(handles_choices)
%Does decoding of the odor concentraiton for a mouse undergoing odor plume navigation 
%following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020


close all

if exist('handles_choices')==0
    clear all

    %Troubleshooting Fabio's files May 14th
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    % dFF_file='20220729_FCM22_withodor_miniscope_sync_L4_ncorre_ext_nonneg.mat';
    % arena_file='20220729_FCM22withodor_odorarena_L4_sync.mat';

    %
    % %First troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    % this_path='/data/SFTP/PreProcessedDR/20220804_FCM22/';
    % dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';


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
    this_path='/data/SFTP/PreProcessedDR/20220804_FCM22/';
    dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';

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

    handles_choices.is_sphgpu=1;
    is_sphgpu=handles_choices.is_sphgpu;

    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

    handles_choices.save_path='/data/SFTP/DecodeOdorConc/';

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


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2; %Time bins for decoding
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
end

if is_sphgpu==1
    addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
    addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
    addpath('/home/restrepd/Documents/MATLAB/drgMaster')
    addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
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


% dFF=readmatrix([this_path dFF_file]); %Timepoints x ROIs

load([this_path arena_file])

%Extract trials
trials=[];
trials.ii_odor=[];
trials.x_odor=[];
trials.y_odor=[];

trials.ii_laneodor1=[];
trials.x_laneodor1=[];
trials.y_laneodor1=[];

trials.ii_laneodor4=[];
trials.x_laneodor4=[];
trials.y_laneodor4=[];

trials.ii_lanewater1=[];
trials.x_lanewater1=[];
trials.y_lanewater1=[];

trials.ii_lanewater4=[];
trials.x_lanewater4=[];
trials.y_lanewater4=[];

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


%Display the location of trial start and water delivery
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on
plot(trials.x_laneodor1,trials.y_laneodor1,'or')
plot(trials.x_laneodor4,trials.y_laneodor4,'ob')
plot([10 10],[385 435],'-r')
plot([10 10],[25 75],'-b')

plot(trials.x_lanewater1,trials.y_lanewater1,'xr')
plot(trials.x_lanewater4,trials.y_lanewater4,'xb')
xlabel('x')
ylabel('y')
set(gca, 'YDir', 'reverse');
title('Trial start (o) and water delivery (x), red lane 1, blue lane 4')


%Bin positions into dt time bins
pos=[];
pos(:,1)=arena.xsync;
pos(:,2)=arena.ysync;
no_time_points=size(pos,1);


dFF_times=[1:no_time_points]*dt_miniscope;

no_neurons=size(dFF,2)-1;
no_time_bins=ceil(dFF_times(end)/dt);
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

%Now calculate the behavioral performance
trials.hit1=0;
trials.miss1=0;
trials.lane1=0;
trials.lane4=0;
trials.hit4=0;
trials.miss4=0;
trials.odor_trNo=0;

trim_factor=no_time_bins/no_time_points;

dii_trial=[];
for trNo=1:length(trials.ii_odor)

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
            trials.hit1_ii_start(trials.hit1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.hit1_ii_end(trials.hit1)=ceil(trim_factor*trials.ii_lanewater1(this_water));
            dii_trial=[dii_trial trials.hit1_ii_end(trials.hit1)-trials.hit1_ii_start(trials.hit1)];
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.hit1_ii_start(trials.hit1);
            trials.odor_ii_end(trials.odor_trNo)=trials.hit1_ii_end(trials.hit1);
            trials.odor_trial_type(trials.odor_trNo)=1;
            trials.odor_lane(trials.odor_trNo)=1;
        else
            trials.miss1=trials.miss1+1;
            trials.miss1_ii_start(trials.miss1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss1_ii_start(trials.miss1);
            trials.odor_trial_type(trials.odor_trNo)=2;
            trials.odor_lane(trials.odor_trNo)=1;
        end
    end

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
            trials.hit4_ii_start(trials.hit4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.hit4_ii_end(trials.hit4)=ceil(trim_factor*trials.ii_lanewater4(this_water));
            dii_trial=[dii_trial trials.hit4_ii_end(trials.hit4)-trials.hit4_ii_start(trials.hit4)];
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.hit4_ii_start(trials.hit4);
            trials.odor_ii_end(trials.odor_trNo)=trials.hit4_ii_end(trials.hit4);
            trials.odor_trial_type(trials.odor_trNo)=3;
            trials.odor_lane(trials.odor_trNo)=4;
        else
            trials.miss4=trials.miss4+1;
            trials.miss4_ii_start(trials.miss4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss4_ii_start(trials.miss4);
            trials.odor_trial_type(trials.odor_trNo)=4;
            trials.odor_lane(trials.odor_trNo)=4;
        end
    end

end

for ii_miss=1:trials.miss1
    trials.miss1_ii_end(ii_miss)=trials.miss1_ii_start(ii_miss)+ceil(mean(dii_trial));
end

for ii_miss=1:trials.miss4
    trials.miss4_ii_end(ii_miss)=trials.miss4_ii_start(ii_miss)+ceil(mean(dii_trial));
end

for ii_odor=1:trials.odor_trNo
    if trials.odor_trial_type(ii_odor)==4
        trials.odor_ii_end(ii_odor)=trials.odor_ii_start(ii_odor)+ceil(mean(dii_trial));
    end
    if trials.odor_trial_type(ii_odor)==2
        trials.odor_ii_end(ii_odor)=trials.odor_ii_start(ii_odor)+ceil(mean(dii_trial));
    end
end

if trials.odor_ii_end(trials.odor_trNo)+handles_choices.trial_end_offset>size(pos_binned,1)
    trials.odor_ii_end(trials.odor_trNo)=size(pos_binned,1)-handles_choices.trial_end_offset;
end
percent_correct=100*(trials.hit4+trials.hit1)/(trials.hit4+trials.hit1+trials.miss4+trials.miss1);
percent_correct1=100*(trials.hit1)/(trials.hit1+trials.miss1);
percent_correct4=100*(trials.hit4)/(trials.hit4+trials.miss4);
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct1 ' num2str(percent_correct1) ' percent correct4 ' num2str(percent_correct4) '\n\n'])


nan_mask=logical(ones(1,no_time_bins));
for ii_neuron=1:no_neurons
    this_nan_mask=ones(1,no_time_bins);
    this_nan_mask(1,:)=~isnan(neural_data(:,ii_neuron));
    nan_mask=nan_mask&this_nan_mask;
end
neural_data_trimmed=neural_data(nan_mask,:);
pos_binned_trimmed=pos_binned(nan_mask,:);
no_time_bins=sum(nan_mask);

%Do z scores
mean_neural_data_trimmed_col=mean(neural_data_trimmed,1);
mean_neural_data_trimmed=repmat(mean_neural_data_trimmed_col,no_time_bins,1);

std_neural_data_trimmed_col=std(neural_data_trimmed,1);
std_neural_data_trimmed=repmat(std_neural_data_trimmed_col,no_time_bins,1);

neural_data_trimmed=(neural_data_trimmed-mean_neural_data_trimmed)./std_neural_data_trimmed;

pffft=1;

% Format for Wiener Filter, Wiener Cascade, XGBoost, and Dense Neural Network
% Put in "flat" format, so each "neuron / time" is a single feature
% i.e. each time point in the before and after window becomes a different
% "neuron"
all_bins_per_window=bins_before+bins_after+bins_current;
no_X_dFF_neurons=no_neurons*all_bins_per_window;
X_dFF=zeros(no_time_bins,no_X_dFF_neurons);

for ii_t=bins_before+1:no_time_bins-bins_after
    ii_n=0;
    for no_win=1:all_bins_per_window
        ii_this_t=ii_t-bins_before+no_win-1;
        X_dFF(ii_t,ii_n+1:ii_n+no_neurons)=neural_data_trimmed(ii_this_t,:);
        ii_n=ii_n+no_neurons;
    end
end

%Now do the neural networks
tic


%Train with trials using a leave one out approach

%training_range_template has all the trials
training_range_template=zeros(1,no_time_bins);
lane_template=zeros(1,no_time_bins);
odor_plume_template=minC*ones(1,no_time_bins);
for trNo=1:trials.odor_trNo
    % y_pred(trNo).data=[];
    op_pred(trNo).data=[];

    op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;

    training_range_template(op_predictedstart:op_predictedend)=1;
    %Now calculate the odor plume concentration
    these_x=[];
    these_y=[];
    these_ca=[];
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

% parfor trNo=1:trials.odor_trNo

for trNo=1:trials.odor_trNo
    start_toc=toc;

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
    this_training_range=logical(training_range_template)&(~logical(this_test_range));

    % ii_valid_range=ceil(valid_range*no_time_bins);
    %     ii_test_range=ceil(test_range*no_time_bins);

    XdFFtrain=X_dFF(this_training_range,:);
    % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
    XdFFtest=X_dFF(logical(this_test_range),:);

    % XYtrain=pos_binned_trimmed(this_training_range,:);

    OPtrain=odor_plume_template(:,this_training_range)';
    % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
    %     XYtest=pos_binned_trimmed(logical(this_test_range),:);

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
            close(fig)
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

    fprintf(1,['Elapsed time ' num2str(toc-start_toc) ' for trial number ' num2str(trNo) ' \n\n'])
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
    % y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

end

%Now do predictions for reversed/permuted training periods
op_predicted_sh=zeros(no_time_bins,n_shuffle);
% y_predicted_sh=zeros(no_time_bins,n_shuffle);


%We will do a reversal and a circular permutation
sh_shift=0;
while sh_shift==0
    sh_shift=floor(rand*n_shuffle);
end

odor_plume_reversed=zeros(size(odor_plume_template,1),size(odor_plume_template,2));
for ii_trl=1:length(odor_plume_template)
    odor_plume_reversed(1,length(odor_plume_template)-ii_trl+1)=odor_plume_template(1,ii_trl);
end

for ii_shuffled=1:n_shuffle

    for trNo=1:trials.odor_trNo
        % y_pred(trNo).data=[];
        op_pred(trNo).data=[];
        MdlY1_pars(trNo).pars=[];
        % MdlY2_pars(trNo).pars=[];
    end

    parfor trNo=1:trials.odor_trNo

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
        this_training_range=logical(training_range_template)&(~logical(this_test_range));

        % ii_valid_range=ceil(valid_range*no_time_bins);
        %     ii_test_range=ceil(test_range*no_time_bins);

        XdFFtrain=X_dFF(this_training_range,:);
        % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
        XdFFtest=X_dFF(logical(this_test_range),:);

        % XYtrain=pos_binned_reversed(this_training_range,:);
        OPtrain=odor_plume_reversed(:,this_training_range)';

        % ob_pred(trNo).MdlY1=MdlY1;
        % y_pred(trNo).MdlY2=MdlY2;

        switch which_training_algorithm
            case 1
                %Decode using neural network
                %Decode using neural network
                bestHyperparameters = op_pred(trNo).MdlY1.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
                activationsCell = cellstr(bestHyperparameters.Activations);
                standardizeCell = cellstr(bestHyperparameters.Standardize);
                layer_sizes=bestHyperparameters.Layer_1_Size;
                if ~isnan(bestHyperparameters.Layer_2_Size)
                    layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
                end
                if ~isnan(bestHyperparameters.Layer_3_Size)
                    layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
                end
                MdlY1 = fitrnet(XdFFtrain,OPtrain,'LayerSizes', layer_sizes, ...
                    'Activations', activationsCell{1}, ...
                    'Lambda', bestHyperparameters.Lambda, ...
                    'Standardize', strcmpi(standardizeCell{1},'true'));

                MdlY1_pars(trNo).pars.activations=activationsCell{1};
                MdlY1_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
                MdlY1_pars(trNo).pars.Standardize=standardizeCell{1};
            case 2
                %Gaussian process regression, inherently Bayesian
                bestHyperparameters = op_pred(trNo).MdlY1.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
                standardizeCell = cellstr(bestHyperparameters.Standardize);

                MdlY1 = fitrgp(XdFFtrain,OPtrain,'Sigma', bestHyperparameters.Sigma,'Standardize', strcmpi(standardizeCell{1},'true'));

                MdlY1_pars(trNo).pars.Sigma=bestHyperparameters.Sigma;
                MdlY1_pars(trNo).pars.Standardize=standardizeCell{1};
            case 3
                %Binary decision tree
                bestHyperparameters = op_pred(trNo).MdlY1.HyperparameterOptimizationResults.XAtMinEstimatedObjective;

                MdlY1 = fitrtree(XdFFtrain,OPtrain,'MinLeafSize', bestHyperparameters.MinLeafSize);

                MdlY1_pars(trNo).pars.MinLeafSize=bestHyperparameters.MinLeafSize;
            case 4
                MdlY1 = fitglm(XdFFtrain,OPtrain);
            case 5
                %SVM
                bestHyperparameters = op_pred(trNo).MdlY1.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
                standardizeCell = cellstr(bestHyperparameters.Standardize);

                MdlY1 = fitrsvm(XdFFtrain,OPtrain,'BoxConstraint', bestHyperparameters.BoxConstraint,...
                    'KernelScale',bestHyperparameters.KernelScale,...
                    'Epsilon',bestHyperparameters.Epsilon,...
                    'Standardize', strcmpi(standardizeCell{1},'true'));

                MdlY1_pars(trNo).pars.BoxConstraint=bestHyperparameters.BoxConstraint;
                MdlY1_pars(trNo).pars.KernelScale=bestHyperparameters.KernelScale;
                MdlY1_pars(trNo).pars.Epsilon=bestHyperparameters.Epsilon;
                MdlY1_pars(trNo).pars.Standardize=standardizeCell{1};

        end

        op_pred(trNo).data=predict(MdlY1,XdFFtest);

        % bestHyperparameters = y_pred(trNo).MdlY2.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
        % activationsCell = cellstr(bestHyperparameters.Activations);
        % standardizeCell = cellstr(bestHyperparameters.Standardize);
        % layer_sizes=bestHyperparameters.Layer_1_Size;
        % if ~isnan(bestHyperparameters.Layer_2_Size)
        %     layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
        % end
        % if ~isnan(bestHyperparameters.Layer_3_Size)
        %     layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
        % end
        % MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2),'LayerSizes', layer_sizes, ...
        %                     'Activations', activationsCell{1}, ...
        %                     'Lambda', bestHyperparameters.Lambda, ...
        %                     'Standardize', strcmpi(standardizeCell{1},'true'));
        %
        % MdlY2_pars(trNo).pars.activations=activationsCell{1};
        % MdlY2_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
        % MdlY2_pars(trNo).pars.Standardize=standardizeCell{1};
        %
        %
        % y_pred(trNo).data=predict(MdlY2,XdFFtest);
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

        op_predicted_sh(logical(this_test_range),ii_shuffled)=op_pred(trNo).data;
        % y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(trNo).data;
    end

end


fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

op_predictedstart=1;
op_predictedend=length(op_predicted(:,1));

XYtest=pos_binned_trimmed;

handles_out.op_predicted_sh=op_predicted_sh;
handles_out.op_predicted=op_predicted;
handles_out.odor_plume_template=odor_plume_template;
handles_out.trials=trials;
handles_out.MdlY1_pars=MdlY1_pars;


no_conv_points=11;
% conv_win=ones(1,no_conv_points)/no_conv_points;
conv_win_gauss = gausswin(no_conv_points);
conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

op_predicted_conv=conv(op_predicted,conv_win_gauss,'same');

%Now limit the x and y to max and min
minop=min(odor_plume_template);
op_predicted_conv(op_predicted_conv<minop)=minop;
maxop=max(odor_plume_template);
op_predicted_conv(op_predicted_conv>maxop)=maxop;


op_predicted_sh_conv=zeros(size(op_predicted_sh,1),size(op_predicted_sh,2));

for ii_sh=1:n_shuffle
    this_op_predicted_sh=zeros(length(op_predicted),1);
    this_op_predicted_sh(:,1)=op_predicted_sh(:,ii_sh);

    this_op_predicted_sh_conv=[];

    this_op_predicted_sh_conv=conv(this_op_predicted_sh,conv_win_gauss,'same');


    %Now limit the x and y to max and min
    minop=min(odor_plume_template);
    this_op_predicted_sh_conv(this_op_predicted_sh_conv<minop)=minop;
    maxop=max(odor_plume_template);
    this_op_predicted_sh_conv(this_op_predicted_sh_conv>maxop)=maxop;

    op_predicted_sh_conv(:,ii_sh)=this_op_predicted_sh_conv;
end

figHandles=get(0,'children');
figureNumbers=[figHandles.Number];
if max(figureNumbers)>figNo
    for fNo=figNo+1:max(figureNumbers)
        if sum(figureNumbers==fNo)>0
            close(fNo)
        end
    end
end

% fprintf(1, 'R2 nn conv x, y entire run: %d %d\n',drgGetR2(XYtest(:,1),op_predicted_conv),drgGetR2(XYtest(:,2),y_predicted_conv));
R1=corrcoef(odor_plume_template,op_predicted_conv);
fprintf(1, 'Correlation coefficient nn conv odor concentration entire run %d\n\n',R1(1,2));

fractionExplained_op = drgMini_calculateVarianceExplained(odor_plume_template', op_predicted_conv);
fprintf(1, 'Fraction of variance  predicted nn conv odor concentration entire run %d\n\n',fractionExplained_op);
fprintf(1, '\n\n');

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on

op_predictedstart=1;
op_predictedend=length(op_predicted(:,1));
plot(odor_plume_template(op_predictedstart:op_predictedend),'Color',[0.6 0.6 1],'LineWidth',3)
plot(op_predicted_conv(op_predictedstart:op_predictedend),'-k','LineWidth',1.5)

maxop=max(odor_plume_template(op_predictedstart:op_predictedend));
minop=min(odor_plume_template(op_predictedstart:op_predictedend));

for trNo=1:trials.odor_trNo
    this_op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    this_op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
    plot([this_op_predictedstart:this_op_predictedend],(maxC+(maxC-minC)*0.1)*ones(this_op_predictedend-this_op_predictedstart+1,1),'-k','LineWidth',2)
end



title('Odor concentration, b:original, k:predicted (trained per trial)')


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(odor_plume_template(op_predictedstart:op_predictedend),op_predicted_conv(op_predictedstart:op_predictedend),'.b')

xlabel('Actual odor concentraiton')
ylabel('Predicted')


title('Odor concentration for nn (trained per trial)')


%Keep track of the per trial decoding
op_all_trials=[];
op_decod_all_trials=[];
op_all_trials_sh=[];
op_all_hits=[];
op_decod_all_hits=[];
op_all_hits_sh=[];
op_decod_all_hits_sh=[];
op_all_miss=[];
op_decod_all_miss=[];
op_all_miss_sh=[];
op_decod_all_miss_sh=[];
op_decod_all_trials_sh=[];
op_between_trials=[];
op_decod_between_trials=[];
op_between_trials_sh=[];
op_decod_between_trials_sh=[];


%Plot the per trial results for x and re-classify the wall hits where the
%mouse chooses to go along a wall with low odor and then move along the x=0
%cm wall towards the other side
spout4_xy=[0 70];
spout1_xy=[0 450];

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;

last_op_predictedend=1;
wall_threshold=50;

for trNo=1:trials.odor_trNo

    op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;

    op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
    op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];

    op_between_trials=[op_between_trials; odor_plume_template(last_op_predictedend:op_predictedstart)'];
    op_decod_between_trials=[op_decod_between_trials; op_predicted_conv(last_op_predictedend:op_predictedstart)];
    last_op_predictedend=op_predictedend;

    ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;

   

    %Okabe_Ito colors
    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits vermillion
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[213/255 94/255 0],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];
        case 2
            %Lane 1 miss orange
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[230/255 159/255 0],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];

        case 3
            %Lane 4 hit blue
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[0 114/255 178/255],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];
        case 4
            %Lane 4 miss sky blue
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[86/255 180/255 233/255],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
    end
    ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
end

ylim([min(op_all_trials)-0.1*(max(op_all_trials)-min(op_all_trials)) max(op_all_trials)+0.1*(max(op_all_trials)-min(op_all_trials))])

title('Odor concentration per trial, verm:hit1, or:mis1, b:hit4, bsky:miss4 k:predicted')


%Plot the per trial results for x with nn trained with permuted input
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;
last_op_predictedend=1;
for trNo=1:trials.odor_trNo

    op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;

    op_all_trials_sh=[op_all_trials_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
    op_decod_all_trials_sh=[op_decod_all_trials_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];

    op_between_trials_sh=[op_between_trials_sh; odor_plume_template(last_op_predictedend:op_predictedstart)'];
    op_decod_between_trials_sh=[op_decod_between_trials_sh; op_predicted_sh_conv(last_op_predictedend:op_predictedstart,1)];
    last_op_predictedend=op_predictedend;

    ii_end=ii_start+length(XYtest(op_predictedstart:op_predictedend,2))-1;

    %Plot accuracy per trial for permuted training control
    CIsp = bootci(1000, @mean, op_predicted_sh_conv(op_predictedstart:op_predictedend,:)');
    meansp=mean(op_predicted_sh_conv(op_predictedstart:op_predictedend,:)',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;

    [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(op_predicted_sh_conv(op_predictedstart:op_predictedend,:)',1)', CIsp', '-k');


    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[213/255 94/255 0],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_hits_sh=[op_all_hits_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_hits_sh=[op_decod_all_hits_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[230/255 159/255 0],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_miss_sh=[op_all_miss_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_miss_sh=[op_decod_all_miss_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[0 114/255 178/255],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_hits_sh=[op_all_hits_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_hits_sh=[op_decod_all_hits_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
        case 4
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[86/255 180/255 233/255],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_miss_sh=[op_all_miss_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_miss_sh=[op_decod_all_miss_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
    end
    ii_start=ii_start+length(XYtest(op_predictedstart:op_predictedend,2))+20;
end

title('Odor concentration per trial permuted, verm:hit1, or:mis1, b:hit4, bsky:miss4 k:predicted')


R1=corrcoef(op_all_trials,op_decod_all_trials);
fprintf(1, 'Correlation coefficient nn conv odor conc within trials %d\n\n',R1(1,2));
handles_out.R1.all_trials=R1(1,2);

R1=corrcoef(op_all_hits,op_decod_all_hits);
fprintf(1, 'Correlation coefficient nn conv odor conc hits %d\n\n',R1(1,2));
handles_out.R1.all_hits=R1(1,2);

R1=corrcoef(op_all_miss,op_decod_all_miss);
fprintf(1, 'Correlation coefficient nn conv odor conc miss %d\n\n',R1(1,2));
handles_out.R1.all_miss=R1(1,2);

R1=corrcoef(op_between_trials,op_decod_between_trials);
fprintf(1, 'Correlation coefficient nn conv odor conc between trials %d\n\n',R1(1,2));
handles_out.R1.between_trials=R1(1,2);

R1=corrcoef(op_all_trials_sh,op_decod_all_trials_sh);
fprintf(1, 'Correlation coefficient nn permuted odor conc per trial %d\n\n',R1(1,2));
handles_out.R1.all_trials_sh=R1(1,2);

R1=corrcoef(op_all_hits_sh,op_decod_all_hits_sh);
fprintf(1, 'Correlation coefficient nn permuted odor conc hits %d\n\n',R1(1,2));
handles_out.R1.hits_sh=R1(1,2);

R1=corrcoef(op_all_miss_sh,op_decod_all_miss_sh);
fprintf(1, 'Correlation coefficient nn permuted odor conc miss %d\n\n',R1(1,2));
handles_out.R1.miss_sh=R1(1,2);

R1=corrcoef(op_between_trials_sh,op_decod_between_trials_sh);
fprintf(1, 'Correlation coefficient nn permuted odor conc between trial %d\n\n',R1(1,2));
handles_out.R1.between_trials_sh=R1(1,2);


%Fraction of variance  explaind
fractionExplained_op = drgMini_calculateVarianceExplained(op_all_trials, op_decod_all_trials);
fprintf(1, 'Fraction of variance  predicted nn conv odor conc within trials %d\n\n',fractionExplained_op);
handles_out.pVar.all_trials=fractionExplained_op;

fractionExplained_op = drgMini_calculateVarianceExplained(op_all_hits, op_decod_all_hits);
fprintf(1, 'Fraction of variance  predicted nn conv odor conc hits %d\n\n',fractionExplained_op);
handles_out.pVar.all_hits=fractionExplained_op;

fractionExplained_op = drgMini_calculateVarianceExplained(op_all_miss, op_decod_all_miss);
fprintf(1, 'Fraction of variance  predicted nn conv odor conc misss %d\n\n',fractionExplained_op);
handles_out.pVar.all_miss=fractionExplained_op;

fractionExplained_op = drgMini_calculateVarianceExplained(op_between_trials, op_decod_between_trials);
fprintf(1, '\n\nFraction of variance  predicted nn conv odor conc between trials %d\n\n',fractionExplained_op);
handles_out.pVar.between_trials=fractionExplained_op;

fractionExplained_op = drgMini_calculateVarianceExplained(op_all_trials_sh, op_decod_all_trials_sh);
fprintf(1, '\n\nFraction of variance  predicted nn permuted odor conc within trials %d\n\n',fractionExplained_op);
handles_out.pVar.all_trials_sh=fractionExplained_op;

fractionExplained_op = drgMini_calculateVarianceExplained(op_all_hits_sh, op_decod_all_hits_sh);
fprintf(1, '\n\nFraction of variance  predicted nn permuted odor conc hit %d\n\n',fractionExplained_op);
handles_out.pVar.all_hits_sh=fractionExplained_op;

fractionExplained_op = drgMini_calculateVarianceExplained(op_all_miss_sh, op_decod_all_miss_sh);
fprintf(1, '\n\nFraction of variance  predicted nn permuted odor conc miss %d\n\n',fractionExplained_op);
handles_out.pVar.all_miss_sh=fractionExplained_op;

fractionExplained_op = drgMini_calculateVarianceExplained(op_between_trials_sh, op_decod_between_trials_sh);
fprintf(1, '\n\nFraction of variance  predicted nn permuted odor conc between trials %d\n\n',fractionExplained_op);
handles_out.pVar.between_trials_sh=fractionExplained_op;



%Plot odor conc vs decoded
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

plot(op_all_trials,op_decod_all_trials,'.b')
xlabel('Odor concentration')
ylabel('Predicted')

title('odor concentration per trial (trained per trial)')

minC=min(odor_plume_template);
maxC=max(odor_plume_template);

figure(4)
shading interp
caxis([minC maxC]);

figure(3)
shading interp
caxis([minC maxC]);

%Plot pseudocolor odor concentrations for each lane 1 trial
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on

this_cmap=fire;

for trNo=1:trials.odor_trNo
    if (trials.odor_trial_type(trNo)==1)||(trials.odor_trial_type(trNo)==2)
        op_predictedstart=trials.odor_ii_start(trNo);
        op_predictedend=trials.odor_ii_end(trNo);


        this_x=pos_binned_trimmed(op_predictedstart:op_predictedend,1);
        this_y=pos_binned_trimmed(op_predictedstart:op_predictedend,2);
        this_op=odor_plume_template(op_predictedstart:op_predictedend)';

        for ii_point=1:length(this_op)-1
            this_color_ii=round(256*(((this_op(ii_point)+this_op(ii_point+1))/2)-minC)/(maxC-minC));
            if this_color_ii>256
                this_color_ii=256;
            end
            if this_color_ii<1
                this_color_ii=1;
            end
            this_color=zeros(1,3);
            this_color(1,:)=this_cmap(this_color_ii,:);
            plot([this_x(ii_point) this_x(ii_point+1)],[this_y(ii_point) this_y(ii_point+1)],'-','LineWidth', 2,'Color',this_color)
        end
    end
end


shading interp
% caxis([minC maxC]);

xlabel('x')
ylabel('y')
set(gca, 'YDir', 'reverse');
title('Odor concentration for each lane 1 trial')




%Plot pseudocolor odor concentration predicted for each lane 1 trial
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on

this_cmap=fire;

for trNo=1:trials.odor_trNo
    if (trials.odor_trial_type(trNo)==1)||(trials.odor_trial_type(trNo)==2)
        op_predictedstart=trials.odor_ii_start(trNo);
        op_predictedend=trials.odor_ii_end(trNo);


        this_x=pos_binned_trimmed(op_predictedstart:op_predictedend,1);
        this_y=pos_binned_trimmed(op_predictedstart:op_predictedend,2);
        this_op=op_predicted_conv(op_predictedstart:op_predictedend)';

        for ii_point=1:length(this_op)-1
            this_color_ii=round(256*(((this_op(ii_point)+this_op(ii_point+1))/2)-minC)/(maxC-minC));
            if this_color_ii>256
                this_color_ii=256;
            end
            if this_color_ii<1
                this_color_ii=1;
            end
            this_color=zeros(1,3);
            this_color(1,:)=this_cmap(this_color_ii,:);
            plot([this_x(ii_point) this_x(ii_point+1)],[this_y(ii_point) this_y(ii_point+1)],'-','LineWidth', 2,'Color',this_color)
        end
    end
end

shading interp
% caxis([minC maxC]);

xlabel('x')
ylabel('y')
set(gca, 'YDir', 'reverse');
title('Predicted odor concentration for each lane 1 trial')


%Plot pseudocolor odor concentrations for each lane 4 trial
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on

this_cmap=fire;

for trNo=1:trials.odor_trNo
    if (trials.odor_trial_type(trNo)==3)||(trials.odor_trial_type(trNo)==4)
        op_predictedstart=trials.odor_ii_start(trNo);
        op_predictedend=trials.odor_ii_end(trNo);


        this_x=pos_binned_trimmed(op_predictedstart:op_predictedend,1);
        this_y=pos_binned_trimmed(op_predictedstart:op_predictedend,2);
        this_op=odor_plume_template(op_predictedstart:op_predictedend)';

        for ii_point=1:length(this_op)-1
            this_color_ii=round(256*(((this_op(ii_point)+this_op(ii_point+1))/2)-minC)/(maxC-minC));
            if this_color_ii>256
                this_color_ii=256;
            end
            if this_color_ii<1
                this_color_ii=1;
            end
            this_color=zeros(1,3);
            this_color(1,:)=this_cmap(this_color_ii,:);
            plot([this_x(ii_point) this_x(ii_point+1)],[this_y(ii_point) this_y(ii_point+1)],'-','LineWidth', 2,'Color',this_color)
        end
    end
end

shading interp
caxis([minC maxC]);

xlim([0 500])
ylim([0 480])

xlabel('x')
ylabel('y')
set(gca, 'YDir', 'reverse');
title('Odor concentration for each lane 4 trial')




%Plot pseudocolor odor concentration predicted for each lane 4 trial
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on

this_cmap=fire;

for trNo=1:trials.odor_trNo
    if (trials.odor_trial_type(trNo)==3)||(trials.odor_trial_type(trNo)==4)
        op_predictedstart=trials.odor_ii_start(trNo);
        op_predictedend=trials.odor_ii_end(trNo);


        this_x=pos_binned_trimmed(op_predictedstart:op_predictedend,1);
        this_y=pos_binned_trimmed(op_predictedstart:op_predictedend,2);
        this_op=op_predicted_conv(op_predictedstart:op_predictedend)';

        for ii_point=1:length(this_op)-1
            this_color_ii=round(256*(((this_op(ii_point)+this_op(ii_point+1))/2)-minC)/(maxC-minC));
            if this_color_ii>256
                this_color_ii=256;
            end
            if this_color_ii<1
                this_color_ii=1;
            end
            this_color=zeros(1,3);
            this_color(1,:)=this_cmap(this_color_ii,:);
            plot([this_x(ii_point) this_x(ii_point+1)],[this_y(ii_point) this_y(ii_point+1)],'-','LineWidth', 2,'Color',this_color)
        end
    end
end

shading interp
caxis([minC maxC]);

xlim([0 500])
ylim([0 480])

xlabel('x')
ylabel('y')
set(gca, 'YDir', 'reverse');
title('Predicted odor concentration for each lane 4 trial')


%Plot pseudocolor odor concentrations for each trial one at a time
figNo=figNo+1;


this_cmap=fire;

for trNo=1:trials.odor_trNo
    op_predictedstart=trials.odor_ii_start(trNo);
    op_predictedend=trials.odor_ii_end(trNo);


    this_x=pos_binned_trimmed(op_predictedstart:op_predictedend,1);
    this_y=pos_binned_trimmed(op_predictedstart:op_predictedend,2);
    this_op=odor_plume_template(op_predictedstart:op_predictedend)';

    handles_out.per_trial.trial(trNo).this_x=this_x;
    handles_out.per_trial.trial(trNo).this_y=this_y;
    handles_out.per_trial.trial(trNo).this_op=this_op;

    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.7 .2 .3 .3])

    xlabel('x')
    ylabel('y')
    set(gca, 'YDir', 'reverse');

    shading interp
    caxis([minC maxC]);

    xlim([0 500])
    ylim([0 480])

    switch trials.odor_trial_type(trNo)
        case 1
            title(['Odor concentration for trial ' num2str(trNo) ' hit lane 1'])
        case 2
            title(['Odor concentration for trial ' num2str(trNo) ' miss lane 1'])
        case 3
            title(['Odor concentration for trial ' num2str(trNo) ' hit lane 4'])
        case 4
            title(['Odor concentration for trial ' num2str(trNo) ' miss lane 4'])
    end


    hold on

    for ii_point=1:length(this_op)-1
        this_color_ii=round(256*(((this_op(ii_point)+this_op(ii_point+1))/2)-minC)/(maxC-minC));
        if this_color_ii>256
            this_color_ii=256;
        end
        if this_color_ii<1
            this_color_ii=1;
        end
        this_color=zeros(1,3);
        this_color(1,:)=this_cmap(this_color_ii,:);
        plot([this_x(ii_point) this_x(ii_point+1)],[this_y(ii_point) this_y(ii_point+1)],'-','LineWidth', 2,'Color',this_color)
    end

    op_predictedstart=trials.odor_ii_start(trNo);
    op_predictedend=trials.odor_ii_end(trNo);


    this_x=pos_binned_trimmed(op_predictedstart:op_predictedend,1);
    this_y=pos_binned_trimmed(op_predictedstart:op_predictedend,2);
    this_op_pred=op_predicted_conv(op_predictedstart:op_predictedend)';

    handles_out.per_trial.trial(trNo).this_op_pred=this_op_pred;


    try
        close(figNo+1)
    catch
    end

    hFig = figure(figNo+1);

    set(hFig, 'units','normalized','position',[.4 .2 .3 .3])

    xlabel('x')
    ylabel('y')
    set(gca, 'YDir', 'reverse');

    shading interp
    caxis([minC maxC]);

    xlim([0 500])
    ylim([0 480])

    switch trials.odor_trial_type(trNo)
        case 1
            title(['Predicted odor concentration for trial ' num2str(trNo) ' hit lane 1'])
        case 2
            title(['Predicted odor concentration for trial ' num2str(trNo) ' miss lane 1'])
        case 3
            title(['Predicted odor concentration for trial ' num2str(trNo) ' hit lane 4'])
        case 4
            title(['Predicted odor concentration for trial ' num2str(trNo) ' miss lane 4'])
    end

    hold on

    for ii_point=1:length(this_op)-1
        this_color_ii=round(256*(((this_op_pred(ii_point)+this_op_pred(ii_point+1))/2)-minC)/(maxC-minC));
        if this_color_ii>256
            this_color_ii=256;
        end
        if this_color_ii<1
            this_color_ii=1;
        end
        this_color=zeros(1,3);
        this_color(1,:)=this_cmap(this_color_ii,:);
        plot([this_x(ii_point) this_x(ii_point+1)],[this_y(ii_point) this_y(ii_point+1)],'-','LineWidth', 2,'Color',this_color)
    end
 
    pffft=1;
end

save([handles_choices.save_path arena_file(1:end-4) handles_choices.save_tag '.mat'],'handles_out','handles_choices','-v7.3')

pffft=1;
