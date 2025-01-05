function handles_out=drgMini_AnalyzePathv2(handles_choices)
%Does decoding of the navigation path for a mouse undergoing odor plume navigation 
%following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020
%Using trials code from drgMini_DecodeOdorConc

close all

if exist('handles_choices')==0
    clear all

    handles_choices.save_results=1;
    handles_choices.is_sphgpu=0; %0 Diego's Mac, 1 sphgpu, 2 Alpine
    is_sphgpu=handles_choices.is_sphgpu;
    %Troubleshooting Fabio's files May 14th
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    % dFF_file='20220729_FCM22_withodor_miniscope_sync_L4_ncorre_ext_nonneg.mat';
    % arena_file='20220729_FCM22withodor_odorarena_L4_sync.mat';

    %
    % %First troubleshooting files
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';

    handles_choices.save_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/Temp/';


    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';

    %     %Second troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';

%sphgpu
%     this_path='/data/SFTP/PreProcessedDR/20220713_FCM6/';
%     dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
%     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';
    % 
    % %Alpine
    % this_path='/scratch/alpine/drestrepo@xsede.org/PreProcessed/20220713_FCM6/';
    % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';
    % 
    % 
    % handles_choices.save_path='/scratch/alpine/drestrepo@xsede.org/PreProcessed/Temp/';
    

    %No odor troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';


    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

    %algorithm
    handles_choices.algo=2;
    %1 fitrnet
    %2 fitrtree

    handles_choices.save_tag='treexy';

    %The user can define what time period to use spikes from (with respect to the output).
    bins_before=16; %How many bins of neural data prior to the output are used for decoding, 4
    bins_current=1; %Whether to use concurrent time bin of neural data, 1
    bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    handles_choices.bins_before=bins_before;
    handles_choices.bins_current=bins_current;
    handles_choices.bins_after=bins_after;


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.1; %Time bins for decoding, this was 0.2
    dt_miniscope=1/30;
    n_shuffle=3; %Note that n_shuffle is changed to a maximum of ii_n_training
        %n_shuffle=1 will yield an infinite loop

    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    handles_choices.n_shuffle=n_shuffle;

    handles_choices.trial_start_offset=-15; %This was -10
    handles_choices.trial_end_offset=15;
else
    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;
%     training_fraction=handles_choices.training_fraction;
    bins_before=handles_choices.bins_before;
    bins_current=handles_choices.bins_current;
    bins_after=handles_choices.bins_after;
    dt=handles_choices.dt;
    dt_miniscope=handles_choices.dt_miniscope;
    n_shuffle=handles_choices.n_shuffle;
    is_sphgpu=handles_choices.is_sphgpu;
end

try
    delete(gcp('nocreate'));
catch
end

setenv('MW_PCT_TRANSPORT_HEARTBEAT_INTERVAL', '700')

switch is_sphgpu
    case 1
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
    case 2
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/drgMiniscope')
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/m new/Chi Squared')
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/drgMaster')
        addpath(genpath('/projects/drestrepo@xsede.org/software/DR_matlab/m new/kakearney-boundedline-pkg-32f2a1f'))
end

if ~exist(handles_choices.save_path(1:end-1),'dir')
    mkdir(handles_choices.save_path(1:end-1))
end
 
fileID = fopen([this_path 'analyze_path_output.txt'],'w');

% 
% fprintf(1,['\nTrained with within trial data\n\n'])
% fprintf(fileID,['\nTrained with within trial data\n\n'])


figNo=0;

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
trials.ii_r_r=[];
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



handles_out.no_neurons=no_neurons;

XYtest=pos_binned_trimmed;

%Now do the neural networks
tic
%Set what part of data should be part of the training/testing/validation sets
%Note that there was a long period of no movement after about 80% of recording, so I did not use this data.

%Train within trials using a leave one out approach

%training_range_template has all the trials
training_range_template=zeros(1,no_time_bins);
for trNo=1:trials.odor_trNo
    y_pred(trNo).data=[];
    x_pred(trNo).data=[];


    x_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    x_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
    if x_predictedend>no_time_bins
        x_predictedend=no_time_bins;
    end
    training_range_template(x_predictedstart:x_predictedend)=1;
end

%Now plot per trial trajectories
for trNo=1:trials.odor_trNo
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    hold on

    ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;

    ii_trialstart=trials.odor_ii_start(trNo);
    ii_trialend=trials.odor_ii_end(trNo);

    %Okabe_Ito colors
    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(XYtest(ii_predictedstart:ii_predictedend,1),XYtest(ii_predictedstart:ii_predictedend,2),'Color',[213/255 94/255 0],'LineWidth',3)
            title(['Lane 1 hit, trial number ' num2str(trNo)])
        case 2
            %Lane 1 miss
            plot(XYtest(ii_predictedstart:ii_predictedend,1),XYtest(ii_predictedstart:ii_predictedend,2),'Color',[230/255 159/255 0],'LineWidth',3)
            title(['Lane 1 miss, trial number ' num2str(trNo)])
        case 3
            %Lane 4 hit
            plot(XYtest(ii_predictedstart:ii_predictedend,1),XYtest(ii_predictedstart:ii_predictedend,2),'Color',[0 114/255 178/255],'LineWidth',3)
            title(['Lane 4 hit, trial number ' num2str(trNo)])
        case 4
            %Lane 4 hit
            plot(XYtest(ii_predictedstart:ii_predictedend,1),XYtest(ii_predictedstart:ii_predictedend,2),'Color',[86/255 180/255 233/255],'LineWidth',3)
            title(['Lane 4 miss, trial number ' num2str(trNo)])
    end

    plot(XYtest(ii_trialstart,1),XYtest(ii_trialstart,2),'ob','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b')
    plot(XYtest(ii_trialend,1),XYtest(ii_trialend,2),'or','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')

    xlabel('x')
    ylabel('y')
    set(gca, 'YDir', 'reverse');
    xlim([0 500])
    ylim([0 480])

    yticks([50 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})
    
    pffft1=1;

end


pfft=1;