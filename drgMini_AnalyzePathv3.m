function handles_out=drgMini_AnalyzePathv3(handles_choices)
%Does decoding of the navigation path for a mouse undergoing odor plume navigation 
%following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020
%Using trials code from drgMini_DecodeOdorConc

%Code for trials modified to take on acount sycnronization between FLIR and
%miniscope

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
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    % dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';

    handles_choices.save_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/Temp/';


    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';

    %     %Second troubleshooting files
        this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220713_FCM6/';
        dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
        arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';

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
            trials.lane_per_trial(trNo)=1;
        end

        if sum(arena.laneodor4(ii-3:ii+3)==1)>0
            %Note: laneodor4 is 1 only for one time point
            trials.lane_per_trial(trNo)=4;
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
    if trials.lane_per_trial(trNo)==1
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
    if trials.lane_per_trial(trNo)==1
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
end


%Bin positions into dt time bins
pos=[];
pos(:,1)=arena.xsync;
pos(:,2)=arena.ysync;
no_time_points=size(pos,1);


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
all_bins_per_window=bins_before+bins_after+bins_current;
no_X_dFF_neurons=no_neurons*all_bins_per_window;
X_dFF=zeros(no_time_bins,no_X_dFF_neurons);

for ii_t=bins_before+1:no_time_bins-bins_after
    ii_n=0;
    for no_win=1:all_bins_per_window
        ii_this_t=ii_t-bins_before+no_win-1;
        X_dFF(ii_t,ii_n+1:ii_n+no_neurons)=neural_data(ii_this_t,:);
        ii_n=ii_n+no_neurons;
    end
end



handles_out.no_neurons=no_neurons;

XYtest=pos_binned;

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

%Calculate the angle
angles=[];
for trNo=1:trials.odor_trNo

    ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;

    no_points=5;
    for ii=ii_predictedstart:ii_predictedend-no_points
        % Define line
        x = XYtest(ii:ii+no_points-1,1); 
        y = XYtest(ii:ii+no_points-1,2); 


        % Fit a linear polynomial (line) to the data
        p = polyfit(x, y, 1);
        slope = p(1); % The slope of the line

        % Calculate the angle in radians
        angle_rad = atan2(x(end)-x(1),slope*(x(end)-x(1)));
        angles.trial(trNo).angle_rad(ii-ii_predictedstart+1) = angle_rad;
        % Convert to degrees
        angles.trial(trNo).angle_deg(ii-ii_predictedstart+1) = rad2deg(angle_rad);
    end
    angles.trial(trNo).ii_for_angle=ii_predictedstart+floor(no_points/2):ii_predictedend-ceil(no_points/2);
    angles.trial(trNo).per_trial_ii_for_angle=[1:length(angles.trial(trNo).angle_deg)];

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    %Plot angle
    subplot(3,1,1)

    plot([ii_predictedstart+floor(no_points/2):ii_predictedend-ceil(no_points/2)],angles.trial(trNo).angle_deg,'-b')
    xlabel('Data point')
    ylabel('Direction (degrees)')

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            title(['Angles for lane 1 hit, trial number ' num2str(trNo)])
        case 2
            %Lane 1 miss
            title(['Angles for lane 1 miss, trial number ' num2str(trNo)])
        case 3
            %Lane 4 hit
            title(['Angles for lane 4 hit, trial number ' num2str(trNo)])
        case 4
            %Lane 4 miss
            title(['Angles for lane 4 miss, trial number ' num2str(trNo)])
    end

    %Plot x
    subplot(3,1,2)

    x = XYtest(ii_predictedstart:ii_predictedend,1);
    plot([ii_predictedstart:ii_predictedend],x,'-b')
    xlabel('Data point')
    ylabel('x')

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            title(['x for lane 1 hit, trial number ' num2str(trNo)])
        case 2
            %Lane 1 miss
            title(['x for lane 1 miss, trial number ' num2str(trNo)])
        case 3
            %Lane 4 hit
            title(['x for lane 4 hit, trial number ' num2str(trNo)])
        case 4
            %Lane 4 miss
            title(['x for lane 4 miss, trial number ' num2str(trNo)])
    end

    %Plot y
    subplot(3,1,3)

    y = XYtest(ii_predictedstart:ii_predictedend,2);
    plot([ii_predictedstart:ii_predictedend],y,'-b')
    xlabel('Data point')
    ylabel('y')

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            title(['y for lane 1 hit, trial number ' num2str(trNo)])
        case 2
            %Lane 1 miss
            title(['y for lane 1 miss, trial number ' num2str(trNo)])
        case 3
            %Lane 4 hit
            title(['y for lane 4 hit, trial number ' num2str(trNo)])
        case 4
            %Lane 4 miss
            title(['y for lane 4 miss, trial number ' num2str(trNo)])
    end
end


%Now plot per trial trajectories
mean_end_angles=[];
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
            %Lane 4 miss
            plot(XYtest(ii_predictedstart:ii_predictedend,1),XYtest(ii_predictedstart:ii_predictedend,2),'Color',[86/255 180/255 233/255],'LineWidth',3)
            title(['Lane 4 miss, trial number ' num2str(trNo)])
    end

    plot(XYtest(ii_trialstart,1),XYtest(ii_trialstart,2),'ob','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b')
    plot(XYtest(ii_trialend,1),XYtest(ii_trialend,2),'or','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')

    angles.trial(trNo).ii_for_angle
    [min_delta_ii,ii_min]=min(abs(angles.trial(trNo).ii_for_angle-ii_trialend));
    angles.trial(trNo).end_angle=angles.trial(trNo).angle_deg(angles.trial(trNo).per_trial_ii_for_angle(ii_min));

    %Find mean angle for the trajectory within 150 mm of the end
    x = XYtest(ii_predictedstart:ii_predictedend,1);
    y = XYtest(ii_predictedstart:ii_predictedend,2);
    these_points=[ii_predictedstart:ii_predictedend];
    this_ii_end=find(these_points==ii_trialend);
    x_end=x(this_ii_end);
    y_end=y(this_ii_end);

    this_point=ii_trialend-1;
    this_ii=find(these_points==this_point);
    these_angles=[];
    at_end=0;
    while (sqrt((x(this_ii)-x_end)^2 +(y(this_ii)-y_end)^2)<150)&(at_end==0)
        [min_delta_ii,ii_min]=min(abs(angles.trial(trNo).ii_for_angle-this_point));
        these_angles=[these_angles angles.trial(trNo).angle_rad(angles.trial(trNo).per_trial_ii_for_angle(ii_min))];
        this_point=this_point-1;
        this_ii=find(these_points==this_point);
        if this_point<ii_predictedstart
            at_end=1;
        end
    end

    angles.trial(trNo).mean_end_angle=rad2deg(circ_mean(these_angles'));
  
    %Find where the mouse makes turns towards the odor spouts
    %Let's see if we can find the decision point using a shift to an angle
    %between -60 and -120 

    max_angle=-60;
    min_angle=-120;

    at_end=0;
    ii=0;
    ii_turns=0;
    while at_end==0
        ii_next=find((angles.trial(trNo).angle_deg(ii+1:end)<max_angle)&(angles.trial(trNo).angle_deg(ii+1:end)>min_angle),1,'first');
        if ~isempty(ii_next)
            ii_turns=ii_turns+1;
            angles.trial(trNo).ii_turns(ii_turns)=ii+ii_next+ceil(no_points/2); %This is where the animal turns towards the spout
            ii=ii+ii_next;
            ii_next=find((angles.trial(trNo).angle_deg(ii+1:end)>max_angle)|(angles.trial(trNo).angle_deg(ii+1:end)<min_angle),1,'first');
            if ~isempty(ii_next)
                % angles.trial(trNo).ii_turns_delta(ii_turns)=ii_next; %This is how many points the angle remains going toward the spout
                angles.trial(trNo).delta_x(ii_turns)=abs(x(ii+floor(no_points/2))-x(ii+ii_next+floor(no_points/2))); %This is how far along the x axis the mouse moves torwards the spout
                ii=ii+ii_next;
            else
                angles.trial(trNo).delta_x(ii_turns)=abs(x(ii+floor(no_points/2))-x(length(x)-floor(no_points/2)));
                at_end=1;
            end
        else
            at_end=1;
        end
    end

        % ii_turn=find((angles.trial(trNo).angle_deg<max_angle)&(angles.trial(trNo).angle_deg>min_angle),1,'first');
    % ii_turn=find((angles.trial(trNo).angle_deg<-60)&(angles.trial(trNo).angle_deg>-120),1,'first');


 
    % %Let's see if we can find the decision point using a shift to an angle
    % %between -60 and -120 and concomitant decrease in x
    % ii_shift=5; %Number of points for shifted x
    % x_mm_decrease=200;
    % x = XYtest(ii_predictedstart:ii_predictedend,1);
    % x_shifted = zeros(length(x),1);
    % x_shifted(1:end-ii_shift)=x(ii_shift+1:end);
    % x_shifted(end-ii_shift+1:end)=x(end);

    % ii_turn=find((angles.trial(trNo).angle_deg<-60)&(angles.trial(trNo).angle_deg>-120)&(x_shifted-x<-x_mm_decrease),1,'first');


    % ii_angle_shift=round(ii_predictedstart+(no_points/2)-1+ii_turn);
    % if ii_angle_shift>size(XYtest,1)
    %     ii_angle_shift=size(XYtest,1);
    % end

    %Plot the point of turns towards the spout
    for ii_t=1:ii_turns
        if angles.trial(trNo).delta_x(ii_t)>100
            this_ii_turn=angles.trial(trNo).ii_turns(ii_t);
            plot(x(this_ii_turn),y(this_ii_turn),'or','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','k')
        end
    end

    text(100,100,['End angle ' num2str(angles.trial(trNo).end_angle)])
    text(100,130,['Mean end angle ' num2str(angles.trial(trNo).mean_end_angle)])
    xlabel('x')
    ylabel('y')
    set(gca, 'YDir', 'reverse');
    xlim([0 500])
    ylim([0 480])

    yticks([50 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})
    
    pffft1=1;
 
    mean_end_angles=[mean_end_angles angles.trial(trNo).mean_end_angle];
end

%Do a histogram of mean end angle
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

edges=[-180:10:180];
histogram(mean_end_angles((trials.odor_trial_type==1)|(trials.odor_trial_type==3)), edges)
title('Histogram of mean end angles for hits')

pfft=1;