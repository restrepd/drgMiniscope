function handles_out=drgMini_DecodeOdorArenav4(handles_choices)
%Does decoding of the navigation path for a mouse undergoing odor plume navigation 
%following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020


close all

if exist('handles_choices')==0
    clear all

    handles_choices.save_results=1;
    handles_choices.is_sphgpu=2; %1 sphgpu, 2 Alpine
    is_sphgpu=handles_choices.is_sphgpu;
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

%sphgpu
%     this_path='/data/SFTP/PreProcessedDR/20220713_FCM6/';
%     dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
%     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';

    %Alpine
    this_path='/scratch/alpine/drestrepo@xsede.org/PreProcessed/20220713_FCM6/';
    dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';


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
    training_fraction=0.9;
    handles_choices.training_fraction=training_fraction;

    %     training_range=[0, 0.5];
    %     valid_range=[0.5,0.65];
    %     test_range=[0.5, 1];

    %The user can define what time period to use spikes from (with respect to the output).
    bins_before=10; %How many bins of neural data prior to the output are used for decoding, 4
    bins_current=1; %Whether to use concurrent time bin of neural data, 1
    bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    handles_choices.bins_before=bins_before;
    handles_choices.bins_current=bins_current;
    handles_choices.bins_after=bins_after;


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2; %Time bins for decoding
    dt_miniscope=1/30;
    n_shuffle=3; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    handles_choices.n_shuffle=n_shuffle;

    which_training_algorithm=2;
    handles_choices.which_training_algorithm=which_training_algorithm;
    %1=entire session is used for training
    %2=only per trial is used for training
    %3=trained with data between trials

else

    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;
    training_fraction=handles_choices.training_fraction;
    bins_before=handles_choices.bins_before;
    bins_current=handles_choices.bins_current;
    bins_after=handles_choices.bins_after;
    dt=handles_choices.dt;
    dt_miniscope=handles_choices.dt_miniscope;
    n_shuffle=handles_choices.n_shuffle;
    which_training_algorithm=handles_choices.which_training_algorithm;
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
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/Chi Squared')
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/drgMaster')
        addpath(genpath('/projects/drestrepo@xsede.org/software/DR_matlab/m new/kakearney-boundedline-pkg-32f2a1f'))
end

fileID = fopen([this_path 'R1R2output.txt'],'w');

switch which_training_algorithm
    case 1
        fprintf(1,['\nTrained with data for the entire session\n\n'])
    case 2
        fprintf(1,['\nTrained with within trial data\n\n'])
    case 3
        fprintf(1,['\nTrained with between trial data\n\n'])
end

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
%
% %Extract odor on using the camera sync
% at_end=0;
% ii=0;
% jj=0;
% while at_end==0
%     next_ii=find(arena.odorsync(ii+1:end)==1,1,'first');
%     if ~isempty(next_ii)
%         jj=jj+1;
%         trials.ii_odor(jj)=ii+next_ii;
%         trials.x_odor(jj)=arena.xsync(ii+next_ii);
%         trials.y_odor(jj)=arena.ysync(ii+next_ii);
%         ii=ii+next_ii;
%         next_ii=find(arena.odorsync(ii+1:end)==0,1,'first');
%         if ~isempty(next_ii)
%             ii=ii+next_ii;
%         else
%             at_end=1;
%         end
%     else
%          at_end=1;
%     end
% end
%
% %Extract lane1
% at_end=0;
% ii_flir=0;
% jj=0;
% while at_end==0
%     next_ii_flir=find(arena.laneodor1(ii_flir+1:end)==1,1,'first');
%     if ~isempty(next_ii_flir)
%         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
%         jj=jj+1;
%         trials.ii_laneodor1(jj)=next_ii;
%         trials.x_laneodor1(jj)=arena.xsync(next_ii);
%         trials.y_laneodor1(jj)=arena.ysync(next_ii);
%         next_ii_flir2=find(arena.laneodor1(ii_flir+next_ii_flir+1:end)==0,1,'first');
% %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
%         if ~isempty(next_ii_flir2)
%             ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
%         else
%             at_end=1;
%         end
%     else
%          at_end=1;
%     end
% end
%
% %Extract lane4
% at_end=0;
% ii_flir=0;
% jj=0;
% while at_end==0
%     next_ii_flir=find(arena.laneodor4(ii_flir+1:end)==1,1,'first');
%     if ~isempty(next_ii_flir)
%         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
%         jj=jj+1;
%         trials.ii_laneodor4(jj)=next_ii;
%         trials.x_laneodor4(jj)=arena.xsync(next_ii);
%         trials.y_laneodor4(jj)=arena.ysync(next_ii);
%         next_ii_flir2=find(arena.laneodor4(ii_flir+next_ii_flir+1:end)==0,1,'first');
% %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
%         if ~isempty(next_ii_flir2)
%             ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
%         else
%             at_end=1;
%         end
%     else
%          at_end=1;
%     end
% end

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

if handles_choices.save_results==1
    savefig([this_path 'figure1.fig'])
end

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
        else
            trials.miss1=trials.miss1+1;
            trials.miss1_ii_start(trials.miss1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss1_ii_start(trials.miss1);
            trials.odor_trial_type(trials.odor_trNo)=2;
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
        else
            trials.miss4=trials.miss4+1;
            trials.miss4_ii_start(trials.miss4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss4_ii_start(trials.miss4);
            trials.odor_trial_type(trials.odor_trNo)=4;
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

percent_correct=100*(trials.hit4+trials.hit1)/(trials.hit4+trials.hit1+trials.miss4+trials.miss1);
percent_correct1=100*(trials.hit1)/(trials.hit1+trials.miss1);
percent_correct4=100*(trials.hit4)/(trials.hit4+trials.miss4);
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct1 ' num2str(percent_correct1) ' percent correct4 ' num2str(percent_correct4) '\n\n'])
fprintf(fileID,['\nPercent correct ' num2str(percent_correct) ' percent correct1 ' num2str(percent_correct1) ' percent correct4 ' num2str(percent_correct4) '\n\n'])


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
%Set what part of data should be part of the training/testing/validation sets
%Note that there was a long period of no movement after about 80% of recording, so I did not use this data.

switch which_training_algorithm
    case 1
        %Train with entire session
        ii_n_training=floor(1/(1-training_fraction));

        for ii_train=1:ii_n_training
            y_pred(ii_train).data=[];
            x_pred(ii_train).data=[];
        end

        % parfor ii_train=1:ii_n_training
        for ii_train=1:ii_n_training
            start_toc=toc;
            this_test_range=zeros(1,no_time_bins);
            ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
            ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
            if ii_train==ii_n_training
                ii_test_range_end=no_time_bins;
            end
            this_test_range(ii_test_range_start:ii_test_range_end)=1;

            % ii_valid_range=ceil(valid_range*no_time_bins);
            %     ii_test_range=ceil(test_range*no_time_bins);

            XdFFtrain=X_dFF(~logical(this_test_range),:);
            % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
            XdFFtest=X_dFF(logical(this_test_range),:);

            XYtrain=pos_binned_trimmed(~logical(this_test_range),:);
            % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
            %     XYtest=pos_binned_trimmed(logical(this_test_range),:);

            %Decode using neural network
            opts = struct('ShowPlots', false, ...
                'Verbose', 0, ...
                'MaxObjectiveEvaluations', 15,...
                'UseParallel',true);

            MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1),'OptimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions', opts);
            MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2),'OptimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions', opts);

            %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);

            %     x_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
            %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

            x_pred(ii_train).data=predict(MdlY1,XdFFtest);
            y_pred(ii_train).data=predict(MdlY2,XdFFtest);
            fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])
        end
        

        %Parse out the parfor loop output
        x_predicted=zeros(no_time_bins,1);
        y_predicted=zeros(no_time_bins,1);
        for ii_train=1:ii_n_training
            this_test_range=zeros(1,no_time_bins);
            ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
            ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
            if ii_train==ii_n_training
                ii_test_range_end=no_time_bins;
            end
            this_test_range(ii_test_range_start:ii_test_range_end)=1;

            x_predicted(logical(this_test_range),1)=x_pred(ii_train).data;
            y_predicted(logical(this_test_range),1)=y_pred(ii_train).data;

        end

        %Now do predictions for reversed/permuted training periods
        if n_shuffle>ii_n_training
            n_shuffle=ii_n_training;
        end
        handles_choices.n_shuffle=n_shuffle;
        x_predicted_sh=zeros(no_time_bins,n_shuffle);
        y_predicted_sh=zeros(no_time_bins,n_shuffle);


        %We will do a reversal and a circular permutation
        sh_shift=0;
        while sh_shift==0
            sh_shift=floor(rand*n_shuffle);
        end

        pos_binned_reversed=zeros(size(pos_binned_trimmed,1),size(pos_binned_trimmed,2));
        for ii_trl=1:size(pos_binned_trimmed,1)
            pos_binned_reversed(size(pos_binned_trimmed,1)-ii_trl+1,:)=pos_binned_trimmed(ii_trl,:);
        end

        for ii_shuffled=1:n_shuffle

            for ii_train=1:ii_n_training
                y_pred(ii_train).data=[];
                x_pred(ii_train).data=[];
            end

            parfor ii_train=1:ii_n_training

                this_test_range=zeros(1,no_time_bins);
                ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
                ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
                if ii_train==ii_n_training
                    ii_test_range_end=no_time_bins;
                end
                this_test_range(ii_test_range_start:ii_test_range_end)=1;

                ii_train_sh=ii_train+sh_shift;
                if ii_train_sh>ii_n_training
                    ii_train_sh=ii_train_sh-ii_n_training;
                end
                this_test_range_sh=zeros(1,no_time_bins);
                ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train_sh-1));
                ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train_sh);
                if ii_train_sh==ii_n_training
                    ii_test_range_end=no_time_bins;
                end
                this_test_range_sh(ii_test_range_start:ii_test_range_end)=1;

                %Make sure that this_test_range_sh is the same size as this_test_range
                if sum(this_test_range_sh)>sum(this_test_range)
                    first_ii=find(this_test_range_sh==1,1,'first');
                    this_test_range_sh(first_ii:first_ii+sum(this_test_range_sh)-sum(this_test_range)-1)=0;
                end

                if sum(this_test_range_sh)<sum(this_test_range)
                    first_ii=find(this_test_range_sh==1,1,'first');
                    if first_ii+(sum(this_test_range_sh)-sum(this_test_range))>0
                        this_test_range_sh(first_ii+(sum(this_test_range_sh)-sum(this_test_range)):first_ii-1)=1;
                    else
                        last_ii=find(this_test_range_sh==1,1,'last');
                        this_test_range_sh(last_ii+1:last_ii+(sum(this_test_range)-sum(this_test_range_sh)))=1;
                    end
                end

                % ii_valid_range=ceil(valid_range*no_time_bins);
                %     ii_test_range=ceil(test_range*no_time_bins);

                XdFFtrain=X_dFF(~logical(this_test_range),:);

                XdFFtest=X_dFF(logical(this_test_range),:);


                XYtrain=pos_binned_reversed(~logical(this_test_range_sh),:);

                %Decode using neural network
                MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
                MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));

                x_pred(ii_train).data=predict(MdlY1,XdFFtest);
                y_pred(ii_train).data=predict(MdlY2,XdFFtest);
            end

            for ii_train=1:ii_n_training
                this_test_range=zeros(1,no_time_bins);
                ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
                ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
                if ii_train==ii_n_training
                    ii_test_range_end=no_time_bins;
                end
                this_test_range(ii_test_range_start:ii_test_range_end)=1;

                x_predicted_sh(logical(this_test_range),ii_shuffled)=x_pred(ii_train).data;
                y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(ii_train).data;
            end

        end

    case 2
        %Train with trials using a leave one out approach

        %training_range_template has all the trials
        training_range_template=zeros(1,no_time_bins);
        for trNo=1:trials.odor_trNo
            y_pred(trNo).data=[];
            x_pred(trNo).data=[];

            x_predictedstart=trials.odor_ii_start(trNo)-10;
            x_predictedend=trials.odor_ii_end(trNo)+15;
            training_range_template(x_predictedstart:x_predictedend)=1;
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

            XYtrain=pos_binned_trimmed(this_training_range,:);
            % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
            %     XYtest=pos_binned_trimmed(logical(this_test_range),:);

            XdFFtrain_gpu=gpuArray(XdFFtrain);
            XdFFtest_gpu=gpuArray(XdFFtest);
            XYtrain_gpu=gpuArray(XYtrain);
% 
%             %Decode using neural network
%             opts = struct('ShowPlots', false, ...
%                 'Verbose', 0, ...
%                 'MaxObjectiveEvaluations', 15,...
%                 'UseParallel',true);
% 
%             try
%                 delete(gcp('nocreate'));
%             catch
%             end
%             parpool;

            MdlY1_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,1),'Standardize',true...
                );

            %             MdlY1_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,1),'Standardize',true,...
            %                 'OptimizeHyperparameters','auto',...
            %                 'HyperparameterOptimizationOptions', opts);

            MdlY1=gather(MdlY1_gpu);
% 
%             try
%                 delete(gcp('nocreate'));
%             catch
%             end
%             parpool;

            MdlY2_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,2),'Standardize',true...
                );

%             MdlY2_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,2),'Standardize',true,...
%                 'OptimizeHyperparameters','auto',...
%                 'HyperparameterOptimizationOptions', opts);

            MdlY2=gather(MdlY2_gpu);
            %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);

            %     x_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
            %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

            x_pred(trNo).data=predict(MdlY1,XdFFtest);
            y_pred(trNo).data=predict(MdlY2,XdFFtest);

            x_pred(trNo).MdlY1=MdlY1;
            y_pred(trNo).MdlY2=MdlY2;
            
            fprintf(1,['Elapsed time ' num2str(toc-start_toc) ' for trial number ' num2str(trNo) ' \n\n'])
        end
        

        %Parse out the parfor loop output
        x_predicted=zeros(no_time_bins,1);
        y_predicted=zeros(no_time_bins,1);
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

            x_predicted(logical(this_test_range),1)=x_pred(trNo).data;
            y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

        end

        %Now do predictions for reversed/permuted training periods
        x_predicted_sh=zeros(no_time_bins,n_shuffle);
        y_predicted_sh=zeros(no_time_bins,n_shuffle);


        %We will do a reversal and a circular permutation
        sh_shift=0;
        while sh_shift==0
            sh_shift=floor(rand*n_shuffle);
        end

        pos_binned_reversed=zeros(size(pos_binned_trimmed,1),size(pos_binned_trimmed,2));
        for ii_trl=1:size(pos_binned_trimmed,1)
            pos_binned_reversed(size(pos_binned_trimmed,1)-ii_trl+1,:)=pos_binned_trimmed(ii_trl,:);
        end

        for ii_shuffled=1:n_shuffle

            for trNo=1:trials.odor_trNo
                y_pred(trNo).data=[];
                x_pred(trNo).data=[];
                MdlY1_pars(trNo).pars=[];
                MdlY2_pars(trNo).pars=[];
            end

            for trNo=1:trials.odor_trNo
                %                 parfor trNo=1:trials.odor_trNo

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

                XYtrain=pos_binned_reversed(this_training_range,:);


                XdFFtrain_gpu=gpuArray(XdFFtrain);
                XdFFtest_gpu=gpuArray(XdFFtest);
                XYtrain_gpu=gpuArray(XYtrain);

%                 try
%                     delete(gcp('nocreate'));
%                 catch
%                 end
%                 parpool;

                MdlY1_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,1),'Standardize',true...
                    );

                %             MdlY1_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,1),'Standardize',true,...
                %                 'OptimizeHyperparameters','auto',...
                %                 'HyperparameterOptimizationOptions', opts);

                MdlY1=gather(MdlY1_gpu);

                % x_pred(trNo).MdlY1=MdlY1;
                % y_pred(trNo).MdlY2=MdlY2;
                %
                %                 %Decode using neural network
                %                 bestHyperparameters = x_pred(trNo).MdlY1.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
                %                 activationsCell = cellstr(bestHyperparameters.Activations);
                %                 standardizeCell = cellstr(bestHyperparameters.Standardize);
                %                 layer_sizes=bestHyperparameters.Layer_1_Size;
                %                 if ~isnan(bestHyperparameters.Layer_2_Size)
                %                     layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
                %                 end
                %                 if ~isnan(bestHyperparameters.Layer_3_Size)
                %                     layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
                %                 end
                %                 MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1),'LayerSizes', layer_sizes, ...
                %                     'Activations', activationsCell{1}, ...
                %                     'Lambda', bestHyperparameters.Lambda, ...
                %                     'Standardize', strcmpi(standardizeCell{1},'true'));
                %
                %                 MdlY1_pars(trNo).pars.activations=activationsCell{1};
                %                 MdlY1_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
                %                 MdlY1_pars(trNo).pars.Standardize=standardizeCell{1};

                x_pred(trNo).data=predict(MdlY1,XdFFtest);

                %                 bestHyperparameters = y_pred(trNo).MdlY2.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
                %                 activationsCell = cellstr(bestHyperparameters.Activations);
                %                 standardizeCell = cellstr(bestHyperparameters.Standardize);
                %                 layer_sizes=bestHyperparameters.Layer_1_Size;
                %                 if ~isnan(bestHyperparameters.Layer_2_Size)
                %                     layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
                %                 end
                %                 if ~isnan(bestHyperparameters.Layer_3_Size)
                %                     layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
                %                 end
                %                 MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2),'LayerSizes', layer_sizes, ...
                %                                     'Activations', activationsCell{1}, ...
                %                                     'Lambda', bestHyperparameters.Lambda, ...
                %                                     'Standardize', strcmpi(standardizeCell{1},'true'));
                %
                %                 MdlY2_pars(trNo).pars.activations=activationsCell{1};
                %                 MdlY2_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
                %                 MdlY2_pars(trNo).pars.Standardize=standardizeCell{1};

%                 try
%                     delete(gcp('nocreate'));
%                 catch
%                 end
%                 parpool;

                MdlY2_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,2),'Standardize',true...
                    );

                %             MdlY2_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,2),'Standardize',true,...
                %                 'OptimizeHyperparameters','auto',...
                %                 'HyperparameterOptimizationOptions', opts);

                MdlY2=gather(MdlY2_gpu);

                y_pred(trNo).data=predict(MdlY2,XdFFtest);

                fprintf(1,['Elapsed time ' num2str(toc-start_toc) ' for trial number ' num2str(trNo) ' shuffle no ' num2str(ii_shuffled) '\n\n'])
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

                x_predicted_sh(logical(this_test_range),ii_shuffled)=x_pred(trNo).data;
                y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(trNo).data;
            end

        end

    case 3
        %Train between trials using a leave one out approach

        %training_range_template has all the trials
        training_range_template=zeros(1,no_time_bins);
        for trNo=1:trials.odor_trNo
            y_pred(trNo).data=[];
            x_pred(trNo).data=[];

            if trNo==1
                last_x_predictedend=1;
            else
                last_x_predictedend=trials.odor_ii_end(trNo-1)+15;
            end


            this_x_predictedstart=trials.odor_ii_start(trNo)-10;


            training_range_template(last_x_predictedend:this_x_predictedstart)=1;
        end

        last_x_predictedend=trials.odor_ii_end(trials.odor_trNo)+15;
        this_x_predictedstart=no_time_bins;
        training_range_template(last_x_predictedend:this_x_predictedstart)=1;

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

            XYtrain=pos_binned_trimmed(this_training_range,:);
            % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
            %     XYtest=pos_binned_trimmed(logical(this_test_range),:);

            %Decode using neural network
            MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
            MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));

            %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);

            %     x_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
            %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

            x_pred(trNo).data=predict(MdlY1,XdFFtest);
            y_pred(trNo).data=predict(MdlY2,XdFFtest);
        end
        fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

        %Parse out the parfor loop output
        x_predicted=zeros(no_time_bins,1);
        y_predicted=zeros(no_time_bins,1);
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

            x_predicted(logical(this_test_range),1)=x_pred(trNo).data;
            y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

        end

        %Now do predictions for reversed/permuted training periods
        x_predicted_sh=zeros(no_time_bins,n_shuffle);
        y_predicted_sh=zeros(no_time_bins,n_shuffle);


        %We will do a reversal and a circular permutation
        sh_shift=0;
        while sh_shift==0
            sh_shift=floor(rand*n_shuffle);
        end

        pos_binned_reversed=zeros(size(pos_binned_trimmed,1),size(pos_binned_trimmed,2));
        for ii_trl=1:size(pos_binned_trimmed,1)
            pos_binned_reversed(size(pos_binned_trimmed,1)-ii_trl+1,:)=pos_binned_trimmed(ii_trl,:);
        end

        for ii_shuffled=1:n_shuffle

            for trNo=1:trials.odor_trNo
                y_pred(trNo).data=[];
                x_pred(trNo).data=[];
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

                XYtrain=pos_binned_reversed(this_training_range,:);



                %Decode using neural network
                MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
                MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));

                x_pred(trNo).data=predict(MdlY1,XdFFtest);
                y_pred(trNo).data=predict(MdlY2,XdFFtest);
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

                x_predicted_sh(logical(this_test_range),ii_shuffled)=x_pred(trNo).data;
                y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(trNo).data;
            end

        end

end

fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

x_predictedstart=1;
x_predictedend=length(x_predicted(:,1));

XYtest=pos_binned_trimmed;

handles_out.x_predicted_sh=x_predicted_sh;
handles_out.y_predicted_sh=y_predicted_sh;
handles_out.x_predicted=x_predicted;
handles_out.y_predicted=y_predicted;
handles_out.XYtest=XYtest;
handles_out.trials=trials;
handles_out.MdlY2_pars=MdlY2_pars;
handles_out.MdlY1_pars=MdlY1_pars;

save([this_path arena_file(1:end-4) 'decxy' num2str(which_training_algorithm) '.mat'],'handles_out','handles_choices','-v7.3')

no_conv_points=11;
% conv_win=ones(1,no_conv_points)/no_conv_points;
conv_win_gauss = gausswin(no_conv_points);
conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

x_predicted_conv=conv(x_predicted,conv_win_gauss,'same');
y_predicted_conv=conv(y_predicted,conv_win_gauss,'same');

%Now limit the x and y to max and min
minY1=min(XYtest(:,1));
x_predicted_conv(x_predicted_conv<minY1)=minY1;
maxY1=max(XYtest(:,1));
x_predicted_conv(x_predicted_conv>maxY1)=maxY1;

minY2=min(XYtest(:,2));
y_predicted_conv(y_predicted_conv<minY2)=minY2;
maxY2=max(XYtest(:,2));
y_predicted_conv(y_predicted_conv>maxY2)=maxY2;

x_predicted_sh_conv=zeros(size(x_predicted_sh,1),size(x_predicted_sh,2));
y_predicted_sh_conv=zeros(size(y_predicted_sh,1),size(y_predicted_sh,2));
for ii_sh=1:n_shuffle
    this_x_predicted_sh=zeros(size(x_predicted_sh,1),1);
    this_x_predicted_sh(:,1)=x_predicted_sh(:,ii_sh);
    this_y_predicted_sh=zeros(size(y_predicted_sh,1),1);
    this_y_predicted_sh(:,1)=y_predicted_sh(:,ii_sh);
    this_x_predicted_sh_conv=[];
    this_y_predicted_sh_conv=[];
    this_x_predicted_sh_conv=conv(this_x_predicted_sh,conv_win_gauss,'same');
    this_y_predicted_sh_conv=conv(this_y_predicted_sh,conv_win_gauss,'same');

    %Now limit the x and y to max and min
    minY1=min(XYtest(:,1));
    this_x_predicted_sh_conv(this_x_predicted_sh_conv<minY1)=minY1;
    maxY1=max(XYtest(:,1));
    this_x_predicted_sh_conv(this_x_predicted_sh_conv>maxY1)=maxY1;

    minY2=min(XYtest(:,2));
    this_y_predicted_sh_conv(this_y_predicted_sh_conv<minY2)=minY2;
    maxY2=max(XYtest(:,2));
    this_y_predicted_sh_conv(this_y_predicted_sh_conv>maxY2)=maxY2;

    x_predicted_sh_conv(:,ii_sh)=this_x_predicted_sh_conv;
    y_predicted_sh_conv(:,ii_sh)=this_y_predicted_sh_conv;
end




% fprintf(1, 'R2 nn conv x, y entire run: %d %d\n',drgGetR2(XYtest(:,1),x_predicted_conv),drgGetR2(XYtest(:,2),y_predicted_conv));
R1=corrcoef(XYtest(:,1),x_predicted_conv);
R2=corrcoef(XYtest(:,2),y_predicted_conv);
fprintf(1, 'Correlation coefficient nn conv x, y entire run %d %d\n\n',R1(1,2),R2(1,2));
fprintf(fileID, 'Correlation coefficient nn conv x, y entire run %d %d\n\n',R1(1,2),R2(1,2));

fractionExplained_x = drgMini_calculateVarianceExplained(XYtest(:,1), x_predicted_conv);
fractionExplained_y = drgMini_calculateVarianceExplained(XYtest(:,2), y_predicted_conv);

fprintf(1, 'Fraction of variance predicted nn conv x, y entire run %d %d\n\n',fractionExplained_x,fractionExplained_y);
fprintf(fileID, 'Fraction of variance predicted nn conv x, y entire run %d %d\n\n',fractionExplained_x,fractionExplained_y);

fprintf(1, 'Fraction of variance predicted nn conv x, y entire run %d %d\n\n',fractionExplained_x,fractionExplained_y);
fprintf(fileID, 'Fraction of variance predicted nn conv x, y entire run %d %d\n\n',fractionExplained_x,fractionExplained_y);

fractionExplainedXY = drgMini_calculateVarianceExplainedXY(XYtest(:,1),XYtest(:,2),x_predicted_conv, y_predicted_conv);

fprintf(1, 'Fraction of variance predicted nn conv XY entire run %d %d\n\n',fractionExplainedXY);
fprintf(fileID, 'Fraction of variance predicted nn conv XY entire run %d %d\n\n',fractionExplainedXY);

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(x_predicted_conv(x_predictedstart:x_predictedend,1),y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
plot(XYtest(x_predictedstart:x_predictedend,1),XYtest(x_predictedstart:x_predictedend,2),'-b','LineWidth',1.5)
switch which_training_algorithm
    case 1
        title('xy for neural network convolved (trained per session)')
    case 2
        title('xy for neural network convolved (trained per trial)')
    case 3
        title('xy for neural network convolved (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure2.fig'])
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on

x_predictedstart=1;
x_predictedend=length(x_predicted(:,1));
plot(XYtest(x_predictedstart:x_predictedend,1),'-b','LineWidth',3)
plot(x_predicted_conv(x_predictedstart:x_predictedend,1),'-r','LineWidth',1)

for trNo=1:trials.odor_trNo
    this_x_predictedstart=trials.odor_ii_start(trNo)-10;
    this_x_predictedend=trials.odor_ii_end(trNo)+15;
    plot([this_x_predictedstart:this_x_predictedend],75*ones(this_x_predictedend-this_x_predictedstart+1,1),'-k','LineWidth',2)
end


switch which_training_algorithm
    case 1
        title('x for nn, b:original, r:predicted (trained per session)')
    case 2
        title('x for nn, b:original, r:predicted (trained per trial)')
    case 3
        title('x for nn, b:original, r:predicted (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure3.fig'])
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(XYtest(x_predictedstart:x_predictedend,1),x_predicted_conv(x_predictedstart:x_predictedend,1),'.b')

xlabel('Actual x')
ylabel('Decoded x')

switch which_training_algorithm
    case 1
        title('x for nn (trained per session)')
    case 2
        title('x for nn (trained per trial)')
    case 3
        title('x for nn (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure4.fig'])
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
plot(XYtest(x_predictedstart:x_predictedend,2),'-b','LineWidth',3)
plot(y_predicted_conv(x_predictedstart:x_predictedend,1),'-r','LineWidth',1)

for trNo=1:trials.odor_trNo
    this_x_predictedstart=trials.odor_ii_start(trNo)-10;
    this_x_predictedend=trials.odor_ii_end(trNo)+15;
    plot([this_x_predictedstart:this_x_predictedend],75*ones(this_x_predictedend-this_x_predictedstart+1,1),'-k','LineWidth',2)
end


switch which_training_algorithm
    case 1
        title('y for nn, b:original, r:predicted (trained per session)')
    case 2
        title('y for nn, b:original, r:predicted (trained per trial)')
    case 3
        title('y for nn, b:original, r:predicted (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure5.fig'])
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(XYtest(x_predictedstart:x_predictedend,2),y_predicted_conv(x_predictedstart:x_predictedend,1),'.b')
xlabel('Actual y')
ylabel('Decoded y')

switch which_training_algorithm
    case 1
        title('y for nn(trained per session)')
    case 2
        title('y for nn (trained per trial)')
    case 3
        title('y for nn (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure6.fig'])
end

%Keep track of the per trial decoding
x_all_trials=[];
y_all_trials=[];
x_decod_all_trials=[];
y_decod_all_trials=[];
x_all_trials_sh=[];
y_all_trials_sh=[];
x_decod_all_trials_sh=[];
y_decod_all_trials_sh=[];

x_between_trials=[];
y_between_trials=[];
x_decod_between_trials=[];
y_decod_between_trials=[];
x_between_trials_sh=[];
y_between_trials_sh=[];
x_decod_between_trials_sh=[];
y_decod_between_trials_sh=[];

%Plot the per trial results for y for permuted input
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;
last_x_predictedend=1;
for trNo=1:trials.odor_trNo

    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    y_all_trials_sh=[y_all_trials_sh; XYtest(x_predictedstart:x_predictedend,2)];
    y_decod_all_trials_sh=[y_decod_all_trials_sh; y_predicted_sh_conv(x_predictedstart:x_predictedend,1)];

    y_between_trials_sh=[y_between_trials_sh; XYtest(last_x_predictedend:x_predictedstart,2)];
    y_decod_between_trials_sh=[y_decod_between_trials_sh; y_predicted_sh_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;

    %Plot accuracy per trial
    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    %Plot accuracy per trial for permuted training control
    CIsp = bootci(1000, @mean, y_predicted_sh_conv(x_predictedstart:x_predictedend,:)');
    meansp=mean(y_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;

    [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(y_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1)', CIsp', '-k');

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'Color',[1 0.6 0.6],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-c','LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'Color',[0.4 0.4 1],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
        case 4
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-m','LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;


end

switch which_training_algorithm
    case 1
        title('y for nn, permuted b:original, r:predicted (trained per session)')
    case 2
        title('y for nn, permuted b:original, r:predicted (trained per trial)')
    case 3
        title('y for nn, permuted b:original, r:predicted (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure7.fig'])
end

%Plot the within trial results for y
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;
last_x_predictedend=1;
for trNo=1:trials.odor_trNo

    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    y_all_trials=[y_all_trials; XYtest(x_predictedstart:x_predictedend,2)];
    y_decod_all_trials=[y_decod_all_trials; y_predicted_conv(x_predictedstart:x_predictedend,1)];

    y_between_trials=[y_between_trials; XYtest(last_x_predictedend:x_predictedstart,2)];
    y_decod_between_trials=[y_decod_between_trials; y_predicted_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;

    %Plot accuracy per trial
    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'Color',[1 0.4 0.4],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-c','LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'Color',[0.6 0.6 1],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
        case 4
            %Lane 4 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-m','LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 480],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 480],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;


end

switch which_training_algorithm
    case 1
        title('y for nn within trial, r:hit1, c:mis1, b:hit4, m:miss4 k:predicted (trained per session)')
    case 2
        title('y for nn within trial, r:hit1, c:mis1, b:hit4, m:miss4 k:predicted (trained per trial)')
    case 3
        title('y for nn within trial, r:hit1, c:mis1, b:hit4, m:miss4 k:predicted (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure8.fig'])
end

%Plot the per trial results for x
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;

last_x_predictedend=1;
for trNo=1:trials.odor_trNo

    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    x_all_trials=[x_all_trials; XYtest(x_predictedstart:x_predictedend,1)];
    x_decod_all_trials=[x_decod_all_trials; x_predicted_conv(x_predictedstart:x_predictedend,1)];

    x_between_trials=[x_between_trials; XYtest(last_x_predictedend:x_predictedstart,1)];
    x_decod_between_trials=[x_decod_between_trials; x_predicted_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;

    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'Color',[1 0.4 0.4],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-c','LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'Color',[0.6 0.6 1],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
        case 4
            %Lane 4 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-m','LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;
end

switch which_training_algorithm
    case 1
        title('x for nn per trial, r:hit1, c:mis1, b:hit4, m:miss4 (trained per session)')
    case 2
        title('x for nn per trial, r:hit1, c:mis1, b:hit4, m:miss4 (trained per trial)')
    case 3
        title('x for nn per trial, r:hit1, c:mis1, b:hit4, m:miss4 (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure9.fig'])
end

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
last_x_predictedend=1;
for trNo=1:trials.odor_trNo

    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    x_all_trials_sh=[x_all_trials_sh; XYtest(x_predictedstart:x_predictedend,1)];
    x_decod_all_trials_sh=[x_decod_all_trials_sh; x_predicted_sh_conv(x_predictedstart:x_predictedend,1)];

    x_between_trials_sh=[x_between_trials_sh; XYtest(last_x_predictedend:x_predictedstart,1)];
    x_decod_between_trials_sh=[x_decod_between_trials_sh; x_predicted_sh_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;

    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    %Plot accuracy per trial for permuted training control
    CIsp = bootci(1000, @mean, x_predicted_sh_conv(x_predictedstart:x_predictedend,:)');
    meansp=mean(x_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;

    [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(x_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1)', CIsp', '-k');


    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'Color',[1 0.6 0.6],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-c','LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'Color',[0.4 0.4 1],'LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
        case 4
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-m','LineWidth',3)
            %             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 500],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 500],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;
end

switch which_training_algorithm
    case 1
        title('x for nn per trial permuted, r:hit1, c:miss1, b:hit4, m:miss4, k:predicted (trained per session)')
    case 2
        title('x for nn per trial permuted, r:hit1, c:miss1, b:hit4, m:miss4, k:predicted  (trained per trial)')
    case 3
        title('x for nn per trial permuted, r:hit1, c:miss1, b:hit4, m:miss4, k:predicted (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure10.fig'])
end



% fprintf(1, 'R2 nn conv x, y per trial run: %d %d\n',drgGetR2(x_all_trials,x_decod_all_trials),drgGetR2(y_all_trials,y_decod_all_trials));
R1=corrcoef(x_all_trials,x_decod_all_trials);
R2=corrcoef(y_all_trials,y_decod_all_trials);
fprintf(1, '\n\nCorrelation coefficient nn conv x, y within trials %d %d\n',R1(1,2),R2(1,2));
fprintf(fileID, '\n\nCorrelation coefficient nn conv x, y within trials %d %d\n',R1(1,2),R2(1,2));

R1=corrcoef(x_between_trials,x_decod_between_trials);
R2=corrcoef(y_between_trials,y_decod_between_trials);
fprintf(1, 'Correlation coefficient nn conv x, y between trials %d %d\n\n',R1(1,2),R2(1,2));
fprintf(fileID, 'Correlation coefficient nn conv x, y betweeen trials %d %d\n',R1(1,2),R2(1,2));

% fprintf(1, 'R2 nn permuted x, y per trial run: %d %d\n',drgGetR2(x_all_trials_sh,x_decod_all_trials_sh),drgGetR2(y_all_trials_sh,y_decod_all_trials_sh));
R1=corrcoef(x_all_trials_sh,x_decod_all_trials_sh);
R2=corrcoef(y_all_trials_sh,y_decod_all_trials_sh);
fprintf(1, 'Correlation coefficient nn permuted x, y within trials %d %d\n',R1(1,2),R2(1,2));
fprintf(fileID, 'Correlation coefficient nn permuted x, y within trials %d %d\n',R1(1,2),R2(1,2));

R1=corrcoef(x_between_trials_sh,x_decod_between_trials_sh);
R2=corrcoef(y_between_trials_sh,y_decod_between_trials_sh);
fprintf(1, 'Correlation coefficient nn permuted x, y between trial %d %d\n\n',R1(1,2),R2(1,2));
fprintf(fileID, 'Correlation coefficient nn permuted x, y between trials %d %d\n',R1(1,2),R2(1,2));


fractionExplained_x = drgMini_calculateVarianceExplained(x_all_trials, x_decod_all_trials);
fractionExplained_y = drgMini_calculateVarianceExplained(y_all_trials, y_decod_all_trials);
fprintf(1, 'Fraction of variance predicted nn conv x, y within trials %d %d\n\n',fractionExplained_x,fractionExplained_y);
fprintf(fileID, 'Fraction of variance predicted nn conv x, y within trials %d %d\n\n',fractionExplained_x,fractionExplained_y);

fractionExplainedXY = drgMini_calculateVarianceExplainedXY(x_all_trials,y_all_trials,x_decod_all_trials, y_decod_all_trials);
fprintf(1, 'Fraction of variance predicted nn conv XY within trials %d %d\n\n',fractionExplainedXY);
fprintf(fileID, 'Fraction of variance predicted nn conv XY within trials %d %d\n\n',fractionExplainedXY);

fractionExplained_x = drgMini_calculateVarianceExplained(x_between_trials, x_decod_between_trials);
fractionExplained_y = drgMini_calculateVarianceExplained(y_between_trials, y_decod_between_trials);
fprintf(1, '\n\nFraction of variance predicted nn conv x, y between trials %d %d\n\n',fractionExplained_x,fractionExplained_y);
fprintf(fileID, '\n\nFraction of variance predicted nn conv x, y between trials %d %d\n\n',fractionExplained_x,fractionExplained_y);

fractionExplainedXY = drgMini_calculateVarianceExplainedXY(x_between_trials,y_between_trials,x_decod_between_trials, y_decod_between_trials);
fprintf(1, 'Fraction of variance predicted nn conv XY between trials %d %d\n\n',fractionExplainedXY);
fprintf(fileID, 'Fraction of variance predicted nn conv XY between trials %d %d\n\n',fractionExplainedXY);

fractionExplained_x = drgMini_calculateVarianceExplained(x_all_trials_sh, x_decod_all_trials_sh);
fractionExplained_y = drgMini_calculateVarianceExplained(y_all_trials_sh, y_decod_all_trials_sh);
fprintf(1, '\n\nFraction of variance predicted nn permuted x, y within trials %d %d\n\n',fractionExplained_x,fractionExplained_y);
fprintf(fileID, '\n\nFraction of variance predicted nn permuted x, y within trials %d %d\n\n',fractionExplained_x,fractionExplained_y);

fractionExplainedXY = drgMini_calculateVarianceExplainedXY(x_all_trials_sh,y_all_trials_sh,x_decod_all_trials_sh, y_decod_all_trials_sh);
fprintf(1, 'Fraction of variance predicted nn permuted XY within trials %d %d\n\n',fractionExplainedXY);
fprintf(fileID, 'Fraction of variance predicted nn permuted XY within trials %d %d\n\n',fractionExplainedXY);

fractionExplained_x = drgMini_calculateVarianceExplained(x_between_trials_sh, x_decod_between_trials_sh);
fractionExplained_y = drgMini_calculateVarianceExplained(y_between_trials_sh, y_decod_between_trials_sh);
fprintf(1, '\n\nFraction of variance predicted nn permuted x, y between trials %d %d\n\n',fractionExplained_x,fractionExplained_y);
fprintf(fileID, '\n\nFraction of variance predicted nn permuted x, y between trials %d %d\n\n',fractionExplained_x,fractionExplained_y);

fractionExplainedXY = drgMini_calculateVarianceExplainedXY(x_between_trials_sh,y_between_trials_sh,x_decod_between_trials_sh, y_decod_between_trials_sh);
fprintf(1, 'Fraction of variance predicted nn permuted XY between trials %d %d\n\n',fractionExplainedXY);
fprintf(fileID, 'Fraction of variance predicted nn cpermuted XY within trials %d %d\n\n',fractionExplainedXY);




%Plot x vs decoded
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

plot(x_all_trials,x_decod_all_trials,'.b')
xlabel('x')
ylabel('x decoded')
switch which_training_algorithm
    case 1
        title('x for nn per trial(trained per session)')
    case 2
        title('x for nn per trial (trained per trial)')
    case 3
        title('x for nn per trial (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure11.fig'])
end

%Plot y vs decoded
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

plot(y_all_trials,y_decod_all_trials,'.b')
xlabel('y')
ylabel('y decoded')


switch which_training_algorithm
    case 1
        title('y for nn per trial (trained per session)')
    case 2
        title('y for nn per trial (trained per trial)')
    case 3
        title('y for nn per trial (trained between)')
end

if handles_choices.save_results==1
    savefig([this_path 'figure12.fig'])
end

fclose(fileID);

pffft=1;