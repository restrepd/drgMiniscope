function handles_outXY=drgMini_DecodeOdorArena_SubROI_gpu(handles_choices)
%Does decoding of the navigation path for a mouse undergoing odor plume navigation 
%following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020


close all

if exist('handles_choices')==0
    clear all

    handles_choices.save_results=1;
    handles_choices.is_sphgpu=1; %1 sphgpu, 2 Alpine
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

    % sphgpu
    this_path='/data2/SFTP/PreProcessed/20220713_FCM6/';
    dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';

    handles_choices.save_path='/data2/SFTP/PreProcessed/Temp/';

    % %Alpine
    % this_path='/scratch/alpine/drestrepo@xsede.org/PreProcessed/20220713_FCM6/';
    % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';


    
    

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
    %3 fitrsvz

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

    % handles_choices.no_sub_ROIs=[1 5 10 20 40 500];
    % handles_choices.no_runs__per_sub_ROI=[1 1 1 1 1];
    handles_choices.no_sub_ROIs=[1 5 10 20 40 500];
    handles_choices.no_runs__per_sub_ROI=[20 10 10 10 10 1];

    handles_choices.save_path_conc='/data2/SFTP/DecodeDynOdorConcSub_1_01082026/';
    handles_choices.generate_idx=0;
    handles_choices.save_tag_conc='dynconctree1';
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
        addpath('/data2/DRMatlab/drgMiniscope')
        addpath('/data2/DRMatlab/m new/Chi Squared')
        addpath('/data2/DRMatlab/drgMaster')
        addpath(genpath('/data2/DRMatlab/m new/kakearney-boundedline-pkg-32f2a1f'))
    case 2
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/drgMiniscope')
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/m new/Chi Squared')
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/drgMaster')
        addpath(genpath('/projects/drestrepo@xsede.org/software/DR_matlab/m new/kakearney-boundedline-pkg-32f2a1f'))
end

if ~exist(handles_choices.save_path(1:end-1),'dir')
    mkdir(handles_choices.save_path(1:end-1))
end
 
fileID = fopen([this_path 'R1R2output.txt'],'w');


fprintf(1,['\nTrained with within trial data\n\n'])
fprintf(fileID,['\nTrained with within trial data\n\n'])


figNo=0;



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

XYtest=pos_binned;

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



handles_outXY.no_neurons=no_neurons;

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

%Perform decoding for each subset of ROIs
% handles_choices.no_sub_ROIs=[1 5 10 20 500];
%  handles_choices.no_runs__per_sub_ROI=[1 1 1 1 1];
ROI_indices=1:no_neurons;

for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    start_toc=toc;
    % parfor trNo=1:trials.odor_trNo
    if handles_choices.no_sub_ROIs(ROIsubset)<=no_neurons
        this_no_neurons=handles_choices.no_sub_ROIs(ROIsubset);
    else
        this_no_neurons=no_neurons;
    end
    % handles_choices.no_sub_ROIs(ROIsubset_ii) subsets (rows), each of length handles_choices.no_runs__per_sub_ROI(ROIsubset_ii) (columns), with replacement
    if handles_choices.generate_idx==1
        idx = randi(numel(ROI_indices), handles_choices.no_runs__per_sub_ROI(ROIsubset), this_no_neurons);
    else
        %Read idx
        load([handles_choices.save_path_conc arena_file(1:end-4) handles_choices.save_tag_conc '_subROI.mat'],'handles_out')
        idx=handles_out.ROIsubset(ROIsubset).idx;
    end

    handles_outXY.ROIsubset(ROIsubset).idx=idx;
    for ii_sub=1:handles_choices.no_runs__per_sub_ROI(ROIsubset)

        X_dFF_trimmed=zeros(size(X_dFF,1),this_no_neurons);
        these_iiROIs=zeros(1,this_no_neurons);
        these_iiROIs(1,:)=idx(ii_sub,:);

        X_dFF_trimmed(:,:)=X_dFF(:,these_iiROIs);

        % parfor trNo=1:trials.odor_trNo
        for trNo=1:trials.odor_trNo
            

            %Please note that I decode the entire time course
            %Below I parse out the within and between trials
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


            XdFFtrain=X_dFF_trimmed(this_training_range,:);
            XdFFtest=X_dFF_trimmed(logical(this_test_range),:);

            XYtrain=pos_binned(this_training_range,:);

            %Now save the start to end
            this_trial_test_range=zeros(1,no_time_bins);
            ii_trial_test_range_start=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
            if ii_trial_test_range_start<1
                ii_trial_test_range_start=1;
            end
            ii_trial_test_range_end=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
            if ii_trial_test_range_end>size(pos_binned,1)
                ii_trial_test_range_end=size(pos_binned,1);
            end
            this_trial_test_range(ii_trial_test_range_start:ii_trial_test_range_end)=1;

            trials.trial(trNo).XYtest=pos_binned(logical(this_trial_test_range),:);
            trials.trial(trNo).XdFFtest=X_dFF(logical(this_trial_test_range),:);

            %Use gpu
            XdFFtrain_gpu=gpuArray(XdFFtrain);
            XdFFtest_gpu=gpuArray(XdFFtest);
            XYtrain_gpu=gpuArray(XYtrain);

            switch handles_choices.algo
                case 1
                    MdlY1_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,1),'Standardize',true);
                    MdlY2_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,2),'Standardize',true);

                    %Gather back the data from the gpu
                    MdlY1=gather(MdlY1_gpu);
                    MdlY2=gather(MdlY2_gpu);
                case 2
                    MdlY1_gpu = fitrtree(XdFFtrain_gpu,XYtrain_gpu(:,1));
                    impx_gpu=predictorImportance(MdlY1_gpu);
                    %Gather back imp from the gpu
                    impx=gather(impx_gpu);
                    MdlY2_gpu = fitrtree(XdFFtrain_gpu,XYtrain_gpu(:,2));
                    impy_gpu=predictorImportance(MdlY2_gpu);
                    %Gather back imp from the gpu
                    impy=gather(impy_gpu);
                    % handles_outXY.imp.trial(trNo).impx=impx;
                    % handles_outXY.imp.trial(trNo).impy=impy;

                    %Gather back the data from the gpu
                    MdlY1=gather(MdlY1_gpu);
                    MdlY2=gather(MdlY2_gpu);
                case 3
                    %SVM
                    opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                        'Verbose', 0, ...
                        'MaxObjectiveEvaluations', 15,...
                        'UseParallel',true);

                    fig = figure('Visible','off');
                    MdlY1 = fitrsvm(XdFFtrain,XYtrain(:,1),...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                    MdlY2 = fitrsvm(XdFFtrain,XYtrain(:,1),...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                case 4
                    %GP
                    opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                        'Verbose', 0, ...
                        'MaxObjectiveEvaluations', 15,...
                        'UseParallel',true);

                    fig = figure('Visible','off');
                    %             MdlY1_gpu = fitrgp(XdFFtrain_gpu,XYtrain_gpu(:,1));
                    %             MdlY2_gpu = fitrgp(XdFFtrain_gpu,XYtrain_gpu(:,1));
                    MdlY1 = fitrgp(XdFFtrain,XYtrain(:,1),'KernelFunction','squaredexponential',...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                    MdlY2 = fitrgp(XdFFtrain,XYtrain(:,1),'KernelFunction','squaredexponential',...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                    %Gather back the data from the gpu
                    %             MdlY1=gather(MdlY1_gpu);
                    %             MdlY2=gather(MdlY2_gpu);
                case 5
                    MdlY1 = fitglm(XdFFtrain,XYtrain(:,1));
                    MdlY2 = fitglm(XdFFtrain,XYtrain(:,2));

            end



            x_pred(trNo).data=predict(MdlY1,XdFFtest);
            y_pred(trNo).data=predict(MdlY2,XdFFtest);

            x_pred(trNo).MdlY1=MdlY1;
            y_pred(trNo).MdlY2=MdlY2;

            %     fprintf(1,['Elapsed time ' num2str(toc-start_toc) ' for trial number ' num2str(trNo) ' \n\n'])
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
        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).x_predicted=x_predicted;
        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).y_predicted=y_predicted;

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

        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).x_predicted_conv=x_predicted_conv;
        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).y_predicted_conv=y_predicted_conv;
        
        fprintf(1,['Done with ROI susbet ' num2str(ROIsubset) ' run number ' num2str(ii_sub) '\n'])
    end
    fprintf(1,['Elapsed time for ROI susbet ' num2str(ROIsubset) ' ' num2str((toc-start_toc)/(60*60)) ' hrs\n\n'])

end


%Now do predictions for reversed/permuted training periods
for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    start_toc=toc;
    % parfor trNo=1:trials.odor_trNo
    if handles_choices.no_sub_ROIs(ROIsubset)<=no_neurons
        this_no_neurons=handles_choices.no_sub_ROIs(ROIsubset);
    else
        this_no_neurons=no_neurons;
    end
    % handles_choices.no_sub_ROIs(ROIsubset_ii) subsets (rows), each of length handles_choices.no_runs__per_sub_ROI(ROIsubset_ii) (columns), with replacement
    if handles_choices.generate_idx==1
        idx = randi(numel(ROI_indices), handles_choices.no_runs__per_sub_ROI(ROIsubset), this_no_neurons);
    else
        %Read idx
        load([handles_choices.save_path_conc arena_file(1:end-4) handles_choices.save_tag_conc '_subROI.mat'],'handles_out')
        idx=handles_out.ROIsubset(ROIsubset).idx;
    end

    handles_outXY.ROIsubset(ROIsubset).idx=idx;
    for ii_sub=1:handles_choices.no_runs__per_sub_ROI(ROIsubset)

        %We will do a reversal and a circular permutation
        pos_binned_reversed=zeros(size(pos_binned,1),size(pos_binned,2));
        offset_ii=(ii_sub-1)*floor(size(pos_binned,1)/handles_choices.no_runs__per_sub_ROI(ROIsubset));
        % offset_ii=(ii_shuffled-1)*floor(size(pos_binned,1)/n_shuffle);
        for ii_trl=1:size(pos_binned,1)
            this_ii_trl=ii_trl+offset_ii;
            if this_ii_trl>size(pos_binned,1)
                offset_ii=-ii_trl+1;
                this_ii_trl=ii_trl+offset_ii;
            end
            pos_binned_reversed(size(pos_binned,1)-ii_trl+1,:)=pos_binned(this_ii_trl,:);
        end

        X_dFF_trimmed=zeros(size(X_dFF,1),this_no_neurons);
        these_iiROIs=zeros(1,this_no_neurons);
        these_iiROIs(1,:)=idx(ii_sub,:);

        X_dFF_trimmed(:,:)=X_dFF(:,these_iiROIs);

        % parfor trNo=1:trials.odor_trNo
        for trNo=1:trials.odor_trNo
            

            %Please note that I decode the entire time course
            %Below I parse out the within and between trials
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


            XdFFtrain=X_dFF_trimmed(this_training_range,:);
            XdFFtest=X_dFF_trimmed(logical(this_test_range),:);

            XYtrain=pos_binned_reversed(this_training_range,:);

            %Now save the start to end
            this_trial_test_range=zeros(1,no_time_bins);
            ii_trial_test_range_start=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
            if ii_trial_test_range_start<1
                ii_trial_test_range_start=1;
            end
            ii_trial_test_range_end=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
            if ii_trial_test_range_end>size(pos_binned,1)
                ii_trial_test_range_end=size(pos_binned,1);
            end
            this_trial_test_range(ii_trial_test_range_start:ii_trial_test_range_end)=1;

            trials.trial(trNo).XYtest=pos_binned(logical(this_trial_test_range),:);
            trials.trial(trNo).XdFFtest=X_dFF(logical(this_trial_test_range),:);

            %Use gpu
            XdFFtrain_gpu=gpuArray(XdFFtrain);
            XdFFtest_gpu=gpuArray(XdFFtest);
            XYtrain_gpu=gpuArray(XYtrain);

            switch handles_choices.algo
                case 1
                    MdlY1_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,1),'Standardize',true);
                    MdlY2_gpu = fitrnet(XdFFtrain_gpu,XYtrain_gpu(:,2),'Standardize',true);

                    %Gather back the data from the gpu
                    MdlY1=gather(MdlY1_gpu);
                    MdlY2=gather(MdlY2_gpu);
                case 2
                    MdlY1_gpu = fitrtree(XdFFtrain_gpu,XYtrain_gpu(:,1));
                    impx_gpu=predictorImportance(MdlY1_gpu);
                    %Gather back imp from the gpu
                    impx=gather(impx_gpu);
                    MdlY2_gpu = fitrtree(XdFFtrain_gpu,XYtrain_gpu(:,2));
                    impy_gpu=predictorImportance(MdlY2_gpu);
                    %Gather back imp from the gpu
                    impy=gather(impy_gpu);
                    % handles_outXY.imp.trial(trNo).impx=impx;
                    % handles_outXY.imp.trial(trNo).impy=impy;

                    %Gather back the data from the gpu
                    MdlY1=gather(MdlY1_gpu);
                    MdlY2=gather(MdlY2_gpu);
                case 3
                    %SVM
                    opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                        'Verbose', 0, ...
                        'MaxObjectiveEvaluations', 15,...
                        'UseParallel',true);

                    fig = figure('Visible','off');
                    MdlY1 = fitrsvm(XdFFtrain,XYtrain(:,1),...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                    MdlY2 = fitrsvm(XdFFtrain,XYtrain(:,1),...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                case 4
                    %GP
                    opts = struct('AcquisitionFunctionName','expected-improvement-plus',...
                        'Verbose', 0, ...
                        'MaxObjectiveEvaluations', 15,...
                        'UseParallel',true);

                    fig = figure('Visible','off');
                    %             MdlY1_gpu = fitrgp(XdFFtrain_gpu,XYtrain_gpu(:,1));
                    %             MdlY2_gpu = fitrgp(XdFFtrain_gpu,XYtrain_gpu(:,1));
                    MdlY1 = fitrgp(XdFFtrain,XYtrain(:,1),'KernelFunction','squaredexponential',...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                    MdlY2 = fitrgp(XdFFtrain,XYtrain(:,1),'KernelFunction','squaredexponential',...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        opts);
                    %Gather back the data from the gpu
                    %             MdlY1=gather(MdlY1_gpu);
                    %             MdlY2=gather(MdlY2_gpu);
                case 5
                    MdlY1 = fitglm(XdFFtrain,XYtrain(:,1));
                    MdlY2 = fitglm(XdFFtrain,XYtrain(:,2));

            end



            x_pred(trNo).data=predict(MdlY1,XdFFtest);
            y_pred(trNo).data=predict(MdlY2,XdFFtest);

            x_pred(trNo).MdlY1=MdlY1;
            y_pred(trNo).MdlY2=MdlY2;

            %     fprintf(1,['Elapsed time ' num2str(toc-start_toc) ' for trial number ' num2str(trNo) ' \n\n'])
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
        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).x_predicted_sh=x_predicted;
        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).y_predicted_sh=y_predicted;

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

        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).x_predicted_conv_sh=x_predicted_conv;
        handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).y_predicted_conv_sh=y_predicted_conv;
        
        fprintf(1,['Done with ROI susbet ' num2str(ROIsubset) ' run number ' num2str(ii_sub) '\n'])
    end
    fprintf(1,['Elapsed time for ROI susbet ' num2str(ROIsubset) ' ' num2str((toc-start_toc)/(60*60)) ' hrs\n\n'])

end
fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

x_predictedstart=1;
x_predictedend=length(x_predicted(:,1));

XYtest=pos_binned;


handles_outXY.XYtest=XYtest;
handles_outXY.trials=trials;


save([handles_choices.save_path arena_file(1:end-4) handles_choices.save_tag '_subROI.mat'],'handles_outXY','handles_choices','-v7.3')



pffft=1;
