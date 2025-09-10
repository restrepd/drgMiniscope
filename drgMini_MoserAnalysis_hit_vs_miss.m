function handles_out=drgMini_MoserAnalysis_hit_vs_miss(handles_choices)
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020

warning('off')

close all

if exist('handles_choices')==0
    clear all
  

    % handles_choices.resume_processing=1; %Make this 1 if you want the program to keep all fits already calculated

  

    handles_choice.is_sphgpu=0;
    is_sphgpu=handles_choice.is_sphgpu;

    switch is_sphgpu
        case 0 %Diego's mac
            handles_choices.save_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/temp/';
        case 1 %sphgpu
            handles_choices.save_path='/data/SFTP/DecodeOdorConc/';
        case 2 %Alpine
    end

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
   
    %Files for dFF data

    %First troubleshooting file
    if is_sphgpu==1
        this_path='/data/SFTP/PreProcessedDR/20220804_FCM22/';
    else
        this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    end
    dFF_file='20220804_FCM22_withodor_miniscope_sync__hitand_miss_ncorre_ext.mat';

    % if handles_choices.distances_in_mm==0
    %     arena_file='20220804_FCM22withodor_odorarena__hitand_miss_sync.mat';
    % else
    arena_file='20220804_FCM22withodor_odorarena__hitand_miss_sync_mm.mat';
    % end
    handles_choices.group=1;


    %Miss only
    % this_path='/data/SFTP/PreProcessedDR/20220526_FCM6_withodor_lane4/';
    % 
    % dFF_file='20220526_FCM6_withodor_miniscope_sync__miss_ncorre_ext_nonneg.mat';
    % 
    % arena_file='20220526_FCM6withodor_odorarena__miss_sync_mm.mat';
    % 
    % handles_choices.group=2;

    %Group 1 is rewarded, odor ISO1 in both Hit and Miss
    %Group 2 is rewarded, with odor Miss, no odor in Hit
    %Group 3 is rewarded, with odor Hit, no odor in Miss
    %Group 4 is rewarded, with no odor in Hit and Miss


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
    % dFF_file='20220713_FCM6_withodor_miniscope_sync__hitand_miss_ncorre_ext.mat';
    % 
    % if handles_choices.distances_in_mm==0
    %     arena_file='20220713_FCM6withodor_odorarena__hitand_miss_syn.mat';
    % else
    %     arena_file='20220713_FCM6withodor_odorarena__hitand_miss_syn_mm.mat'; %Distances converted to mm
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
    % dFF_file='20220727_FCM19_withodor_miniscope_sync__hitand_miss_ncorre_ext_nonneg.mat';
    % arena_file='20220727_FCM19withodor_odorarena__hitand_miss_sync_mm.mat'; %Distances converted to mm
   
    % handles_choices.dt_decoding_op=2; %All bins from -dt_decoding to 0 will be included in the decoder
    % %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    % %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    % handles_choices.dt_decoding_xy=0.6;


    %No odor troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    % dFF_file='20220824_FCM6_withoutodor_miniscope_sync__hitand_miss_ncorre_ext.mat';
    % arena_file='20220824_FCM6withoutodor_odorarena__hitand_miss_sync.mat';
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




   
    % no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
    % fprintf(1,['Number of time bins for op decoding ' num2str(no_dec_time_bins_op) '\n\n'])
    % 
    % 
    % 
    % no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
    % fprintf(1,['Number of time bins for xy decoding ' num2str(no_dec_time_bins_xy) '\n\n'])


    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    % handles_choices.n_shuffle=n_shuffle;


      handles_choices.trial_start_offset=-15; %This was -10
    handles_choices.trial_end_offset=15;

    handles_choices.save_tag='_moser';

else

    %Note: The data brought into the Kording lab jupyter notebook seems to be
    %binned in 200 msec bins
    dt=handles_choices.dt;
    dt_miniscope=handles_choices.dt_miniscope;
    % n_shuffle=handles_choices.n_shuffle;

    % distances_in_mm=handles_choices.distances_in_mm;

    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;



    % which_odor_plume=handles_choices.which_odor_plume;   %1=use odor plume from the publication, 2=use odor plume simulated by Aaron
    % cm_from_floor=handles_choices.cm_from_floor;
    % 
    % weber_fechner=handles_choices.weber_fechner;
    % alpha=handles_choices.alpha;
    % multiplier=handles_choices.multiplier;
    % %0 is Stevens Law, R proportional to C^alpha
    % %1 is Weber-Flechner law R proportional to log(C)
    % %See Copelli et al DOI: 10.1103/PhysRevE.65.060901
    % 
    % algo=handles_choices.algo;
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




    % dt_decoding_op=handles_choices.dt_decoding_op; %All bins from -dt_decoding to 0 will be included in the decoder
    %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    % no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
    % fprintf(1,['Number of time bins for op decoding ' num2str(no_dec_time_bins_op) '\n\n'])


    % dt_decoding_xy=handles_choices.dt_decoding_xy;
    % no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
    % fprintf(1,['Number of time bins for xy decoding ' num2str(no_dec_time_bins_xy) '\n\n'])
    is_sphgpu=handles_choices.is_sphgpu;
    n_shuffle_SI=handles_choices.n_shuffle_SI;

end

if is_sphgpu==1
    addpath('/data2/DRMatlab/drgMiniscope')
    addpath('/data2/DRMatlab/m new/Chi Squared')
    addpath('/data2/DRMatlab/drgMaster')
end

if ~exist(handles_choices.save_path(1:end-1),'dir')
    mkdir(handles_choices.save_path(1:end-1))
end

max_overlap=handles_choices.max_overlap;

figNo=0;



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

%Extract lanes using FLIR data
at_end=0;
ii=0;
trNo=0;
trNo__hit=0;
trNo__miss=0;
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
trNo__hit=0;
trNo__miss=0;
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
title('Trial start (o) and water delivery (x), red Hit, blue Miss')



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
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct for Hit ' num2str(percent_correct1) ' percent correct for Miss ' num2str(percent_correct4) '\n\n'])

%Do z scores
mean_neural_data_col=mean(neural_data,1);
mean_neural_data=repmat(mean_neural_data_col,no_time_bins,1);

std_neural_data_col=std(neural_data,1);
std_neural_data=repmat(std_neural_data_col,no_time_bins,1);

neural_data=(neural_data-mean_neural_data)./std_neural_data;

no_time_bins=size(neural_data,1);



pffft=1;


for this_ROI=1:no_neurons

    for trNo=1:trials.odor_trNo

        this_test_range=zeros(1,no_time_bins);

        ii_test_range_start=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
        ii_test_range_end=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;

        if ii_test_range_start<1
            ii_test_range_start=1;
        end

        if ii_test_range_end>no_time_bins
            ii_test_range_end=no_time_bins;
        end

        this_test_range(ii_test_range_start:ii_test_range_end)=1;

        this_dFFtest=zeros(sum(this_test_range),1);
        this_dFFtest(:,1)=neural_data(logical(this_test_range),this_ROI);
        dFF_pred(trNo).dFF=this_dFFtest;

        dFF_pred(trNo).this_XY=pos_binned(logical(this_test_range),:);

    end

    all_dFF=[];
    all_hit_vs_miss=[];
    all_XY=[];
    all_dFFHit=[];
    all_dFFMiss=[];

    all_XYHit=[];
    all_XYMiss=[];
    all_tr=[];
    all_trHit=[];
    all_trMiss=[];
    ii_trHit=0;
    ii_trMiss=0;

    for trNo=1:trials.odor_trNo

        all_dFF=[all_dFF dFF_pred(trNo).dFF'];
        all_tr(trNo)=size(all_XY,1)+1;
        all_XY=[all_XY; dFF_pred(trNo).this_XY];

        if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
            ii_trHit=ii_trHit+1;
            all_trHit(ii_trHit)=size(all_XYHit,1)+1;
            all_XYHit=[all_XYHit; dFF_pred(trNo).this_XY];
            all_dFFHit=[all_dFFHit dFF_pred(trNo).dFF'];
            all_hit_vs_miss=[all_hit_vs_miss ones(1,length(dFF_pred(trNo).dFF))];
        else
            ii_trMiss=ii_trMiss+1;
            all_trMiss(ii_trMiss)=size(all_XYMiss,1)+1;
            all_XYMiss=[all_XYMiss; dFF_pred(trNo).this_XY];
            all_dFFMiss=[all_dFFMiss dFF_pred(trNo).dFF'];
            all_hit_vs_miss=[all_hit_vs_miss 4*ones(1,length(dFF_pred(trNo).dFF))];
        end
    end

    per_ROI(this_ROI).results.all_hit_vs_miss=all_hit_vs_miss;

    per_ROI(this_ROI).results.all_XY=all_XY;
    per_ROI(this_ROI).results.all_XYHit=all_XYHit;
    per_ROI(this_ROI).results.all_XYMiss=all_XYMiss;
    per_ROI(this_ROI).results.all_tr=all_tr;
    per_ROI(this_ROI).results.all_trHit=all_trHit;
    per_ROI(this_ROI).results.all_trMiss=all_trMiss;
    per_ROI(this_ROI).results.all_dFF=all_dFF;
    per_ROI(this_ROI).results.all_dFFHit=all_dFFHit;
    per_ROI(this_ROI).results.all_dFFMiss=all_dFFMiss;

end


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



%Calculate activity maps and information content
figNo=figNo+1;
information_content=zeros(no_neurons,1);
sparsity=zeros(no_neurons,1);
information_content_hit=zeros(no_neurons,1);
sparsity_hit=zeros(no_neurons,1);
information_content_miss=zeros(no_neurons,1);
sparsity_miss=zeros(no_neurons,1);
information_content_hit_miss=zeros(no_neurons,1);
sparsity_miss=zeros(no_neurons,1);

spatial_rho_hit_vs_miss=zeros(no_neurons,1);
delta_center_of_mass=zeros(no_neurons,1);

sh_information_content=zeros(no_neurons,n_shuffle_SI);
sh_sparsity=zeros(no_neurons,n_shuffle_SI);
sh_information_content_hit=zeros(no_neurons,n_shuffle_SI);
sh_sparsity_hit=zeros(no_neurons,n_shuffle_SI);
sh_information_content_miss=zeros(no_neurons,n_shuffle_SI);
sh_sparsity_miss=zeros(no_neurons,n_shuffle_SI);
sh_information_content_hit_miss=zeros(no_neurons,n_shuffle_SI);
sh_sparsity_miss=zeros(no_neurons,n_shuffle_SI);

sh_spatial_rho_hit_vs_miss=zeros(no_neurons,n_shuffle_SI);
sh_delta_center_of_mass=zeros(no_neurons,n_shuffle_SI);

%Diego's mutual information calculation
dr_information_content_hit=zeros(no_neurons,1); %Mutual information between binary dFF and xy
dr_information_content_miss=zeros(no_neurons,1);
dr_information_contentl_mutual_1_miss=zeros(no_neurons,1);

sh_dr_information_content_hit=zeros(no_neurons,n_shuffle_SI); %Mutual information between binary dFF and xy
sh_dr_information_content_miss=zeros(no_neurons,n_shuffle_SI);
sh_dr_information_contentl_mutual_1_miss=zeros(no_neurons,n_shuffle_SI);

glm_pvalues=[];
for this_ROI=1:no_neurons

    %Initialize variables
    this_dFF_activity=zeros(10,10);
    this_dFF_activity_n=zeros(10,10);
    sum_dFF_activity=0;

    this_dFF_hit_activity=zeros(10,10);
    this_dFF_hit_activity_n=zeros(10,10);
    sum_dFF_hit_activity=0;

    this_dFF_miss_activity=zeros(10,10);
    this_dFF_miss_activity_n=zeros(10,10);
    sum_dFF_miss_activity=0;

    ii_repeat=1;
    these_all_dFF=(per_ROI(this_ROI).results.all_dFF)-min(per_ROI(this_ROI).results.all_dFF);
    these_all_hit_vs_miss=per_ROI(this_ROI).results.all_hit_vs_miss;

    %GLM
    glm_dFF=[];
    glm_dFF_ii=0;

    glm_dFF_hit=[];
    glm_dFF_hit_ii=0;

    glm_dFF_miss=[];
    glm_dFF_miss_ii=0;

    glm_dFF_both=[];
    glm_dFF_both_ii=0;

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
            glm_dFF.lane(glm_dFF_ii)=0; %This is the Miss area in the arena
        else
            if this_y_ii>=7
                glm_dFF.lane(glm_dFF_ii)=1; %This is Hit area in the areana
            else
                glm_dFF.lane(glm_dFF_ii)=2; %This is the center
            end
        end

        if these_all_hit_vs_miss(ii_t)==4
            glm_dFF.hit_vs_miss_trial(glm_dFF_ii)=0; %Miss trial
            glm_dFF_both_ii=glm_dFF_both_ii+1;
            glm_dFF_both.data(glm_dFF_both_ii)=these_all_dFF(ii_t);
            glm_dFF_both.hit_vs_miss_trial(glm_dFF_both_ii)=0; %Miss trial
            glm_dFF_miss_ii=glm_dFF_miss_ii+1;
            glm_dFF_miss.data(glm_dFF_miss_ii)=these_all_dFF(ii_t);
            glm_dFF_miss.xy(glm_dFF_miss_ii)=(this_y_ii-1)*10+this_x_ii;
        else
            glm_dFF.hit_vs_miss_trial(glm_dFF_ii)=1; %Hit trial
            glm_dFF_both_ii=glm_dFF_both_ii+1;
            glm_dFF_both.data(glm_dFF_both_ii)=these_all_dFF(ii_t);
            glm_dFF_both.hit_vs_miss_trial(glm_dFF_both_ii)=1; %Hit trial
            glm_dFF_hit_ii=glm_dFF_hit_ii+1;
            glm_dFF_hit.data(glm_dFF_hit_ii)=these_all_dFF(ii_t);
            glm_dFF_hit.xy(glm_dFF_hit_ii)=(this_y_ii-1)*10+this_x_ii;
        end

        if these_all_hit_vs_miss(ii_t)==1
            this_dFF_hit_activity(this_x_ii,this_y_ii)=this_dFF_hit_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
            this_dFF_hit_activity_n(this_x_ii,this_y_ii)=this_dFF_hit_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFF_hit_activity=sum_dFF_hit_activity+these_all_dFF(ii_t);
        else
            this_dFF_miss_activity(this_x_ii,this_y_ii)=this_dFF_miss_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
            this_dFF_miss_activity_n(this_x_ii,this_y_ii)=this_dFF_miss_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFF_miss_activity=sum_dFF_miss_activity+these_all_dFF(ii_t);
        end
    end




    %Perform the glm
    fprintf(1, ['glm for hit vs miss trials with xy\n'])

    fprintf(1, ['\n\nglm for dFF\n'])
    tbl = table(glm_dFF.data',glm_dFF.xy',glm_dFF.hit_vs_miss_trial',...
        'VariableNames',{'dFF','xy','trial'});
    mdl = fitglm(tbl,'dFF~xy+trial'...
        ,'CategoricalVars',[3])

    glm_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;

    %Perform the glm
    fprintf(1, ['glm for hit vs miss trials without xy\n'])

    fprintf(1, ['\n\nglm for dFF hit vs miss, no xy\n'])
    tbl = table(glm_dFF_both.data',glm_dFF_both.hit_vs_miss_trial',...
        'VariableNames',{'dFF','trial'});
    mdl = fitglm(tbl,'dFF~trial'...
        ,'CategoricalVars',[2])

    glm__both_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;

    %Perform the Hit glm
    fprintf(1, ['glm for hit trials with xy\n'])

    fprintf(1, ['\n\nglm for dFF hit\n'])
    if ~isempty(glm_dFF_hit)
        tbl = table(glm_dFF_hit.data',glm_dFF_hit.xy',...
            'VariableNames',{'dFF','xy'});
        mdl = fitglm(tbl,'dFF~xy')

        glm__hit_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;
    end

    %Perform the Miss glm
    if ~isempty(glm_dFF_miss)
        fprintf(1, ['glm for miss trials with xy\n'])

        fprintf(1, ['\n\nglm for dFF miss\n'])
        tbl = table(glm_dFF_miss.data',glm_dFF_miss.xy',...
            'VariableNames',{'dFF','xy'});
        mdl = fitglm(tbl,'dFF~xy')

        glm__miss_pvalues.ROI(this_ROI).pValues=mdl.Coefficients.pValue;
    end


    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)~=0
                this_dFF_activity(ii_x,ii_y)=this_dFF_activity(ii_x,ii_y)/this_dFF_activity_n(ii_x,ii_y);
                if this_dFF_hit_activity_n(ii_x,ii_y)>0
                    this_dFF_hit_activity(ii_x,ii_y)=this_dFF_hit_activity(ii_x,ii_y)/this_dFF_hit_activity_n(ii_x,ii_y);
                end
                if this_dFF_miss_activity_n(ii_x,ii_y)>0
                    this_dFF_miss_activity(ii_x,ii_y)=this_dFF_miss_activity(ii_x,ii_y)/this_dFF_miss_activity_n(ii_x,ii_y);
                end
            end
        end
    end




    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .8 .6])

    %Plot the fits
    no_time_bins1=length(per_ROI(this_ROI).results.all_dFFHit);
    no_time_bins4=length(per_ROI(this_ROI).results.all_dFFMiss);

    y_shift=max(per_ROI(this_ROI).results.all_dFF);
    if y_shift==0
        y_shift=1;
    end

    %Plot dFF timecourses
    subplot(2,3,1:3)
    hold on
    ii_plot=0;
    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFHit+y_shift*ii_plot,'-k','LineWidth',1.5)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFMiss+y_shift*ii_plot,'-k','LineWidth',1.5)


    %Hit
    for ii_tr=1:length(per_ROI(this_ROI).results.all_trHit)
        this_x=per_ROI(this_ROI).results.all_trHit(ii_tr);
        plot([this_x this_x], [-0.2*y_shift y_shift*(ii_plot+1)],'-','Color',[0.6 0.6 0.6])
    end

    %Miss
    for ii_tr=1:length(per_ROI(this_ROI).results.all_trMiss)
        this_x=no_time_bins1+100+per_ROI(this_ROI).results.all_trMiss(ii_tr);
        plot([this_x this_x], [-0.2*y_shift y_shift*(ii_plot+1)],'-','Color',[0.6 0.6 0.6])
    end

    try
        ylim([-0.2*y_shift 1.2*y_shift*(ii_plot+1)])
    catch
    end


    %Entire arena
    subplot(2,3,4)
    max_this_dFF_activity=max(this_dFF_activity(:));
    max_this_dFF_miss_activity=max(this_dFF_miss_activity(:));
    max_this_dFF_hit_activity=max(this_dFF_hit_activity(:));
    max_activity=max([max_this_dFF_activity max_this_dFF_hit_activity max_this_dFF_miss_activity]);
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


    yticks(0:48:480)
    xticks(0:50:500)
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('All trials')

    %Hit
    subplot(2,3,5)
    this_masked_dFF_hit_activity=this_dFF_hit_activity;
    this_masked_dFF_hit_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
    % drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFF_hit_activity')
    drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFF_hit_activity)
    colormap(this_cmap)
    clim([-1.5*delta_ac max_activity])
    shading interp
    set(gca, 'YDir', 'reverse');

    yticks(0:48:480)
    xticks(0:50:500)
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Hit')

    %Miss
    subplot(2,3,6)

    this_masked_dFF_miss_activity=this_dFF_miss_activity;
    this_masked_dFF_miss_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
    drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFF_miss_activity)
    colormap(this_cmap)
    clim([-1.5*delta_ac max_activity])
    shading interp
    set(gca, 'YDir', 'reverse');

    yticks(0:48:480)
    xticks(0:50:500)
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Miss')


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

        sh_maps(ii_sh).this_dFF_hit_activity=zeros(10,10);
        sh_maps(ii_sh).this_dFF_hit_activity_n=zeros(10,10);
        sh_maps(ii_sh).sum_dFF_hit_activity=0;

        sh_maps(ii_sh).this_dFF_miss_activity=zeros(10,10);
        sh_maps(ii_sh).this_dFF_miss_activity_n=zeros(10,10);
        sh_maps(ii_sh).sum_dFF_miss_activity=0;

        these_all_dFF=(per_ROI(this_ROI).results.all_dFF)-min(per_ROI(this_ROI).results.all_dFF);
        these_all_hit_vs_miss=per_ROI(this_ROI).results.all_hit_vs_miss;

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

            if these_all_hit_vs_miss(ii_t)==1
                sh_maps(ii_sh).this_dFF_hit_activity(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFF_hit_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
                sh_maps(ii_sh).this_dFF_hit_activity_n(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFF_hit_activity_n(this_x_ii,this_y_ii)+1;
                sh_maps(ii_sh).sum_dFF_hit_activity=sum_dFF_hit_activity+these_all_dFF(ii_t);
            else
                sh_maps(ii_sh).this_dFF_miss_activity(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFF_miss_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
                sh_maps(ii_sh).this_dFF_miss_activity_n(this_x_ii,this_y_ii)=sh_maps(ii_sh).this_dFF_miss_activity_n(this_x_ii,this_y_ii)+1;
                sh_maps(ii_sh).sum_dFF_miss_activity=sum_dFF_miss_activity+these_all_dFF(ii_t);
            end
        end

        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)~=0
                    sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)=sh_maps(ii_sh).this_dFF_activity(ii_x,ii_y)/sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y);
                    if sh_maps(ii_sh).this_dFF_hit_activity_n(ii_x,ii_y)>0
                        sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)=sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)/sh_maps(ii_sh).this_dFF_hit_activity_n(ii_x,ii_y);
                    end
                    if sh_maps(ii_sh).this_dFF_miss_activity_n(ii_x,ii_y)>0
                        sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)=sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)/sh_maps(ii_sh).this_dFF_miss_activity_n(ii_x,ii_y);
                    end
                end
            end
        end


    end

    %Calculate information content
    %I will use information theory Markus et al 1994 https://doi.org/10.1002/hipo.450040404


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

    %Hit
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    this_information_content=0;
    this_sparsity=0;
    overall_mean_activity=mean(this_dFF_hit_activity(:));
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFF_hit_activity_n(ii_x,ii_y)/sum(this_dFF_hit_activity_n(:));
                if this_dFF_hit_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFF_hit_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end
    information_content_hit(this_ROI)=this_information_content;
    sparsity_hit(this_ROI)=this_sparsity;

    %Miss
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    this_information_content=0;
    this_sparsity=0;
    overall_mean_activity=mean(this_dFF_miss_activity(:));
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFF_miss_activity_n(ii_x,ii_y)/sum(this_dFF_miss_activity_n(:));
                if this_dFF_miss_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFF_miss_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end
    information_content_miss(this_ROI)=this_information_content;
    sparsity_miss(this_ROI)=this_sparsity;

    %Information content Hit Miss
    %Hit
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    this_information_content=0;
    this_sparsity=0;
    overall_mean_activity=mean([this_dFF_hit_activity(:); this_dFF_miss_activity(:)]);
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFF_hit_activity_n(ii_x,ii_y)/sum(this_dFF_activity_n(:));
                if this_dFF_hit_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFF_hit_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end


    %Miss
    %Please note I am doing this calculations in all space bins within the
    %mouse's trajectory
    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)>0
                this_prob=this_dFF_miss_activity_n(ii_x,ii_y)/sum(this_dFF_activity_n(:));
                if this_dFF_miss_activity(ii_x,ii_y)>0
                    this_information_content=this_information_content+this_prob*...
                        (this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity)*log2(this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity);
                end
                this_sparsity=this_sparsity+this_prob*...
                    (this_dFF_miss_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
            end
        end
    end
    information_content_hit_miss(this_ROI)=this_information_content;
    sparsity_miss(this_ROI)=this_sparsity;


    %Spatial correlation between hit and miss
    this_rho=corrcoef([this_dFF_hit_activity(this_dFF_activity_n>0) this_dFF_miss_activity(this_dFF_activity_n>0)]);
    spatial_rho_hit_vs_miss(this_ROI)=this_rho(1,2);

    %Center of mass
    cmx1=sum(this_dFF_hit_activity(this_dFF_activity_n>0).*X_square(this_dFF_activity_n>0))/sum(this_dFF_hit_activity(this_dFF_activity_n>0));
    cmx4=sum(this_dFF_miss_activity(this_dFF_activity_n>0).*X_square(this_dFF_activity_n>0))/sum(this_dFF_miss_activity(this_dFF_activity_n>0));

    cmy1=sum(this_dFF_hit_activity(this_dFF_activity_n>0).*Y_square(this_dFF_activity_n>0))/sum(this_dFF_hit_activity(this_dFF_activity_n>0));
    cmy4=sum(this_dFF_miss_activity(this_dFF_activity_n>0).*Y_square(this_dFF_activity_n>0))/sum(this_dFF_miss_activity(this_dFF_activity_n>0));

    delta_center_of_mass(this_ROI)=sqrt((cmx1-cmx4)^2 +(cmy1-cmy4)^2);

    %Now calculate the shuffled values
    these_sh_information_content=zeros(n_shuffle_SI,1);
    these_sh_sparsity=zeros(n_shuffle_SI,1);
    these_sh_information_content_hit=zeros(n_shuffle_SI,1);
    these_sh_sparsity_hit=zeros(n_shuffle_SI,1);
    these_sh_information_content_miss=zeros(n_shuffle_SI,1);
    these_sh_sparsity_miss=zeros(n_shuffle_SI,1);
    these_sh_spatial_rho_hit_vs_miss=zeros(n_shuffle_SI,1);
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

        %Hit
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory
        this_information_content=0;
        this_sparsity=0;
        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFF_hit_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end
        these_sh_information_content_hit(ii_sh,1)=this_information_content;
        these_sh_sparsity_hit(ii_sh,1)=this_sparsity;

        %Miss
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory
        this_information_content=0;
        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFF_miss_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end
        these_sh_information_content_miss(ii_sh,1)=this_information_content;
        these_sh_sparsity_miss(ii_sh,1)=this_sparsity;

        %Hit and Miss
        %Hit
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory
        this_information_content=0;
        this_sparsity=0;
        overall_mean_activity=mean([sh_maps(ii_sh).this_dFF_hit_activity(:); sh_maps(ii_sh).this_dFF_miss_activity(:)]);
        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFF_hit_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFF_hit_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end


        %Miss
        %Please note I am doing this calculations in all space bins within the
        %mouse's trajectory

        for ii_x=1:10
            for ii_y=1:10
                if sh_maps(ii_sh).this_dFF_activity_n(ii_x,ii_y)>0
                    this_prob=sh_maps(ii_sh).this_dFF_miss_activity_n(ii_x,ii_y)/sum(sh_maps(ii_sh).this_dFF_activity_n(:));
                    if sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)>0
                        this_information_content=this_information_content+this_prob*...
                            (sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity)*log2(sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)/overall_mean_activity);
                    end
                    this_sparsity=this_sparsity+this_prob*...
                        (sh_maps(ii_sh).this_dFF_miss_activity(ii_x,ii_y)^2)/(overall_mean_activity^2);
                end
            end
        end
        these_sh_information_content_miss(ii_sh,1)=this_information_content;
        these_sh_sparsity_miss(ii_sh,1)=this_sparsity;


        %Spatial correlation between hit and miss
        this_rho=corrcoef([sh_maps(ii_sh).this_dFF_hit_activity(sh_maps(ii_sh).this_dFF_activity_n>0) sh_maps(ii_sh).this_dFF_miss_activity(sh_maps(ii_sh).this_dFF_activity_n>0)]);
        these_sh_spatial_rho_hit_vs_miss(ii_sh)=this_rho(1,2);

        %Delta center of mass
        cmx1=sum(sh_maps(ii_sh).this_dFF_hit_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*X_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFF_hit_activity(sh_maps(ii_sh).this_dFF_activity_n>0));
        cmx4=sum(sh_maps(ii_sh).this_dFF_miss_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*X_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFF_miss_activity(this_dFF_activity_n>0));

        cmy1=sum(sh_maps(ii_sh).this_dFF_hit_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*Y_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFF_hit_activity(sh_maps(ii_sh).this_dFF_activity_n>0));
        cmy4=sum(sh_maps(ii_sh).this_dFF_miss_activity(sh_maps(ii_sh).this_dFF_activity_n>0).*Y_square(sh_maps(ii_sh).this_dFF_activity_n>0))/sum(sh_maps(ii_sh).this_dFF_miss_activity(sh_maps(ii_sh).this_dFF_activity_n>0));

        these_sh_delta_center_of_mass(ii_sh)=sqrt((cmx1-cmx4)^2 +(cmy1-cmy4)^2);
    end

    sh_information_content(this_ROI,:)=these_sh_information_content;
    sh_sparsity(this_ROI,:)=these_sh_sparsity;
    sh_information_content_hit(this_ROI,:)=these_sh_information_content_hit;
    sh_sparsity_hit(this_ROI,:)=these_sh_sparsity_hit;
    sh_information_content_miss(this_ROI,:)=these_sh_information_content_miss;
    sh_sparsity_miss(this_ROI,:)=these_sh_sparsity_miss;
    sh_spatial_rho_hit_vs_miss(this_ROI,:)=these_sh_spatial_rho_hit_vs_miss;
    sh_delta_center_of_mass(this_ROI,:)=these_sh_delta_center_of_mass;

    %Now calcualte information content for the binary dFF
    %Cumulative counts for info
    cum_xy=zeros(1,10*10);
    cum_hit_vs_miss=zeros(1,2);
    cum_bindFF=zeros(1,2);

    cum_xy_bindFF_hit=zeros(2,10*10);
    cum_xy_bindFF_miss=zeros(2,10*10);
    cum_xy_bindFF_both=zeros(2,10*10);

    cum_xy_bindFF_both=zeros(2,10*10,2);

    cum_xy_hit=zeros(1,10*10);
    cum_xy_miss=zeros(1,10*10);
    cum_xy_both=zeros(1,10*10);

    cum_bindFF_hit=zeros(1,2);
    cum_bindFF_miss=zeros(1,2);
    cum_bindFF_both=zeros(1,2);

    these_x=zeros(1,length(per_ROI(this_ROI).results.all_dFF));
    these_y=zeros(1,length(per_ROI(this_ROI).results.all_dFF));
    these_dFF=zeros(1,length(per_ROI(this_ROI).results.all_dFF));

    these_x(1,:)=per_ROI(this_ROI).results.all_XY(:,1);
    these_y(1,:)=per_ROI(this_ROI).results.all_XY(:,2);
    these_dFF(1,:)=per_ROI(this_ROI).results.all_dFF;

    for ii_t=1:length(per_ROI(this_ROI).results.all_dFF)
        this_x_ii=ceil(these_x(ii_t)/50);
        if this_x_ii==11
            this_x_ii=10;
        end

        this_y_ii=ceil(these_y(ii_t)/48);
        if this_y_ii==11
            this_y_ii=10;
        end

        % this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
        % this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
        % sum_dFF_activity=sum_dFF_activity+these_dFF(ii_t);

        %Tally info
        this_bin_dFF=(these_dFF(ii_t)>0)+1;
        cum_bindFF(this_bin_dFF)=cum_bindFF(this_bin_dFF)+1;
        cum_bindFF_both(this_bin_dFF)=cum_bindFF_both(this_bin_dFF)+1;
        this_xy_ii=this_x_ii+10*(this_y_ii-1);
        cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
        cum_xy_bindFF_both(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_both(this_bin_dFF,this_xy_ii)+1;
        cum_xy_both(this_xy_ii)=cum_xy_both(this_xy_ii)+1;

        if trials.lane_per_trial(trNo)==1
            this_dFF_hit_activity(this_x_ii,this_y_ii)=this_dFF_hit_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
            this_dFF_hit_activity_n(this_x_ii,this_y_ii)=this_dFF_hit_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFF_hit_activity=sum_dFF_hit_activity+these_dFF(ii_t);
            cum_xy_bindFF_hit(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_hit(this_bin_dFF,this_xy_ii)+1;
            cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,1)=cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,1)+1;
            cum_xy_hit(this_xy_ii)=cum_xy_hit(this_xy_ii)+1;
            cum_bindFF_hit(this_bin_dFF)=cum_bindFF_hit(this_bin_dFF)+1;
            cum_hit_vs_miss(1)=cum_hit_vs_miss(1)+1;
        else
            this_dFF_miss_activity(this_x_ii,this_y_ii)=this_dFF_miss_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
            this_dFF_miss_activity_n(this_x_ii,this_y_ii)=this_dFF_miss_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFF_miss_activity=sum_dFF_miss_activity+these_dFF(ii_t);
            cum_xy_bindFF_miss(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_miss(this_bin_dFF,this_xy_ii)+1;
            cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,2)=cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,2)+1;
            cum_xy_miss(this_xy_ii)=cum_xy_miss(this_xy_ii)+1;
            cum_bindFF_miss(this_bin_dFF)=cum_bindFF_miss(this_bin_dFF)+1;
            cum_hit_vs_miss(2)=cum_hit_vs_miss(2)+1;
        end

    end

    %Calculate information

    p_xy=cum_xy/sum(cum_xy(:));
    p_both=cum_hit_vs_miss/sum(cum_hit_vs_miss(:));
    p_bindFF=cum_bindFF/sum(cum_bindFF(:));

    p_xy_bindFF_hit=cum_xy_bindFF_hit/sum(cum_xy_bindFF_hit(:));
    p_xy_bindFF_miss=cum_xy_bindFF_miss/sum(cum_xy_bindFF_miss(:));
    p_xy_bindFF_both=cum_xy_bindFF_both/sum(cum_xy_bindFF_both(:));

    p_xy_bindFF_lane=cum_xy_bindFF_both/sum(cum_xy_bindFF_both(:));

    p_xy_Lane1=cum_xy_hit/sum(cum_xy_hit(:));
    p_xy_Lane4=cum_xy_miss/sum(cum_xy_miss(:));
    p_xy_BothLanes=cum_xy_both/sum(cum_xy_both(:));

    p_bindFF_hit=cum_bindFF_hit/sum(cum_bindFF_hit(:));
    p_bindFF_miss=cum_bindFF_miss/sum(cum_bindFF_miss(:));
    p_bindFF_both=cum_bindFF_both/sum(cum_bindFF_both(:));

    %Calculate info for Hit
    for ii_xy=1:10*10
        if sum(p_xy_bindFF_hit(:,ii_xy))>0
            for ii_bin_dFF=1:2
                if p_xy_bindFF_hit(ii_bin_dFF,ii_xy)~=0
                    dr_information_content_hit(this_ROI)=dr_information_content_hit(this_ROI)+p_xy_bindFF_hit(ii_bin_dFF,ii_xy)*...
                        log2(p_xy_bindFF_hit(ii_bin_dFF,ii_xy)/(p_xy_Lane1(ii_xy)*p_bindFF_hit(ii_bin_dFF)));
                end
            end
        end
    end

    %Calculate info for Miss
    for ii_xy=1:10*10
        if sum(p_xy_bindFF_miss(:,ii_xy))>0
            for ii_bin_dFF=1:2
                if p_xy_bindFF_miss(ii_bin_dFF,ii_xy)~=0
                    dr_information_content_miss(this_ROI)=dr_information_content_miss(this_ROI)+p_xy_bindFF_miss(ii_bin_dFF,ii_xy)*...
                        log2(p_xy_bindFF_miss(ii_bin_dFF,ii_xy)/(p_xy_Lane4(ii_xy)*p_bindFF_miss(ii_bin_dFF)));
                end
            end
        end
    end

    %Calculate mutual info between hit and miss
    for ii_xy=1:10*10
        if sum(p_xy_bindFF_miss(:,ii_xy))>0
            for ii_both=1:2
                for ii_bin_dFF=1:2
                    if p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_both)~=0
                        dr_information_contentl_mutual_1_miss(this_ROI)=dr_information_contentl_mutual_1_miss(this_ROI)+p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_both)*...
                            log2(p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_both)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_both(ii_both)));
                    end
                end
            end
        end
    end

    %Now calculate information contend with the shuffled dFF
    rng(this_ROI)
    these_rnd=rand(1,n_shuffle_SI);

    for ii_sh=1:n_shuffle_SI

        %We will do a reversal and a circular permutation of XY
        these_x_reversed=zeros(size(these_x,1),size(these_x,2));
        these_y_reversed=zeros(size(these_y,1),size(these_y,2));

        %Flip and roll Stefanini et al 2020 https://doi.org/10.1016/j.neuron.2020.05.022
        offset_ii=floor(these_rnd(ii_sh)*length(these_x));
        for ii_trl=1:length(these_x)
            this_ii_trl=ii_trl+offset_ii;
            if this_ii_trl>length(these_x)
                offset_ii=-ii_trl+1;
                this_ii_trl=ii_trl+offset_ii;
            end
            these_x_reversed(1,length(these_x)-ii_trl+1,:)=these_x(1,this_ii_trl);
            these_y_reversed(1,length(these_x)-ii_trl+1,:)=these_y(1,this_ii_trl);
        end

        %Cumulative counts for info
        cum_xy=zeros(1,10*10);
        cum_hit_vs_miss=zeros(1,2);
        cum_bindFF=zeros(1,2);

        cum_xy_bindFF_hit=zeros(2,10*10);
        cum_xy_bindFF_miss=zeros(2,10*10);
        cum_xy_bindFF_both=zeros(2,10*10);

        cum_xy_bindFF_both=zeros(2,10*10,2);

        cum_xy_hit=zeros(1,10*10);
        cum_xy_miss=zeros(1,10*10);
        cum_xy_both=zeros(1,10*10);

        cum_bindFF_hit=zeros(1,2);
        cum_bindFF_miss=zeros(1,2);
        cum_bindFF_both=zeros(1,2);


        for ii_t=1:length(per_ROI(this_ROI).results.all_dFF)
            this_x_ii=ceil(these_x_reversed(ii_t)/50);
            if this_x_ii==11
                this_x_ii=10;
            end

            this_y_ii=ceil(these_y_reversed(ii_t)/48);
            if this_y_ii==11
                this_y_ii=10;
            end

            % this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
            % this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
            % sum_dFF_activity=sum_dFF_activity+these_dFF(ii_t);

            %Tally info
            this_bin_dFF=(these_dFF(ii_t)>0)+1;
            cum_bindFF(this_bin_dFF)=cum_bindFF(this_bin_dFF)+1;
            cum_bindFF_both(this_bin_dFF)=cum_bindFF_both(this_bin_dFF)+1;
            this_xy_ii=this_x_ii+10*(this_y_ii-1);
            cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
            cum_xy_bindFF_both(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_both(this_bin_dFF,this_xy_ii)+1;
            cum_xy_both(this_xy_ii)=cum_xy_both(this_xy_ii)+1;

            if trials.lane_per_trial(trNo)==1
                this_dFF_hit_activity(this_x_ii,this_y_ii)=this_dFF_hit_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                this_dFF_hit_activity_n(this_x_ii,this_y_ii)=this_dFF_hit_activity_n(this_x_ii,this_y_ii)+1;
                sum_dFF_hit_activity=sum_dFF_hit_activity+these_dFF(ii_t);
                cum_xy_bindFF_hit(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_hit(this_bin_dFF,this_xy_ii)+1;
                cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,1)=cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,1)+1;
                cum_xy_hit(this_xy_ii)=cum_xy_hit(this_xy_ii)+1;
                cum_bindFF_hit(this_bin_dFF)=cum_bindFF_hit(this_bin_dFF)+1;
                cum_hit_vs_miss(1)=cum_hit_vs_miss(1)+1;
            else
                this_dFF_miss_activity(this_x_ii,this_y_ii)=this_dFF_miss_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                this_dFF_miss_activity_n(this_x_ii,this_y_ii)=this_dFF_miss_activity_n(this_x_ii,this_y_ii)+1;
                sum_dFF_miss_activity=sum_dFF_miss_activity+these_dFF(ii_t);
                cum_xy_bindFF_miss(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_miss(this_bin_dFF,this_xy_ii)+1;
                cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,2)=cum_xy_bindFF_both(this_bin_dFF,this_xy_ii,2)+1;
                cum_xy_miss(this_xy_ii)=cum_xy_miss(this_xy_ii)+1;
                cum_bindFF_miss(this_bin_dFF)=cum_bindFF_miss(this_bin_dFF)+1;
                cum_hit_vs_miss(2)=cum_hit_vs_miss(2)+1;
            end

        end

        %Calculate information

        p_xy=cum_xy/sum(cum_xy(:));
        p_both=cum_hit_vs_miss/sum(cum_hit_vs_miss(:));
        p_bindFF=cum_bindFF/sum(cum_bindFF(:));

        p_xy_bindFF_hit=cum_xy_bindFF_hit/sum(cum_xy_bindFF_hit(:));
        p_xy_bindFF_miss=cum_xy_bindFF_miss/sum(cum_xy_bindFF_miss(:));
        p_xy_bindFF_both=cum_xy_bindFF_both/sum(cum_xy_bindFF_both(:));

        p_xy_bindFF_lane=cum_xy_bindFF_both/sum(cum_xy_bindFF_both(:));

        p_xy_Lane1=cum_xy_hit/sum(cum_xy_hit(:));
        p_xy_Lane4=cum_xy_miss/sum(cum_xy_miss(:));
        p_xy_BothLanes=cum_xy_both/sum(cum_xy_both(:));

        p_bindFF_hit=cum_bindFF_hit/sum(cum_bindFF_hit(:));
        p_bindFF_miss=cum_bindFF_miss/sum(cum_bindFF_miss(:));
        p_bindFF_both=cum_bindFF_both/sum(cum_bindFF_both(:));

        %Calculate info for Hit
        for ii_xy=1:10*10
            if sum(p_xy_bindFF_hit(:,ii_xy))>0
                for ii_bin_dFF=1:2
                    if p_xy_bindFF_hit(ii_bin_dFF,ii_xy)~=0
                        sh_dr_information_content_hit(this_ROI,ii_sh)=sh_dr_information_content_hit(this_ROI,ii_sh)+p_xy_bindFF_hit(ii_bin_dFF,ii_xy)*...
                            log2(p_xy_bindFF_hit(ii_bin_dFF,ii_xy)/(p_xy_Lane1(ii_xy)*p_bindFF_hit(ii_bin_dFF)));
                    end
                end
            end
        end

        %Calculate info for Miss
        for ii_xy=1:10*10
            if sum(p_xy_bindFF_miss(:,ii_xy))>0
                for ii_bin_dFF=1:2
                    if p_xy_bindFF_miss(ii_bin_dFF,ii_xy)~=0
                        sh_dr_information_content_miss(this_ROI,ii_sh)=sh_dr_information_content_miss(this_ROI,ii_sh)+p_xy_bindFF_miss(ii_bin_dFF,ii_xy)*...
                            log2(p_xy_bindFF_miss(ii_bin_dFF,ii_xy)/(p_xy_Lane4(ii_xy)*p_bindFF_miss(ii_bin_dFF)));
                    end
                end
            end
        end

        %Calculate mutual info between hit and miss
        for ii_xy=1:10*10
            if sum(p_xy_bindFF_miss(:,ii_xy))>0
                for ii_both=1:2
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_both)~=0
                            sh_dr_information_contentl_mutual_1_miss(this_ROI,ii_sh)=sh_dr_information_contentl_mutual_1_miss(this_ROI,ii_sh)+p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_both)*...
                                log2(p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_both)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_both(ii_both)));
                        end
                    end
                end
            end
        end

    end
end



handles_out.information_content=information_content;
handles_out.sparsity=sparsity;
handles_out.information_content_hit=information_content_hit;
handles_out.sparsity_hit=sparsity_hit;
handles_out.information_content_miss=information_content_miss;
handles_out.sparsity_miss=sparsity_miss;
handles_out.spatial_rho_hit_vs_miss=spatial_rho_hit_vs_miss;
handles_out.delta_center_of_mass=delta_center_of_mass;
handles_out.sh_information_content=sh_information_content;
handles_out.sh_sparsity=sh_sparsity;
handles_out.sh_information_content_hit=sh_information_content_hit;
handles_out.sh_sparsity_hit=sh_sparsity_hit;
handles_out.sh_information_content_miss=sh_information_content_miss;
handles_out.sh_sparsity_miss=sh_sparsity_miss;
handles_out.sh_information_content_hit_miss=sh_information_content_hit_miss;
handles_out.sh_sparsity_miss=sh_sparsity_miss;
handles_out.sh_spatial_rho_hit_vs_miss=sh_spatial_rho_hit_vs_miss;
handles_out.sh_delta_center_of_mass=sh_delta_center_of_mass;
handles_out.dr_information_content_hit=dr_information_content_hit;
handles_out.dr_information_content_miss=dr_information_content_miss;
handles_out.dr_information_contentl_mutual_1_miss=dr_information_contentl_mutual_1_miss;
handles_out.sh_dr_information_content_hit=sh_dr_information_content_hit;
handles_out.sh_dr_information_content_miss=sh_dr_information_content_miss;
handles_out.sh_dr_information_contentl_mutual_1_miss=sh_dr_information_contentl_mutual_1_miss;

handles_out.glm_pvalues=glm_pvalues;
handles_out.glm__both_pvalues=glm__both_pvalues;
if ~isempty(glm_dFF_hit)
    handles_out.glm__hit_pvalues=glm__hit_pvalues;
end
if ~isempty(glm_dFF_miss)
    handles_out.glm__miss_pvalues=glm__miss_pvalues;
end

handles_out.trials=trials;

save([handles_choices.save_path arena_file(1:end-4) handles_choices.save_tag '.mat'],'handles_out','handles_choices','-v7.3')
   


pffft=1;