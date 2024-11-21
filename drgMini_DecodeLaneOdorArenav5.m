function handles_out=drgMini_DecodeLaneOdorArenav5(handles_choices2)
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020
close all

if exist('handles_choices2')==0
    clear all


    %First troubleshooting files
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';
    pred_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm_deczdFFopt2_711103.mat';

    %Second troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220713_FCM6/';
    % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';
    % pred_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm_deczdFFopt2_711103.mat';

    %No odor troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';


    handles_choices2.this_path=this_path;
    handles_choices2.dFF_file=dFF_file;
    handles_choices2.arena_file=arena_file;

    %     isKording=0;

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    %     dt=0.2;

    which_ROIs=1; %1 Use all ROIs, 2 Use place cells, 3 use odor cells, 4 use ROIs in
    %handles_choices2.process_these_ROIs
    handles_choices2.which_ROIs=which_ROIs;

    handles_choices2.process_these_ROIs=8; %This is only used if which_ROIs=4

    no_repeats=15;
    handles_choices.no_repeats=no_repeats; %Number of times to repeat fitrnet

    %Define the different ranges (training, valid and testing)
    % training_fraction=0.9;
    % handles_choices2.training_fraction=training_fraction;

    %     training_range=[0, 0.5];
    %     valid_range=[0.5,0.65];
    %     test_range=[0.5, 1];

    %The user can define what time period to use spikes from (with respect to the output).
    bins_before=5; %How many bins of neural data prior to the output are used for decoding, 5
    bins_current=1; %Whether to use concurrent time bin of neural data, 1
    bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    handles_choices2.bins_before=bins_before;
    handles_choices2.bins_current=bins_current;
    handles_choices2.bins_after=bins_after;


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2;
    dt_miniscope=1/30;
    n_shuffle=10; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_choices2.dt=dt;
    handles_choices2.dt_miniscope=dt_miniscope;
    handles_choices2.n_shuffle=n_shuffle;

    % which_training_algorithm=2;
    % handles_choices2.which_training_algorithm=which_training_algorithm;
    %1=entire session is used for training
    %2=only per trial is used for training
    %3=trained with data between trials

    which_ml_algo=7;
    %1 SVZ
    %2 nn
    %3 tree
    %4 bayesian
    %5 glm
    %6 linear
    %7 NN parameter optimized

    % algo_name{1}='SVZ';
    % algo_name{2}='NN';
    % algo_name{3}='Tree';
    % algo_name{4}='Bayes';
    % algo_name{5}='glm';
    % algo_name{6}='LD';
    % algo_name{7}='NNopt';


    ii_cost=3;

    this_cost=[0 ii_cost;ii_cost 0];

    %start and end of training period
    align_training_start=0; %0= aligned to end, 1=aligned to start
    % align_training_end=1; %0= odor start, 1=odor end
    dt_training_start=-2; %seconds from start alignment -5
    dt_training_end=5; %seconds from end alignment 5
    ii_dt_training_start=fix(dt_training_start/dt); %samples from start alignment
    ii_dt_training_end=fix(dt_training_end/dt); %samples from end alignment

    %start and end of display period
    % align_display=0; %0= odor start, 1=odor end
    dt_display_start=-5; %seconds from start alignment
    dt_display_end=5; %seconds from end alignment (or mean end)
    ii_dt_display_start=fix(dt_display_start/dt); %samples from start alignment
    ii_dt_display_end=fix(dt_display_end/dt); %samples from end alignment

else

    this_path=handles_choices2.this_path;
    dFF_file=handles_choices2.dFF_file;
    arena_file=handles_choices2.arena_file;
    pred_file=handles_choices2.pred_file;

    which_ROIs=handles_choices2.which_ROIs; %when this function is called the user has to specify handles_choices2.process_these_ROIs
    bins_before=handles_choices2.bins_before;
    bins_current=handles_choices2.bins_current;
    bins_after=handles_choices2.bins_after;
    no_repeats=handles_choices2.no_repeats;
    dt=handles_choices2.dt;
    dt_miniscope=handles_choices2.dt_miniscope;
    n_shuffle=handles_choices2.n_shuffle;
    which_ml_algo=handles_choices2.which_ml_algo;
    align_training_start=handles_choices2.align_training_start;
    dt_training_start=handles_choices2.dt_training_start; %seconds from start alignment -5
    dt_training_end=handles_choices2.dt_training_end; %seconds from end alignment 5


    ii_dt_training_start=fix(dt_training_start/dt); %samples from start alignment
    ii_dt_training_end=fix(dt_training_end/dt); %samples from end alignment

    %start and end of display period
    % align_display=0; %0= odor start, 1=odor end
    dt_display_start=-5; %seconds from start alignment
    dt_display_end=5; %seconds from end alignment (or mean end)
    ii_dt_display_start=fix(dt_display_start/dt); %samples from start alignment
    ii_dt_display_end=fix(dt_display_end/dt); %samples from end alignment
end

 %1 SVZ
    %2 nn
    %3 tree
    %4 bayesian
    %5 glm
    %6 linear
    %7 NN parameter optimized

    algo_name{1}='SVZ';
    algo_name{2}='NN';
    algo_name{3}='Tree';
    algo_name{4}='Bayes';
    algo_name{5}='glm';
    algo_name{6}='LD';
    algo_name{7}='NNopt';

if isempty(gcp('nocreate'))
    parpool;  % This will use the default number of workers
else
    delete(gcp('nocreate'));
    parpool;
end

%Load prediction file and define place and odor cells
load([this_path pred_file])

spatial_rhol1l4=handles_out.spatial_rhol1l4;
delta_center_of_mass=handles_out.delta_center_of_mass;

no_neurons=length(spatial_rhol1l4);

switch which_ROIs
    case 1
        process_these_ROIs=[1:no_neurons];
    case {2,3}
        p_value_odor=[];
        p_value_lanexy_trial=[];
        p_value_xy=[];
        p_value_xyl1=[];
        p_value_xyl4=[];
        for ii_ROI=1:no_neurons
            p_value_odor=[p_value_odor handles_out.glm_l14_pvalues.ROI(ii_ROI).pValues(2)];
            p_value_lanexy_trial=[p_value_lanexy_trial handles_out.glm_pvalues.ROI(ii_ROI).pValues(3)];
            p_value_xy=[p_value_xy handles_out.glm_pvalues.ROI(ii_ROI).pValues(2)];
            p_value_xyl1=[p_value_xy handles_out.glm_l1_pvalues.ROI(ii_ROI).pValues(2)];
            p_value_xyl4=[p_value_xy handles_out.glm_l4_pvalues.ROI(ii_ROI).pValues(2)];
        end


        no_place_cells=0;
        place_cells=[];
        no_odor_cells=0;
        odor_cells=[];

        for this_ROI=1:no_neurons
            %Now show the place cells
            if (abs(spatial_rhol1l4(this_ROI))>=0.4)&(delta_center_of_mass(this_ROI)<=100)...
                    &(p_value_xy(ii_ROI)<drsFDRpval(p_value_xy))&(p_value_xyl1(ii_ROI)<drsFDRpval(p_value_xyl1))&(p_value_xyl4(ii_ROI)<drsFDRpval(p_value_xyl4))
                no_place_cells=no_place_cells+1;
                place_cells=[place_cells this_ROI];
                pffft=1;
            end

            %Now show the odor cells
            if (abs(spatial_rhol1l4(this_ROI))<=0.1)&(delta_center_of_mass(this_ROI)>=300)...
                    &(p_value_odor(ii_ROI)<drsFDRpval(p_value_odor))
                no_odor_cells=no_odor_cells+1;
                odor_cells=[odor_cells this_ROI];
                pffft=1;
            end
        end
        switch which_ROIs
            case 2
                process_these_ROIs=place_cells;
            case 3
                process_these_ROIs=odor_cells;
        end
    case 4
        process_these_ROIs=handles_choices2.process_these_ROIs;

end

% switch which_training_algorithm
%     case 1
%         fprintf(1,['\nTrained with data for the entire session\n\n'])
%     case 2
fprintf(1,['\nTrained with within trial data\n\n'])
%     case 3
%         fprintf(1,['\nTrained with between trial data\n\n'])
% end

figNo=0;

%Restart random seeds
rng('shuffle');

% try
%     delete(gcp('nocreate'))
% catch
% end
% % gcp;
% parpool;

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

dFF=dFF(:,process_these_ROIs);
% dFF=readmatrix([this_path dFF_file]); %Timepoints x ROIs

load([this_path arena_file])

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
%         %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
%         if ~isempty(next_ii_flir2)
%             ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
%         else
%             at_end=1;
%         end
%     else
%         at_end=1;
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
%         %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
%         if ~isempty(next_ii_flir2)
%             ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
%         else
%             at_end=1;
%         end
%     else
%         at_end=1;
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
plot([95 95],[200 250],'-r')
plot([95 95],[0 50],'-b')

plot(trials.x_lanewater1,trials.y_lanewater1,'xr')
plot(trials.x_lanewater4,trials.y_lanewater4,'xb')

set(gca, 'YDir', 'reverse');
xlabel('x')
ylabel('y')
title('Location of trial start (o) and water delivery (x)')


%Bin positions into dt time bins
pos=[];
pos(:,1)=arena.xsync;
pos(:,2)=arena.ysync;
no_time_points=size(pos,1);


dFF_times=[1:no_time_points]*dt_miniscope;

no_neurons=size(dFF,2);
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
trials.hit4=0;
trials.miss4=0;
trials.odor_trNo=0;

trim_factor=no_time_bins/no_time_points;

dii_trial=[];
all_lane_identity=[];
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
            trials.lane_identity(trials.odor_trNo)=0;
        else
            trials.miss1=trials.miss1+1;
            trials.miss1_ii_start(trials.miss1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss1_ii_start(trials.miss1);
            trials.odor_trial_type(trials.odor_trNo)=2;
            trials.lane_identity(trials.odor_trNo)=0;
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
            trials.lane_identity(trials.odor_trNo)=1;
        else
            trials.miss4=trials.miss4+1;
            trials.miss4_ii_start(trials.miss4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss4_ii_start(trials.miss4);
            trials.odor_trial_type(trials.odor_trNo)=4;
            trials.lane_identity(trials.odor_trNo)=1;
        end
    end

end

%Now generate shuffled lane identity
ii_shuffled=0;
these_perms=[];
for ii_sh=1:n_shuffle
    got_perm=0;
    while got_perm==0
        this_perm=randperm(length(trials.lane_identity));
        if ii_sh==1
            got_perm=1;
        else
            got_perm=1;
            for ii_p=1:size(these_perms,1)
                this_these_perms=zeros(1,length(trials.lane_identity));
                this_these_perms(1,:)=these_perms(ii_p,:);
                if all(this_these_perms==this_perm)
                    got_perm=0;
                end
            end
        end
    end
    these_perms=[these_perms; this_perm];
    for ii_tr=1:length(this_perm)
        trials.shuffled(ii_sh).lane_identity(ii_tr)=trials.lane_identity(this_perm(ii_tr));
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

%Report the average time period from start to end
delta_t_per_trial=zeros(1,trials.odor_trNo);
delta_t_per_trial=trials.odor_ii_end-trials.odor_ii_start;
fprintf(1,['\nMean trial period (sec) ' num2str(mean(delta_t_per_trial)*dt) '\n\n'])


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

% neural_data_trimmed=(neural_data_trimmed-mean_neural_data_trimmed)./std_neural_data_trimmed;
neural_data_trimmed=(neural_data_trimmed)./std_neural_data_trimmed;

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

%Now do decoding
tic
%Set what part of data should be part of the training/testing/validation sets
%Note that there was a long period of no movement after about 80% of recording, so I did not use this data.
% switch which_training_algorithm
%     case 1
%         %Train with entire session
%         ii_n_training=floor(1/(1-training_fraction));
%
%         for ii_train=1:ii_n_training
%             y_pred(ii_train).data=[];
%             label_pred(ii_train).data=[];
%         end
%
%         parfor ii_train=1:ii_n_training
%             this_test_range=zeros(1,no_time_bins);
%             ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
%             ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
%             if ii_train==ii_n_training
%                 ii_test_range_end=no_time_bins;
%             end
%             this_test_range(ii_test_range_start:ii_test_range_end)=1;
%
%             % ii_valid_range=ceil(valid_range*no_time_bins);
%             %     ii_test_range=ceil(test_range*no_time_bins);
%
%             XdFFtrain=X_dFF(~logical(this_test_range),:);
%             % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
%             XdFFtest=X_dFF(logical(this_test_range),:);
%
%             XYtrain=pos_binned_trimmed(~logical(this_test_range),:);
%             % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
%             %     XYtest=pos_binned_trimmed(logical(this_test_range),:);
%
%             %Decode using neural network
%             MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
%             MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));
%
%             %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);
%
%             %     label_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
%             %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);
%
%             label_pred(ii_train).data=predict(MdlY1,XdFFtest);
%             y_pred(ii_train).data=predict(MdlY2,XdFFtest);
%         end
%         fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])
%
%         %Parse out the parfor loop output
%         label_predicted=zeros(no_time_bins,1);
%         y_predicted=zeros(no_time_bins,1);
%         for ii_train=1:ii_n_training
%             this_test_range=zeros(1,no_time_bins);
%             ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
%             ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
%             if ii_train==ii_n_training
%                 ii_test_range_end=no_time_bins;
%             end
%             this_test_range(ii_test_range_start:ii_test_range_end)=1;
%
%             label_predicted(logical(this_test_range),1)=label_pred(ii_train).data;
%             y_predicted(logical(this_test_range),1)=y_pred(ii_train).data;
%
%         end
%
%         %Now do predictions for reversed/permuted training periods
%         if n_shuffle>ii_n_training
%             n_shuffle=ii_n_training;
%         end
%         handles_choices2.n_shuffle=n_shuffle;
%         label_predicted_sh=zeros(no_time_bins,n_shuffle);
%         y_predicted_sh=zeros(no_time_bins,n_shuffle);
%
%
%         %We will do a reversal and a circular permutation
%         sh_shift=0;
%         while sh_shift==0
%             sh_shift=floor(rand*n_shuffle);
%         end
%
%         pos_binned_reversed=zeros(size(pos_binned_trimmed,1),size(pos_binned_trimmed,2));
%         for ii_trl=1:size(pos_binned_trimmed,1)
%             pos_binned_reversed(size(pos_binned_trimmed,1)-ii_trl+1,:)=pos_binned_trimmed(ii_trl,:);
%         end
%
%         for ii_shuffled=1:n_shuffle
%
%             for ii_train=1:ii_n_training
%                 y_pred(ii_train).data=[];
%                 label_pred(ii_train).data=[];
%             end
%
%             parfor ii_train=1:ii_n_training
%
%                 this_test_range=zeros(1,no_time_bins);
%                 ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
%                 ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
%                 if ii_train==ii_n_training
%                     ii_test_range_end=no_time_bins;
%                 end
%                 this_test_range(ii_test_range_start:ii_test_range_end)=1;
%
%                 ii_train_sh=ii_train+sh_shift;
%                 if ii_train_sh>ii_n_training
%                     ii_train_sh=ii_train_sh-ii_n_training;
%                 end
%                 this_test_range_sh=zeros(1,no_time_bins);
%                 ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train_sh-1));
%                 ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train_sh);
%                 if ii_train_sh==ii_n_training
%                     ii_test_range_end=no_time_bins;
%                 end
%                 this_test_range_sh(ii_test_range_start:ii_test_range_end)=1;
%
%                 %Make sure that this_test_range_sh is the same size as this_test_range
%                 if sum(this_test_range_sh)>sum(this_test_range)
%                     first_ii=find(this_test_range_sh==1,1,'first');
%                     this_test_range_sh(first_ii:first_ii+sum(this_test_range_sh)-sum(this_test_range)-1)=0;
%                 end
%
%                 if sum(this_test_range_sh)<sum(this_test_range)
%                     first_ii=find(this_test_range_sh==1,1,'first');
%                     if first_ii+(sum(this_test_range_sh)-sum(this_test_range))>0
%                         this_test_range_sh(first_ii+(sum(this_test_range_sh)-sum(this_test_range)):first_ii-1)=1;
%                     else
%                         last_ii=find(this_test_range_sh==1,1,'last');
%                         this_test_range_sh(last_ii+1:last_ii+(sum(this_test_range)-sum(this_test_range_sh)))=1;
%                     end
%                 end
%
%                 % ii_valid_range=ceil(valid_range*no_time_bins);
%                 %     ii_test_range=ceil(test_range*no_time_bins);
%
%                 XdFFtrain=X_dFF(~logical(this_test_range),:);
%
%                 XdFFtest=X_dFF(logical(this_test_range),:);
%
%
%                 XYtrain=pos_binned_reversed(~logical(this_test_range_sh),:);
%
%                 %Decode using neural network
%                 MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
%                 MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));
%
%                 label_pred(ii_train).data=predict(MdlY1,XdFFtest);
%                 y_pred(ii_train).data=predict(MdlY2,XdFFtest);
%             end
%
%             for ii_train=1:ii_n_training
%                 this_test_range=zeros(1,no_time_bins);
%                 ii_test_range_start=floor(1+((1-training_fraction)*no_time_bins)*(ii_train-1));
%                 ii_test_range_end=ceil(((1-training_fraction)*no_time_bins)*ii_train);
%                 if ii_train==ii_n_training
%                     ii_test_range_end=no_time_bins;
%                 end
%                 this_test_range(ii_test_range_start:ii_test_range_end)=1;
%
%                 label_predicted_sh(logical(this_test_range),ii_shuffled)=label_pred(ii_train).data;
%                 y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(ii_train).data;
%             end
%
%         end
%
%     case 2
%Train with trials using a leave one out approach


%training_range_template has all the trials
training_range_template=zeros(1,no_time_bins);
lane_labels=0.5*ones(1,no_time_bins);
for trNo=1:trials.odor_trNo
    y_pred(trNo).data=[];
    label_pred(trNo).data=[];

    switch align_training_start
        case 0
            %Training period aligned to end
            label_predictedstart=trials.odor_ii_end(trNo)+ii_dt_training_start;
            label_predictedend=trials.odor_ii_end(trNo)+ii_dt_training_end;
        case 1
            %Training period aligned to start
            label_predictedstart=trials.odor_ii_start(trNo)+ii_dt_training_start;
            label_predictedend=trials.odor_ii_start(trNo)+ii_dt_training_end;
    end

    display_start=trials.odor_ii_end(trNo)+ii_dt_training_start;
    display_end=trials.odor_ii_end(trNo)+ii_dt_training_end;

    training_range_template(label_predictedstart:label_predictedend)=1;
    lane_labels(label_predictedstart:label_predictedend)=trials.lane_identity(trNo);
end

for trNo=1:trials.odor_trNo
    % parfor trNo=1:trials.odor_trNo

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
    tblTrn=[];
    tblTrn = array2table(XdFFtrain);
    Y=lane_labels(this_training_range);

    switch which_ml_algo
        case 1
            Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
        case 2
            Mdl = fitcnet(tblTrn,Y);
        case 3
            Mdl = fitctree(tblTrn,Y);
        case 4
            Mdl=fitcnb(tblTrn,Y,'Cost',this_cost);
        case 5
            Mdl = fitglm(XdFFtrain,Y','Distribution','binomial');
        case 6
            Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
        case 7
            opts = struct('ShowPlots', false, ...
                'Verbose', 0, ...
                'MaxObjectiveEvaluations', 15,...
                'UseParallel',true);
            Mdl = fitcnet(tblTrn,Y,'OptimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions', opts);
    end

    label_pred(trNo).Mdl=Mdl;
    %
    %             MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
    %             MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));

    %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);

    %     label_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
    %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

    label_pred(trNo).data=predict(Mdl,XdFFtest);
    %             y_pred(trNo).data=predict(MdlY2,XdFFtest);
    pfffft=1;

end
fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

%Parse out the parfor loop output
label_predicted=zeros(no_time_bins,no_repeats);
accuracy=zeros(no_time_bins,no_repeats);
%         y_predicted=zeros(no_time_bins,1);
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

    label_predicted(logical(this_test_range),1)=label_pred(trNo).data;
    %             y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

    %Compute accuracy
    this_accuracy=zeros(size(label_pred(trNo).data,1),size(label_pred(trNo).data,2));
    if trials.lane_identity(trNo)==1
        these_times=logical(label_pred(trNo).data==1);
        this_accuracy(these_times)=1;
        this_accuracy(~these_times)=0;
    else
        these_times=logical(label_pred(trNo).data==1);
        this_accuracy(these_times)=0;
        this_accuracy(~these_times)=1;
    end
    accuracy(logical(this_test_range),1)=this_accuracy;
end

%Now do the rest of the repeats
for ii_repeats=2:no_repeats
    ii_repeats

    label_pred_r=[];
    for trNo=1:trials.odor_trNo
        % parfor trNo=1:trials.odor_trNo

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
        tblTrn=[];
        tblTrn = array2table(XdFFtrain);
        Y=lane_labels(this_training_range);

        switch which_ml_algo
            case 1
                Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
            case 2
                Mdl = fitcnet(tblTrn,Y);
            case 3
                Mdl = fitctree(tblTrn,Y);
            case 4
                Mdl=fitcnb(tblTrn,Y,'Cost',this_cost);
            case 5
                Mdl = fitglm(XdFFtrain,Y','Distribution','binomial');
            case 6
                Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
            case 7


                bestHyperparameters = label_pred(trNo).Mdl.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
                activationsCell = cellstr(bestHyperparameters.Activations);
                standardizeCell = cellstr(bestHyperparameters.Standardize);
                layer_sizes=bestHyperparameters.Layer_1_Size;
                if ~isnan(bestHyperparameters.Layer_2_Size)
                    layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
                end
                if ~isnan(bestHyperparameters.Layer_3_Size)
                    layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
                end
                Mdl = fitcnet(tblTrn,Y,'LayerSizes', layer_sizes, ...
                    'Activations', activationsCell{1}, ...
                    'Lambda', bestHyperparameters.Lambda, ...
                    'Standardize', strcmpi(standardizeCell{1},'true'),...
                    'ClassNames', unique(Y),...
                    'Prior', 'uniform');

                handles_out.decode.Mdl_pars(trNo).pars.activations=activationsCell{1};
                handles_out.decode.Mdl_pars(trNo).pars.Lambda=bestHyperparameters.Lambda;
                handles_out.decode.Mdl_pars(trNo).pars.Standardize=standardizeCell{1};

                % opts = struct('ShowPlots', false, ...
                %     'Verbose', 0, ...
                %     'MaxObjectiveEvaluations', 15,...
                %     'UseParallel',true);
                % Mdl = fitcnet(tblTrn,Y,'OptimizeHyperparameters','auto',...
                %     'HyperparameterOptimizationOptions', opts);
        end

        %
        %             MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
        %             MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));

        %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);

        %     label_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
        %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

        label_predr(trNo).data=predict(Mdl,XdFFtest);
        %             y_pred(trNo).data=predict(MdlY2,XdFFtest);
        pfffft=1;

    end
    fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

   
    %         y_predicted=zeros(no_time_bins,1);
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

        label_predicted(logical(this_test_range),ii_repeats)=label_predr(trNo).data;
        %             y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

        %Compute accuracy
        this_accuracy=zeros(size(label_predr(trNo).data,1),size(label_predr(trNo).data,2));
        if trials.lane_identity(trNo)==1
            these_times=logical(label_predr(trNo).data==1);
            this_accuracy(these_times)=1;
            this_accuracy(~these_times)=0;
        else
            these_times=logical(label_predr(trNo).data==1);
            this_accuracy(these_times)=0;
            this_accuracy(~these_times)=1;
        end
        accuracy(logical(this_test_range),ii_repeats)=this_accuracy;
    end


end

%Now do shuffled decoding
label_predicted_sh=zeros(no_time_bins,n_shuffle);
accuracy_sh=zeros(no_time_bins,n_shuffle);

for ii_shuffled=1:n_shuffle

    for trNo=1:trials.odor_trNo
        y_pred(trNo).data=[];
        label_pred(trNo).data=[];
    end

    %get the lane_labels for this shuffling run
    sh_lane_labels=zeros(1,no_time_bins);
    for trNo=1:trials.odor_trNo
        switch align_training_start
            case 0
                %Training period aligned to end
                label_predictedstart=trials.odor_ii_end(trNo)+ii_dt_training_start;
                label_predictedend=trials.odor_ii_end(trNo)+ii_dt_training_end;
            case 1
                %Training period aligned to start
                label_predictedstart=trials.odor_ii_start(trNo)+ii_dt_training_start;
                label_predictedend=trials.odor_ii_start(trNo)+ii_dt_training_end;
        end
        sh_lane_labels(label_predictedstart:label_predictedend)=trials.shuffled(ii_shuffled).lane_identity(trNo);
    end

    % parfor trNo=1:trials.odor_trNo
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
        tblTrn=[];
        tblTrn = array2table(XdFFtrain);
        Y=sh_lane_labels(this_training_range);

        switch which_ml_algo
            case 1
                Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
            case 2
                Mdl = fitcnet(tblTrn,Y);
            case 3
                Mdl = fitctree(tblTrn,Y);
            case 4
                Mdl=fitcnb(tblTrn,Y,'Cost',this_cost);
            case 5
                Mdl = fitglm(XdFFtrain,Y','Distribution','binomial');
            case 6
                Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
            case 7
                %Note that we do not optimize the shuffled runs

                  bestHyperparameters = label_pred(trNo).Mdl.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
                activationsCell = cellstr(bestHyperparameters.Activations);
                standardizeCell = cellstr(bestHyperparameters.Standardize);
                layer_sizes=bestHyperparameters.Layer_1_Size;
                if ~isnan(bestHyperparameters.Layer_2_Size)
                    layer_sizes=[layer_sizes bestHyperparameters.Layer_2_Size];
                end
                if ~isnan(bestHyperparameters.Layer_3_Size)
                    layer_sizes=[layer_sizes bestHyperparameters.Layer_3_Size];
                end
                Mdl = fitcnet(tblTrn,Y,'LayerSizes', layer_sizes, ...
                    'Activations', activationsCell{1}, ...
                    'Lambda', bestHyperparameters.Lambda, ...
                    'Standardize', strcmpi(standardizeCell{1},'true'),...
                    'ClassNames', unique(Y),...
                    'Prior', 'uniform');

                % Mdl = fitcnet(tblTrn,Y);
        end
        %
        %             MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
        %             MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));

        %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);

        %     label_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
        %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);

        label_pred(trNo).data=predict(Mdl,XdFFtest);
        %             y_pred(trNo).data=predict(MdlY2,XdFFtest);


        %
        %                 this_test_range=zeros(1,no_time_bins);
        %                 if trNo==1
        %                     ii_test_range_start=1;
        %                 else
        %                     ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
        %                 end
        %
        %                 if trNo==trials.odor_trNo
        %                     ii_test_range_end=no_time_bins;
        %                 else
        %                     ii_test_range_end=trials.odor_ii_end(trNo)+15;
        %                 end
        %
        %                 this_test_range(ii_test_range_start:ii_test_range_end)=1;
        %                 this_training_range=logical(training_range_template)&(~logical(this_test_range));
        %
        %                 % ii_valid_range=ceil(valid_range*no_time_bins);
        %                 %     ii_test_range=ceil(test_range*no_time_bins);
        %
        %                 XdFFtrain=X_dFF(this_training_range,:);
        %                 % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
        %                 XdFFtest=X_dFF(logical(this_test_range),:);
        %
        %                 XYtrain=pos_binned_reversed(this_training_range,:);
        %
        %
        %
        %                 %Decode using neural network
        %                 MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
        %                 MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));
        %
        %                 label_pred(trNo).data=predict(MdlY1,XdFFtest);
        %                 y_pred(trNo).data=predict(MdlY2,XdFFtest);
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

        label_predicted_sh(logical(this_test_range),ii_shuffled)=label_pred(trNo).data;
        %Compute accuracy
        this_accuracy=zeros(size(label_pred(trNo).data,1),size(label_pred(trNo).data,2));
        if trials.lane_identity(trNo)==1
            these_times=logical(label_pred(trNo).data==1);
            this_accuracy(these_times)=1;
            this_accuracy(~these_times)=0;
        else
            these_times=logical(label_pred(trNo).data==1);
            this_accuracy(these_times)=0;
            this_accuracy(~these_times)=1;
        end
        accuracy_sh(logical(this_test_range),ii_shuffled)=this_accuracy;
    end

end
%
%         case 3
%         %Train between trials using a leave one out approach
%
%         %training_range_template has all the trials
%         training_range_template=zeros(1,no_time_bins);
%         for trNo=1:trials.odor_trNo
%             y_pred(trNo).data=[];
%             label_pred(trNo).data=[];
%
%             if trNo==1
%                last_label_predictedend=1;
%             else
%             last_label_predictedend=trials.odor_ii_end(trNo-1)+15;
%             end
%
%
%             this_label_predictedstart=trials.odor_ii_start(trNo)-10;
%
%
%             training_range_template(last_label_predictedend:this_label_predictedstart)=1;
%         end
%
%         last_label_predictedend=trials.odor_ii_end(trials.odor_trNo)+15;
%         this_label_predictedstart=no_time_bins;
%         training_range_template(last_label_predictedend:this_label_predictedstart)=1;
%
%         parfor trNo=1:trials.odor_trNo
%
%             this_test_range=zeros(1,no_time_bins);
%             if trNo==1
%                 ii_test_range_start=1;
%             else
%                 ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
%             end
%
%             if trNo==trials.odor_trNo
%                 ii_test_range_end=no_time_bins;
%             else
%                 ii_test_range_end=trials.odor_ii_end(trNo)+15;
%             end
%
%             this_test_range(ii_test_range_start:ii_test_range_end)=1;
%             this_training_range=logical(training_range_template)&(~logical(this_test_range));
%
%             % ii_valid_range=ceil(valid_range*no_time_bins);
%             %     ii_test_range=ceil(test_range*no_time_bins);
%
%             XdFFtrain=X_dFF(this_training_range,:);
%             % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
%             XdFFtest=X_dFF(logical(this_test_range),:);
%
%             XYtrain=pos_binned_trimmed(this_training_range,:);
%             % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
%             %     XYtest=pos_binned_trimmed(logical(this_test_range),:);
%
%             %Decode using neural network
%             MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
%             MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));
%
%             %     MdlY1and2 = fitrnet(XdFFtrain,XYtrain','LayerSizes',[2 10]);
%
%             %     label_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
%             %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);
%
%             label_pred(trNo).data=predict(MdlY1,XdFFtest);
%             y_pred(trNo).data=predict(MdlY2,XdFFtest);
%         end
%         fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])
%
%         %Parse out the parfor loop output
%         label_predicted=zeros(no_time_bins,1);
%         y_predicted=zeros(no_time_bins,1);
%         for trNo=1:trials.odor_trNo
%             this_test_range=zeros(1,no_time_bins);
%
%             if trNo==1
%                 ii_test_range_start=1;
%             else
%                 ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
%             end
%
%             if trNo==trials.odor_trNo
%                 ii_test_range_end=no_time_bins;
%             else
%                 ii_test_range_end=trials.odor_ii_end(trNo)+15;
%             end
%
%             this_test_range(ii_test_range_start:ii_test_range_end)=1;
%
%             label_predicted(logical(this_test_range),1)=label_pred(trNo).data;
%             y_predicted(logical(this_test_range),1)=y_pred(trNo).data;
%
%         end
%
%         %Now do predictions for reversed/permuted training periods
%         label_predicted_sh=zeros(no_time_bins,n_shuffle);
%         y_predicted_sh=zeros(no_time_bins,n_shuffle);
%
%
%         %We will do a reversal and a circular permutation
%         sh_shift=0;
%         while sh_shift==0
%             sh_shift=floor(rand*n_shuffle);
%         end
%
%         pos_binned_reversed=zeros(size(pos_binned_trimmed,1),size(pos_binned_trimmed,2));
%         for ii_trl=1:size(pos_binned_trimmed,1)
%             pos_binned_reversed(size(pos_binned_trimmed,1)-ii_trl+1,:)=pos_binned_trimmed(ii_trl,:);
%         end
%
%         for ii_shuffled=1:n_shuffle
%
%             for trNo=1:trials.odor_trNo
%                 y_pred(trNo).data=[];
%                 label_pred(trNo).data=[];
%             end
%
%             parfor trNo=1:trials.odor_trNo
%
%                 this_test_range=zeros(1,no_time_bins);
%                 if trNo==1
%                     ii_test_range_start=1;
%                 else
%                     ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
%                 end
%
%                 if trNo==trials.odor_trNo
%                     ii_test_range_end=no_time_bins;
%                 else
%                     ii_test_range_end=trials.odor_ii_end(trNo)+15;
%                 end
%
%                 this_test_range(ii_test_range_start:ii_test_range_end)=1;
%                 this_training_range=logical(training_range_template)&(~logical(this_test_range));
%
%                 % ii_valid_range=ceil(valid_range*no_time_bins);
%                 %     ii_test_range=ceil(test_range*no_time_bins);
%
%                 XdFFtrain=X_dFF(this_training_range,:);
%                 % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
%                 XdFFtest=X_dFF(logical(this_test_range),:);
%
%                 XYtrain=pos_binned_reversed(this_training_range,:);
%
%
%
%                 %Decode using neural network
%                 MdlY1 = fitrnet(XdFFtrain,XYtrain(:,1));
%                 MdlY2 = fitrnet(XdFFtrain,XYtrain(:,2));
%
%                 label_pred(trNo).data=predict(MdlY1,XdFFtest);
%                 y_pred(trNo).data=predict(MdlY2,XdFFtest);
%             end
%
%             for trNo=1:trials.odor_trNo
%                 this_test_range=zeros(1,no_time_bins);
%
%                 if trNo==1
%                     ii_test_range_start=1;
%                 else
%                     ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
%                 end
%
%                 if trNo==trials.odor_trNo
%                     ii_test_range_end=no_time_bins;
%                 else
%                     ii_test_range_end=trials.odor_ii_end(trNo)+15;
%                 end
%
%                 this_test_range(ii_test_range_start:ii_test_range_end)=1;
%
%                 label_predicted_sh(logical(this_test_range),ii_shuffled)=label_pred(trNo).data;
%                 y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(trNo).data;
%             end
%
%         end
%
% end

fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

% label_predictedstart=1;
% label_predictedend=length(label_predicted(:,1));

XYtest=pos_binned_trimmed;

handles_out.label_predicted_sh=label_predicted_sh;
handles_out.accuracy_sh=accuracy_sh;
handles_out.label_predicted=label_predicted;
handles_out.accuracy=accuracy;
handles_out.XYtest=XYtest;
handles_out.trials=trials;



% no_conv_points=11;
% % conv_win=ones(1,no_conv_points)/no_conv_points;
% conv_win_gauss = gausswin(no_conv_points);
% conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);
% 
% label_predicted_conv=conv(label_predicted,conv_win_gauss,'same');
% accuracy_conv=conv(accuracy,conv_win_gauss,'same');

% y_predicted_conv=conv(y_predicted,conv_win_gauss,'same');

% %Now limit the x and y to max and min
% minY1=min(XYtest(:,1));
% label_predicted_conv(label_predicted_conv<minY1)=minY1;
% maxY1=max(XYtest(:,1));
% label_predicted_conv(label_predicted_conv>maxY1)=maxY1;
%
% minY2=min(XYtest(:,2));
% y_predicted_conv(y_predicted_conv<minY2)=minY2;
% maxY2=max(XYtest(:,2));
% y_predicted_conv(y_predicted_conv>maxY2)=maxY2;

% label_predicted_sh_conv=zeros(size(label_predicted_sh,1),size(label_predicted_sh,2));
% accuracy_sh_conv=zeros(size(label_predicted_sh,1),size(label_predicted_sh,2));
% for ii_sh=1:n_shuffle
%     this_label_predicted_sh=zeros(size(label_predicted_sh,1),1);
%     this_label_predicted_sh(:,1)=label_predicted_sh(:,ii_sh);
%     this_label_predicted_sh_conv=[];
%     this_label_predicted_sh_conv=conv(this_label_predicted_sh,conv_win_gauss,'same');
%     label_predicted_sh_conv(:,ii_sh)=this_label_predicted_sh_conv;
% 
%     this_accuracy_sh=zeros(size(accuracy_sh,1),1);
%     this_accuracy_sh(:,1)=accuracy_sh(:,ii_sh);
%     this_accuracy_sh_conv=[];
%     this_accuracy_sh_conv=conv(this_accuracy_sh,conv_win_gauss,'same');
%     accuracy_sh_conv(:,ii_sh)=this_accuracy_sh_conv;
% end




% fprintf(1, 'R2 nn conv x, y entire run: %d %d\n',drgGetR2(XYtest(:,1),label_predicted_conv),drgGetR2(XYtest(:,2),y_predicted_conv));
% R1=corrcoef(XYtest(:,1),label_predicted_conv);
% % R2=corrcoef(XYtest(:,2),y_predicted_conv);
% fprintf(1, 'Correlation coefficient nn conv x, y entire run %d\n\n',R1(1,2));

%
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
%
%
% hold on
% plot(label_predicted_conv(label_predictedstart:label_predictedend,1),y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1.5)
% plot(XYtest(label_predictedstart:label_predictedend,1),XYtest(label_predictedstart:label_predictedend,2),'-b','LineWidth',1.5)
% switch which_training_algorithm
%     case 1
%         title('xy for neural network convolved (trained per session)')
%     case 2
%         title('xy for neural network convolved (trained per trial)')
%     case 3
%         title('xy for neural network convolved (trained between)')
% end

% 
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% 
% hFig = figure(figNo);
% 
% set(hFig, 'units','normalized','position',[.1 .1 .7 .3])
% 
% 
% hold on
% 
% label_predictedstart=1;
% label_predictedend=length(label_predicted(:,1));
% % plot(lane_labels(1,label_predictedstart:label_predictedend),'-b','LineWidth',3)
% plot(mean(label_predicted(label_predictedstart:label_predictedend,:),2),'-k','LineWidth',1)
% 
% for trNo=1:trials.odor_trNo
% 
%     this_label_predictedstart=trials.odor_ii_start(trNo)-10;
%     this_label_predictedend=trials.odor_ii_end(trNo)+15;
%     if trials.lane_identity(trNo)==0
%         plot([this_label_predictedstart:this_label_predictedend],...
%             0.5*ones(this_label_predictedend-this_label_predictedstart+1,1),'-r','LineWidth',2)
%     else
%         plot([this_label_predictedstart:this_label_predictedend],...
%             0.5*ones(this_label_predictedend-this_label_predictedstart+1,1),'-b','LineWidth',2)
%     end
% end
% 
% 
% % switch which_training_algorithm
% %     case 1
% %         title('lane prediciton for nn, b:original, r:predicted (trained per session)')
% %     case 2
% title('lane prediction for nn, b:lane 4, r:lane 1 (trained per trial)')
%     case 3
%         title('lane prediciton for nn, b:original, r:predicted (trained between)')
% end
%
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
%
%
% hold on
% plot(XYtest(label_predictedstart:label_predictedend,1),label_predicted_conv(label_predictedstart:label_predictedend,1),'.b')
%
% xlabel('Actual x')
% ylabel('Decoded x')

% switch which_training_algorithm
%     case 1
%         title('x for nn (trained per session)')
%     case 2
% title('x for nn (trained per trial)')
%     case 3
%         title('x for nn (trained between)')
% end
%
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .7 .3])
%
%
% hold on
% plot(XYtest(label_predictedstart:label_predictedend,2),'-b','LineWidth',3)
% plot(y_predicted_conv(label_predictedstart:label_predictedend,1),'-r','LineWidth',1)
%
% for trNo=1:trials.odor_trNo
%     this_label_predictedstart=trials.odor_ii_start(trNo)-10;
%     this_label_predictedend=trials.odor_ii_end(trNo)+15;
%     plot([this_label_predictedstart:this_label_predictedend],75*ones(this_label_predictedend-this_label_predictedstart+1,1),'-k','LineWidth',2)
% end
%
%
% switch which_training_algorithm
%     case 1
%         title('y for nn, b:original, r:predicted (trained per session)')
%     case 2
%         title('y for nn, b:original, r:predicted (trained per trial)')
%     case 3
%         title('y for nn, b:original, r:predicted (trained between)')
% end


% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
%
%
% hold on
% plot(XYtest(label_predictedstart:label_predictedend,2),y_predicted_conv(label_predictedstart:label_predictedend,1),'.b')
% xlabel('Actual y')
% ylabel('Decoded y')
%
% switch which_training_algorithm
%     case 1
%         title('y for nn(trained per session)')
%     case 2
%         title('y for nn (trained per trial)')
%     case 3
%         title('y for nn (trained between)')
% end


%Keep track of the per trial decoding
lane_labels_all_trials=[];
% y_all_trials=[];
lane_labels_decod_all_trials=[];
% y_decod_all_trials=[];
lane_labels_all_trials_sh=[];
% y_all_trials_sh=[];
lane_labels_decod_all_trials_sh=[];
% y_decod_all_trials_sh=[];

lane_labels_between_trials=[];
% y_between_trials=[];
lane_labels_decod_between_trials=[];
% y_decod_between_trials=[];
lane_labels_between_trials_sh=[];
% y_between_trials_sh=[];
lane_labels_decod_between_trials_sh=[];
% y_decod_between_trials_sh=[];
% %
% %Plot the per trial results for y for permuted input
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .7 .3])
%
%
% hold on
% ii_start=0;
% last_label_predictedend=1;
% for trNo=1:trials.odor_trNo
%
%     label_predictedstart=trials.odor_ii_start(trNo)-10;
%     label_predictedend=trials.odor_ii_end(trNo)+15;
%
%     y_all_trials_sh=[y_all_trials_sh; XYtest(label_predictedstart:label_predictedend,2)];
%     y_decod_all_trials_sh=[y_decod_all_trials_sh; y_predicted_sh_conv(label_predictedstart:label_predictedend,1)];
%
%     y_between_trials_sh=[y_between_trials_sh; XYtest(last_label_predictedend:label_predictedstart,2)];
%     y_decod_between_trials_sh=[y_decod_between_trials_sh; y_predicted_sh_conv(last_label_predictedend:label_predictedstart,1)];
%     last_label_predictedend=label_predictedend;
%
%     %Plot accuracy per trial
%     ii_end=ii_start+length(XYtest(label_predictedstart:label_predictedend,2))-1;
%
%     %Plot accuracy per trial for permuted training control
%     CIsp = bootci(1000, @mean, y_predicted_sh_conv(label_predictedstart:label_predictedend,:)');
%     meansp=mean(y_predicted_sh_conv(label_predictedstart:label_predictedend,:)',1);
%     CIsp(1,:)=meansp-CIsp(1,:);
%     CIsp(2,:)=CIsp(2,:)-meansp;
%
%     [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(y_predicted_sh_conv(label_predictedstart:label_predictedend,:)',1)', CIsp', '-k');
%
%     switch trials.odor_trial_type(trNo)
%         case 1
%             %Lane 1 hits
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-r','LineWidth',3)
% %             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%         case 2
%             %Lane 1 miss
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-c','LineWidth',3)
% %             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%         case 3
%             %Lane 4 hit
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-b','LineWidth',3)
% %             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%         case 4
%            %Lane 4 hit
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-m','LineWidth',3)
% %             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%     end
%     ii_start=ii_start+length(XYtest(label_predictedstart:label_predictedend,2))+20;
%
%
% end
%
% switch which_training_algorithm
%     case 1
%         title('y for nn, permuted b:original, r:predicted (trained per session)')
%     case 2
%         title('y for nn, permuted b:original, r:predicted (trained per trial)')
%     case 3
%         title('y for nn, permuted b:original, r:predicted (trained between)')
% end
%
%
% %Plot the per trial results for y
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .7 .3])
%
%
% hold on
% ii_start=0;
% last_label_predictedend=1;
% for trNo=1:trials.odor_trNo
%
%     label_predictedstart=trials.odor_ii_start(trNo)-10;
%     label_predictedend=trials.odor_ii_end(trNo)+15;
%
%     y_all_trials=[y_all_trials; XYtest(label_predictedstart:label_predictedend,2)];
%     y_decod_all_trials=[y_decod_all_trials; y_predicted_conv(label_predictedstart:label_predictedend,1)];
%
%     y_between_trials=[y_between_trials; XYtest(last_label_predictedend:label_predictedstart,2)];
%     y_decod_between_trials=[y_decod_between_trials; y_predicted_conv(last_label_predictedend:label_predictedstart,1)];
%     last_label_predictedend=label_predictedend;
%
%     %Plot accuracy per trial
%     ii_end=ii_start+length(XYtest(label_predictedstart:label_predictedend,2))-1;
%
%     switch trials.odor_trial_type(trNo)
%         case 1
%             %Lane 1 hits
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-r','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%         case 2
%             %Lane 1 miss
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-c','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%         case 3
%             %Lane 4 hit
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-b','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%         case 4
%            %Lane 4 hit
%             plot(dt*[ii_start:ii_end]',XYtest(label_predictedstart:label_predictedend,2),'-m','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
%     end
%     ii_start=ii_start+length(XYtest(label_predictedstart:label_predictedend,2))+20;
%
%
% end
%
% switch which_training_algorithm
%     case 1
%         title('y for nn per trial, b:original, r:predicted (trained per session)')
%     case 2
%         title('y for nn per trial, b:original, r:predicted (trained per trial)')
%     case 3
%         title('y for nn per trial, b:original, r:predicted (trained between)')
% end

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

last_label_predictedend=1;
for trNo=1:trials.odor_trNo

    label_predictedstart=trials.odor_ii_start(trNo)-10;
    label_predictedend=trials.odor_ii_end(trNo)+15;
    % 
    % lane_labels_all_trials=[lane_labels_all_trials mean(lane_labels(label_predictedstart:label_predictedend,2))];
    % % lane_labels_decod_all_trials=[lane_labels_decod_all_trials label_predicted_conv(label_predictedstart:label_predictedend)'];
    % 
    % % lane_labels_between_trials=[lane_labels_between_trials lane_labels(last_label_predictedend:label_predictedstart)];
    % % lane_labels_decod_between_trials=[lane_labels_decod_between_trials label_predicted_conv(last_label_predictedend:label_predictedstart)'];
    % last_label_predictedend=label_predictedend;

    ii_end=ii_start+length(lane_labels(label_predictedstart:label_predictedend))-1;

    switch trials.lane_identity(trNo)
        case 0
            %Lane 1
            plot(dt*[ii_start:ii_end]',ones(1,length(dt*[ii_start:ii_end]))+0.05,'-r','LineWidth',3)

            if size(label_predicted,2)>=3
                CIsp = bootci(1000, @mean, label_predicted(label_predictedstart:label_predictedend,:)');
                meansp=mean(label_predicted(label_predictedstart:label_predictedend,:)');
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;

                [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(label_predicted(label_predictedstart:label_predictedend,:)')', CIsp', '-k');
            else
                plot(dt*[ii_start:ii_end]',mean(label_predicted(label_predictedstart:label_predictedend,:)')', '-k');
            end

            % plot(dt*[ii_start:ii_end]',mean(label_predicted(label_predictedstart:label_predictedend,:),2),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 1.1],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 1.1],'-k')
        case 1
            %Lane 4
            plot(dt*[ii_start:ii_end]',zeros(1,length(dt*[ii_start:ii_end]))-0.05,'-b','LineWidth',3)

            if size(label_predicted,2)>=3
                CIsp = bootci(1000, @mean, label_predicted(label_predictedstart:label_predictedend,:)');
                meansp=mean(label_predicted(label_predictedstart:label_predictedend,:)');
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;

                [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(label_predicted(label_predictedstart:label_predictedend,:)')', CIsp', '-k');
            else
                plot(dt*[ii_start:ii_end]',mean(label_predicted(label_predictedstart:label_predictedend,:)')', '-k');
            end

            % plot(dt*[ii_start:ii_end]',mean(label_predicted(label_predictedstart:label_predictedend),:),2),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 1.1],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 1.1],'-k')
    end
    ii_start=ii_start+length(lane_labels(label_predictedstart:label_predictedend))+20;
end

% switch which_training_algorithm
%     case 1
%         title('lane labels for nn per trial, r: hit1, c: mis1, b: hit4, m: miss4 (trained per session)')
%     case 2
title('Predicted lane labels per trial, r: lane 1, b: lane 4')
%     case 3
%         title('lane labels for nn per trial, r: hit1, c: mis1, b: hit4, m: miss4 (trained between)')
% end

ylim([-0.1 1.1])


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

last_label_predictedend=1;
for trNo=1:trials.odor_trNo

    label_predictedstart=trials.odor_ii_start(trNo)-10;
    label_predictedend=trials.odor_ii_end(trNo)+15;

    % lane_labels_all_trials=[lane_labels_all_trials lane_labels(label_predictedstart:label_predictedend)];
    % % lane_labels_decod_all_trials=[lane_labels_decod_all_trials label_predicted_conv(label_predictedstart:label_predictedend)'];
    % 
    % lane_labels_between_trials=[lane_labels_between_trials lane_labels(last_label_predictedend:label_predictedstart)];
    % lane_labels_decod_between_trials=[lane_labels_decod_between_trials label_predicted_conv(last_label_predictedend:label_predictedstart)'];
    % last_label_predictedend=label_predictedend;

    ii_end=ii_start+length(lane_labels(label_predictedstart:label_predictedend))-1;

    switch trials.lane_identity(trNo)
        case 0
            %Lane 1
            plot(dt*[ii_start:ii_end]',ones(1,length(dt*[ii_start:ii_end]))+0.05,'-r','LineWidth',3)

            
            CIsp = bootci(1000, @mean, label_predicted_sh(label_predictedstart:label_predictedend,:)');
            meansp=mean(label_predicted_sh(label_predictedstart:label_predictedend,:)');
            CIsp(1,:)=meansp-CIsp(1,:);
            CIsp(2,:)=CIsp(2,:)-meansp;

            [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(label_predicted_sh(label_predictedstart:label_predictedend,:)')', CIsp', '-k');


            % plot(dt*[ii_start:ii_end]',mean(label_predicted_sh(label_predictedstart:label_predictedend,:),2),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 1],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 1],'-k')
        case 1
            %Lane 4
            plot(dt*[ii_start:ii_end]',zeros(1,length(dt*[ii_start:ii_end]))-0.05,'-b','LineWidth',3)


            CIsp = bootci(1000, @mean, label_predicted_sh(label_predictedstart:label_predictedend,:)');
            meansp=mean(label_predicted_sh(label_predictedstart:label_predictedend,:)');
            CIsp(1,:)=meansp-CIsp(1,:);
            CIsp(2,:)=CIsp(2,:)-meansp;

            [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(label_predicted_sh(label_predictedstart:label_predictedend,:)')', CIsp', '-k');

            % plot(dt*[ii_start:ii_end]',mean(label_predicted_sh(label_predictedstart:label_predictedend,:),2),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 1],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 1],'-k')
    end
    ii_start=ii_start+length(lane_labels(label_predictedstart:label_predictedend))+20;
end

% switch which_training_algorithm
%     case 1
%         title('lane labels for nn per trial, r: hit1, c: mis1, b: hit4, m: miss4 (trained per session)')
%     case 2
title('Predicted mean lane labels per trial, shuffled decon, r: lane 1, b: lane 4')
%     case 3
%         title('lane labels for nn per trial, r: hit1, c: mis1, b: hit4, m: miss4 (trained between)')
% end

ylim([-0.1 1.1])



% 
% %Plot the per trial results for x with nn trained with permuted input
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% 
% hFig = figure(figNo);
% 
% set(hFig, 'units','normalized','position',[.1 .1 .7 .3])
% 
% 
% hold on
% ii_start=0;
% last_label_predictedend=1;
% for trNo=1:trials.odor_trNo
% 
%     label_predictedstart=trials.odor_ii_start(trNo)-10;
%     label_predictedend=trials.odor_ii_end(trNo)+15;
% 
%     lane_labels_all_trials_sh=[lane_labels_all_trials_sh lane_labels(label_predictedstart:label_predictedend)];
%     lane_labels_decod_all_trials_sh=[lane_labels_decod_all_trials_sh label_predicted_sh_conv(label_predictedstart:label_predictedend)];
% 
%     lane_labels_between_trials_sh=[lane_labels_between_trials_sh lane_labels(last_label_predictedend:label_predictedstart)];
%     lane_labels_decod_between_trials_sh=[lane_labels_decod_between_trials_sh label_predicted_sh_conv(last_label_predictedend:label_predictedstart)];
%     last_label_predictedend=label_predictedend;
% 
%     ii_end=ii_start+length(lane_labels(label_predictedstart:label_predictedend))-1;
% 
    % %Plot label predicted per trial for permuted training control
    % CIsp = bootci(1000, @mean, label_predicted_sh_conv(label_predictedstart:label_predictedend,:)');
    % meansp=mean(label_predicted_sh_conv(label_predictedstart:label_predictedend,:)');
    % CIsp(1,:)=meansp-CIsp(1,:);
    % CIsp(2,:)=CIsp(2,:)-meansp;
    % 
    % [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(label_predicted_sh_conv(label_predictedstart:label_predictedend,:)')', CIsp', '-k');
    % 

%     switch trials.odor_trial_type(trNo)
%         case 1
%             %Lane 1 hits
%             plot(dt*[ii_start:ii_end]',lane_labels(label_predictedstart:label_predictedend),'-r','LineWidth',3)
%             %             plot(dt*[ii_start:ii_end]',label_predicted_conv(label_predictedstart:label_predictedend),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 4.5],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 4.5],'-k')
%         case 2
%             %Lane 1 miss
%             plot(dt*[ii_start:ii_end]',lane_labels(label_predictedstart:label_predictedend),'-c','LineWidth',3)
%             %             plot(dt*[ii_start:ii_end]',label_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 4.5],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 4.5],'-k')
%         case 3
%             %Lane 4 hit
%             plot(dt*[ii_start:ii_end]',lane_labels(label_predictedstart:label_predictedend),'-b','LineWidth',3)
%             %             plot(dt*[ii_start:ii_end]',label_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 4.5],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 4.5],'-k')
%         case 4
%             %Lane 4 hit
%             plot(dt*[ii_start:ii_end]',lane_labels(label_predictedstart:label_predictedend),'-m','LineWidth',3)
%             %             plot(dt*[ii_start:ii_end]',label_predicted_conv(label_predictedstart:label_predictedend,1),'-k','LineWidth',1)
%             plot(dt*[ii_start+10 ii_start+10],[0 4.5],'-k')
%             plot(dt*[ii_end-15 ii_end-15],[0 4.5],'-k')
%     end
%     ii_start=ii_start+length(lane_labels(label_predictedstart:label_predictedend))+20;
% end
% 
% % switch which_training_algorithm
% %     case 1
% %         title('x for nn per trial permuted, r: hit1, c: mis1, b: hit4, m: miss4 (trained per session)')
% %     case 2
% title('x for nn per trial permuted, r: hit1, c: mis1, b: hit4, m: miss4 (trained per trial)')
% %     case 3
% %         title('x for nn per trial permuted, r: hit1, c: mis1, b: hit4, m: miss4 (trained between)')
% % end
% 
% ylim([-0.1 1.1])

% % fprintf(1, 'R2 nn conv x, y per trial run: %d %d\n',drgGetR2(lane_labels_all_trials,lane_labels_decod_all_trials),drgGetR2(y_all_trials,y_decod_all_trials));
% R1=corrcoef(lane_labels_all_trials,lane_labels_decod_all_trials);
% % R2=corrcoef(y_all_trials,y_decod_all_trials);
% fprintf(1, 'Correlation coefficient nn conv x, y per trial %d\n',R1(1,2));
%
% R1=corrcoef(lane_labels_between_trials,lane_labels_decod_between_trials);
% % R2=corrcoef(y_between_trials,y_decod_between_trials);
% fprintf(1, 'Correlation coefficient nn conv x, y between trials %d\n\n',R1(1,2));
%
% % fprintf(1, 'R2 nn permuted x, y per trial run: %d %d\n',drgGetR2(lane_labels_all_trials_sh,lane_labels_decod_all_trials_sh),drgGetR2(y_all_trials_sh,y_decod_all_trials_sh));
% R1=corrcoef(lane_labels_all_trials_sh,lane_labels_decod_all_trials_sh);
% % R2=corrcoef(y_all_trials_sh,y_decod_all_trials_sh);
% fprintf(1, 'Correlation coefficient nn permuted x, y per trial %d\n',R1(1,2));
%
% R1=corrcoef(lane_labels_between_trials_sh,lane_labels_decod_between_trials_sh);
% % R2=corrcoef(y_between_trials_sh,y_decod_between_trials_sh);
% fprintf(1, 'Correlation coefficient nn permuted x, y between trial %d\n',R1(1,2));
%
% %Plot x vs decoded
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
%
%
% hold on
%
% plot(lane_labels_all_trials,lane_labels_decod_all_trials,'.b')
% xlabel('x')
% ylabel('x decoded')
% % switch which_training_algorithm
% %     case 1
% %         title('x for nn per trial(trained per session)')
% %     case 2
%         title('x for nn per trial (trained per trial)')
%     case 3
%         title('x for nn per trial (trained between)')
% end

%Now plot average convolved accuracy per trial


% %start and end of display period
% align_display=0; %0= odor start, 1=odor end
% dt_display_start=-5; %seconds from start alignment
% dt_display_end=5; %seconds from end alignment
% ii_dt_display_start=fix(dt_display_start/dt); %samples from start alignment
% ii_dt_display_end=fix(dt_display_end/dt); %samples from end alignment

for trNo=1:trials.odor_trNo
    label_predictedstart(trNo)=trials.odor_ii_start(trNo)-10;
    label_predictedend(trNo)=trials.odor_ii_end(trNo)+15;
end

trial_dt=delta_t_per_trial*dt;
mean_trial_dt=mean(delta_t_per_trial*dt);


%Show accuracy aligned to the start of the trial
align_display=0;
accuracy_all_trials=[];
accuracy_conv_all_trials=[];
accuracy_sh_conv_all_trials=[];
accuracy_sh_all_trials=[];
delta_odor_start_to_end=ceil(mean(trials.odor_ii_end-trials.odor_ii_start));
for trNo=1:trials.odor_trNo
    switch align_display
        case 0
            %Aligned to start
            this_display_start=trials.odor_ii_start(trNo)+ii_dt_display_start;
            this_display_end=trials.odor_ii_start(trNo)+delta_odor_start_to_end+ii_dt_display_end;
        case 1
            this_display_start=trials.odor_ii_end(trNo)-delta_odor_start_to_end+ii_dt_display_start;
            this_display_end=trials.odor_ii_end(trNo)+ii_dt_display_end;
    end
    for ii_rp=1:no_repeats
        accuracy_all_trials=[accuracy_all_trials accuracy(this_display_start:this_display_end,ii_rp)];
    end
    % accuracy_conv_all_trials=[accuracy_conv_all_trials accuracy_conv(this_display_start:this_display_end)];

    for ii_sh=1:n_shuffle
        accuracy_sh_all_trials=[accuracy_sh_all_trials accuracy_sh(this_display_start:this_display_end,ii_sh)];
        % accuracy_sh_conv_all_trials=[accuracy_sh_conv_all_trials accuracy_sh_conv(this_display_start:this_display_end,ii_sh)];
    end
    %
    %     lane_labels_between_trials=[lane_labels_between_trials lane_labels(last_label_predictedend:label_predictedstart)];
    %     lane_labels_decod_between_trials=[lane_labels_decod_between_trials label_predicted_conv(last_label_predictedend:label_predictedstart)'];
    %     last_label_predictedend=label_predictedend;

    %     ii_end=ii_start+length(lane_labels(label_predictedstart:label_predictedend))-1;
end

time_span=dt*[1:size(accuracy_all_trials,1)]+dt_display_start;


%Plot accuracy
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

CIpvsm = bootci(1000, @mean, accuracy_sh_all_trials');
meanpvsm=mean(accuracy_sh_all_trials',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(time_span,mean(accuracy_sh_all_trials'), CIpvsm','r');

CIpvsm = bootci(1000, @mean, accuracy_all_trials');
meanpvsm=mean(accuracy_all_trials',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(time_span,mean(accuracy_all_trials'), CIpvsm','k');

plot(time_span,mean(accuracy_sh_all_trials'), '-r');
plot(time_span,mean(accuracy_all_trials'), '-k');

ylim([0 1.1])

plot([0 0],[0 1.1],'-k')

plot([mean_trial_dt mean_trial_dt],[0 1.1],'-k')

xlabel('Time (sec)')
ylabel('Accuracy')
if align_training_start==0
    title(['Accuracy aligned to trial start, trained with trial end'])
else
    title(['Accuracy aligned to trial start, trained with trial start'])
end

pfft=1;


%Plot accuracy separately for each event
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

%hit1 red
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==1),2), '-r');

%miss1 cyan
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==2),2), '-c');

%hit4 blue
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==3),2), '-b');

%miss4 magenta
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==4),2), '-m');

plot([0 0],[0 1.1],'-k')

plot([mean_trial_dt mean_trial_dt],[0 1.1],'-k')

xlabel('Time (sec)')
ylabel('Accuracy')
title(['Accuracy trial start, hit1=r, miss1=c, hit4=b, miss4=m'])



%Now show accuracy time course aligned to the end
align_display=1;
accuracy_all_trials=[];
accuracy_conv_all_trials=[];
accuracy_sh_conv_all_trials=[];
accuracy_sh_all_trials=[];
for trNo=1:trials.odor_trNo
    switch align_display
        case 0
            %Aligned to start
            this_display_start=trials.odor_ii_start(trNo)+ii_dt_display_start;
            this_display_end=trials.odor_ii_start(trNo)+delta_odor_start_to_end+ii_dt_display_end;
        case 1
            this_display_start=trials.odor_ii_end(trNo)-delta_odor_start_to_end+ii_dt_display_start;
            this_display_end=trials.odor_ii_end(trNo)+ii_dt_display_end;
    end

    for ii_rp=1:no_repeats
        accuracy_all_trials=[accuracy_all_trials accuracy(this_display_start:this_display_end,ii_rp)];
        % accuracy_conv_all_trials=[accuracy_conv_all_trials accuracy_conv(this_display_start:this_display_end)];
    end

    for ii_sh=1:n_shuffle
        accuracy_sh_all_trials=[accuracy_sh_all_trials accuracy_sh(this_display_start:this_display_end,ii_sh)];
        % accuracy_sh_conv_all_trials=[accuracy_sh_conv_all_trials accuracy_sh_conv(this_display_start:this_display_end,ii_sh)];
    end
    %
    %     lane_labels_between_trials=[lane_labels_between_trials lane_labels(last_label_predictedend:label_predictedstart)];
    %     lane_labels_decod_between_trials=[lane_labels_decod_between_trials label_predicted_conv(last_label_predictedend:label_predictedstart)'];
    %     last_label_predictedend=label_predictedend;

    %     ii_end=ii_start+length(lane_labels(label_predictedstart:label_predictedend))-1;
end




%Plot accuracy
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

CIpvsm = bootci(1000, @mean, accuracy_sh_all_trials');
meanpvsm=mean(accuracy_sh_all_trials',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(time_span,mean(accuracy_sh_all_trials'), CIpvsm','r');

CIpvsm = bootci(1000, @mean, accuracy_all_trials');
meanpvsm=mean(accuracy_all_trials',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(time_span,mean(accuracy_all_trials'), CIpvsm','k');

plot(time_span,mean(accuracy_sh_all_trials'), '-r');

ylim([0 1.1])


plot([0 0],[0 1.1],'-k')

plot([mean_trial_dt mean_trial_dt],[0 1.1],'-k')

% plot([mean_trial_dt mean_trial_dt],[0 1.1],'-k')

xlabel('Time (sec)')
ylabel('Accuracy')
title(['Accuracy aligned to trial end'])

pfft=1;


%Plot accuracy separately for each event
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

%hit1 red
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==1),2), '-r');
plot(time_span,mean(accuracy_sh_all_trials(:,trials.odor_trial_type==1),2), '-','Color',[1 0.7 0.7]);

%miss1 cyan
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==2),2), '-c');
plot(time_span,mean(accuracy_sh_all_trials(:,trials.odor_trial_type==2),2), '-','Color',[0 0.7 0.7]);

%hit4 blue
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==3),2), '-b');
plot(time_span,mean(accuracy_sh_all_trials(:,trials.odor_trial_type==3),2), '-','Color',[0.7 0.7 1]);

%miss4 magenta
plot(time_span,mean(accuracy_all_trials(:,trials.odor_trial_type==4),2), '-m');
plot(time_span,mean(accuracy_sh_all_trials(:,trials.odor_trial_type==4),2), '-','Color',[1 0.7 1]);



plot([0 0],[0 1.1],'-k')

plot([mean_trial_dt mean_trial_dt],[0 1.1],'-k')

xlabel('Time (sec)')
ylabel('Accuracy')
title(['Accuracy trial end, hit1=r, miss1=c, hit4=b, miss4=m '])
% %Plot y vs decoded
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
%
%
% hold on
%
% plot(y_all_trials,y_decod_all_trials,'.b')
% xlabel('y')
% ylabel('y decoded')
%
%
% switch which_training_algorithm
%     case 1
%         title('y for nn per trial (trained per session)')
%     case 2
%         title('y for nn per trial (trained per trial)')
%     case 3
%         title('y for nn per trial (trained between)')
% end
handles_out.accuracy_conv_all_trials=accuracy_conv_all_trials;
handles_out.time_span=time_span;
handles_out.accuracy_sh_conv_all_trials=accuracy_sh_conv_all_trials;
handles_out.trial_dt=trial_dt;
try
    delete(gcp('nocreate'));
catch
end
% save([this_path arena_file(1:end-4) 'lane_algo' num2str(which_ml_algo) '.mat'],'handles_out','handles_choices2','-v7.3')

pffft=1;
