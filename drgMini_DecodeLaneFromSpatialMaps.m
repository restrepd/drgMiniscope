function handles_out=drgMini_DecodeLaneFromSpatialMaps(handles_choices)
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020
close all

if exist('handles_choices')==0
    clear all


    %First troubleshooting files
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';
    pred_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm_deczdFFopt2_711103.mat';

    %     %Second troubleshooting files
        % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
        % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
        % arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';

    %No odor troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';


    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

    handles_choices.z=1; %Convert dFF to z
    handles_choices.displayFigures=1;

    no_steps=10;
    handles_choices.no_steps=no_steps; %The activity map will be a square divided in no_steps x no_steps small squares

    no_repeats=20;
    handles_choices.no_repeats=no_repeats; %Number of times to repeat fitrnet

    %     isKording=0;

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    %     dt=0.2;

    which_ROIs=3; %1 Use all ROIs, 2 Use place cells, 3 use odor cells

    %Define the different ranges (training, valid and testing)
    % training_fraction=0.9;
    % handles_choices.training_fraction=training_fraction;

    %     training_range=[0, 0.5];
    %     valid_range=[0.5,0.65];
    %     test_range=[0.5, 1];

    % %The user can define what time period to use spikes from (with respect to the output).
    % bins_before=5; %How many bins of neural data prior to the output are used for decoding, 10
    % bins_current=1; %Whether to use concurrent time bin of neural data, 1
    % bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    % handles_choices.bins_before=bins_before;
    % handles_choices.bins_current=bins_current;
    % handles_choices.bins_after=bins_after;


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2;
    dt_miniscope=1/30;
    n_shuffle=10; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    handles_choices.n_shuffle=n_shuffle;

    % which_training_algorithm=2;
    % handles_choices.which_training_algorithm=which_training_algorithm;
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

    algo_name{1}='SVZ';
    algo_name{2}='NN';
    algo_name{3}='Tree';
    algo_name{4}='Bayes';
    algo_name{5}='glm';
    algo_name{6}='LD';
    algo_name{7}='NNopt';

    fprintf(1,['\n' num2str(no_steps) 'x' num2str(no_steps) ' square. Algorithm is ' algo_name{which_ml_algo} '\n\n'])

    ii_cost=3;

    this_cost=[0 ii_cost;ii_cost 0];

    %start and end of training period
    align_training_start=0; %0= odor start, 1=odor end
    align_training_end=1; %0= odor start, 1=odor end
    dt_training_start=-5; %seconds from start alignment
    dt_training_end=5; %seconds from end alignment
    ii_dt_training_start=fix(dt_training_start/dt); %samples from start alignment
    ii_dt_training_end=fix(dt_training_end/dt); %samples from end alignment

    %start and end of display period
    align_display=0; %0= odor start, 1=odor end
    dt_display_start=-5; %seconds from start alignment
    dt_display_end=5; %seconds from end alignment (or mean end)
    ii_dt_display_start=fix(dt_display_start/dt); %samples from start alignment
    ii_dt_display_end=fix(dt_display_end/dt); %samples from end alignment

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

    %Load prediction file and define place and odor cells
    load([this_path pred_file])

    spatial_rhol1l4=handles_out.spatial_rhol1l4;
    delta_center_of_mass=handles_out.delta_center_of_mass;

    no_neurons=length(spatial_rhol1l4);

if which_ROIs~=1
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
else
    process_these_ROIs=[1:no_neurons];
end

% switch which_training_algorithm
%     case 1
%         fprintf(1,['\nTrained with data for the entire session\n\n'])
%     case 2
%         fprintf(1,['\nTrained with within trial data\n\n'])
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
 
%Use a subset of the ROIs
dFF=dFF(:,process_these_ROIs);


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

load([this_path arena_file])
% 
% %Extract trials
% trials=[];
% 
% %Extract odor on using the camera sync
% at_end=0;
% ii=0;
% jj=0;
% jj_l1=0;
% jj_l4=0;
% while at_end==0
%     next_ii=find(arena.odorsync(ii+1:end)==1,1,'first');
%     if ~isempty(next_ii)
%         jj=jj+1;
%         trials.ii_odor(jj)=ii+next_ii;
%         trials.x_odor(jj)=arena.xsync(ii+next_ii);
%         trials.y_odor(jj)=arena.ysync(ii+next_ii);
% 
%         ii=ii+next_ii;
%         ii_mini=arena.index_flirsynctominiscope(ii);
% 
%         if sum(arena.laneodor1(ii_mini-5:ii_mini+5)==1)>0
%             %Note: laneodor4 is 1 only for one time point
%             jj_l1=jj_l1+1;
%             trials.ii_laneodor1(jj_l1)=ii;
%             trials.x_laneodor1(jj_l1)=trials.x_odor(jj);
%             trials.y_laneodor1(jj_l1)=trials.y_odor(jj);
%         end
% 
%         if sum(arena.laneodor4(ii_mini-5:ii_mini+5)==1)>0
%             %Note: laneodor4 is 1 only for one time point
%             jj_l4=jj_l4+1;
%             trials.ii_laneodor4(jj_l4)=ii;
%             trials.x_laneodor4(jj_l4)=trials.x_odor(jj);
%             trials.y_laneodor4(jj_l4)=trials.y_odor(jj);
%         end
% 
%         next_ii=find(arena.odorsync(ii+1:end)==0,1,'first');
%         if ~isempty(next_ii)
%             ii=ii+next_ii;
%         else
%             at_end=1;
%         end
%     else
%         at_end=1;
%     end
% end
% % 
% % %Extract lane1
% % at_end=0;
% % ii_flir=0;
% % jj=0;
% % while at_end==0
% %     next_ii_flir=find(arena.laneodor1(ii_flir+1:end)==1,1,'first');
% %     if ~isempty(next_ii_flir)
% %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
% %         jj=jj+1;
% %         trials.ii_laneodor1(jj)=next_ii;
% %         trials.x_laneodor1(jj)=arena.xsync(next_ii);
% %         trials.y_laneodor1(jj)=arena.ysync(next_ii);
% %         next_ii_flir2=find(arena.laneodor1(ii_flir+next_ii_flir+1:end)==0,1,'first');
% %         %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
% %         if ~isempty(next_ii_flir2)
% %             ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
% %         else
% %             at_end=1;
% %         end
% %     else
% %         at_end=1;
% %     end
% % end
% % 
% % %Extract lane4
% % at_end=0;
% % ii_flir=0;
% % jj=0;
% % while at_end==0
% %     next_ii_flir=find(arena.laneodor4(ii_flir+1:end)==1,1,'first');
% %     if ~isempty(next_ii_flir)
% %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
% %         jj=jj+1;
% %         trials.ii_laneodor4(jj)=next_ii;
% %         trials.x_laneodor4(jj)=arena.xsync(next_ii);
% %         trials.y_laneodor4(jj)=arena.ysync(next_ii);
% %         next_ii_flir2=find(arena.laneodor4(ii_flir+next_ii_flir+1:end)==0,1,'first');
% %         %         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
% %         if ~isempty(next_ii_flir2)
% %             ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
% %         else
% %             at_end=1;
% %         end
% %     else
% %         at_end=1;
% %     end
% % end
% 
% %Extract lanewater1
% at_end=0;
% ii_flir=0;
% jj=0;
% while at_end==0
%     next_ii_flir=find(arena.lanewater1(ii_flir+1:end)==1,1,'first');
%     if ~isempty(next_ii_flir)
%         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
%         jj=jj+1;
%         trials.ii_lanewater1(jj)=next_ii;
%         trials.x_lanewater1(jj)=arena.xsync(next_ii);
%         trials.y_lanewater1(jj)=arena.ysync(next_ii);
%         next_ii_flir2=find(arena.lanewater1(ii_flir+next_ii_flir+1:end)==0,1,'first');
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
% 
% %Extract lanewater4
% at_end=0;
% ii_flir=0;
% jj=0;
% while at_end==0
%     next_ii_flir=find(arena.lanewater4(ii_flir+1:end)==1,1,'first');
%     if ~isempty(next_ii_flir)
%         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
%         jj=jj+1;
%         trials.ii_lanewater4(jj)=next_ii;
%         trials.x_lanewater4(jj)=arena.xsync(next_ii);
%         trials.y_lanewater4(jj)=arena.ysync(next_ii);
%         next_ii_flir2=find(arena.lanewater4(ii_flir+next_ii_flir+1:end)==0,1,'first');
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





% dFF_times=[1:no_time_points]*dt_miniscope;

no_neurons=size(dFF,2);


%Now do space activity maps and calculate information content
%Now show the pseudocolor activity plots
dx=500/no_steps;
x=dx/2:dx:500-dx/2;
dy=480/no_steps;
y=dy/2:dy:480-dy/2;

X_square=zeros(no_steps,no_steps);
Y_square=zeros(no_steps,no_steps);
for ii_x=1:no_steps
    for ii_y=1:no_steps
        X_square(ii_x,ii_y)=x(ii_x);
        Y_square(ii_x,ii_y)=y(ii_y);
    end
end


%Calculate activity maps and information content
figNo=figNo+1;

activity_maps=[];

for this_ROI=1:no_neurons

    %Initialize variables
    activity_maps.ROI(this_ROI).all_trial_map=zeros(no_steps,no_steps);
    activity_maps.ROI(this_ROI).all_trial_map_n=zeros(no_steps,no_steps);

    activity_maps.ROI(this_ROI).l1_trial_map=zeros(no_steps,no_steps);
    activity_maps.ROI(this_ROI).l1_trial_map_n=zeros(no_steps,no_steps);

    activity_maps.ROI(this_ROI).l4_trial_map=zeros(no_steps,no_steps);
    activity_maps.ROI(this_ROI).l4_trial_map_n=zeros(no_steps,no_steps);

    %Calculate the activity maps per trial
    for trNo=1:trials.odor_trNo

        activity_maps.ROI(this_ROI).trial(trNo).this_trial_map=zeros(no_steps,no_steps);
        activity_maps.ROI(this_ROI).trial(trNo).this_trial_map_n=zeros(no_steps,no_steps);

        for ii=trials.odor_ii_start(trNo):trials.odor_ii_end(trNo)

            this_x_ii=ceil(arena.xsync(ii)/dx);
            if this_x_ii==no_steps+1
                this_x_ii=no_steps;
            end

            this_y_ii=ceil(arena.ysync(ii)/dy);
            if this_y_ii==no_steps+1
                this_y_ii=no_steps;
            end

            activity_maps.ROI(this_ROI).trial(trNo).this_trial_map(this_x_ii,this_y_ii)=...
                activity_maps.ROI(this_ROI).trial(trNo).this_trial_map(this_x_ii,this_y_ii)+dFF(ii,this_ROI);

            activity_maps.ROI(this_ROI).trial(trNo).this_trial_map_n(this_x_ii,this_y_ii)=...
                activity_maps.ROI(this_ROI).trial(trNo).this_trial_map_n(this_x_ii,this_y_ii)+1;

            activity_maps.ROI(this_ROI).all_trial_map(this_x_ii,this_y_ii)=...
                activity_maps.ROI(this_ROI).all_trial_map(this_x_ii,this_y_ii)+dFF(ii,this_ROI);

            activity_maps.ROI(this_ROI).all_trial_map_n(this_x_ii,this_y_ii)=...
                activity_maps.ROI(this_ROI).all_trial_map_n(this_x_ii,this_y_ii)+1;

            if trials.odor_lane(trNo)==1
                activity_maps.ROI(this_ROI).l1_trial_map(this_x_ii,this_y_ii)=...
                    activity_maps.ROI(this_ROI).l1_trial_map(this_x_ii,this_y_ii)+dFF(ii,this_ROI);

                activity_maps.ROI(this_ROI).l1_trial_map_n(this_x_ii,this_y_ii)=...
                    activity_maps.ROI(this_ROI).l1_trial_map_n(this_x_ii,this_y_ii)+1;
            else
                activity_maps.ROI(this_ROI).l4_trial_map(this_x_ii,this_y_ii)=...
                    activity_maps.ROI(this_ROI).l4_trial_map(this_x_ii,this_y_ii)+dFF(ii,this_ROI);

                activity_maps.ROI(this_ROI).l4_trial_map_n(this_x_ii,this_y_ii)=...
                    activity_maps.ROI(this_ROI).l4_trial_map_n(this_x_ii,this_y_ii)+1;
            end



        end

        %Now normalize for this trial
        activity_maps.ROI(this_ROI).trial(trNo).this_trial_map=...
            activity_maps.ROI(this_ROI).trial(trNo).this_trial_map./...
            activity_maps.ROI(this_ROI).trial(trNo).this_trial_map_n;

    end

    %Now normalize for all trials
    activity_maps.ROI(this_ROI).all_trial_map=...
        activity_maps.ROI(this_ROI).all_trial_map./...
        activity_maps.ROI(this_ROI).all_trial_map_n;

    activity_maps.ROI(this_ROI).l1_trial_map=...
        activity_maps.ROI(this_ROI).l1_trial_map./...
        activity_maps.ROI(this_ROI).l1_trial_map_n;

    activity_maps.ROI(this_ROI).l4_trial_map=...
        activity_maps.ROI(this_ROI).l4_trial_map./...
        activity_maps.ROI(this_ROI).l4_trial_map_n;

    %
end



%Now do the decoding for the first run
label_pred=[];
Mdl_pred=[];
for trNo1=1:trials.odor_trNo
    % parfor trNo=1:trials.odor_trNo

    %Setup the training data set
    SpActrain=[];

    these_trNos=[];
    for trNo2=1:trials.odor_trNo
        if trNo2~=trNo1
            these_trNos=[these_trNos trNo2];
            this_row=[];
            for this_ROI=1:no_neurons
                for x_ii=1:no_steps
                    for y_ii=1:no_steps
                        if ~isnan(activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii))
                            pffft=1;
                            if ~isnan(activity_maps.ROI(this_ROI).trial(trNo2).this_trial_map(x_ii,y_ii))
                                this_row=[this_row activity_maps.ROI(this_ROI).trial(trNo2).this_trial_map(x_ii,y_ii)];
                            else
                                this_row=[this_row activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii)];
                            end
                        end
                    end
                end
            end
    
            SpActrain=[SpActrain; this_row];
        end
    end

    this_trNo=0;
    Y=[];
    for trNo3=these_trNos
        this_trNo=this_trNo+1;
        if trials.odor_lane(trNo3)==1
            Y(this_trNo)=0;
        else
            Y(this_trNo)=1;
        end
    end
    
    %Decode using neural network
    tblTrn=[];
    tblTrn = array2table(SpActrain);

  
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
            Mdl = fitglm(tblTrn,Y','Distribution','binomial');
        case 6
            Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
        case 7
            opts = struct('ShowPlots', false, ...
                'Verbose', 0, ...
                'MaxObjectiveEvaluations', 15,...
                'UseParallel',true);
            Mdl = fitcnet(tblTrn,Y,'OptimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions', opts,'ClassNames', unique(Y),...
                'Prior', 'uniform');
    end

    Mdl_pred(trNo1).Mdl=Mdl;

    %Calculate the prediction
    this_row=[];
    for this_ROI=1:no_neurons
        for x_ii=1:no_steps
            for y_ii=1:no_steps
                if ~isnan(activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii))
                    pffft=1;
                    if ~isnan(activity_maps.ROI(this_ROI).trial(trNo1).this_trial_map(x_ii,y_ii))
                        this_row=[this_row activity_maps.ROI(this_ROI).trial(trNo1).this_trial_map(x_ii,y_ii)];
                    else
                        this_row=[this_row activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii)];
                    end
                end
            end
        end
    end
    label_pred(trNo1).data=predict(Mdl,this_row);
    
end
fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

%Parse out the parfor loop output
label_predicted=zeros(trials.odor_trNo,1);
accuracy=zeros(trials.odor_trNo,1);
accuracy_lane1=[];
accuracy_lane4=[];
sh_accuracy=zeros(trials.odor_trNo,1);
sh_accuracy_lane1=[];
sh_accuracy_lane4=[];
odor_lanes=[];
%         y_predicted=zeros(no_time_bins,1);
perm_trials=randperm(trials.odor_trNo);
for trNo=1:trials.odor_trNo
    label_predicted(trNo)=label_pred(trNo).data;
    %             y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

    %Compute accuracy
    this_accuracy=0;
    this_lane1_accuracy=0;

    odor_lanes=[odor_lanes trials.odor_lane(trNo)];
    if (trials.odor_lane(trNo)==1)
        if (label_predicted(trNo)==0)
            this_accuracy=1;
            accuracy_lane1=[accuracy_lane1 1];
        else
            accuracy_lane1=[accuracy_lane1 0];
        end
    else
        if (trials.odor_lane(trNo)==4)
            if (label_predicted(trNo)==1)
                this_accuracy=1;
                accuracy_lane4=[accuracy_lane4 1];
            else
                accuracy_lane4=[accuracy_lane4 0];
            end
        end
    end
    accuracy(trNo)=this_accuracy;

     %Compute shifted accuracy
    this_accuracy=0;
    this_lane1_accuracy=0;

    odor_lanes=[odor_lanes trials.odor_lane(perm_trials(trNo))];
    if (trials.odor_lane(perm_trials(trNo))==1)
        if (label_predicted(perm_trials(trNo))==0)
            this_accuracy=1;
            sh_accuracy_lane1=[sh_accuracy_lane1 1];
        else
            sh_accuracy_lane1=[sh_accuracy_lane1 0];
        end
    else
        if (trials.odor_lane(perm_trials(trNo))==4)
            if (label_predicted(perm_trials(trNo))==1)
                this_accuracy=1;
                sh_accuracy_lane4=[sh_accuracy_lane4 1];
            else
                sh_accuracy_lane4=[sh_accuracy_lane4 0];
            end
        end
    end
    sh_accuracy(perm_trials(trNo))=this_accuracy;
end

handles_out.decode.repeat(1).label_predicted=label_predicted;
handles_out.decode.repeat(1).accuracy=accuracy;
handles_out.decode.repeat(1).accuracy_lane1=accuracy_lane1;
handles_out.decode.repeat(1).accuracy_lane4=accuracy_lane4;

handles_out.decode.repeat(1).sh_accuracy=accuracy;
handles_out.decode.repeat(1).sh_accuracy_lane1=accuracy_lane1;
handles_out.decode.repeat(1).sh_accuracy_lane4=accuracy_lane4;

%Now decode for the rest of the repeats
for ii_repeat=2:no_repeats
    ii_repeat
    label_pred=[];
    for trNo1=1:trials.odor_trNo
        % parfor trNo=1:trials.odor_trNo

        %Setup the training data set
        SpActrain=[];

        these_trNos=[];
        for trNo2=1:trials.odor_trNo
            if trNo2~=trNo1
                these_trNos=[these_trNos trNo2];
                this_row=[];
                for this_ROI=1:no_neurons
                    for x_ii=1:no_steps
                        for y_ii=1:no_steps
                            if ~isnan(activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii))
                                pffft=1;
                                if ~isnan(activity_maps.ROI(this_ROI).trial(trNo2).this_trial_map(x_ii,y_ii))
                                    this_row=[this_row activity_maps.ROI(this_ROI).trial(trNo2).this_trial_map(x_ii,y_ii)];
                                else
                                    this_row=[this_row activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii)];
                                end
                            end
                        end
                    end
                end

                SpActrain=[SpActrain; this_row];
            end
        end

        this_trNo=0;
        Y=[];
        for trNo3=these_trNos
            this_trNo=this_trNo+1;
            if trials.odor_lane(trNo3)==1
                Y(this_trNo)=0;
            else
                Y(this_trNo)=1;
            end
        end

        %Decode using neural network
        tblTrn=[];
        tblTrn = array2table(SpActrain);


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
                Mdl = fitglm(tblTrn,Y','Distribution','binomial');
            case 6
                Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
            case 7

                bestHyperparameters = Mdl_pred(trNo).Mdl.HyperparameterOptimizationResults.XAtMinEstimatedObjective;
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

        end


        %Calculate the prediction
        this_row=[];
        for this_ROI=1:no_neurons
            for x_ii=1:no_steps
                for y_ii=1:no_steps
                    if ~isnan(activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii))
                        pffft=1;
                        if ~isnan(activity_maps.ROI(this_ROI).trial(trNo1).this_trial_map(x_ii,y_ii))
                            this_row=[this_row activity_maps.ROI(this_ROI).trial(trNo1).this_trial_map(x_ii,y_ii)];
                        else
                            this_row=[this_row activity_maps.ROI(this_ROI).all_trial_map(x_ii,y_ii)];
                        end
                    end
                end
            end
        end
        label_pred(trNo1).data=predict(Mdl,this_row);

    end

    %Parse out the parfor loop output
    label_predicted=zeros(trials.odor_trNo,1);
    accuracy=zeros(trials.odor_trNo,1);
    accuracy_lane1=[];
    accuracy_lane4=[];
    odor_lanes=[];
    %         y_predicted=zeros(no_time_bins,1);
    for trNo=1:trials.odor_trNo
        label_predicted(trNo)=label_pred(trNo).data;
        %             y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

        %Compute accuracy
        this_accuracy=0;
        this_lane1_accuracy=0;

        odor_lanes=[odor_lanes trials.odor_lane(trNo)];
        if (trials.odor_lane(trNo)==1)
            if (label_predicted(trNo)==0)
                this_accuracy=1;
                accuracy_lane1=[accuracy_lane1 1];
            else
                accuracy_lane1=[accuracy_lane1 0];
            end
        else
            if (trials.odor_lane(trNo)==4)
                if (label_predicted(trNo)==1)
                    this_accuracy=1;
                    accuracy_lane4=[accuracy_lane4 1];
                else
                    accuracy_lane4=[accuracy_lane4 0];
                end
            end
        end
        accuracy(trNo)=this_accuracy;

        %Compute shifted accuracy
        this_accuracy=0;
        this_lane1_accuracy=0;

        odor_lanes=[odor_lanes trials.odor_lane(perm_trials(trNo))];
        if (trials.odor_lane(perm_trials(trNo))==1)
            if (label_predicted(perm_trials(trNo))==0)
                this_accuracy=1;
                sh_accuracy_lane1=[sh_accuracy_lane1 1];
            else
                sh_accuracy_lane1=[sh_accuracy_lane1 0];
            end
        else
            if (trials.odor_lane(perm_trials(trNo))==4)
                if (label_predicted(perm_trials(trNo))==1)
                    this_accuracy=1;
                    sh_accuracy_lane4=[sh_accuracy_lane4 1];
                else
                    sh_accuracy_lane4=[sh_accuracy_lane4 0];
                end
            end
        end
        sh_accuracy(perm_trials(trNo))=this_accuracy;
    end

    handles_out.decode.repeat(ii_repeat).label_predicted=label_predicted;
    handles_out.decode.repeat(ii_repeat).accuracy=accuracy;
    handles_out.decode.repeat(ii_repeat).accuracy_lane1=accuracy_lane1;
    handles_out.decode.repeat(ii_repeat).accuracy_lane4=accuracy_lane4;

    handles_out.decode.repeat(ii_repeat).sh_accuracy=sh_accuracy;
    handles_out.decode.repeat(ii_repeat).sh_accuracy_lane1=sh_accuracy_lane1;
    handles_out.decode.repeat(ii_repeat).sh_accuracy_lane4=sh_accuracy_lane4;

end


switch which_ROIs
    case 1
            fprintf(1,['\nDecoding performed with all ROIs\n\n'])
    case 2
            fprintf(1,['\nDecoding performed with place cells\n\n'])
    case 3
            fprintf(1,['\nDecoding performed with odor cells\n\n'])
end


percent_correct=100*(trials.hit4+trials.hit1)/(trials.hit4+trials.hit1+trials.miss4+trials.miss1);
percent_correct1=100*(trials.hit1)/(trials.hit1+trials.miss1);
percent_correct4=100*(trials.hit4)/(trials.hit4+trials.miss4);
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct1 ' num2str(percent_correct1) ' percent correct4 ' num2str(percent_correct4) '\n\n'])

accuracy=[];
accuracy_lane1=[];
accuracy_lane4=[];

for ii_repeat=1:no_repeats
    accuracy(ii_repeat)=mean(handles_out.decode.repeat(ii_repeat).accuracy);
    accuracy_lane1(ii_repeat)=mean(handles_out.decode.repeat(ii_repeat).accuracy_lane1);
    accuracy_lane4(ii_repeat)=mean(handles_out.decode.repeat(ii_repeat).accuracy_lane4);
end

mean_accuracy=mean(accuracy);
mean_accuracy_lane1=mean(accuracy_lane1);
mean_accuracy_lane4=mean(accuracy_lane4);
fprintf(1,['Overall accuracy ' num2str(mean_accuracy) ' lane 1 accuracy ' num2str(mean_accuracy_lane1)...
    ' lane 4 accuracy ' num2str(mean_accuracy_lane4) '\n\n'])

sh_accuracy=[];
sh_accuracy_lane1=[];
sh_accuracy_lane4=[];

for ii_repeat=1:no_repeats
    sh_accuracy(ii_repeat)=mean(handles_out.decode.repeat(ii_repeat).sh_accuracy);
    sh_accuracy_lane1(ii_repeat)=mean(handles_out.decode.repeat(ii_repeat).sh_accuracy_lane1);
    sh_accuracy_lane4(ii_repeat)=mean(handles_out.decode.repeat(ii_repeat).sh_accuracy_lane4);
end

mean_sh_accuracy=mean(sh_accuracy);
mean_sh_accuracy_lane1=mean(sh_accuracy_lane1);
mean_sh_accuracy_lane4=mean(sh_accuracy_lane4);
fprintf(1,['Overall shifted accuracy ' num2str(mean_sh_accuracy) ' lane 1 shifted accuracy ' num2str(mean_sh_accuracy_lane1)...
    ' lane 4 shifted accuracy ' num2str(mean_sh_accuracy_lane4) '\n\n'])
 

pffft=1;
