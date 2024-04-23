function handles_out=drgMini_dFFprediction(handles_choices)
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020

warning('off')

if exist('handles_choices')==0
    clear all
    close all


    %Files for dFF data
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';



    %File for odor plume
    %11282017_20cms_unbounded__README
    op_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Data_Plume/';
    op_file='dataset.mat';


%     %Second troubleshooting files
%     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
%     dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
%     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';

    %No odor troubleshooting files
%     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
%     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
%     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';


    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

    handles_choices.group=1; %1=odor plume, 2=spatial
    handles_choices.weber_fechner=1; %0 is Stevens Law, 1 is Weber-Flechner law
    handles_choices.alpha=1;
    handles_choices.repeats=1;
    handles_choices.sh_repeats=1;
    handles_choices.algo=2; %1 ann, 2 glm
    handles_choices.max_overlap=3; %Maximum overlap of the shuffled segments
    handles_choices.displayFigures=1;



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
    dt=0.2;
    dt_miniscope=1/30;
    n_shuffle=5; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    handles_choices.n_shuffle=n_shuffle;

    % which_training_algorithm=2; 
    % handles_choices.which_training_algorithm=which_training_algorithm;
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
    
    % which_training_algorithm=handles_choices.which_training_algorithm;

    op_path=handles_choices.op_path;
    op_file=handles_choices.op_file;

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

%Plot odor plume
load([op_path op_file])
x_range=[0 500]; %in mm
y_range=[0 480];
y_lane1=430;
y_lane4=50;

mean_plume=mean_plume-min(mean_plume(:));
mean_plume(mean_plume==0)=min(mean_plume(mean_plume~=0));
if handles_choices.weber_fechner==1
    mean_plume=log10(mean_plume);
end
 

if handles_choices.displayFigures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the timecourse
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.07 .5 .4 .4])

    drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),mean_plume')
    colormap fire
    shading interp
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';

    title('Mean odor plume')
end

%Now make mean_plume matrices for lanes 1 and 4
delta_x=x(2)-x(1); %please note this is the same as delta_y
x_for_plume=0:delta_x:500;
y_for_plume=0:delta_x:480;


[minabsy,minabsy_ii]=min(abs(y));

%lane 4
% mean_plume_lane4=prctile(mean_plume(:),1)*ones(length(y_for_plume),length(x_for_plume));
if handles_choices.weber_fechner==0
    mean_plume_lane4=zeros(length(y_for_plume),length(x_for_plume));
else
    mean_plume_lane4=prctile(mean_plume(:),1)*ones(length(y_for_plume),length(x_for_plume));
end


[minabsy4,minabsy4_ii]=min(abs(y_for_plume-50));

iiy4_to=minabsy4_ii+(length(y)-minabsy_ii);
delta_iiy4_to=(length(y)-minabsy_ii);
if iiy4_to>size(mean_plume_lane4,1)
    delta_iiy4_to=size(mean_plume_lane4,1)-minabsy4_ii;
    iiy4_to=size(mean_plume_lane4,1);
end

iiy4_from=minabsy4_ii-(length(y)-minabsy_ii);
delta_iiy4_from=-(length(y)-minabsy_ii);
if iiy4_from<1
    iiy4_from=1;
    delta_iiy4_from=-(minabsy4_ii-1);
end

mean_plume_lane4(iiy4_from:iiy4_to,1:size(mean_plume,2))=(mean_plume(minabsy_ii+delta_iiy4_from:minabsy_ii+delta_iiy4_to,:)).^handles_choices.alpha;


%lane 1
if handles_choices.weber_fechner==0
    mean_plume_lane1=zeros(length(y_for_plume),length(x_for_plume));
else
    mean_plume_lane1=prctile(mean_plume(:),1)*ones(length(y_for_plume),length(x_for_plume));
end


[minabsy1,minabsy1_ii]=min(abs(y_for_plume-430));

iiy1_to=minabsy1_ii+(length(y)-minabsy_ii);
delta_iiy1_to=(length(y)-minabsy_ii);
if iiy1_to>size(mean_plume_lane1,1)
    delta_iiy1_to=size(mean_plume_lane1,1)-minabsy1_ii;
    iiy1_to=size(mean_plume_lane1,1);
end

iiy1_from=minabsy1_ii-(length(y)-minabsy_ii);
delta_iiy1_from=-(length(y)-minabsy_ii);
if iiy1_from<1
    iiy1_from=1;
    delta_iiy1_from=minabsy1_ii-1;
end

mean_plume_lane1(iiy1_from:iiy1_to,1:size(mean_plume,2))=(mean_plume(minabsy_ii+delta_iiy1_from:minabsy_ii+delta_iiy1_to,:)).^handles_choices.alpha;

minC=min(mean_plume(:));
maxC=max(mean_plume(:));

if handles_choices.displayFigures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the timecourse
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.07 .5 .4 .4])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_lane4')
    colormap fire
    shading interp
    caxis([minC maxC]);
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';

    title('Mean odor plume shifted to lane 4')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the timecourse
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.07 .5 .4 .4])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_lane1')
    colormap fire
    shading interp
    caxis([minC maxC]);
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';

    title('Mean odor plume shifted to lane 1')
end

%Restart random seeds
rng('shuffle');

gcp;

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
    plot(trials.x_laneodor1,trials.y_laneodor1,'or')
    plot(trials.x_laneodor4,trials.y_laneodor4,'ob')
    plot([95 95],[200 250],'-r')
    plot([95 95],[0 50],'-b')

    plot(trials.x_lanewater1,trials.y_lanewater1,'xr')
    plot(trials.x_lanewater4,trials.y_lanewater4,'xb')
    xlabel('x')
    ylabel('y')
    title('Location of trial start (o) and water delivery (x)')
end

%The physical odor arena is 50 cm along the air flow and 48 cm wide,
%with lanes 1 and 4 5 cme from the wall
%We need to morph the video camera locations to the physical locaiton
%For now I will do this isometrically
y_mean_video_lane1=mean(trials.y_lanewater1);
y_mean_video_lane4=mean(trials.y_lanewater4);
y_mean_video_lanes14=mean([y_mean_video_lane1 y_mean_video_lane4]);
delta_y_physical_lane14=380;
y_half_physical=480/2;


%Bin positions into dt time bins
pos=[];
pos(:,1)=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(arena.xsync-95);
pos(:,2)=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(arena.ysync-y_mean_video_lanes14)+y_half_physical;
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

percent_correct=100*(trials.hit4+trials.hit1)/(trials.hit4+trials.hit1+trials.miss4+trials.miss1);
percent_correct1=100*(trials.hit1)/(trials.hit1+trials.miss1);
percent_correct4=100*(trials.hit4)/(trials.hit4+trials.miss4);
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct lane 1 ' num2str(percent_correct1) ' percent correct lane 4 ' num2str(percent_correct4) '\n\n'])

lane1_trials=trials.hit1+trials.miss1;
lane4_trials=trials.hit4+trials.miss4;

[p_hit_miss, Q]= chi2test([trials.hit1, trials.miss1; trials.hit4, trials.miss4]);
fprintf(1,['Chi squared testing difference in hit vs miss between lane 1 and 4 ' num2str(p_hit_miss) '\n\n'])

[p_hit_miss_lane1, Q]= chi2test([trials.hit1, trials.miss1; (trials.hit1+trials.miss1)/2, (trials.hit1+trials.miss1)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for lane 1 ' num2str(p_hit_miss_lane1) '\n\n'])

[p_hit_miss_lane4, Q]= chi2test([trials.hit4, trials.miss4; (trials.hit4+trials.miss4)/2, (trials.hit4+trials.miss4)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for lane 4 ' num2str(p_hit_miss_lane4) '\n\n'])

[p_hit_miss_both_lanes, Q]= chi2test([trials.hit1+trials.hit4, trials.miss1+trials.miss4; (trials.hit4+trials.hit1+trials.miss4+trials.miss1)/2, (trials.hit4+trials.hit1+trials.miss4+trials.miss1)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for both lanes ' num2str(p_hit_miss_both_lanes) '\n\n'])

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
% switch which_training_algorithm
%     case 1
%         %Train with entire session
%         ii_n_training=floor(1/(1-training_fraction));
%
%         for ii_train=1:ii_n_training
%             y_pred(ii_train).data=[];
%             x_pred(ii_train).data=[];
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
%             %     x_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
%             %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);
%
%             x_pred(ii_train).data=predict(MdlY1,XdFFtest);
%             y_pred(ii_train).data=predict(MdlY2,XdFFtest);
%         end
%         fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])
%
%         %Parse out the parfor loop output
%         x_predicted=zeros(no_time_bins,1);
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
%             x_predicted(logical(this_test_range),1)=x_pred(ii_train).data;
%             y_predicted(logical(this_test_range),1)=y_pred(ii_train).data;
%
%         end
%
%         %Now do predictions for reversed/permuted training periods
%         if n_shuffle>ii_n_training
%             n_shuffle=ii_n_training;
%         end
%         handles_choices.n_shuffle=n_shuffle;
%         x_predicted_sh=zeros(no_time_bins,n_shuffle);
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
%                 x_pred(ii_train).data=[];
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
%                 x_pred(ii_train).data=predict(MdlY1,XdFFtest);
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
%                 x_predicted_sh(logical(this_test_range),ii_shuffled)=x_pred(ii_train).data;
%                 y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(ii_train).data;
%             end
%
%         end

% case 2

%Train with per trial data using a leave one out approach
per_ROI=[];
for ii_repeat=1:handles_choices.repeats
    for this_ROI=1:no_neurons
        %training_range_template has all the trials
        training_range_template=zeros(1,no_time_bins);
        lane_template=zeros(1,no_time_bins);
        odor_plume_template=zeros(1,no_time_bins);
        for trNo=1:trials.odor_trNo
            dFF_pred(trNo).dFFpred_xy=[];
            % dFF_pred(trNo).dFFpred_xyl=[];
            % dFF_pred(trNo).dFFpred_xyop=[];
            dFF_pred(trNo).dFFpred_op=[];
            dFF_pred(trNo).dFF=[];
            % x_pred(trNo).data=[];

            x_predictedstart=trials.odor_ii_start(trNo)-10;
            x_predictedend=trials.odor_ii_end(trNo)+15;
            training_range_template(x_predictedstart:x_predictedend)=1;
            if trials.odor_lane(trNo)==1
                lane_template(x_predictedstart:x_predictedend)=1;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned_trimmed(x_ii,1);
                    this_y=pos_binned_trimmed(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane1(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            else
                lane_template(x_predictedstart:x_predictedend)=4;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned_trimmed(x_ii,1);
                    this_y=pos_binned_trimmed(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane4(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            end
        end

        % parfor trNo=1:trials.odor_trNo
        start_toc=toc;
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
            this_dFFtrain=zeros(sum(this_training_range),1);
            this_dFFtrain=neural_data_trimmed(this_training_range,this_ROI);
            % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
            this_dFFtest=zeros(sum(this_test_range),1);
            this_dFFtest(:,1)=neural_data_trimmed(logical(this_test_range),this_ROI);
            dFF_pred(trNo).dFF=this_dFFtest;

            XYtrain=pos_binned_trimmed(this_training_range,:);
            % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
            XYtest=pos_binned_trimmed(logical(this_test_range),:);

            lane_train=lane_template(this_training_range);
            lane_test=lane_template(logical(this_test_range));

            odor_plume_train=odor_plume_template(this_training_range);
            odor_plume_test=odor_plume_template(logical(this_test_range));

            %Decode using neural network

            %Decode using x and y only
            switch handles_choices.algo
                case 1
                    Mdl_dFFxy = fitrnet(XYtrain,this_dFFtrain);
                case 2
                    Mdl_dFFxy = fitglm(XYtrain,this_dFFtrain);
            end
            dFF_pred(trNo).dFFpred_xy=predict(Mdl_dFFxy,XYtest);

            % %Decode using x, y and odor lane
            % switch handles_choices.algo
            %     case 1
            %         Mdl_dFFxyl = fitrnet([XYtrain lane_train'],this_dFFtrain);
            %     case 2
            %         Mdl_dFFxyl = fitglm([XYtrain lane_train'],this_dFFtrain);
            % end
            % dFF_pred(trNo).dFFpred_xyl=predict(Mdl_dFFxyl,[XYtest lane_test']);
            % 
            % %Decode using x, y and odor plume
            % switch handles_choices.algo
            %     case 1
            %         Mdl_dFFxyop = fitrnet([XYtrain odor_plume_train'],this_dFFtrain);
            %     case 2
            %         Mdl_dFFxyop = fitglm([XYtrain odor_plume_train'],this_dFFtrain);
            % end
            % dFF_pred(trNo).dFFpred_xyop=predict(Mdl_dFFxyop,[XYtest odor_plume_test']);

            %Decode using odor plume only
            switch handles_choices.algo
                case 1
                    Mdl_dFFop = fitrnet([odor_plume_train'],this_dFFtrain);
                case 2
                    Mdl_dFFop = fitglm([odor_plume_train'],this_dFFtrain);
            end
            dFF_pred(trNo).dFFpred_op=predict(Mdl_dFFop,[odor_plume_test']);

            pffft=1;
            % y_pred(trNo).data=predict(MdlY2,XdFFtest);
        end
        fprintf(1,['Elapsed time ' num2str((toc-start_toc)/(60)) ' mins\n\n'])

        %Parse out the parfor loop output
        all_dFFpred_xy=[];
        % all_dFFpred_xyl=[];
        % all_dFFpred_xyop=[];
        all_dFFpred_op=[];
        all_dFF=[];
        for trNo=1:trials.odor_trNo
            per_ROI(this_ROI).repeats(ii_repeat).trial(trNo).dFFpred_xy=dFF_pred(trNo).dFFpred_xy;
            all_dFFpred_xy=[all_dFFpred_xy dFF_pred(trNo).dFFpred_xy'];
            % per_ROI(this_ROI).repeats(ii_repeat).trial(trNo).dFFpred_xyl=dFF_pred(trNo).dFFpred_xyl;
            % all_dFFpred_xyl=[all_dFFpred_xyl dFF_pred(trNo).dFFpred_xyl'];
            % all_dFFpred_xyop=[all_dFFpred_xyop dFF_pred(trNo).dFFpred_xyop'];
            all_dFFpred_op=[all_dFFpred_op dFF_pred(trNo).dFFpred_op'];
            per_ROI(this_ROI).repeats(ii_repeat).trial(trNo).dFF=dFF_pred(trNo).dFF;
            all_dFF=[all_dFF dFF_pred(trNo).dFF'];
        end

        per_ROI(this_ROI).repeats(ii_repeat).all_dFF=all_dFF;
        per_ROI(this_ROI).repeats(ii_repeat).all_dFFpred_xy=all_dFFpred_xy;
        % per_ROI(this_ROI).repeats(ii_repeat).all_dFFpred_xyl=all_dFFpred_xyl;
        % per_ROI(this_ROI).repeats(ii_repeat).all_dFFpred_xyop=all_dFFpred_xyop;
        per_ROI(this_ROI).repeats(ii_repeat).all_dFFpred_op=all_dFFpred_op;

        thisRxy=corrcoef([all_dFF' all_dFFpred_xy']);
        per_ROI(this_ROI).repeats(ii_repeat).Rxy=thisRxy(1,2);
        % thisRxyl=corrcoef([all_dFF' all_dFFpred_xyl']);
        % per_ROI(this_ROI).repeats(ii_repeat).Rxyl=thisRxyl(1,2);
        % thisRxyop=corrcoef([all_dFF' all_dFFpred_xyop']);
        % per_ROI(this_ROI).repeats(ii_repeat).Rxyop=thisRxyop(1,2);
        thisRop=corrcoef([all_dFF' all_dFFpred_op']);
        per_ROI(this_ROI).repeats(ii_repeat).Rop=thisRop(1,2);
        % fprintf(1, 'Repeat %d ROI No %d, rho for prediction dFFxy, dFFxyl, dFFxyop, dFFop %d %d %d %d\n\n'...
        %     ,ii_repeat,this_ROI, per_ROI(this_ROI).repeats(ii_repeat).Rxy,per_ROI(this_ROI).repeats(ii_repeat).Rxyl,per_ROI(this_ROI).repeats(ii_repeat).Rxyop,per_ROI(this_ROI).repeats(ii_repeat).Rop);
        fprintf(1, 'Repeat %d ROI No %d, rho for prediction dFFxy, dFFop %d %d\n\n'...
            ,ii_repeat,this_ROI, per_ROI(this_ROI).repeats(ii_repeat).Rxy,per_ROI(this_ROI).repeats(ii_repeat).Rop);

        % if (thisRxy(1,2)>0.5)||(thisRxyl(1,2)>0.5)
        %     pffft=1;
        % end

    end
end

%Train with per trial with shuffled data using a leave one out approach
per_ROI_sh=[];
for ii_repeat=1:handles_choices.sh_repeats
    for this_ROI=1:no_neurons
        %training_range_template has all the trials
        training_range_template=zeros(1,no_time_bins);
        lane_template=zeros(1,no_time_bins);
        odor_plume_template=zeros(1,no_time_bins);
        for trNo=1:trials.odor_trNo
            dFF_pred(trNo).dFFpred_xy=[];
            % dFF_pred(trNo).dFFpred_xyl=[];
            dFF_pred(trNo).dFF=[];
            % x_pred(trNo).data=[];

            x_predictedstart=trials.odor_ii_start(trNo)-10;
            x_predictedend=trials.odor_ii_end(trNo)+15;
            training_range_template(x_predictedstart:x_predictedend)=1;
            if trials.odor_lane(trNo)==1
                lane_template(x_predictedstart:x_predictedend)=1;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned_trimmed(x_ii,1);
                    this_y=pos_binned_trimmed(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane1(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            else
                lane_template(x_predictedstart:x_predictedend)=4;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned_trimmed(x_ii,1);
                    this_y=pos_binned_trimmed(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane4(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            end
        end

        % parfor trNo=1:trials.odor_trNo
        start_toc=toc;
        max_overlap=handles_choices.max_overlap;
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
            this_dFFtrain=zeros(sum(this_training_range),1);
            this_dFFtrain=neural_data_trimmed(this_training_range,this_ROI);
            % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
            this_dFFtest=zeros(sum(this_test_range),1);
            this_dFFtest(:,1)=neural_data_trimmed(logical(this_test_range),this_ROI);
            dFF_pred(trNo).dFF=this_dFFtest;

            %Now shuffle
            XYtrain_pre=pos_binned_trimmed(this_training_range,:);
            shuffle_length=floor(size(XYtrain_pre,1)/(trials.odor_trNo*2));

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

            % Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
            XYtest=pos_binned_trimmed(logical(this_test_range),:);

            lane_train_pre=lane_template(this_training_range);
            lane_train_pre=lane_train_pre';
            lane_test=lane_template(logical(this_test_range));

            odor_plume_train_pre=odor_plume_template(this_training_range);
            odor_plume_train_pre=odor_plume_train_pre';
            odor_plume_test=odor_plume_template(logical(this_test_range));

            XYtrain=zeros(size(XYtrain_pre,1),size(XYtrain_pre,2));
            lane_train=zeros(size(XYtrain_pre,1),1);
            odor_plume_train=zeros(size(XYtrain_pre,1),1);
            for sh_ii=1:trials.odor_trNo*2
                XYtrain((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
                    XYtrain_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
                lane_train((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
                    lane_train_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
                odor_plume_train((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
                    odor_plume_train_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
            end

            %Make the last few equal to the first few
            if trials.odor_trNo*2*shuffle_length<size(XYtrain,1)
                delta_ii=size(XYtrain,1)-trials.odor_trNo*2*shuffle_length;
                XYtrain(trials.odor_trNo*2*shuffle_length+1:end,:)=...
                    XYtrain_pre(1:delta_ii,:);
                lane_train(trials.odor_trNo*2*shuffle_length+1:end,:)=...
                    lane_train_pre(1:delta_ii,:);
                odor_plume_train(trials.odor_trNo*2*shuffle_length+1:end,:)=...
                    odor_plume_train_pre(1:delta_ii,:);
            end


            %Decode using neural network

            %Decode using x and y only
            switch handles_choices.algo
                case 1
                    Mdl_dFFxy = fitrnet(XYtrain,this_dFFtrain');
                case 2
                    Mdl_dFFxy = fitglm(XYtrain,this_dFFtrain');
            end
            dFF_pred(trNo).dFFpred_xy=predict(Mdl_dFFxy,XYtest);

            % %Decode using x, y and odor lane
            % switch handles_choices.algo
            %     case 1
            %         Mdl_dFFxyl = fitrnet([XYtrain lane_train],this_dFFtrain');
            %     case 2
            %         Mdl_dFFxyl = fitglm([XYtrain lane_train],this_dFFtrain');
            % end
            % dFF_pred(trNo).dFFpred_xyl=predict(Mdl_dFFxyl,[XYtest lane_test']);
            % 
            % %Decode using x, y and odor plume
            % switch handles_choices.algo
            %     case 1
            %         Mdl_dFFxyop = fitrnet([XYtrain odor_plume_train],this_dFFtrain);
            %     case 2
            %         Mdl_dFFxyop = fitglm([XYtrain odor_plume_train],this_dFFtrain);
            % end
            % dFF_pred(trNo).dFFpred_xyop=predict(Mdl_dFFxyop,[XYtest odor_plume_test']);

            %Decode using odor plume
            switch handles_choices.algo
                case 1
                    Mdl_dFFop = fitrnet([odor_plume_train],this_dFFtrain);
                case 2
                    Mdl_dFFop = fitglm([odor_plume_train],this_dFFtrain);
            end
            dFF_pred(trNo).dFFpred_op=predict(Mdl_dFFop,[odor_plume_test']);


            pffft=1;
            % y_pred(trNo).data=predict(MdlY2,XdFFtest);
        end
        fprintf(1,['Elapsed time ' num2str((toc-start_toc)/(60)) ' mins\n\n'])

        %Parse out the parfor loop output
        all_dFFpred_xy=[];
        all_dFFpred_xyl=[];
        all_dFFpred_xyop=[];
        all_dFFpred_op=[];
        all_dFF=[];
        for trNo=1:trials.odor_trNo
            per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).dFFpred_xy=dFF_pred(trNo).dFFpred_xy;
            all_dFFpred_xy=[all_dFFpred_xy dFF_pred(trNo).dFFpred_xy'];
            % per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).dFFpred_xyl=dFF_pred(trNo).dFFpred_xyl;
            % all_dFFpred_xyl=[all_dFFpred_xyl dFF_pred(trNo).dFFpred_xyl'];
            % all_dFFpred_xyop=[all_dFFpred_xyop dFF_pred(trNo).dFFpred_xyop'];
            all_dFFpred_op=[all_dFFpred_op dFF_pred(trNo).dFFpred_op'];
            per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).dFF=dFF_pred(trNo).dFF;
            all_dFF=[all_dFF dFF_pred(trNo).dFF'];
        end

        per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFF=all_dFF;
        per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_xy=all_dFFpred_xy;
        % per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_xyl=all_dFFpred_xyl;
        % per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_xyop=all_dFFpred_xyop;
        per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_op=all_dFFpred_op;

        thisRxy=corrcoef([all_dFF' all_dFFpred_xy']);
        per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy=thisRxy(1,2);
        % thisRxyl=corrcoef([all_dFF' all_dFFpred_xyl']);
        % per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyl=thisRxyl(1,2);
        % thisRxyop=corrcoef([all_dFF' all_dFFpred_xyop']);
        % per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop=thisRxyop(1,2);
        thisRop=corrcoef([all_dFF' all_dFFpred_op']);
        per_ROI_sh(this_ROI).repeats(ii_repeat).Rop=thisRop(1,2);
        % fprintf(1, 'Repeat %d ROI No %d, rho for shuffled prediction dFFxy, dFFxyl dFFxyop, dFFop %d %d %d %d\n\n',ii_repeat,...
        %     this_ROI,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyl,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop,per_ROI_sh(this_ROI).repeats(ii_repeat).Rop);
                fprintf(1, 'Repeat %d ROI No %d, rho for shuffled prediction dFFxy, dFFop %d %d \n\n',ii_repeat,...
            this_ROI,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rop);
        % if (thisRxy(1,2)>0.5)||(thisRxyl(1,2)>0.5)
        %     pffft=1;
        % end

    end
end

%         pffft=1;
%     case 3
%         %Train between trials using a leave one out approach
%
%         %training_range_template has all the trials
%         training_range_template=zeros(1,no_time_bins);
%         for trNo=1:trials.odor_trNo
%             y_pred(trNo).data=[];
%             x_pred(trNo).data=[];
%
%             if trNo==1
%                last_x_predictedend=1;
%             else
%             last_x_predictedend=trials.odor_ii_end(trNo-1)+15;
%             end
%
%
%             this_x_predictedstart=trials.odor_ii_start(trNo)-10;
%
%
%             training_range_template(last_x_predictedend:this_x_predictedstart)=1;
%         end
%
%         last_x_predictedend=trials.odor_ii_end(trials.odor_trNo)+15;
%         this_x_predictedstart=no_time_bins;
%         training_range_template(last_x_predictedend:this_x_predictedstart)=1;
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
%             %     x_predicted(logical(this_test_range),1)=predict(MdlY1,XdFFtest);
%             %     y_predicted(logical(this_test_range),1)=predict(MdlY2,XdFFtest);
%
%             x_pred(trNo).data=predict(MdlY1,XdFFtest);
%             y_pred(trNo).data=predict(MdlY2,XdFFtest);
%         end
%         fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])
%
%         %Parse out the parfor loop output
%         x_predicted=zeros(no_time_bins,1);
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
%             x_predicted(logical(this_test_range),1)=x_pred(trNo).data;
%             y_predicted(logical(this_test_range),1)=y_pred(trNo).data;
%
%         end
%
%         %Now do predictions for reversed/permuted training periods
%         x_predicted_sh=zeros(no_time_bins,n_shuffle);
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
%                 x_pred(trNo).data=[];
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
%                 x_pred(trNo).data=predict(MdlY1,XdFFtest);
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
%                 x_predicted_sh(logical(this_test_range),ii_shuffled)=x_pred(trNo).data;
%                 y_predicted_sh(logical(this_test_range),ii_shuffled)=y_pred(trNo).data;
%             end
%
%         end
%
% end

fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

if handles_choices.displayFigures==1
    % figNo=figNo+1;
    % try
    %     close(figNo)
    % catch
    % end
    % 
    % hFig = figure(figNo);
    % 
    % set(hFig, 'units','normalized','position',[.1 .1 .4 .4])
    % 
    % 
    % hold on
    % 
    % %Plot the xy vs xy rhos
    % for ii_repeat1=1:handles_choices.repeats
    %     for ii_repeat2=ii_repeat1+1:handles_choices.repeats
    %         for this_ROI=1:no_neurons
    %             thisRxy1=corrcoef([per_ROI(this_ROI).repeats(ii_repeat1).all_dFF' per_ROI(this_ROI).repeats(ii_repeat1).all_dFFpred_xy']);
    %             thisRxy2=corrcoef([per_ROI(this_ROI).repeats(ii_repeat2).all_dFF' per_ROI(this_ROI).repeats(ii_repeat2).all_dFFpred_xy']);
    %             plot(thisRxy1(1,2),thisRxy2(1,2),'.k','MarkerSize',8)
    %         end
    %     end
    % end
    % 
    % %Plot the xy vs xyop rhos
    % for ii_repeat=1:handles_choices.repeats
    %     for this_ROI=1:no_neurons
    %         plot(per_ROI(this_ROI).repeats(ii_repeat).Rxy,per_ROI(this_ROI).repeats(ii_repeat).Rxyl,'.b','MarkerSize',8)
    %         plot(per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyl,'.r','MarkerSize',8)
    %     end
    % end
    % 
    % 
    % 
    % this_ylim=ylim;
    % this_xlim=xlim;
    % 
    % plot([-0.4 0.8],[-0.4 0.8],'-k')
    % plot([0 0],[-0.4 0.8],'-k')
    % plot([-0.4 0.8],[0 0],'-k')
    % xlim([-0.4 0.8])
    % ylim([-0.4 0.8])
    % title(['rho ANN decoding, xy vs. xy + lane, log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha)])
    % xlabel('rho for predicted dFF computed with x,y')
    % ylabel('rho for predicted dFF computed with x,y and odor lane')
    % 
    % 
    % 
    % figNo=figNo+1;
    % try
    %     close(figNo)
    % catch
    % end
    % 
    % hFig = figure(figNo);
    % 
    % set(hFig, 'units','normalized','position',[.1 .1 .4 .4])
    % 
    % 
    % hold on
    % 
    % %Plot the xy vs xy rhos
    % for ii_repeat1=1:handles_choices.repeats
    %     for ii_repeat2=ii_repeat1+1:handles_choices.repeats
    %         for this_ROI=1:no_neurons
    %             thisRxy1=corrcoef([per_ROI(this_ROI).repeats(ii_repeat1).all_dFF' per_ROI(this_ROI).repeats(ii_repeat1).all_dFFpred_xy']);
    %             thisRxy2=corrcoef([per_ROI(this_ROI).repeats(ii_repeat2).all_dFF' per_ROI(this_ROI).repeats(ii_repeat2).all_dFFpred_xy']);
    %             plot(thisRxy1(1,2),thisRxy2(1,2),'.k','MarkerSize',8)
    %         end
    %     end
    % end
    % 
    % %Plot the xy vs xyop rhos
    % for ii_repeat=1:handles_choices.sh_repeats
    %     for this_ROI=1:no_neurons
    %         plot(per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop,'.r','MarkerSize',8)
    %     end
    % end
    % 
    % 
    % for ii_repeat=1:handles_choices.repeats
    %     for this_ROI=1:no_neurons
    %         plot(per_ROI(this_ROI).repeats(ii_repeat).Rxy,per_ROI(this_ROI).repeats(ii_repeat).Rxyop,'.b','MarkerSize',8)
    %     end
    % end
    % 
    % 
    % 
    % 
    % 
    % this_ylim=ylim;
    % this_xlim=xlim;
    % 
    % plot([-0.4 0.8],[-0.4 0.8],'-k')
    % plot([0 0],[-0.4 0.8],'-k')
    % plot([-0.4 0.8],[0 0],'-k')
    % xlim([-0.4 0.8])
    % ylim([-0.4 0.8])
    % title(['rho ANN decoding, xy vs. xy + op, log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha)])
    % xlabel('rho for predicted dFF computed with x,y')
    % ylabel('rho for predicted dFF computed with x,y and odor plume')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .4 .4])


    hold on

    for ii_repeat1=1:handles_choices.sh_repeats
        for this_ROI=1:no_neurons
            plot(per_ROI_sh(this_ROI).repeats(ii_repeat1).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat1).Rop,'.r','MarkerSize',8)
        end
    end

    for ii_repeat1=1:handles_choices.repeats
        for this_ROI=1:no_neurons
            plot(per_ROI(this_ROI).repeats(ii_repeat1).Rxy,per_ROI(this_ROI).repeats(ii_repeat1).Rop,'.b','MarkerSize',8)
        end
    end


    this_ylim=ylim;
    this_xlim=xlim;

    plot([-0.4 0.8],[-0.4 0.8],'-k')
    plot([0 0],[-0.4 0.8],'-k')
    plot([-0.4 0.8],[0 0],'-k')
    xlim([-0.4 0.8])
    ylim([-0.4 0.8])
    title(['rho ANN decoding, xy vs. op, log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha)])
    xlabel('rho for predicted dFF computed with x,y')
    ylabel('rho for predicted dFF computed with odor plume')
end
    
handles_out.per_ROI=per_ROI;
handles_out.per_ROI_sh=per_ROI_sh;
handles_out.trials=trials;
save([this_path arena_file(1:end-4) 'decdFFa' num2str(handles_choices.algo) ...
    'wf' num2str(handles_choices.weber_fechner) 'g' num2str(handles_choices.group) ...
    'a' num2str(handles_choices.alpha) '.mat'],'handles_out','handles_choices','-v7.3')
 
%Now plot pseudocolor maps of activity


pffft=1;