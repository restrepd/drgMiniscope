%drgMini_plot_path_from_angle_file
%This file is for Kira to use to find the start and end points for
%each trial

close all
clear all


%20220727_FCM19 file 3 to troubleshoot
% this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220727_FCM19/';
% % dFF_file='20221117_FCM22_withodor_nearfloor_miniscope_sync_L1andL4_ncorre_fix_ext.mat';
% arena_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm.mat';
% angle_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm_aangle.mat';


%20221117_FCM22_lanes_1_4 for first figure in manuscript
this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20221117_FCM22_lanes_1_4/';
% dFF_file='20221117_FCM22_withodor_nearfloor_miniscope_sync_L1andL4_ncorre_fix_ext.mat';
arena_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm.mat';
angle_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm_aangle.mat';

% %20221117_FCM22_lanes_1_4 for first figure in manuscript
% this_path='/data2/SFTP/20221117_FCM22_lanes_1_4/';
% % dFF_file='20221117_FCM22_withodor_nearfloor_miniscope_sync_L1andL4_ncorre_fix_ext.mat';
% arena_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm.mat';
% angle_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm_aangle.mat';

%definition of variables
dt=0.1; %Time bins for decoding, this was 0.2
air_flow_speed=50; %mm/sec
dt_miniscope=1/30;
figNo=0;
trial_start_offset=-15;
trial_end_offset=15;
 
%This loads the start and end times
load([this_path angle_file])

trials=handles_out.trials;
angles=handles_out.angles;


load([this_path arena_file])


%Bin positions into dt time bins
pos=[];

%x and y corrected to mm 
pos(:,1)=arena.xsync;
pos(:,2)=arena.ysync;

%Note: I think arena.xsync_orig and arena.ysync_orig are in pixels

no_time_points=size(pos,1);

%I take the mean of x and y for dt (0.1 sec)
dFF_times=[1:no_time_points]*dt_miniscope;
no_time_bins=round(dFF_times(end)/dt);
time_binned=[1:no_time_bins]*dt-dt/2;
pos_binned=zeros(no_time_bins,2);

for ii_time_bin=1:no_time_bins
    time_from=time_binned(ii_time_bin)-dt/2;
    time_to=time_binned(ii_time_bin)+dt/2;
    pos_binned(ii_time_bin,:)=mean(pos((dFF_times>=time_from)&(dFF_times<time_to),:),1);
end

trim_factor=no_time_bins/no_time_points;

%Now plot per trial trajectories
mean_end_angles=[];
ii_first_odor_encounter=[];
dt_first_odor_encounter=[];
dt_turn_to_odor_encounter=[];
for trNo=1:trials.odor_trNo
    
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


        hold on
    

    ii_predictedstart=trials.odor_ii_start(trNo)+trial_start_offset;
    ii_predictedend=trials.odor_ii_end(trNo)+trial_end_offset;

    if ii_predictedend>no_time_bins
        ii_predictedend=no_time_bins;
    end

    ii_trialstart=trials.odor_ii_start(trNo);
    ii_trialend=trials.odor_ii_end(trNo);

    x = pos_binned(ii_predictedstart:ii_predictedend,1);
    y = pos_binned(ii_predictedstart:ii_predictedend,2);

    %find first odor encounter
    ii_delta_start=1-trial_start_offset;
    odor_detected=0;
    ii_first_odor_encounter(trNo)=ii_delta_start;
    dt_first_odor_encounter(trNo)=ii_delta_start*dt;
    for ii=ii_delta_start+1:length(x)
        x_on=ii*dt*air_flow_speed;
        if (odor_detected==0)&(x(ii)<=x_on)
            odor_detected=1;
            ii_first_odor_encounter(trNo)=ii;
            dt_first_odor_encounter(trNo)=ii*dt;
        end
    end


    % plot(arena.xsync(ii_predictedstart*(53707/17902):ii_predictedend*(53707/17902)),arena.ysync(ii_predictedstart*(53707/17902):ii_predictedend*(53707/17902)),'-k')
    % set(gca, 'YDir', 'reverse');

    %Okabe_Ito colors
    
        switch trials.odor_trial_type(trNo)
            case 1
                %Lane 1 hits
                plot(pos_binned(ii_predictedstart:ii_predictedend,1),pos_binned(ii_predictedstart:ii_predictedend,2),'Color',[213/255 94/255 0],'LineWidth',3)
                title(['Lane 1 hit, trial number ' num2str(trNo)])
            case 2
                %Lane 1 miss
                plot(pos_binned(ii_predictedstart:ii_predictedend,1),pos_binned(ii_predictedstart:ii_predictedend,2),'Color',[230/255 159/255 0],'LineWidth',3)
                title(['Lane 1 miss, trial number ' num2str(trNo)])
            case 3
                %Lane 4 hit
                plot(pos_binned(ii_predictedstart:ii_predictedend,1),pos_binned(ii_predictedstart:ii_predictedend,2),'Color',[0 114/255 178/255],'LineWidth',3)
                title(['Lane 4 hit, trial number ' num2str(trNo)])
            case 4
                %Lane 4 miss
                plot(pos_binned(ii_predictedstart:ii_predictedend,1),pos_binned(ii_predictedstart:ii_predictedend,2),'Color',[86/255 180/255 233/255],'LineWidth',3)
                title(['Lane 4 miss, trial number ' num2str(trNo)])
        end

        plot(pos_binned(ii_trialstart,1),pos_binned(ii_trialstart,2),'ob','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b')
        plot(pos_binned(ii_trialend,1),pos_binned(ii_trialend,2),'or','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')
    


    %Find the last turn
    ii_turns=length(angles.trial(trNo).delta_x);
    last_turn_ii=1;
    last_turn_found=0;
    for ii_t=1:ii_turns
        if angles.trial(trNo).delta_x(ii_t)>100
            this_ii_turn=angles.trial(trNo).ii_turns(ii_t);
            last_turn_ii=this_ii_turn;
            last_turn_found=1;
        end
    end

    if last_turn_found==1
        plot(x(last_turn_ii),y(last_turn_ii),'or','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
    end
    dt_turn_to_odor_encounter(trNo)=(ii_first_odor_encounter(trNo)-last_turn_ii)*dt;

    plot(x(ii_first_odor_encounter(trNo)),y(ii_first_odor_encounter(trNo)),'oy','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')

    text(100,100,['End angle ' num2str(angles.trial(trNo).end_angle)])
    text(100,130,['Mean end angle ' num2str(angles.trial(trNo).mean_end_angle)])
    xlabel('x')
    ylabel('y')
    set(gca, 'YDir', 'reverse');
    xlim([0 500])
    ylim([0 480])

    yticks([50 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})
 

    if (trNo==28)||(trNo==30)
        fprintf(1,['\n\nTrial number ' num2str(trNo) '\n'])
        fprintf(1,['Frame number for start of movie (5 sec before trial start)  ' num2str((ii_trialstart-5*10)*3) '\n'])
        fprintf(1,['Frame number for trial start  ' num2str(ii_trialstart*3) '\n'])
        fprintf(1,['Frame number for last turn  ' num2str((ii_predictedstart+last_turn_ii-1)*3) '\n'])
        fprintf(1,['Frame number for odor encounter  ' num2str((ii_predictedstart+ii_first_odor_encounter(trNo)-1)*3) '\n'])
        fprintf(1,['Frame number end  ' num2str(ii_trialend*3) '\n'])
    end

    pffft1=1;

end

pffft=1;