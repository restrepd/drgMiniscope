%drgMini_analyze_batch_PathAngle
close all
clear all

is_sphgpu=0; %0=Mac

switch is_sphgpu
    case 0
        %Mac

        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';

        %Trained with hits only
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

        %Trained with hits only taking on account when mouse detects the odor
        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc04192024/';
        choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_04192024.m'

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        %This was used initially
        % save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        % choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        %This one has the odor encounter
        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle05152025/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_05102025.m';

        % choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
        fileID = fopen([choiceBatchPathName 'pathAngle_stats.txt'],'w');

        addpath('/Users/restrepd/Dropbox/m new/Chi Squared')
    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end

%Speed for background air flow in mm/sec
air_flow_speed=50; %5 cm/sec = 50 mm/sec

display_dt=[-1.5 1.5]; %The speed will be displayed in this span (seconds)

%These are used to classify angle approach
low_angle=-130;
high_angle=-50;

y_lane1=430;
y_lane4=70;
hit_radius=100;

figNo=0;

%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4, 2 cm from floor
%Group 2 is rewarded, with odor lane 4, no odor in lane 1
%Group 3 is rewarded, with odor lane 1, no odor in lane 4
%Group 4 is rewarded, with no odor in lane 1 and lane 4
%Group 5 is rewarded, with ISO1 in both lane 1 and lane 4, 1 cm from floor

group_label{1}='Odor both lanes 2 cm';
group_label{2}='Odor lane 4';
group_label{3}='Odor lane 1';
group_label{4}='No odor';
group_label{5}='Odor both lanes 1 cm';

run_label{1}='0 bins before';
run_label{5}='16 bins before';

% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])

dt=handles_conc.dt;

figNo=0;


ii_run=1;

%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4, 2 cm from floor
%Group 2 is rewarded, with odor lane 4, no odor in lane 1
%Group 3 is rewarded, with odor lane 1, no odor in lane 4
%Group 4 is rewarded, with no odor in lane 1 and lane 4
%Group 5 is rewarded, with ISO1 in both lane 1 and lane 4, 1 cm from floor
these_groups=[1 5];
% these_groups=[4];
% these_groups=[2 3];

%Determine which files will be included
%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc,1);

%Calculate percent correct behavior
behavior_analysis=[];
behavior_analysis.files_included=[];
percent_correct=[];
time_set=0;
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
        if time_set==0
            display_t=[display_dt(1):handles_XY.dt:display_dt(2)];
            time_set=1;
        end
        behavior_analysis.files_included=[behavior_analysis.files_included fileNo];
        %Load conc data
        arena_file=handles_conc.arena_file{fileNo};
        %load the ouptut file
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;

        behavior_analysis.file(fileNo).percent_correct=100*(sum(trials.hit4)+sum(trials.hit1))/length(trials.hit4);
        percent_correct=[percent_correct behavior_analysis.file(fileNo).percent_correct];
        behavior_analysis.file(fileNo).percent_correct1=100*(sum(trials.hit1))/(sum(trials.hit1)+sum(trials.miss1));
        behavior_analysis.file(fileNo).percent_correct4=100*(sum(trials.hit4))/(sum(trials.hit4)+sum(trials.miss4));
        fprintf(1,['\nFor file No ' num2str(fileNo) '\n'])
        fprintf(1,['percent correct ' num2str(behavior_analysis.file(fileNo).percent_correct) '\n'])
        fprintf(1,['percent correct lane 1 ' num2str(behavior_analysis.file(fileNo).percent_correct1) '\n'])
        fprintf(1,['percent correct lane 4 ' num2str(behavior_analysis.file(fileNo).percent_correct4) '\n'])


        [p_hit_miss_between_lanes, Q]= chi2test([sum(trials.hit1), sum(trials.miss1); sum(trials.hit4), sum(trials.miss4)]);
        fprintf(1,['Chi squared testing difference in hit vs miss between lane 1 and 4 ' num2str(p_hit_miss_between_lanes) '\n'])

        [p_hit_miss_lane1_vs_random, Q]= chi2test([sum(trials.hit1), sum(trials.miss1); (sum(trials.hit1)+sum(trials.miss1))/2, (sum(trials.hit1)+sum(trials.miss1))/2]);
        fprintf(1,['Chi squared testing hit/miss vs random for lane 1 ' num2str(p_hit_miss_lane1_vs_random) '\n'])

        [p_hit_miss_lane4_vs_random, Q]= chi2test([sum(trials.hit4), sum(trials.miss4); (sum(trials.hit4)+sum(trials.miss4))/2, (sum(trials.hit4)+sum(trials.miss4))/2]);
        fprintf(1,['Chi squared testing hit/miss vs random for lane 4 ' num2str(p_hit_miss_lane4_vs_random) '\n'])

        [p_hit_miss_both_lanes_vs_random, Q]= chi2test([sum(trials.hit1)+sum(trials.hit4), sum(trials.miss1)+sum(trials.miss4); length(trials.hit4)/2, length(trials.hit4)/2]);
        fprintf(1,['Chi squared testing hit/miss vs random for both lanes ' num2str(p_hit_miss_both_lanes_vs_random) '\n'])

        fprintf(1,['\n\n'])

    end
end


%Now plot pseudocolor locations for last turn angles, start and end points,
%etc

%First make point maps
ii_for_corr=0;

y_length=480;
x_length=500;

y_values=(y_length/20):y_length/10:y_length-(y_length/20);
x_values=(x_length/20):x_length/10:x_length-(x_length/20);

dt_turn_to_odor_lane1=[];
dt_turn_to_odor_lane4=[];

distance_turn_to_odor_lane1=[];
distance_turn_to_odor_lane4=[];

lane1_hit_turn_angle_positions=zeros(10,10);
lane1_miss_turn_angle_positions=zeros(10,10);

lane4_hit_turn_angle_positions=zeros(10,10);
lane4_miss_turn_angle_positions=zeros(10,10);

lane1_hit_start_positions=zeros(10,10);
lane1_miss_start_positions=zeros(10,10);

lane4_hit_start_positions=zeros(10,10);
lane4_miss_start_positions=zeros(10,10);

lane1_hit_end_positions=zeros(10,10);
lane1_miss_end_positions=zeros(10,10);

lane4_hit_end_positions=zeros(10,10);
lane4_miss_end_positions=zeros(10,10);

dt_spout_hits=[];
dt_spout_miss=[];
per_session_dt_spout_hits=[];
per_session_dt_spout_miss=[];

speed_before_turn_hits=[];
speed_after_turn_hits=[];
speed_before_turn_misses=[];
speed_after_turn_misses=[];
average_speed=[];

normalized_speed_before_turn_hits=[];
normalized_speed_after_turn_hits=[];
normalized_speed_before_turn_misses=[];
normalized_speed_after_turn_misses=[];

per_session_normalized_speed_before_turn_hits=[];
per_session_normalized_speed_after_turn_hits=[];
per_session_normalized_speed_before_turn_misses=[];
per_session_normalized_speed_after_turn_misses=[];

per_session_speed_timecourse_aligned_to_odor_hits=zeros(length(display_t),length(handles_conc.arena_file));
per_session_speed_timecourse_aligned_to_odor_miss=zeros(length(display_t),length(handles_conc.arena_file));

per_session_normalized_speed_timecourse_aligned_to_odor_hits=zeros(length(display_t),length(handles_conc.arena_file));
per_session_normalized_speed_timecourse_aligned_to_odor_miss=zeros(length(display_t),length(handles_conc.arena_file));

ii_files_included=0;

glm_speed=[];
glm_speed_ii=0;

id_speed_ii=0;
input_speed_data=[];
for ii=1:4
    input_speed_data(ii).data=[];
end
input_speed_data(1).description='Speed before turn, hits';
input_speed_data(2).description='Speed after turn, hits';
input_speed_data(3).description='Speed before turn, misses';
input_speed_data(4).description='Speed after turn, misses';


input_norm_speed_data(1).description='Speed before turn, hits';
input_norm_speed_data(2).description='Speed after turn, hits';
input_norm_speed_data(3).description='Speed before turn, misses';
input_norm_speed_data(4).description='Speed after turn, misses';

for fileNo=1:length(handles_conc.arena_file)
    % if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
      if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
     
        %Load angle file
        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        meanAngles=[];
        for trNo=1:length(handles_out.angles.trial)
            if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
                meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
            else
                meanAngles(trNo)=0; %I need to fix this
            end
        end

        handles_out_angle=handles_out;


        %Load conc data
        arena_file=handles_conc.arena_file{fileNo};
        %load the ouptut file
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        odor_plume_template=handles_out.odor_plume_template;
        op_predicted=handles_out.op_predicted;


         %Load conc data
        arena_file=handles_XY.arena_file{fileNo};
        %load the ouptut file
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        x_predicted=handles_out.x_predicted;
        y_predicted=handles_out.y_predicted;
        x_predicted_sh=handles_out.x_predicted_sh(:,1);
        y_predicted_sh=handles_out.y_predicted_sh(:,1); %Note: all sh are identical
        XYtest=handles_out.XYtest;

        this_file_speed_before_turn_hits=[];
        this_file_speed_after_turn_hits=[];
        this_file_speed_before_turn_misses=[];
        this_file_speed_after_turn_misses=[];

        this_file_speed_timecourse_aligned_to_turn_hits=zeros(length(display_t),trials.odor_trNo);
        this_file_speed_timecourse_aligned_to_odor_hits=zeros(length(display_t),trials.odor_trNo);
        no_trials_odor_hits=0;

        this_file_speed_timecourse_aligned_to_turn_miss=zeros(length(display_t),trials.odor_trNo);
        this_file_speed_timecourse_aligned_to_odor_miss=zeros(length(display_t),trials.odor_trNo);
        no_trials_odor_miss=0;

        this_file_dt_spout_hits=[];
        this_file_dt_spout_miss=[];

        for trNo=1:trials.odor_trNo

            these_x=trials.trial(trNo).XYtest(:,1);
            these_y=trials.trial(trNo).XYtest(:,2);
            ii_t_odor_on=1-handles_choices.trial_start_offset;

            these_distances=zeros(1,length(these_x));
            for ii_x=2:length(these_x)
                these_distances(ii_x)=these_distances(ii_x-1)+sqrt( (these_x(ii_x)-these_x(ii_x-1))^2 + (these_y(ii_x)-these_y(ii_x-1))^2);
            end
            

            %Find the last turn
            this_ii_turn=find(handles_out_angle.angles.trial(trNo).delta_x>100,1,'last');
            if ~isempty(this_ii_turn)

                ii_turns=handles_out_angle.angles.trial(trNo).ii_turns(this_ii_turn);
                ii_odor=handles_out_angle.angles.trial(trNo).ii_first_odor_encounter;
                these_distances=these_distances-these_distances(ii_turns);

                ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if ii_predictedend>size(XYtest,1)
                    ii_predictedend=size(XYtest,1);
                end
                ii_end=trials.odor_ii_end(trNo);
                if ii_end>size(XYtest,1)
                    ii_end=size(XYtest,1);
                end

                ii_start=trials.odor_ii_start(trNo);

                %Find positions for start points
                this_start_x=XYtest(trials.odor_ii_start(trNo),1);
                this_start_y=XYtest(trials.odor_ii_start(trNo),2);

                ii_start_x=ceil(this_start_x/(x_length/10));
                if ii_start_x==0
                    ii_start_x=1;
                end
                if ii_start_x>10
                    ii_start_x=10;
                end
 
                ii_start_y=ceil(this_start_y/(y_length/10));
                if ii_start_y==0
                    ii_start_y=1;
                end
                if ii_start_y>10
                    ii_start_y=10;
                end

                %Find positions for end points
                this_end_x=XYtest(ii_end,1);
                this_end_y=XYtest(ii_end,2);

                ii_end_x=ceil(this_end_x/(x_length/10));
                if ii_end_x==0
                    ii_end_x=1;
                end
                if ii_end_x>10
                    ii_end_x=10;
                end
 
                ii_end_y=ceil(this_end_y/(y_length/10));
                if ii_end_y==0
                    ii_end_y=1;
                end
                if ii_end_y>10
                    ii_end_y=10;
                end

                %Find the position for this last turn
                this_odor_ii=ii_predictedstart+ii_odor-1;

                %Find the position for this odor encounter
                this_ap_ii=ii_predictedstart+ii_turns-1;

                this_ap_x=XYtest(this_ap_ii,1);
                this_ap_y=XYtest(this_ap_ii,2);

                ii_ap_x=ceil(this_ap_x/(x_length/10));
                if ii_ap_x==0
                    ii_ap_x=1;
                end
                if ii_ap_x>10
                    ii_ap_x=10;
                end
 
                ii_ap_y=ceil(this_ap_y/(y_length/10));
                if ii_ap_y==0
                    ii_ap_y=1;
                end
                if ii_ap_y>10
                    ii_ap_y=10;
                end

                %Find the position when the mouse encounters the odor
                found_it=0;
                for ii_t=1:length(these_x)
                    %Odor truns on at ii_t_odor_on
                     if ii_t>=ii_t_odor_on
                        x_on=(ii_t-ii_t_odor_on)*dt*air_flow_speed;
                        if (these_x(ii_t)<=x_on)&(found_it==0)
                            found_it=1;
                            this_odor_on_x=these_x(ii_t);
                            this_odor_on_y=these_y(ii_t);
                            this_odor_on_ii_t=ii_t;
                            this_dt_odor_to_turn=dt*(this_odor_on_ii_t-ii_turns);
                            this_distance_from_turn=these_distances(ii_t);
                        end
                     end
                end

               

                %How long does the animal linger around the water spout after the angle
                %turn?
                this_ii_dt_spout=0;
                for ii=this_ap_ii:ii_predictedend
                    if trials.lane_per_trial(trNo)==1
                        %lane 1
                        spout_d=sqrt((XYtest(ii,1))^2+(XYtest(ii,2)-y_lane1)^2);
                    else
                        spout_d=sqrt((XYtest(ii,1))^2+(XYtest(ii,2)-y_lane4)^2);
                    end
                    if spout_d<=hit_radius
                        this_ii_dt_spout=this_ii_dt_spout+1;
                    end
                end
                behavior_analysis.file(fileNo).trial(trNo).dt_spout=this_ii_dt_spout*handles_XY.dt;
                behavior_analysis.file(fileNo).trial(trNo).odor_trial_type=trials.odor_trial_type(trNo);

                %Calculate speed before and after turn
                these_before_speeds=[];
                for ii=ii_start:this_ap_ii-1
                    delta_d=sqrt((XYtest(ii,1)-XYtest(ii+1,1))^2+(XYtest(ii,2)-XYtest(ii+1,2))^2);
                    these_before_speeds=[these_before_speeds delta_d/handles_XY.dt];
                end
                
                these_after_speeds=[];
                for ii=this_ap_ii:ii_end-1
                    delta_d=sqrt((XYtest(ii,1)-XYtest(ii+1,1))^2+(XYtest(ii,2)-XYtest(ii+1,2))^2);
                    these_after_speeds=[these_after_speeds delta_d/handles_XY.dt];
                end

                %Calculate speed aligned to odor encounter
                these_speeds_aligned_to_odor=[];
                for ii=this_odor_ii-floor(length(display_t)/2):this_odor_ii+floor(length(display_t)/2)
                    if ii+1<=size(XYtest,1)
                        delta_d=sqrt((XYtest(ii,1)-XYtest(ii+1,1))^2+(XYtest(ii,2)-XYtest(ii+1,2))^2);
                        last_speed=delta_d/handles_XY.dt;
                    end
                    these_speeds_aligned_to_odor=[these_speeds_aligned_to_odor last_speed];
                end
                
                %Calculate speed aligned to turn
                these_speeds_aligned_to_turn=[];
                for ii=this_ap_ii-floor(length(display_t)/2):this_ap_ii+floor(length(display_t)/2)
                    if ii+1<=size(XYtest,1)
                        delta_d=sqrt((XYtest(ii,1)-XYtest(ii+1,1))^2+(XYtest(ii,2)-XYtest(ii+1,2))^2);
                        last_speed=delta_d/handles_XY.dt;
                    end
                    these_speeds_aligned_to_turn=[these_speeds_aligned_to_turn last_speed];
                end

                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
 
                        %Time course for speed 
                        no_trials_odor_hits=no_trials_odor_hits+1;
                        this_file_speed_timecourse_aligned_to_odor_hits(:,no_trials_odor_hits)=these_speeds_aligned_to_odor;
                        this_file_speed_timecourse_aligned_to_turn_hits(:,no_trials_odor_hits)=these_speeds_aligned_to_turn;
        
                        %Speed before and after turn
                        dt_turn_to_odor_lane1=[dt_turn_to_odor_lane1 this_dt_odor_to_turn];
                        distance_turn_to_odor_lane1=[distance_turn_to_odor_lane1 this_distance_from_turn];
                        lane1_hit_turn_angle_positions(ii_ap_x,ii_ap_y)=lane1_hit_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane1_hit_start_positions(ii_start_x,ii_start_y)=lane1_hit_start_positions(ii_start_x,ii_start_y)+1;
                        lane1_hit_end_positions(ii_end_x,ii_end_y)=lane1_hit_end_positions(ii_end_x,ii_end_y)+1;
                        dt_spout_hits=[dt_spout_hits this_ii_dt_spout*handles_XY.dt];
                        this_file_dt_spout_hits=[this_file_dt_spout_hits this_ii_dt_spout*handles_XY.dt];
                        
                        speed_before_turn_hits=[speed_before_turn_hits mean(these_before_speeds)];
                        speed_after_turn_hits=[speed_after_turn_hits mean(these_after_speeds)];

                        this_file_speed_before_turn_hits=[this_file_speed_before_turn_hits mean(these_before_speeds)];
                        this_file_speed_after_turn_hits=[this_file_speed_after_turn_hits mean(these_after_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_before_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=0;
                        glm_speed.hit_miss(glm_speed_ii+1)=0;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(1).data=[input_speed_data(1).data mean(these_before_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_after_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=1;
                        glm_speed.hit_miss(glm_speed_ii+1)=0;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(2).data=[input_speed_data(2).data mean(these_after_speeds)];
                      
                    case 2
                        %Lane 1 miss orange

                        %Time course for speed 
                        no_trials_odor_miss=no_trials_odor_miss+1;
                        this_file_speed_timecourse_aligned_to_odor_miss(:,no_trials_odor_miss)=these_speeds_aligned_to_odor;
                        this_file_speed_timecourse_aligned_to_turn_miss(:,no_trials_odor_miss)=these_speeds_aligned_to_turn;
        
                        %Speed before and after turn
                        lane1_miss_turn_angle_positions(ii_ap_x,ii_ap_y)=lane1_miss_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane1_miss_start_positions(ii_start_x,ii_start_y)=lane1_miss_start_positions(ii_start_x,ii_start_y)+1;
                        lane1_miss_end_positions(ii_end_x,ii_end_y)=lane1_miss_end_positions(ii_end_x,ii_end_y)+1;
                        dt_spout_miss=[dt_spout_miss this_ii_dt_spout*handles_XY.dt];
                        this_file_dt_spout_miss=[this_file_dt_spout_miss this_ii_dt_spout*handles_XY.dt];

                        speed_before_turn_misses=[speed_before_turn_misses mean(these_before_speeds)];
                        speed_after_turn_misses=[speed_after_turn_misses mean(these_after_speeds)];

                        this_file_speed_before_turn_misses=[this_file_speed_before_turn_misses mean(these_before_speeds)];
                        this_file_speed_after_turn_misses=[this_file_speed_after_turn_misses mean(these_after_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_before_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=0;
                        glm_speed.hit_miss(glm_speed_ii+1)=1;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(3).data=[input_speed_data(3).data mean(these_before_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_after_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=1;
                        glm_speed.hit_miss(glm_speed_ii+1)=1;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(4).data=[input_speed_data(4).data mean(these_after_speeds)];
                    case 3
                        %Lane 4 hit blue


                        %Time course for speed 
                        no_trials_odor_hits=no_trials_odor_hits+1;
                        this_file_speed_timecourse_aligned_to_odor_hits(:,no_trials_odor_hits)=these_speeds_aligned_to_odor;
                        this_file_speed_timecourse_aligned_to_turn_hits(:,no_trials_odor_hits)=these_speeds_aligned_to_turn;
        
        
                        %Speed before and after turn
                        dt_turn_to_odor_lane4=[dt_turn_to_odor_lane4 this_dt_odor_to_turn];
                        distance_turn_to_odor_lane4=[distance_turn_to_odor_lane4 this_distance_from_turn];
                        lane4_hit_turn_angle_positions(ii_ap_x,ii_ap_y)=lane4_hit_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane4_hit_start_positions(ii_start_x,ii_start_y)=lane4_hit_start_positions(ii_start_x,ii_start_y)+1;
                        lane4_hit_end_positions(ii_end_x,ii_end_y)=lane4_hit_end_positions(ii_end_x,ii_end_y)+1;
                        dt_spout_hits=[dt_spout_hits this_ii_dt_spout*handles_XY.dt];
                        this_file_dt_spout_hits=[this_file_dt_spout_hits this_ii_dt_spout*handles_XY.dt];

                        speed_before_turn_hits=[speed_before_turn_hits mean(these_before_speeds)];
                        speed_after_turn_hits=[speed_after_turn_hits mean(these_after_speeds)];

                        this_file_speed_before_turn_hits=[this_file_speed_before_turn_hits mean(these_before_speeds)];
                        this_file_speed_after_turn_hits=[this_file_speed_after_turn_hits mean(these_after_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_before_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=0;
                        glm_speed.hit_miss(glm_speed_ii+1)=0;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(1).data=[input_speed_data(1).data mean(these_before_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_after_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=1;
                        glm_speed.hit_miss(glm_speed_ii+1)=0;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(2).data=[input_speed_data(2).data mean(these_after_speeds)];
                    case 4
                        %Lane 4 miss sky blue

                        %Time course for speed 
                        no_trials_odor_miss=no_trials_odor_miss+1;
                        this_file_speed_timecourse_aligned_to_odor_miss(:,no_trials_odor_miss)=these_speeds_aligned_to_odor;
                        this_file_speed_timecourse_aligned_to_turn_miss(:,no_trials_odor_miss)=these_speeds_aligned_to_turn;
        
                        %Speed before and after turn
                        lane4_miss_turn_angle_positions(ii_ap_x,ii_ap_y)=lane4_miss_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane4_miss_start_positions(ii_start_x,ii_start_y)=lane4_miss_start_positions(ii_start_x,ii_start_y)+1;
                        lane4_miss_end_positions(ii_end_x,ii_end_y)=lane4_miss_end_positions(ii_end_x,ii_end_y)+1;
                        dt_spout_miss=[dt_spout_miss this_ii_dt_spout*handles_XY.dt];
                        this_file_dt_spout_miss=[this_file_dt_spout_miss this_ii_dt_spout*handles_XY.dt];

                        speed_before_turn_misses=[speed_before_turn_misses mean(these_before_speeds)];
                        speed_after_turn_misses=[speed_after_turn_misses mean(these_after_speeds)];

                        this_file_speed_before_turn_misses=[this_file_speed_before_turn_misses mean(these_before_speeds)];
                        this_file_speed_after_turn_misses=[this_file_speed_after_turn_misses mean(these_after_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_before_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=0;
                        glm_speed.hit_miss(glm_speed_ii+1)=1;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(3).data=[input_speed_data(3).data mean(these_before_speeds)];

                        glm_speed.data(glm_speed_ii+1)=mean(these_after_speeds);
                        glm_speed.before_after(glm_speed_ii+1)=1;
                        glm_speed.hit_miss(glm_speed_ii+1)=1;
                        glm_speed_ii=glm_speed_ii+1;
                        input_speed_data(4).data=[input_speed_data(4).data mean(these_after_speeds)];
                end
            end

        end
        this_file_speed_before_turn_misses=this_file_speed_before_turn_misses(~isnan(this_file_speed_before_turn_misses));
        this_file_speed_after_turn_misses=this_file_speed_after_turn_misses(~isnan(this_file_speed_after_turn_misses));
        this_file_speed_before_turn_hits=this_file_speed_before_turn_hits(~isnan(this_file_speed_before_turn_hits));
        this_file_speed_after_turn_hits=this_file_speed_after_turn_hits(~isnan(this_file_speed_after_turn_hits));

        % average_speed=[average_speed mean([this_file_speed_before_turn_misses this_file_speed_after_turn_misses...
        %     this_file_speed_before_turn_hits])];
        average_speed=[average_speed mean([this_file_speed_before_turn_misses...
            this_file_speed_before_turn_hits])];
        normalized_speed_before_turn_hits=[normalized_speed_before_turn_hits this_file_speed_before_turn_hits/average_speed(end)];
        normalized_speed_after_turn_hits=[normalized_speed_after_turn_hits this_file_speed_after_turn_hits/average_speed(end)];
        normalized_speed_before_turn_misses=[normalized_speed_before_turn_misses this_file_speed_before_turn_misses/average_speed(end)];
        normalized_speed_after_turn_misses=[normalized_speed_after_turn_misses this_file_speed_after_turn_misses/average_speed(end)];

        per_session_normalized_speed_before_turn_hits=[per_session_normalized_speed_before_turn_hits mean(this_file_speed_before_turn_hits/average_speed(end))];
        per_session_normalized_speed_after_turn_hits=[per_session_normalized_speed_after_turn_hits mean(this_file_speed_after_turn_hits/average_speed(end))];
        per_session_normalized_speed_before_turn_misses=[per_session_normalized_speed_before_turn_misses mean(this_file_speed_before_turn_misses/average_speed(end))];
        per_session_normalized_speed_after_turn_misses=[per_session_normalized_speed_after_turn_misses mean(this_file_speed_after_turn_misses/average_speed(end))];

        per_session_dt_spout_hits=[per_session_dt_spout_hits mean(this_file_dt_spout_hits)];
        per_session_dt_spout_miss=[per_session_dt_spout_miss mean(this_file_dt_spout_miss)];

      
        
        this_file_speed_timecourse_aligned_to_odor_hits=this_file_speed_timecourse_aligned_to_odor_hits(:,1:no_trials_odor_hits);
        this_file_speed_timecourse_aligned_to_odor_miss=this_file_speed_timecourse_aligned_to_odor_miss(:,1:no_trials_odor_miss);
        this_file_speed_timecourse_aligned_to_turn_hits=this_file_speed_timecourse_aligned_to_turn_hits(:,1:no_trials_odor_hits);
        this_file_speed_timecourse_aligned_to_turn_miss=this_file_speed_timecourse_aligned_to_turn_miss(:,1:no_trials_odor_miss);
        
        ii_files_included=ii_files_included+1;
        per_session_speed_timecourse_aligned_to_odor_hits(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_odor_hits,2);
        per_session_speed_timecourse_aligned_to_odor_miss(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_odor_miss,2);
        per_session_speed_timecourse_aligned_to_turn_hits(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_turn_hits,2);
        per_session_speed_timecourse_aligned_to_turn_miss(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_turn_miss,2);

        per_session_normalized_speed_timecourse_aligned_to_odor_hits(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_odor_hits,2)/average_speed(end);
        per_session_normalized_speed_timecourse_aligned_to_odor_miss(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_odor_miss,2)/average_speed(end);
        per_session_normalized_speed_timecourse_aligned_to_turn_hits(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_turn_hits,2)/average_speed(end);
        per_session_normalized_speed_timecourse_aligned_to_turn_miss(:,ii_files_included)=mean(this_file_speed_timecourse_aligned_to_turn_miss,2)/average_speed(end);


      end
end


%Perform the glm mouse speed
fprintf(1, ['\nglm for mouse speed\n'])
fprintf(fileID, ['\nglm for mouse speed\n']);


tbl = table(glm_speed.data',glm_speed.before_after',glm_speed.hit_miss',...
    'VariableNames',{'speed','before_after','hit_miss'});
mdl = fitglm(tbl,'speed~before_after+hit_miss+before_after*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for mouse speed\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for mouse speed\n']);


[output_data] = drgMutiRanksumorTtest(input_speed_data, fileID,0);

%Normalize
lane1_hit_turn_angle_positions=lane1_hit_turn_angle_positions/sum(lane1_hit_turn_angle_positions(:));
lane1_miss_turn_angle_positions=lane1_miss_turn_angle_positions/sum(lane1_miss_turn_angle_positions(:));
lane4_hit_turn_angle_positions=lane4_hit_turn_angle_positions/sum(lane4_hit_turn_angle_positions(:));
lane4_miss_turn_angle_positions=lane4_miss_turn_angle_positions/sum(lane4_miss_turn_angle_positions(:));

lane1_hit_start_positions=lane1_hit_start_positions/sum(lane1_hit_start_positions(:));
lane1_miss_start_positions=lane1_miss_start_positions/sum(lane1_miss_start_positions(:));
lane4_hit_start_positions=lane4_hit_start_positions/sum(lane4_hit_start_positions(:));
lane4_miss_start_positions=lane4_miss_start_positions/sum(lane4_miss_start_positions(:));

lane1_hit_end_positions=lane1_hit_end_positions/sum(lane1_hit_end_positions(:));
lane1_miss_end_positions=lane1_miss_end_positions/sum(lane1_miss_end_positions(:));
lane4_hit_end_positions=lane4_hit_end_positions/sum(lane4_hit_end_positions(:));
lane4_miss_end_positions=lane4_miss_end_positions/sum(lane4_miss_end_positions(:));

% colormap fire
colormap gray
this_cmap=colormap;
this_cmap(1,:)=[0.2 0.2 0.2];
mult_cmax=0.35;


%Plot turn points
minC=0;
maxC=max([max(lane1_hit_turn_angle_positions(:)) max(lane1_miss_turn_angle_positions(:)) max(lane4_hit_turn_angle_positions(:)) max(lane4_miss_turn_angle_positions(:))]);


%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_turn_angle_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_hit_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit last turn positions')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_turn_angle_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_miss_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss last turn positions')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_turn_angle_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_hit_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit last turn positions')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_turn_angle_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_miss_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss last turn positions')

%Plot start points
minC=0;
maxC=max([max(lane1_hit_start_positions(:)) max(lane1_miss_start_positions(:)) max(lane4_hit_start_positions(:)) max(lane4_miss_start_positions(:))]);

%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_start_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_hit_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit start positions')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_start_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_miss_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss start positions')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_start_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_hit_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit start positions')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_start_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_miss_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss start positions')

%Plot end points
minC=0;
maxC=max([max(lane1_hit_end_positions(:)) max(lane1_miss_end_positions(:)) max(lane4_hit_end_positions(:)) max(lane4_miss_end_positions(:))]);

%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_end_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_hit_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit end positions')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_end_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_miss_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss end positions')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_end_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_hit_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit end positions')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_end_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_miss_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss end positions')

%Now plot pseudocolor map of trajectoreis before and after turn


%First make point maps
ii_run=1;
these_groups=[1 5];
ii_for_corr=0;

y_length=480;
x_length=500;

y_values=(y_length/20):y_length/10:y_length-(y_length/20);
x_values=(x_length/20):x_length/10:x_length-(x_length/20);

lane1_hit_before_positions=zeros(10,10);
lane1_miss_before_positions=zeros(10,10);

lane4_hit_before_positions=zeros(10,10);
lane4_miss_before_positions=zeros(10,10);

lane1_hit_after_positions=zeros(10,10);
lane1_miss_after_positions=zeros(10,10);

lane4_hit_after_positions=zeros(10,10);
lane4_miss_after_positions=zeros(10,10);


for fileNo=1:length(handles_conc.arena_file)
    % if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Load angle file
        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        meanAngles=[];
        for trNo=1:length(handles_out.angles.trial)
            if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
                meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
            else
                meanAngles(trNo)=0; %I need to fix this
            end
        end

        handles_out_angle=handles_out;


        %Load conc data
        arena_file=handles_conc.arena_file{fileNo};
        %load the ouptut file
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        odor_plume_template=handles_out.odor_plume_template;
        op_predicted=handles_out.op_predicted;


        %Load conc data
        arena_file=handles_XY.arena_file{fileNo};
        %load the ouptut file
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        x_predicted=handles_out.x_predicted;
        y_predicted=handles_out.y_predicted;
        x_predicted_sh=handles_out.x_predicted_sh(:,1);
        y_predicted_sh=handles_out.y_predicted_sh(:,1); %Note: all sh are identical
        XYtest=handles_out.XYtest;

        for trNo=1:trials.odor_trNo



            %Find the last turn
            this_ii_turn=find(handles_out_angle.angles.trial(trNo).delta_x>100,1,'last');
            if ~isempty(this_ii_turn)

                ii_turns=handles_out_angle.angles.trial(trNo).ii_turns(this_ii_turn);

                ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if ii_predictedend>size(XYtest,1)
                    ii_predictedend=size(XYtest,1);
                end
                ii_end=trials.odor_ii_end(trNo);
                if ii_end>size(XYtest,1)
                    ii_end=size(XYtest,1);
                end



                %Find the position for this last turn
                this_ap_ii=ii_predictedstart+ii_turns-1;

                %Before
                for ii=trials.odor_ii_start(trNo):this_ap_ii

                    %Find position for this point
                    this_x=XYtest(ii,1);
                    this_y=XYtest(ii,2);

                    ii_x=ceil(this_x/(x_length/10));
                    if ii_x==0
                        ii_x=1;
                    end
                    if ii_x>10
                        ii_x=10;
                    end

                    ii_y=ceil(this_y/(y_length/10));
                    if ii_y==0
                        ii_y=1;
                    end
                    if ii_y>10
                        ii_y=10;
                    end
                    %Okabe_Ito colors
                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            lane1_hit_before_positions(ii_x,ii_y)=lane1_hit_before_positions(ii_x,ii_y)+1;

                        case 2
                            %Lane 1 miss orange
                            lane1_miss_before_positions(ii_x,ii_y)=lane1_miss_before_positions(ii_x,ii_y)+1;

                        case 3
                            %Lane 4 hit blue
                            lane4_hit_before_positions(ii_x,ii_y)=lane4_hit_before_positions(ii_x,ii_y)+1;

                        case 4
                            %Lane 4 miss sky blue
                            lane4_miss_before_positions(ii_x,ii_y)=lane4_miss_before_positions(ii_x,ii_y)+1;
                    end
                end


                %After
                for ii=this_ap_ii+1:ii_end

                    %Find position for this point
                    this_x=XYtest(ii,1);
                    this_y=XYtest(ii,2);

                    ii_x=ceil(this_x/(x_length/10));
                    if ii_x==0
                        ii_x=1;
                    end
                    if ii_x>10
                        ii_x=10;
                    end

                    ii_y=ceil(this_y/(y_length/10));
                    if ii_y==0
                        ii_y=1;
                    end
                    if ii_y>10
                        ii_y=10;
                    end
                    %Okabe_Ito colors
                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            lane1_hit_after_positions(ii_x,ii_y)=lane1_hit_after_positions(ii_x,ii_y)+1;

                        case 2
                            %Lane 1 miss orange
                            lane1_miss_after_positions(ii_x,ii_y)=lane1_miss_after_positions(ii_x,ii_y)+1;

                        case 3
                            %Lane 4 hit blue
                            lane4_hit_after_positions(ii_x,ii_y)=lane4_hit_after_positions(ii_x,ii_y)+1;

                        case 4
                            %Lane 4 miss sky blue
                            lane4_miss_after_positions(ii_x,ii_y)=lane4_miss_after_positions(ii_x,ii_y)+1;
                    end
                end

            end
            pfft=1;
        end

    end
end

%Normalize
lane1_hit_before_positions=lane1_hit_before_positions/sum(lane1_hit_before_positions(:));
lane1_miss_before_positions=lane1_miss_before_positions/sum(lane1_miss_before_positions(:));
lane4_hit_before_positions=lane4_hit_before_positions/sum(lane4_hit_before_positions(:));
lane4_miss_before_positions=lane4_miss_before_positions/sum(lane4_miss_before_positions(:));

lane1_hit_after_positions=lane1_hit_after_positions/sum(lane1_hit_after_positions(:));
lane1_miss_after_positions=lane1_miss_after_positions/sum(lane1_miss_after_positions(:));
lane4_hit_after_positions=lane4_hit_after_positions/sum(lane4_hit_after_positions(:));
lane4_miss_after_positions=lane4_miss_after_positions/sum(lane4_miss_after_positions(:));


%Plot before points
%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_before_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_hit_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit trajectories before turn')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_before_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_miss_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss trajectories before turn')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_before_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_hit_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit trajectories before turn')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_before_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_miss_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss trajectories before turn')

%Plot after points
minC=0;
maxC=max([max(lane4_hit_after_positions(:)) max(lane1_miss_after_positions(:)) max(lane4_miss_after_positions(:)) max(lane1_hit_after_positions(:))]);

%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_after_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_hit_after_positions(:));
caxis([minC mult_cmax*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit trajectories after turn')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_after_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane1_miss_after_positions(:));
caxis([minC mult_cmax*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss trejectories after turn')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_after_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_hit_after_positions(:));
caxis([minC mult_cmax*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit trajectories after turn')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_after_positions)
colormap(this_cmap)
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

% minC=0;
% maxC=max(lane4_miss_after_positions(:));
caxis([minC mult_cmax*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss trejectories after turn')


%Plot the dt_spout graph
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[0:0.25:6];
rand_offset=0.5;

%Plot the different dt_spouts
bar(bar_offset,mean(dt_spout_hits),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%Violin plot
[mean_out, CIout]=drgViolinPoint(dt_spout_hits...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;


bar(bar_offset,mean(dt_spout_miss),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])


%Violin plot
[mean_out, CIout]=drgViolinPoint(dt_spout_miss...
    ,edges,bar_offset,rand_offset,'k','k',4);

for ii=1:length(per_session_dt_spout_hits)
    plot([bar_offset-1 bar_offset],[per_session_dt_spout_hits(ii),per_session_dt_spout_miss(ii)],...
        'o-','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'markerSize',8,'LineWidth',2)
end


% x_pos=-0.5;
% text(x_pos,5,'Hits','Color',[230/255 159/255 0/255])
% text(x_pos,4.5,'Misses','Color',[86/255 180/255 233/255])
% 

xticks([0 1])
xticklabels({'Hits','Misses'})


title(['Time at the water spout'])
ylabel('Time (sec)')
ylim([-0.5 6])
xlim([-0.75 1.75])

[h,p]=ttest2(dt_spout_hits,dt_spout_miss);

 fprintf(1,['\np value for time near spout ' num2str(p) '\n'])



%Plot the speed graph
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[0:20:400];
rand_offset=0.5;



speed_before_turn_hits=speed_before_turn_hits(~isnan(speed_before_turn_hits))/10;
speed_after_turn_hits=speed_after_turn_hits(~isnan(speed_after_turn_hits))/10;
speed_before_turn_misses=speed_before_turn_misses(~isnan(speed_before_turn_misses))/10;
speed_after_turn_misses=speed_after_turn_misses(~isnan(speed_after_turn_misses))/10;

%Speed before turn hits
bar(bar_offset,mean(speed_before_turn_hits),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(speed_before_turn_hits...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;

%Speed after turn hits
bar(bar_offset,mean(speed_after_turn_hits),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(speed_after_turn_hits...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+2;

%Speed before turn miss
bar(bar_offset,mean(speed_before_turn_misses),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(speed_before_turn_misses...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;

%Speed after turn miss
bar(bar_offset,mean(speed_after_turn_misses),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(speed_after_turn_misses...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;


x_pos=-0.5;
text(x_pos,45,'Hits','Color',[230/255 159/255 0/255])
text(x_pos,40,'Misses','Color',[86/255 180/255 233/255])


xticks([0 1 3 4])
xticklabels({'Before','After','Before','After'})


title(['Mouse speed'])
ylabel('Speed (cm/sec)')
ylim([0 50])
xlim([-1 5])



%Plot the normalized speed graph
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[0:20:400];
rand_offset=0.5;

normalized_speed_before_turn_hits=mean(average_speed)*normalized_speed_before_turn_hits/10;
normalized_speed_after_turn_hits=mean(average_speed)*normalized_speed_after_turn_hits/10;
normalized_speed_before_turn_misses=mean(average_speed)*normalized_speed_before_turn_misses/10;
normalized_speed_after_turn_misses=mean(average_speed)*normalized_speed_after_turn_misses/10;

per_session_normalized_speed_before_turn_hits=mean(average_speed)*per_session_normalized_speed_before_turn_hits/10;
per_session_normalized_speed_after_turn_hits=mean(average_speed)*per_session_normalized_speed_after_turn_hits/10;
per_session_normalized_speed_before_turn_misses=mean(average_speed)*per_session_normalized_speed_before_turn_misses/10;
per_session_normalized_speed_after_turn_misses=mean(average_speed)*per_session_normalized_speed_after_turn_misses/10;

speed_before_turn_hits=speed_before_turn_hits(~isnan(speed_before_turn_hits))/10;
speed_after_turn_hits=speed_after_turn_hits(~isnan(speed_after_turn_hits))/10;
speed_before_turn_misses=speed_before_turn_misses(~isnan(speed_before_turn_misses))/10;
speed_after_turn_misses=speed_after_turn_misses(~isnan(speed_after_turn_misses))/10;

%Speed before turn hits
bar(bar_offset,mean(normalized_speed_before_turn_hits),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(normalized_speed_before_turn_hits...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;

%Speed after turn hits
bar(bar_offset,mean(normalized_speed_after_turn_hits),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(normalized_speed_after_turn_hits...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+2;

%Speed before turn miss
bar(bar_offset,mean(normalized_speed_before_turn_misses),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(normalized_speed_before_turn_misses...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;

%Speed after turn miss
bar(bar_offset,mean(normalized_speed_after_turn_misses),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(normalized_speed_after_turn_misses...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;

for ii=1:length(per_session_normalized_speed_before_turn_hits)
    plot([0 1],[per_session_normalized_speed_before_turn_hits(ii),per_session_normalized_speed_after_turn_hits(ii)],...
        'o-','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'markerSize',8,'LineWidth',2)

    plot([3 4],[per_session_normalized_speed_before_turn_misses(ii),per_session_normalized_speed_after_turn_misses(ii)],...
        'o-','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'markerSize',8,'LineWidth',2)
end

x_pos=-0.5;
text(x_pos,45,'Hits','Color',[230/255 159/255 0/255])
text(x_pos,40,'Misses','Color',[86/255 180/255 233/255])


xticks([0 1 3 4])
xticklabels({'Before','After','Before','After'})

 
title(['Mouse speed'])
ylabel('Speed (cm/sec)')
ylim([0 80])
xlim([-1 5])

glm_norm_speed=[];
glm_norm_speed_ii=0;


glm_norm_speed.data(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_before_turn_hits))=normalized_speed_before_turn_hits;
glm_norm_speed.before_after(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_before_turn_hits))=0;
glm_norm_speed.hit_miss(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_before_turn_hits))=0;
glm_norm_speed_ii=glm_norm_speed_ii+length(normalized_speed_before_turn_hits);

input_norm_speed_data(1).data=normalized_speed_before_turn_hits;

glm_norm_speed.data(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_after_turn_hits))=normalized_speed_after_turn_hits;
glm_norm_speed.before_after(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_after_turn_hits))=1;
glm_norm_speed.hit_miss(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_after_turn_hits))=0;
glm_norm_speed_ii=glm_norm_speed_ii+length(normalized_speed_after_turn_hits);

input_norm_speed_data(2).data=normalized_speed_after_turn_hits;


glm_norm_speed.data(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_before_turn_misses))=normalized_speed_before_turn_misses;
glm_norm_speed.before_after(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_before_turn_misses))=0;
glm_norm_speed.hit_miss(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_before_turn_misses))=1;
glm_norm_speed_ii=glm_norm_speed_ii+length(normalized_speed_before_turn_misses);

input_norm_speed_data(3).data=normalized_speed_before_turn_misses;

glm_norm_speed.data(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_after_turn_misses))=normalized_speed_after_turn_misses;
glm_norm_speed.before_after(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_after_turn_misses))=1;
glm_norm_speed.hit_miss(glm_norm_speed_ii+1:glm_norm_speed_ii+length(normalized_speed_after_turn_misses))=1;
glm_norm_speed_ii=glm_norm_speed_ii+length(normalized_speed_after_turn_misses);

input_norm_speed_data(4).data=normalized_speed_after_turn_misses;



%Perform the glm mouse normalized speed
fprintf(1, ['\nglm for mouse normalized speed\n'])
fprintf(fileID, ['\nglm for mouse normalized speed\n']);


tbl = table(glm_norm_speed.data',glm_norm_speed.before_after',glm_norm_speed.hit_miss',...
    'VariableNames',{'speed','before_after','hit_miss'});
mdl = fitglm(tbl,'speed~before_after+hit_miss+before_after*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for mouse speed\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for mouse speed\n']);


[output_data] = drgMutiRanksumorTtest(input_norm_speed_data, fileID,0);


%Plot the percent correct graph
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[0:2:100];
rand_offset=0.5;

%Plot the different dt_spouts
bar(bar_offset,mean(percent_correct),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%Violin plot
[mean_out, CIout]=drgViolinPoint(percent_correct...
    ,edges,bar_offset,rand_offset,'k','k',8);
bar_offset=bar_offset+1;


title(['Percent correct'])
ylabel('Percent correct')
ylim([0 100])

xticks([0])
xticklabels({''})

[h,p]=ttest(percent_correct-50);

fprintf(1,['\np value for percent corect ' num2str(p) '\n'])

%Plot the histogram for the time for odor encounter with respect to the
%final turn
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
edges=[-10:0.5:3];
histogram([dt_turn_to_odor_lane1 dt_turn_to_odor_lane4],edges)

xlabel('Time (sec)')
title('Time from last turn to odor encounter')
ylabel('Count')

%Plot the histogram for the time for odor encounter with respect to the
%final turn
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
edges=[-700:100:500];
histogram([distance_turn_to_odor_lane1 distance_turn_to_odor_lane4],edges)

xlabel('Distance (mm)')
title('Distance from last turn to odor encounter')
ylabel('Count')


%Trim the per session speed files
per_session_speed_timecourse_aligned_to_odor_hits=per_session_speed_timecourse_aligned_to_odor_hits(:,1:ii_files_included)/10;
per_session_speed_timecourse_aligned_to_odor_miss=per_session_speed_timecourse_aligned_to_odor_miss(:,1:ii_files_included)/10;

per_session_speed_timecourse_aligned_to_turn_hits=per_session_speed_timecourse_aligned_to_turn_hits(:,1:ii_files_included)/10;
per_session_speed_timecourse_aligned_to_turn_miss=per_session_speed_timecourse_aligned_to_turn_miss(:,1:ii_files_included)/10;





%Show the time course for speed aligned to odor encounter
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
%Miss
CIpvsm = bootci(1000, @mean, per_session_speed_timecourse_aligned_to_odor_miss');
meanpvsm=mean(per_session_speed_timecourse_aligned_to_odor_miss',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_speed_timecourse_aligned_to_odor_miss'), CIpvsm','cmap',[86/255 180/255 233/255]);

%Hits
CIpvsm = bootci(1000, @mean, per_session_speed_timecourse_aligned_to_odor_hits');
meanpvsm=mean(per_session_speed_timecourse_aligned_to_odor_hits',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_speed_timecourse_aligned_to_odor_hits'), CIpvsm','cmap',[230/255 159/255 0/255]);

text(-1,80,'Hits','Color',[230/255 159/255 0/255])
text(-1,70,'Misses','Color',[86/255 180/255 233/255])

xlabel('Time (sec)')
ylabel('Speed (cm/sec)')
title('Mouse speed aligned to odor onset')


%Show the time course for speed aligned to turn
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
%Miss
CIpvsm = bootci(1000, @mean, per_session_speed_timecourse_aligned_to_turn_miss');
meanpvsm=mean(per_session_speed_timecourse_aligned_to_turn_miss',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_speed_timecourse_aligned_to_turn_miss'), CIpvsm','cmap',[86/255 180/255 233/255]);

%Hits
CIpvsm = bootci(1000, @mean, per_session_speed_timecourse_aligned_to_turn_hits');
meanpvsm=mean(per_session_speed_timecourse_aligned_to_turn_hits',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_speed_timecourse_aligned_to_turn_hits'), CIpvsm','cmap',[230/255 159/255 0/255]);

text(-1,80,'Hits','Color',[230/255 159/255 0/255])
text(-1,70,'Misses','Color',[86/255 180/255 233/255])

xlabel('Time (sec)')
ylabel('Speed (cm/sec)')
title('Mouse speed aligned to turn')


%Trim the per session speed files
per_session_normalized_speed_timecourse_aligned_to_odor_hits=mean(average_speed)*per_session_normalized_speed_timecourse_aligned_to_odor_hits(:,1:ii_files_included)/10;
per_session_normalized_speed_timecourse_aligned_to_odor_miss=mean(average_speed)*per_session_normalized_speed_timecourse_aligned_to_odor_miss(:,1:ii_files_included)/10;

per_session_normalized_speed_timecourse_aligned_to_turn_hits=mean(average_speed)*per_session_normalized_speed_timecourse_aligned_to_turn_hits(:,1:ii_files_included)/10;
per_session_normalized_speed_timecourse_aligned_to_turn_miss=mean(average_speed)*per_session_normalized_speed_timecourse_aligned_to_turn_miss(:,1:ii_files_included)/10;


%Show the time course for speed aligned to odor encounter
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
%Miss
CIpvsm = bootci(1000, @mean, per_session_normalized_speed_timecourse_aligned_to_odor_miss');
meanpvsm=mean(per_session_normalized_speed_timecourse_aligned_to_odor_miss',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_normalized_speed_timecourse_aligned_to_odor_miss'), CIpvsm','cmap',[86/255 180/255 233/255]);

%Hits
CIpvsm = bootci(1000, @mean, per_session_normalized_speed_timecourse_aligned_to_odor_hits');
meanpvsm=mean(per_session_normalized_speed_timecourse_aligned_to_odor_hits',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_normalized_speed_timecourse_aligned_to_odor_hits'), CIpvsm','cmap',[230/255 159/255 0/255]);

text(-1,80,'Hits','Color',[230/255 159/255 0/255])
text(-1,70,'Misses','Color',[86/255 180/255 233/255])

xlabel('Time (sec)')
ylabel('Speed (cm/sec)')
title('Mouse speed aligned to odor onset (normalized)')


%Show the time course for speed aligned to turn
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
%Miss
CIpvsm = bootci(1000, @mean, per_session_normalized_speed_timecourse_aligned_to_turn_miss');
meanpvsm=mean(per_session_normalized_speed_timecourse_aligned_to_turn_miss',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_normalized_speed_timecourse_aligned_to_turn_miss'), CIpvsm','cmap',[86/255 180/255 233/255]);

%Hits
CIpvsm = bootci(1000, @mean, per_session_normalized_speed_timecourse_aligned_to_turn_hits');
meanpvsm=mean(per_session_normalized_speed_timecourse_aligned_to_turn_hits',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(display_t,mean(per_session_normalized_speed_timecourse_aligned_to_turn_hits'), CIpvsm','cmap',[230/255 159/255 0/255]);

text(-1,80,'Hits','Color',[230/255 159/255 0/255])
text(-1,70,'Misses','Color',[86/255 180/255 233/255])

xlabel('Time (sec)')
ylabel('Speed (cm/sec)')
title('Mouse speed aligned to turn (normalized)')

fclose(fileID);

pffft=1;