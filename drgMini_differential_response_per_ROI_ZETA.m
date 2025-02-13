%drgMini_differential_response_perROIv2
close all
clear all

is_sphgpu=0;

switch is_sphgpu
    case 0

        % %There was a bug and the shuffled runs were the same for all shuffled runs
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';
        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        %For the files below the shuffled runs should be different
        %Trained with all trials
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01062025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01062025.m';

         %Trained with hits only
         save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
         choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        %This is not used here, all the Moser infor is input through PredImp
        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';


        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        save_PathPredImp='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FilePredImp='outputPredictionImportance.mat';

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end


percentile_stringency=95; %mean_imps at this percentile stringency are included as cells of high prediction importance

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


% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])
eval(['handles_Moser=' choiceMoserFileName(1:end-2) ';'])

figureNo=0;

colormap fire
this_cmap=colormap;
this_cmap(1,:)=[0.3 0.3 0.3];

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

these_groups=[1 5];
ii_run=1;

n_shuffle_SI=100;


x=25:50:475;
y=24:48:456;

trial_type_labels{1}='Hit1';
trial_type_labels{2}='Miss1';
trial_type_labels{3}='Hit4';
trial_type_labels{4}='Miss4';

these_adiv_names{1}='all_div_hit1';
these_adiv_names{2}='all_div_miss1';
these_adiv_names{3}='all_div_hit4';
these_adiv_names{4}='all_div_miss4';

%Labels for comparisons
ii_comp=0;
for ii_type1=1:4
    for ii_type2=ii_type1+1:4
        ii_comp=ii_comp+1;
        trial_type_comp_labels{ii_comp}=[trial_type_labels{ii_type1} ' vs ' trial_type_labels{ii_type2}];
        trial_type_comp_ii_type1(ii_comp)=ii_type1;
        trial_type_comp_ii_type2(ii_comp)=ii_type2;
        these_types=[1:4];
        these_types=these_types((these_types~=ii_type1)&(these_types~=ii_type2))
        trial_type_comp_ii_type3(ii_comp)=these_types(1);
        trial_type_comp_ii_type4(ii_comp)=these_types(2);
        trial_type_no_comp_labels{ii_comp}=[trial_type_labels{these_types(1)} ' vs ' trial_type_labels{these_types(2)}];
    end
end

p_values=[]; %dimension 1 is ii_all_pvalues, dimension 2 is ii_comp 1 to 6
ii_all_pvalues=0;
all_ROI_files=[];
all_iiROIs=[];
pFDRs_per_file=[];

%Load information on prediction importance to do the Fusi analysis
load([save_PathPredImp save_FilePredImp])

handles_out2=[];

%Show a subset of the spatial maps
figureNo=figureNo+1;
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Get dFF
        arena_file=handles_XY.arena_file{fileNo};
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        no_neurons=handles_out.no_neurons;

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

        %Get XY and dFF per trial



        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        angles=handles_out.angles;

        % these_important_ROIs=unique([imps.file(fileNo).ROIs_conc imps.file(fileNo).ROIs_x imps.file(fileNo).ROIs_y]);

        %Now show the spatial activity maps
        no_place_cells=0;
        place_cells=[];
        no_lane_trial_cells=0;
        lane_trial_cells=[];
        these_p_values_per_ROI=zeros(no_neurons,6);

        for ii_ROI=1:no_neurons

            % ii_ROI_all=find((all_info_fileNo==fileNo)&(all_info_ii_ROI==ii_ROI));

            % this_ROI=these_important_ROIs(ii_ROI);
            % spatial_rhol1l4=imps.file(fileNo).spatial_rhol1l4;
            % delta_center_of_mass=imps.file(fileNo).delta_center_of_mass;
            % SSI=imps.file(fileNo).information_content;
            % SSIl1=imps.file(fileNo).information_contentl1;
            % SSIl4=imps.file(fileNo).information_contentl4;
            % sparsity=imps.file(fileNo).sparsity;

            %Initialize variables
            this_dFF_activity=zeros(10,10);
            this_dFF_activity_n=zeros(10,10);
            sum_dFF_activity=0;

            this_dFFl1_activity=zeros(10,10);
            this_dFFl1_activity_n=zeros(10,10);
            sum_dFFl1_activity=0;

            this_dFFl4_activity=zeros(10,10);
            this_dFFl4_activity_n=zeros(10,10);
            sum_dFFl4_activity=0;



            % delta_below_zero_ii=max(all_ii_turns(include_trial==1));
            % delta_above_zero_ii=max(all_ii_ends(include_trial==1)-all_ii_turns(include_trial==1));

            
            %These encompass times from the earliest start time to the latest end times in all trials 
            hit1_dFF=[];
            ii_hit1=0;
            miss1_dFF=[];
            ii_miss1=0;
            hit4_dFF=[];
            ii_hit4=0;
            miss4_dFF=[];
            ii_miss4=0;

            %Let's do the glm and ZETA stats from -3 to 3 sec in 1 sec bins
            trimmed_time_range=[-3 3];
            trimmed_dt=1;
            trimmed_time_bins=[trimmed_time_range(1)+(trimmed_dt/2):trimmed_dt:trimmed_time_range(2)-(trimmed_dt/2)];

            

            glm_div_ii=0;
            glm_div=[];


            %Find turn points and time vectors
            all_ii_turns=zeros(1,trials.odor_trNo);
            all_ii_ends=zeros(1,trials.odor_trNo);
            include_trial=zeros(1,trials.odor_trNo);
            for trNo=1:trials.odor_trNo
                this_ii_last_turn=find(angles.trial(trNo).delta_x>=100,1,'last');
                if ~isempty(this_ii_last_turn)
                    all_ii_turns(trNo)=angles.trial(trNo).ii_turns(this_ii_last_turn);
                    all_ii_ends(trNo)= size(trials.trial(trNo).XYtest,1);
                    include_trial(trNo)=1;
                else
                    [maxnum, maxii]=max(angles.trial(trNo).delta_x);
                    all_ii_turns(trNo)=angles.trial(trNo).ii_turns(maxii);
                    all_ii_ends(trNo)= size(trials.trial(trNo).XYtest,1);
                end
            end

            delta_below_zero_ii=max(all_ii_turns);
            delta_above_zero_ii=max(all_ii_ends-all_ii_turns);

            time_bins=handles_XY.dt*([1:delta_below_zero_ii+delta_above_zero_ii]-delta_below_zero_ii);

            %These encompass times from start to the end of trimmed
            %trials
            turn_time_bins=time_bins((time_bins>=trimmed_time_range(1))&(time_bins<=trimmed_time_range(2)));
            hit1_dFF_turn=zeros(sum(trials.hit1),sum((time_bins>=trimmed_time_range(1))&(time_bins<=trimmed_time_range(2))));
            miss1_dFF_turn=zeros(sum(trials.miss1),sum((time_bins>=trimmed_time_range(1))&(time_bins<=trimmed_time_range(2))));
            hit4_dFF_turn=zeros(sum(trials.hit4),sum((time_bins>=trimmed_time_range(1))&(time_bins<=trimmed_time_range(2))));
            miss4_dFF_turn=zeros(sum(trials.miss4),sum((time_bins>=trimmed_time_range(1))&(time_bins<=trimmed_time_range(2))));

            for trNo=1:trials.odor_trNo
                these_x=trials.trial(trNo).XYtest(:,1);
                these_y=trials.trial(trNo).XYtest(:,2);
                these_dFF=trials.trial(trNo).XdFFtest(:,ii_ROI);

                for ii_t=1:length(these_x)
                    this_x_ii=ceil(these_x(ii_t)/50);
                    if this_x_ii==11
                        this_x_ii=10;
                    end

                    this_y_ii=ceil(these_y(ii_t)/48);
                    if this_y_ii==11
                        this_y_ii=10;
                    end

                    this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                    this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
                    sum_dFF_activity=sum_dFF_activity+these_dFF(ii_t);

                    %Tally info
                    this_bin_dFF=(these_dFF(ii_t)>0)+1;
                    % cum_bindFF(this_bin_dFF)=cum_bindFF(this_bin_dFF)+1;
                    % cum_bindFF_BothLanes(this_bin_dFF)=cum_bindFF_BothLanes(this_bin_dFF)+1;
                    this_xy_ii=this_x_ii+10*(this_y_ii-1);
                    % cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
                    % cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)+1;
                    % cum_xy_BothLanes(this_xy_ii)=cum_xy_BothLanes(this_xy_ii)+1;

                    if trials.lane_per_trial(trNo)==1
                        this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
                        sum_dFFl1_activity=sum_dFFl1_activity+these_dFF(ii_t);
                        % cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)+1;
                        % cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)+1;
                        % cum_xy_Lane1(this_xy_ii)=cum_xy_Lane1(this_xy_ii)+1;
                        % cum_bindFF_Lane1(this_bin_dFF)=cum_bindFF_Lane1(this_bin_dFF)+1;
                        % cum_lane(1)=cum_lane(1)+1;
                    else
                        this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
                        sum_dFFl4_activity=sum_dFFl4_activity+these_dFF(ii_t);
                        % cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)+1;
                        % cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)+1;
                        % cum_xy_Lane4(this_xy_ii)=cum_xy_Lane4(this_xy_ii)+1;
                        % cum_bindFF_Lane4(this_bin_dFF)=cum_bindFF_Lane4(this_bin_dFF)+1;
                        % cum_lane(2)=cum_lane(2)+1;
                    end

                end
                % end
                %
                %
                %
                %
                %
                % % %Plot dFF timecourses
                % % subplot(2,3,1:3)
                % % hold on
                % % ii_plot=0;
                %
                %
                % for trNo=1:trials.odor_trNo
                % if include_trial(trNo)==1
                % these_dFF=trials.trial(trNo).XdFFtest(:,ii_ROI);
                these_time_bins=time_bins(delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)));
                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
                        ii_hit1=ii_hit1+1;
                        hit1_dFF(ii_hit1,1:length(time_bins))=min(these_dFF);
                        hit1_dFF(ii_hit1,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                        hit1_dFF_turn(ii_hit1,:)=min(these_dFF);
                        for ii_tur=1:length(turn_time_bins)
                            if sum(these_time_bins==turn_time_bins(ii_tur))
                                hit1_dFF_turn(ii_hit1,ii_tur)=these_dFF(these_time_bins==turn_time_bins(ii_tur));
                            end
                        end
                        for ii_tr=1:length(trimmed_time_bins)
                            these_tr_time_bins=(these_time_bins>=(trimmed_time_range(1)+(ii_tr-1)*trimmed_dt))&(these_time_bins<(trimmed_time_range(1)+ii_tr*trimmed_dt));
                            if sum(these_tr_time_bins)>0
                                glm_div.data(glm_div_ii+1)=mean(these_dFF(these_tr_time_bins));
                                glm_div.trial_type(glm_div_ii+1)=1;
                                glm_div.time(glm_div_ii+1)=trimmed_time_bins(ii_tr);
                                glm_div_ii=glm_div_ii+1;
                            end
                        end
                    case 2
                        %Lane 1 miss orange
                        ii_miss1=ii_miss1+1;
                        miss1_dFF(ii_miss1,1:length(time_bins))=min(these_dFF);
                        miss1_dFF(ii_miss1,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                        miss1_dFF_turn(ii_miss1,:)=min(these_dFF);
                        for ii_tur=1:length(turn_time_bins)
                            if sum(these_time_bins==turn_time_bins(ii_tur))
                                miss1_dFF_turn(ii_miss1,ii_tur)=these_dFF(these_time_bins==turn_time_bins(ii_tur));
                            end
                        end
                        for ii_tr=1:length(trimmed_time_bins)
                            these_tr_time_bins=(these_time_bins>=(trimmed_time_range(1)+(ii_tr-1)*trimmed_dt))&(these_time_bins<(trimmed_time_range(1)+ii_tr*trimmed_dt));
                            if sum(these_tr_time_bins)>0
                                glm_div.data(glm_div_ii+1)=mean(these_dFF(these_tr_time_bins));
                                glm_div.trial_type(glm_div_ii+1)=2;
                                glm_div.time(glm_div_ii+1)=trimmed_time_bins(ii_tr);
                                glm_div_ii=glm_div_ii+1;
                            end
                        end
                    case 3
                        %Lane 4 hit blue
                        ii_hit4=ii_hit4+1;
                        hit4_dFF(ii_hit4,1:length(time_bins))=min(these_dFF);
                        hit4_dFF(ii_hit4,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                        hit4_dFF_turn(ii_hit4,:)=min(these_dFF);
                        for ii_tur=1:length(turn_time_bins)
                            if sum(these_time_bins==turn_time_bins(ii_tur))
                                hit4_dFF_turn(ii_hit4,ii_tur)=these_dFF(these_time_bins==turn_time_bins(ii_tur));
                            end
                        end
                        for ii_tr=1:length(trimmed_time_bins)
                            these_tr_time_bins=(these_time_bins>=(trimmed_time_range(1)+(ii_tr-1)*trimmed_dt))&(these_time_bins<(trimmed_time_range(1)+ii_tr*trimmed_dt));
                            if sum(these_tr_time_bins)>0
                                glm_div.data(glm_div_ii+1)=mean(these_dFF(these_tr_time_bins));
                                glm_div.trial_type(glm_div_ii+1)=3;
                                glm_div.time(glm_div_ii+1)=trimmed_time_bins(ii_tr);
                                glm_div_ii=glm_div_ii+1;
                            end
                        end
                    case 4
                        %Lane 4 miss sky blue
                        ii_miss4=ii_miss4+1;
                        miss4_dFF(ii_miss4,1:length(time_bins))=min(these_dFF);
                        miss4_dFF(ii_miss4,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                        miss4_dFF_turn(ii_miss4,:)=min(these_dFF);
                        for ii_tur=1:length(turn_time_bins)
                            if sum(these_time_bins==turn_time_bins(ii_tur))
                                miss4_dFF_turn(ii_miss4,ii_tur)=these_dFF(these_time_bins==turn_time_bins(ii_tur));
                            end
                        end
                        for ii_tr=1:length(trimmed_time_bins)
                            these_tr_time_bins=(these_time_bins>=(trimmed_time_range(1)+(ii_tr-1)*trimmed_dt))&(these_time_bins<(trimmed_time_range(1)+ii_tr*trimmed_dt));
                            if sum(these_tr_time_bins)>0
                                glm_div.data(glm_div_ii+1)=mean(these_dFF(these_tr_time_bins));
                                glm_div.trial_type(glm_div_ii+1)=4;
                                glm_div.time(glm_div_ii+1)=trimmed_time_bins(ii_tr);
                                glm_div_ii=glm_div_ii+1;
                            end
                        end
                end
                pffft=1;
            end

            handles_out2.file(fileNo).ROI(ii_ROI).hit1_dFF_turn=mean(hit1_dFF_turn,1);
            handles_out2.file(fileNo).ROI(ii_ROI).miss1_dFF_turn=mean(miss1_dFF_turn,1);
            handles_out2.file(fileNo).ROI(ii_ROI).hit4_dFF_turn=mean(hit4_dFF_turn,1);
            handles_out2.file(fileNo).ROI(ii_ROI).miss4_dFF_turn=mean(miss4_dFF_turn,1);
            handles_out2.file(fileNo).ROI(ii_ROI).turn_time_bins=turn_time_bins;

            for ii_x=1:10
                for ii_y=1:10
                    if this_dFF_activity_n(ii_x,ii_y)~=0
                        this_dFF_activity(ii_x,ii_y)=this_dFF_activity(ii_x,ii_y)/this_dFF_activity_n(ii_x,ii_y);
                    end
                end
            end

            for ii_x=1:10
                for ii_y=1:10
                    if this_dFFl1_activity_n(ii_x,ii_y)~=0
                        this_dFFl1_activity(ii_x,ii_y)=this_dFFl1_activity(ii_x,ii_y)/this_dFFl1_activity_n(ii_x,ii_y);
                    end
                end
            end

            for ii_x=1:10
                for ii_y=1:10
                    if this_dFFl4_activity_n(ii_x,ii_y)~=0
                        this_dFFl4_activity(ii_x,ii_y)=this_dFFl4_activity(ii_x,ii_y)/this_dFFl4_activity_n(ii_x,ii_y);
                    end
                end
            end

            %Plot the space activity maps
            try
                close(figureNo)
            catch
            end


            hFig = figure(figureNo);
            set(hFig, 'units','normalized','position',[.1 .1 .75 .75])

            %Plot dFFs for lane 1 and 4
            y_gap=2;
            hit_miss_1=zeros(size(hit1_dFF,1)+size(miss1_dFF,1)+y_gap,delta_below_zero_ii+delta_above_zero_ii);
            hit_miss_1(size(miss1_dFF,1)+y_gap+1:size(hit1_dFF,1)+size(miss1_dFF,1)+y_gap,:)=hit1_dFF;
            hit_miss_1(1:size(miss1_dFF,1),:)=miss1_dFF;

            y_trials1=[1:size(hit_miss_1,1)];

            hit_miss_4=zeros(size(hit4_dFF,1)+size(miss4_dFF,1)+y_gap,delta_below_zero_ii+delta_above_zero_ii);
            hit_miss_4(size(miss4_dFF,1)+y_gap+1:size(hit4_dFF,1)+size(miss4_dFF,1)+y_gap,:)=hit4_dFF;
            hit_miss_4(1:size(miss4_dFF,1),:)=miss4_dFF;

            y_trials4=[1:size(hit_miss_4,1)];

            y_trials_end=max([y_trials1(end) y_trials4(end)]);

            all_hit_miss_4=hit_miss_4(:);
            all_hit_miss_1=hit_miss_1(:);
            all_hit_miss=[hit_miss_4(:); hit_miss_1(:)];
            c_percentile=prctile(all_hit_miss(all_hit_miss>0),99);

            %Plot dFF for lane 1
            subplot(2, 6, [1 2 3]);
            hold on

            drg_pcolor(repmat(time_bins,length(y_trials1),1),repmat(y_trials1,length(time_bins),1)',hit_miss_1)
            colormap(this_cmap)
            if max(max([max(hit_miss_4(:)) max(hit_miss_1(:))]))>0
                clim([0 c_percentile]);
            end

            shading flat
            plot([0 0],[y_trials1(1) y_trials1(end)+1],'-w','LineWidth',3)

            %y_trials1=1 to 2 is the first miss
            %size(miss1_dFF,1) to size(miss1_dFF,1) + 1 is the last miss,
            %size(miss1_dFF,1)+1 to size(miss1_dFF,1)+1+y_gap is the in between miss
            %and hit
            %size(miss1_dFF,1)+1+y_gap to size(miss1_dFF,1)+1+y_gap+size(hit1_dFF,1) is
            %hit

            rectangle('Position', [time_bins(1), size(miss1_dFF,1)+1, 1.03*(time_bins(end)-time_bins(1)), y_gap], ... % [x, y, width, height]
                'FaceColor', 'white', ...    % Fill color (white)
                'EdgeColor', 'none');        % No border for the rectangle

            ylim([0 y_trials_end+2])
            xlim([time_bins(1)-0.05*(time_bins(end)-time_bins(1)) time_bins(end)+0.05*(time_bins(end)-time_bins(1)) ])
            xlabel('Time (sec)')
            yticks([1+(size(miss1_dFF,1)/2) size(miss1_dFF,1)+1+y_gap+(size(hit1_dFF,1)/2)])
            yticklabels({'Misses','Hits'})
            title(['Lane 1 dFF'])

            %Plot dFF for lane 4
            subplot(2, 6, [4 5 6]);
            hold on

            drg_pcolor(repmat(time_bins,length(y_trials4),1),repmat(y_trials4,length(time_bins),1)',hit_miss_4)
            colormap(this_cmap)
            if max([max(hit_miss_4(:)) max(hit_miss_1(:))])>0
                clim([0 c_percentile])
            end
            shading flat
            plot([0 0],[y_trials4(1) y_trials4(end)+1],'-w','LineWidth',3)

            rectangle('Position', [time_bins(1), size(miss4_dFF,1)+1, 1.03*(time_bins(end)-time_bins(1)), y_gap], ... % [x, y, width, height]
                'FaceColor', 'white', ...    % Fill color (white)
                'EdgeColor', 'none');        % No border for the rectangle
            ylim([0 y_trials_end+2])
            xlim([time_bins(1)-0.05*(time_bins(end)-time_bins(1)) time_bins(end)+0.05*(time_bins(end)-time_bins(1)) ])
            xlabel('Time (sec)')
            yticks([1+(size(miss4_dFF,1)/2) size(miss4_dFF,1)+1+y_gap+(size(hit4_dFF,1)/2)])
            yticklabels({'Misses','Hits'})
            title(['Lane 4 dFF'])

            %Plot space activity maps as in Moser https://doi.org/10.1126/science.1114037
            max_this_dFF_activity=max(this_dFF_activity(:));
            max_this_dFFl4_activity=max(this_dFFl4_activity(:));
            max_this_dFFl1_activity=max(this_dFFl1_activity(:));

            colormap fire
            this_cmap=colormap;
            this_cmap(1,:)=[0.3 0.3 0.3];

            if max_this_dFFl1_activity>max_this_dFFl4_activity
                %Lane 1
                subplot(2, 6, [7 8]);
                max_activity=max_this_dFFl1_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_dFFl1_activity=this_dFFl1_activity;
                this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([70 100 200 300 400 430])
                yticklabels({'Lane 4','100','200','300','400','Lane 1'})

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')


                title_legend=['Lane 1'];
                % if sig_all_info_lane1(ii_ROI_all)==1
                %     title_legend=[title_legend ' S'];
                % end

                title(title_legend)

                %Lane 4 normalized to lane 1
                subplot(2, 6, [9 10]);
                max_activity=max_this_dFFl1_activity;
                this_masked_dFFl4_activity=this_dFFl4_activity;
                this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 4 normalized to Lane 1'])

                %Lane 4 normalized to lane 4
                subplot(2, 6, [11 12]);
                max_activity=max_this_dFFl4_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_dFFl4_activity=this_dFFl4_activity;
                this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')


                title_legend=['Lane 4'];
                % if sig_all_info_lane4(ii_ROI_all)==1
                %     title_legend=[title_legend ' S'];
                % end

                title(title_legend)
            else


                %Lane 4 normalized to lane 4
                subplot(2, 6, [7 8]);
                max_activity=max_this_dFFl4_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_dFFl4_activity=this_dFFl4_activity;
                this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                xticks(0:50:500)

                yticks([70 100 200 300 400 430])
                yticklabels({'Lane 4','100','200','300','400','Lane 1'})

                xlabel('x (mm)')
                ylabel('y (mm)')

                title_legend=['Lane 4'];
                % if sig_all_info_lane4(ii_ROI_all)==1
                %     title_legend=[title_legend ' S'];
                % end
                title(title_legend)

                %Lane 1 normalized to lane 4
                subplot(2, 6, [9 10]);
                max_activity=max_this_dFFl4_activity;
                this_masked_dFFl1_activity=this_dFFl1_activity;
                this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 1 normalized to Lane 4'])

                %Lane 1 normalized to lane 1
                subplot(2, 6, [11 12]);
                max_activity=max_this_dFFl1_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.00001;
                end
                this_masked_dFFl1_activity=this_dFFl1_activity;
                this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')

                title_legend=['Lane 1'];
                % if sig_all_info_lane1(ii_ROI_all)==1
                %     title_legend=[title_legend ' S'];
                % end
                title(title_legend)
            end


            sgt_legend=['dFF map file  No ' num2str(fileNo) ' ROI No ' num2str(ii_ROI) ];



            % 
            % 
            % if sig_all_info_mutual_info14(ii_ROI_all)==1
            %     sgt_legend=[sgt_legend ' S']; %MI is significant
            % end
            % 
            if sum(ii_ROI==imps.file(fileNo).ROIs_conc)>0
                sgt_legend=[sgt_legend ' odor '];
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_x)>0
                sgt_legend=[sgt_legend ' x '];
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_y)>0
                sgt_legend=[sgt_legend ' y '];
            end
            % 
            % switch idx(ii_ROI_all)
            %     case 1
            %         sgt_legend=[sgt_legend ' cl1 '];
            %     case 2
            %         sgt_legend=[sgt_legend ' cl2 '];
            %     case 3
            %         sgt_legend=[sgt_legend ' cl3 '];
            % end

            sgtitle(sgt_legend)

            %Now do glm
            tbl = table(glm_div.data',glm_div.trial_type',glm_div.time',...
                'VariableNames',{'dFF','trial_type','time'});
            mdl = fitglm(tbl,'dFF~trial_type+time'...
                ,'CategoricalVars',[2])
            
            %Save the p values
            ii_all_pvalues=ii_all_pvalues+1;
            all_ROI_files(ii_all_pvalues)=fileNo;
            all_iiROIs(ii_all_pvalues)=ii_ROI;

            p_values(ii_all_pvalues,1)=coefTest(mdl,[0 1 0 0 0]);
            p_values(ii_all_pvalues,2)=coefTest(mdl,[0 0 1 0 0]);
            p_values(ii_all_pvalues,3)=coefTest(mdl,[0 0 0 1 0]);
            p_values(ii_all_pvalues,4)=coefTest(mdl,[0 -1 1 0 0]);
            p_values(ii_all_pvalues,5)=coefTest(mdl,[0 -1 0 1 0]);
            p_values(ii_all_pvalues,6)=coefTest(mdl,[0 0 -1 1 0]);

            these_p_values_per_ROI(ii_ROI,1)=coefTest(mdl,[0 1 0 0 0]);
            these_p_values_per_ROI(ii_ROI,2)=coefTest(mdl,[0 0 1 0 0]);
            these_p_values_per_ROI(ii_ROI,3)=coefTest(mdl,[0 0 0 1 0]);
            these_p_values_per_ROI(ii_ROI,4)=coefTest(mdl,[0 1 -1 0 0]);
            these_p_values_per_ROI(ii_ROI,5)=coefTest(mdl,[0 1 0 -1 0]);
            these_p_values_per_ROI(ii_ROI,6)=coefTest(mdl,[0 0 -1 1 0]);


            fprintf(1, '\n')
            for ii_comp=1:6
                fprintf(1, ['\nFor ' trial_type_comp_labels{ii_comp} ' p value: '...
                    num2str(p_values(ii_all_pvalues,ii_comp)) '\n'])
            end

 
            % 
            % if sig_all_info_lane1(ii_ROI_all)==1
            %     pffft=1;
            % end
            % 
            % if sig_all_info_lane4(ii_ROI_all)==1
            %     pffft=1;
            % end
            % 
            % if sig_all_info_mutual_info14(ii_ROI_all)==1
            %     pfft=1; %MI is significant
            % end
            % 
            % if sum(ii_ROI==imps.file(fileNo).ROIs_conc)>0
            %     pfft=1;
            % end
            % 
            % if sum(ii_ROI==imps.file(fileNo).ROIs_x)>0
            %     pfft=1;
            % end
            % 
            % if sum(ii_ROI==imps.file(fileNo).ROIs_y)>0
            %     pffft=1;
            % end

            if (fileNo==7)&(ii_ROI==154)
                pffft=1;
            end

            if (fileNo==13)&(ii_ROI==124)
                pffft=1;
            end

            if (fileNo==14)&(ii_ROI==8)
                pffft=1;
            end

        end

        %Get FDR p values
        fprintf(1, ['\n\n\nFor file No ' num2str(fileNo) ' the total number of neurons is ' num2str(no_neurons) '\n\n'])
        pFDRs_per_file(fileNo)=drsFDRpval(these_p_values_per_ROI(:));
        handles_out2.file(fileNo).pFDR=pFDRs_per_file(fileNo);
        fprintf(1, ['\n p FDR: ' num2str(pFDRs_per_file(fileNo)) '\n'])

        for ii_comp=1:6
            these_p_values=zeros(no_neurons,1);
            these_p_values(:,1)=these_p_values_per_ROI(:,ii_comp);
            handles_out2.file(fileNo).ii_comp(ii_comp).these_p_values=these_p_values;
            handles_out2.file(fileNo).ii_comp(ii_comp).significant=these_p_values<=pFDRs_per_file(fileNo);
            fprintf(1, ['\nFor ' trial_type_comp_labels{ii_comp}  ' number significant ' num2str(sum(these_p_values<=pFDRs_per_file(fileNo))) '\n'])
        end

    end
end


%Now show the divergent responses
sorted_handles_out=[];
no_clusters=[3 2 2 2 2 2];
for ii_comp=1:6

    %Get all divergent ROIs
    all_div_hit1=[];
    all_div_miss1=[];
    all_div_hit4=[];
    all_div_miss4=[];
    all_div_to_sort=[];

  

    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
             these_significant=handles_out2.file(fileNo).ii_comp(ii_comp).significant;
             no_neurons=length(these_significant);
             for ii_ROI=1:no_neurons
                 if these_significant(ii_ROI)==1
                     all_div_hit1=[all_div_hit1; handles_out2.file(fileNo).ROI(ii_ROI).hit1_dFF_turn];
                     all_div_miss1=[all_div_miss1; handles_out2.file(fileNo).ROI(ii_ROI).miss1_dFF_turn];
                     all_div_hit4=[all_div_hit4; handles_out2.file(fileNo).ROI(ii_ROI).hit4_dFF_turn];
                     all_div_miss4=[all_div_miss4; handles_out2.file(fileNo).ROI(ii_ROI).miss4_dFF_turn];
                     
                 end
             end
             %Sort using the two ii_types that were used in the comparison
             eval(['all_div_to_sort=[' these_adiv_names{trial_type_comp_ii_type1(ii_comp)} ' '...
                 these_adiv_names{trial_type_comp_ii_type2(ii_comp)} '];'])
        end
    end
    
    all_div=[all_div_hit1;all_div_miss1;all_div_hit4;all_div_miss4];

    % %Calculate the crosscorrelations
    % these_all_div_dFFspm=zeros(size(handles_out.all_div_dFFspm,1),sum(ROIs_included));
    % these_all_div_dFFspm(:,:)=handles_out.all_div_dFFspm(:,logical(ROIs_included));
    % 
    % these_all_div_dFFsplus=zeros(size(handles_out.all_div_dFFsplus,2),sum(ROIs_included));
    % these_all_div_dFFsplus(:,:)=handles_out.all_div_dFFsplus(logical(ROIs_included),:)';
    % 
    % these_all_div_dFFsminus=zeros(size(handles_out.all_div_dFFsminus,2),sum(ROIs_included));
    % these_all_div_dFFsminus(:,:)=handles_out.all_div_dFFsminus(logical(ROIs_included),:)';
    % 
    % these_all_delta_dFFsplus=zeros(1,sum(ROIs_included));
    % these_all_delta_dFFsplus(1,:)=handles_out.all_div_delta_dFFsplus(1,logical(ROIs_included));
    % 
    % these_all_delta_dFFsminus=zeros(1,sum(ROIs_included));
    % these_all_delta_dFFsminus(1,:)=handles_out.all_div_delta_dFFsminus(1,logical(ROIs_included));
    % 
    % these_all_div_t=zeros(1,sum(ROIs_included));
    % these_all_div_t(1,:)=handles_out.all_div_t(1,logical(ROIs_included));


    %             croscorr_traces=corrcoef(all_div_dFFspm);
    croscorr_traces=corrcoef(all_div_to_sort');

    Z = linkage(croscorr_traces,'complete','correlation');

    
    handles_out2.ii_clus(ii_comp).clusters = cluster(Z,'Maxclust',no_clusters(ii_comp));
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    %Do cutoff for no_clusters
    cutoff = median([Z(end-(no_clusters(ii_comp)-1),3) Z(end-(no_clusters(ii_comp)-2),3)]);
    [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
    set(H,'LineWidth',2)
    hFig=figure(figureNo);
    set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

    %re-sort the matrix
    for ii=1:size(croscorr_traces,1)
        for jj=1:size(croscorr_traces,1)
            perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
        end
    end



    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
    hold on
    pcolor(perm_croscorr_traces)
    colormap fire
    shading flat

    % caxis([-1  1])
    clim([-1 1])
    eval(['title([''Cross correlations for all ROIs ' trial_type_comp_labels{ii_comp}  '''])'])

    ROIs_included=size(all_div_to_sort,1);
    xlim([1 ROIs_included])
    ylim([1 ROIs_included])

    for ii_type=1:4
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
        hold on


        eval(['sorted_handles_out.ii_comp(ii_comp).' these_adiv_names{ii_type} '=[];'])
        eval(['sorted_handles_out.ii_comp(ii_comp).' these_adiv_names{ii_type} '= ' these_adiv_names{ii_type} '(outperm,:);'])

        % sorted_handles_out.ii_comp(ii_comp).all_div_hit1=[];
        % sorted_handles_out.ii_comp(ii_comp).all_div_hit1=all_div_hit1(outperm,:);
        
        eval(['ii_included=size(sorted_handles_out.ii_comp(ii_comp).' these_adiv_names{ii_type} ',1);'])
        
        % ii_included=size(sorted_handles_out.ii_comp(ii_comp).all_div_hit1,1);

        time_span_mat=repmat(turn_time_bins,ii_included,1);
        ROI_mat=repmat(1:ii_included,length(turn_time_bins),1)';

        eval(['pcolor(time_span_mat,ROI_mat,sorted_handles_out.ii_comp(ii_comp).' these_adiv_names{ii_type} ')'])
        % pcolor(time_span_mat,ROI_mat,sorted_handles_out.ii_comp(ii_comp).all_div_hit1)
        colormap fire
        shading flat

        % caxis([prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),1) prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),99.9)])
 
        caxis([prctile(all_div(:),1) prctile(all_div(:),99.9)])

        plot([0 0],[0 ii_included],'-r')


        xlim([-3 3])
        ylim([1 ii_included])
        eval(['title([''' trial_type_comp_labels{ii_comp} ': ' trial_type_labels{ii_type} '''])'])
        xlabel('Time (sec)')
        ylabel('ROI number')
    end
    



    %Plot the average timecourses per cluster
    do_std=0;
    these_ylim=[];
    for clus=1:no_clusters(ii_comp)
        %Report the number of ROIs per cluster
        fprintf(1, ['\nFor '  trial_type_labels{trial_type_comp_ii_type1(ii_comp)}...
            ' vs ' trial_type_labels{trial_type_comp_ii_type2(ii_comp)} ' cluster No '...
            num2str(clus) ' has ' num2str(sum(handles_out2.ii_clus(ii_comp).clusters==clus)) ' ROIs\n'])

        clusters=handles_out2.ii_clus(ii_comp).clusters;

        %Plot the average of this cluster for the type of trials compared
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.3 .2 .18 .18])


        hold on

        %get the dF/F

        % these_adiv_names{1}='all_div_hit1';
        % these_adiv_names{2}='all_div_miss1';
        % these_adiv_names{3}='all_div_hit4';
        % these_adiv_names{4}='all_div_miss1';

        this_adiv_name1=these_adiv_names{trial_type_comp_ii_type1(ii_comp)};
        this_adiv_name2=these_adiv_names{trial_type_comp_ii_type2(ii_comp)};


        %S-
        this_cluster_adiv_name2=[];
        ii_included=0;
        eval(['no_ROIs_this_clus= size(' this_adiv_name2 ',1);'])
        for ii=1:no_ROIs_this_clus
            if clusters(ii)==clus
                ii_included=ii_included+1;
                eval(['this_cluster_adiv_name2(ii_included,:)=' this_adiv_name2 '(ii,:);'])
            end
        end

        % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

        if do_std==0
            CIpv = bootci(1000, @mean, this_cluster_adiv_name2);
            meanpv=mean(this_cluster_adiv_name2,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
        else
            STDpv = std(this_cluster_adiv_name2);
            meanpv=mean(this_cluster_adiv_name2,1);
            CIpv(1,:)=STDpv;
            CIpv(2,:)=STDpv;
        end


        [hlpvl, hppvl] = boundedline(turn_time_bins,mean(this_cluster_adiv_name2), CIpv','cmap',[158/255 31/255 99/255]);

        %S+
        this_cluster_adiv_name1=[];
        ii_included=0;
        eval(['no_ROIs_this_clus= size(' this_adiv_name1 ',1);'])
        for ii=1:no_ROIs_this_clus
            if clusters(ii)==clus
                ii_included=ii_included+1;
                eval(['this_cluster_adiv_name1(ii_included,:)=' this_adiv_name1 '(ii,:);'])
            end
        end

        % dFF_timecourse_per_clus.group(grNo).cluster(clus).splus_timecourses=this_cluster_dFFsplus;

        if do_std==0
            CIpv = bootci(1000, @mean, this_cluster_adiv_name1);
            meanpv=mean(this_cluster_adiv_name1,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
        else
            STDpv = std(this_cluster_adiv_name1);
            meanpv=mean(this_adiv_name1,1);
            CIpv(1,:)=STDpv;
            CIpv(2,:)=STDpv;
        end


        [hlpvl, hppvl] = boundedline(turn_time_bins, mean(this_cluster_adiv_name1), CIpv','cmap',[0 114/255 178/255]);


        plot(turn_time_bins',mean(this_cluster_adiv_name2)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
        plot(turn_time_bins',mean(this_cluster_adiv_name1)','Color',[0 114/255 178/255],'LineWidth',1.5);


        % ylim([-0.5 1.5])
        xlim([-3 3])
        these_ylim=[these_ylim; ylim];

    
       
       
        xlabel('Time(sec)')
        ylabel('dFF')

       eval(['title([''' trial_type_comp_labels{ii_comp} ''' '' significantly diff, cluster no ''  num2str(clus) ])'])

       %Plot the average of this cluster for the other two types of trials
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.3 .2 .18 .18])


        hold on

        %get the dF/F

        % these_adiv_names{1}='all_div_hit1';
        % these_adiv_names{2}='all_div_miss1';
        % these_adiv_names{3}='all_div_hit4';
        % these_adiv_names{4}='all_div_miss1';

        this_adiv_name1=these_adiv_names{trial_type_comp_ii_type3(ii_comp)};
        this_adiv_name2=these_adiv_names{trial_type_comp_ii_type4(ii_comp)};


        %S-
        this_cluster_adiv_name1=[];
        ii_included=0;
        eval(['no_ROIs_this_clus= size(' this_adiv_name1 ',1);'])
        for ii=1:no_ROIs_this_clus
            if clusters(ii)==clus
                ii_included=ii_included+1;
                eval(['this_cluster_adiv_name1(ii_included,:)=' this_adiv_name1 '(ii,:);'])
            end
        end

        % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

        if do_std==0
            CIpv = bootci(1000, @mean, this_cluster_adiv_name1);
            meanpv=mean(this_cluster_adiv_name1,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
        else
            STDpv = std(this_cluster_adiv_name1);
            meanpv=mean(this_cluster_adiv_name1,1);
            CIpv(1,:)=STDpv;
            CIpv(2,:)=STDpv;
        end


        [hlpvl, hppvl] = boundedline(turn_time_bins,mean(this_cluster_adiv_name1), CIpv','cmap',[158/255 31/255 99/255]);

        %S+
        this_cluster_adiv_name2=[];
        ii_included=0;
        eval(['no_ROIs_this_clus= size(' this_adiv_name2 ',1);'])
        for ii=1:no_ROIs_this_clus
            if clusters(ii)==clus
                ii_included=ii_included+1;
                eval(['this_cluster_adiv_name2(ii_included,:)=' this_adiv_name2 '(ii,:);'])
            end
        end

        % dFF_timecourse_per_clus.group(grNo).cluster(clus).splus_timecourses=this_cluster_dFFsplus;

        if do_std==0
            CIpv = bootci(1000, @mean, this_cluster_adiv_name2);
            meanpv=mean(this_cluster_adiv_name2,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
        else
            STDpv = std(this_cluster_adiv_name2);
            meanpv=mean(this_adiv_name2,1);
            CIpv(1,:)=STDpv;
            CIpv(2,:)=STDpv;
        end


        [hlpvl, hppvl] = boundedline(turn_time_bins, mean(this_cluster_adiv_name2), CIpv','cmap',[0 114/255 178/255]);


        plot(turn_time_bins',mean(this_cluster_adiv_name1)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
        plot(turn_time_bins',mean(this_cluster_adiv_name2)','Color',[0 114/255 178/255],'LineWidth',1.5);
 

        % ylim([-0.5 1.5])
        xlim([-3 3])
        

        these_ylim=[these_ylim; ylim];

        
       
        xlabel('Time(sec)')
        ylabel('dFF')

       eval(['title([''' trial_type_no_comp_labels{ii_comp} ''' '' cluster no ''  num2str(clus) ])'])
          
    end

    
    fig_minus=2*no_clusters(ii_comp)-1;
    for clus=1:no_clusters(ii_comp)

        figNo=figureNo-fig_minus;
        figure(figNo)
        ylim([min(these_ylim(:)) max(these_ylim(:))])
        this_ylim=ylim;

         %Odor on markers
        plot([0 0],this_ylim,'-k')

        text(-2.5,this_ylim(1)+0.9*(this_ylim(2)-this_ylim(1)),trial_type_labels{trial_type_comp_ii_type2(ii_comp)},'Color',[158/255 31/255 99/255])
        text(-2.5,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),trial_type_labels{trial_type_comp_ii_type1(ii_comp)},'Color',[0 114/255 178/255])

        fig_minus=fig_minus-1;

        figNo=figureNo-fig_minus;
        figure(figNo)
        ylim([min(these_ylim(:)) max(these_ylim(:))])
        this_ylim=ylim;

        text(-2.5,this_ylim(1)+0.9*(this_ylim(2)-this_ylim(1)),trial_type_labels{trial_type_comp_ii_type3(ii_comp)},'Color',[158/255 31/255 99/255])
        text(-2.5,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),trial_type_labels{trial_type_comp_ii_type4(ii_comp)},'Color',[0 114/255 178/255])

         %Odor on markers
        plot([0 0],this_ylim,'-k')

        fig_minus=fig_minus-1;

    end

    pffft=1;



end

fclose(fileID);

pffft=1;

