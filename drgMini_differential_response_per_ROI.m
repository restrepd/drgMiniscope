%drgMini_differential_response_per_ROI
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

        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

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

%Now plot the space activity maps for each of the high prediction
%importance ROIs
figNo=figureNo+1;
x=25:50:475;
y=24:48:456;

all_info_ii=0;
all_info_fileNo=[];
all_info_ii_ROI=[];
all_info_lane1=[];
all_info_lane4=[];
all_info_mutual_info14=[];

all_info_lane1_sh=[];
all_info_lane4_sh=[];
all_info_mutual_info14_sh=[];

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Get XY and dFF per trial
        arena_file=handles_XY.arena_file{fileNo};
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        no_neurons=handles_out.no_neurons;


        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        angles=handles_out.angles;

 

        %Now show the spatial activity maps
        no_place_cells=0;
        place_cells=[];
        no_lane_trial_cells=0;
        lane_trial_cells=[];
        for ii_ROI=1:no_neurons

 

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

            %Cumulative counts for info
            cum_xy=zeros(1,10*10);
            cum_lane=zeros(1,2);
            cum_bindFF=zeros(1,2);

            cum_xy_bindFF_Lane1=zeros(2,10*10);
            cum_xy_bindFF_Lane4=zeros(2,10*10);
            cum_xy_bindFF_BothLanes=zeros(2,10*10);

            cum_xy_bindFF_lane=zeros(2,10*10,2);

            cum_xy_Lane1=zeros(1,10*10);
            cum_xy_Lane4=zeros(1,10*10);
            cum_xy_BothLanes=zeros(1,10*10);

            cum_bindFF_Lane1=zeros(1,2);
            cum_bindFF_Lane4=zeros(1,2);
            cum_bindFF_BothLanes=zeros(1,2);

            all_these_x=[];
            all_these_y=[];
            all_these_dFF=[];
            all_lanes=[];

 
            for trNo=1:trials.odor_trNo
                these_x=trials.trial(trNo).XYtest(:,1);
                these_y=trials.trial(trNo).XYtest(:,2);
                these_dFF=trials.trial(trNo).XdFFtest(:,ii_ROI);

                all_these_x=[all_these_x; these_x];
                all_these_y=[all_these_y; these_y];
                all_these_dFF=[all_these_dFF; these_dFF];
                all_lanes=[all_lanes; trials.lane_per_trial(trNo)*ones(length(these_x),1)];

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
                    cum_bindFF(this_bin_dFF)=cum_bindFF(this_bin_dFF)+1;
                    cum_bindFF_BothLanes(this_bin_dFF)=cum_bindFF_BothLanes(this_bin_dFF)+1;
                    this_xy_ii=this_x_ii+10*(this_y_ii-1);
                    cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
                    cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)+1;
                    cum_xy_BothLanes(this_xy_ii)=cum_xy_BothLanes(this_xy_ii)+1;

                    if trials.lane_per_trial(trNo)==1
                        % this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl1_activity=sum_dFFl1_activity+these_dFF(ii_t);
                        cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)+1;
                        cum_xy_Lane1(this_xy_ii)=cum_xy_Lane1(this_xy_ii)+1;
                        cum_bindFF_Lane1(this_bin_dFF)=cum_bindFF_Lane1(this_bin_dFF)+1;
                        cum_lane(1)=cum_lane(1)+1;
                    else
                        % this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl4_activity=sum_dFFl4_activity+these_dFF(ii_t);
                        cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)+1;
                        cum_xy_Lane4(this_xy_ii)=cum_xy_Lane4(this_xy_ii)+1;
                        cum_bindFF_Lane4(this_bin_dFF)=cum_bindFF_Lane4(this_bin_dFF)+1;
                        cum_lane(2)=cum_lane(2)+1;
                    end

                end
            end

            %Calculate information

            p_xy=cum_xy/sum(cum_xy(:));
            p_lane=cum_lane/sum(cum_lane(:));
            p_bindFF=cum_bindFF/sum(cum_bindFF(:));

            p_xy_bindFF_Lane1=cum_xy_bindFF_Lane1/sum(cum_xy_bindFF_Lane1(:));
            p_xy_bindFF_Lane4=cum_xy_bindFF_Lane4/sum(cum_xy_bindFF_Lane4(:));
            p_xy_bindFF_BothLanes=cum_xy_bindFF_BothLanes/sum(cum_xy_bindFF_BothLanes(:));

            p_xy_bindFF_lane=cum_xy_bindFF_lane/sum(cum_xy_bindFF_lane(:));

            p_xy_Lane1=cum_xy_Lane1/sum(cum_xy_Lane1(:));
            p_xy_Lane4=cum_xy_Lane4/sum(cum_xy_Lane4(:));
            p_xy_BothLanes=cum_xy_BothLanes/sum(cum_xy_BothLanes(:));

            p_bindFF_Lane1=cum_bindFF_Lane1/sum(cum_bindFF_Lane1(:));
            p_bindFF_Lane4=cum_bindFF_Lane4/sum(cum_bindFF_Lane4(:));
            p_bindFF_BothLanes=cum_bindFF_BothLanes/sum(cum_bindFF_BothLanes(:));

            all_info_ii=all_info_ii+1;
            all_info_fileNo(all_info_ii)=fileNo;
            all_info_ii_ROI(all_info_ii)=ii_ROI;

            all_info_lane1(all_info_ii)=0;
            all_info_lane4(all_info_ii)=0;
            all_info_mutual_info14(all_info_ii)=0;

            %Calculate info for lane 1
            for ii_xy=1:10*10
                if sum(p_xy_bindFF_Lane1(:,ii_xy))>0
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)~=0
                            all_info_lane1(all_info_ii)=all_info_lane1(all_info_ii)+p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)*...
                                log2(p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)/(p_xy_Lane1(ii_xy)*p_bindFF_Lane1(ii_bin_dFF)));
                        end
                    end
                end
            end

            %Calculate info for lane 4
            for ii_xy=1:10*10
                if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)~=0
                            all_info_lane4(all_info_ii)=all_info_lane4(all_info_ii)+p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)*...
                                log2(p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)/(p_xy_Lane4(ii_xy)*p_bindFF_Lane4(ii_bin_dFF)));
                        end
                    end
                end
            end

            %Calculate mutual info between lanes 1 and 4
            for ii_xy=1:10*10
                if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_lane=1:2
                        for ii_bin_dFF=1:2
                            if p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)~=0
                                all_info_mutual_info14(all_info_ii)=all_info_mutual_info14(all_info_ii)+p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)*...
                                    log2(p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_lane(ii_lane)));
                            end
                        end
                    end
                end
            end

            %Now calculate information content with the shuffled dFF
            rng(ii_ROI)
            these_rnd=rand(1,n_shuffle_SI);

            for ii_sh=1:n_shuffle_SI

                %We will do a reversal and a circular permutation of XY
                these_x_reversed=zeros(1,length(all_these_x));
                these_y_reversed=zeros(1,length(all_these_y));

                %Flip and roll Stefanini et al 2020 https://doi.org/10.1016/j.neuron.2020.05.022
                offset_ii=floor(these_rnd(ii_sh)*length(all_these_x));
                for ii_trl=1:length(all_these_x)
                    this_ii_trl=ii_trl+offset_ii;
                    if this_ii_trl>length(all_these_x)
                        offset_ii=-ii_trl+1;
                        this_ii_trl=ii_trl+offset_ii;
                    end
                    these_x_reversed(1,length(all_these_x)-ii_trl+1)=all_these_x(this_ii_trl);
                    these_y_reversed(1,length(all_these_y)-ii_trl+1)=all_these_y(this_ii_trl);
                end

                %Cumulative counts for info
                cum_xy=zeros(1,10*10);
                cum_lane=zeros(1,2);
                cum_bindFF=zeros(1,2);

                cum_xy_bindFF_Lane1=zeros(2,10*10);
                cum_xy_bindFF_Lane4=zeros(2,10*10);
                cum_xy_bindFF_BothLanes=zeros(2,10*10);

                cum_xy_bindFF_lane=zeros(2,10*10,2);

                cum_xy_Lane1=zeros(1,10*10);
                cum_xy_Lane4=zeros(1,10*10);
                cum_xy_BothLanes=zeros(1,10*10);

                cum_bindFF_Lane1=zeros(1,2);
                cum_bindFF_Lane4=zeros(1,2);
                cum_bindFF_BothLanes=zeros(1,2);


                for ii_t=1:length(these_x_reversed)
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
                    this_bin_dFF=(all_these_dFF(ii_t)>0)+1;
                    cum_bindFF(this_bin_dFF)=cum_bindFF(this_bin_dFF)+1;
                    cum_bindFF_BothLanes(this_bin_dFF)=cum_bindFF_BothLanes(this_bin_dFF)+1;
                    this_xy_ii=this_x_ii+10*(this_y_ii-1);
                    cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
                    cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)+1;
                    cum_xy_BothLanes(this_xy_ii)=cum_xy_BothLanes(this_xy_ii)+1;

                    if all_lanes(ii_t)==1
                        % this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl1_activity=sum_dFFl1_activity+these_dFF(ii_t);
                        cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)+1;
                        cum_xy_Lane1(this_xy_ii)=cum_xy_Lane1(this_xy_ii)+1;
                        cum_bindFF_Lane1(this_bin_dFF)=cum_bindFF_Lane1(this_bin_dFF)+1;
                        cum_lane(1)=cum_lane(1)+1;
                    else
                        % this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl4_activity=sum_dFFl4_activity+these_dFF(ii_t);
                        cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)+1;
                        cum_xy_Lane4(this_xy_ii)=cum_xy_Lane4(this_xy_ii)+1;
                        cum_bindFF_Lane4(this_bin_dFF)=cum_bindFF_Lane4(this_bin_dFF)+1;
                        cum_lane(2)=cum_lane(2)+1;
                    end

                end

                %Calculate information

                p_xy=cum_xy/sum(cum_xy(:));
                p_lane=cum_lane/sum(cum_lane(:));
                p_bindFF=cum_bindFF/sum(cum_bindFF(:));

                p_xy_bindFF_Lane1=cum_xy_bindFF_Lane1/sum(cum_xy_bindFF_Lane1(:));
                p_xy_bindFF_Lane4=cum_xy_bindFF_Lane4/sum(cum_xy_bindFF_Lane4(:));
                p_xy_bindFF_BothLanes=cum_xy_bindFF_BothLanes/sum(cum_xy_bindFF_BothLanes(:));

                p_xy_bindFF_lane=cum_xy_bindFF_lane/sum(cum_xy_bindFF_lane(:));

                p_xy_Lane1=cum_xy_Lane1/sum(cum_xy_Lane1(:));
                p_xy_Lane4=cum_xy_Lane4/sum(cum_xy_Lane4(:));
                p_xy_BothLanes=cum_xy_BothLanes/sum(cum_xy_BothLanes(:));

                p_bindFF_Lane1=cum_bindFF_Lane1/sum(cum_bindFF_Lane1(:));
                p_bindFF_Lane4=cum_bindFF_Lane4/sum(cum_bindFF_Lane4(:));
                p_bindFF_BothLanes=cum_bindFF_BothLanes/sum(cum_bindFF_BothLanes(:));

                all_info_lane1_sh(all_info_ii,ii_sh)=0;
                all_info_lane4_sh(all_info_ii,ii_sh)=0;
                all_info_mutual_info14_sh(all_info_ii,ii_sh)=0;

                %Calculate info for lane 1
                for ii_xy=1:10*10
                    if sum(p_xy_bindFF_Lane1(:,ii_xy))>0
                        for ii_bin_dFF=1:2
                            if p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)~=0
                                all_info_lane1_sh(all_info_ii,ii_sh)=all_info_lane1_sh(all_info_ii,ii_sh)+p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)*...
                                    log2(p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)/(p_xy_Lane1(ii_xy)*p_bindFF_Lane1(ii_bin_dFF)));
                            end
                        end
                    end
                end

                %Calculate info for lane 4
                for ii_xy=1:10*10
                    if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                        for ii_bin_dFF=1:2
                            if p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)~=0
                                all_info_lane4_sh(all_info_ii,ii_sh)=all_info_lane4_sh(all_info_ii,ii_sh)+p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)*...
                                    log2(p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)/(p_xy_Lane4(ii_xy)*p_bindFF_Lane4(ii_bin_dFF)));
                            end
                        end
                    end
                end

                %Calculate mutual info between lanes 1 and 4
                for ii_xy=1:10*10
                    if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                        for ii_lane=1:2
                            for ii_bin_dFF=1:2
                                if p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)~=0
                                    all_info_mutual_info14_sh(all_info_ii,ii_sh)=all_info_mutual_info14_sh(all_info_ii,ii_sh)+p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)*...
                                        log2(p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_lane(ii_lane)));
                                end
                            end
                        end
                    end
                end

            end

            % %Calculate teh space maps
            % for ii_x=1:10
            %     for ii_y=1:10
            %         if this_dFF_activity_n(ii_x,ii_y)~=0
            %             this_dFF_activity(ii_x,ii_y)=this_dFF_activity(ii_x,ii_y)/this_dFF_activity_n(ii_x,ii_y);
            %         end
            %     end
            % end
            %
            % for ii_x=1:10
            %     for ii_y=1:10
            %         if this_dFFl1_activity_n(ii_x,ii_y)~=0
            %             this_dFFl1_activity(ii_x,ii_y)=this_dFFl1_activity(ii_x,ii_y)/this_dFFl1_activity_n(ii_x,ii_y);
            %         end
            %     end
            % end
            %
            % for ii_x=1:10
            %     for ii_y=1:10
            %         if this_dFFl4_activity_n(ii_x,ii_y)~=0
            %             this_dFFl4_activity(ii_x,ii_y)=this_dFFl4_activity(ii_x,ii_y)/this_dFFl4_activity_n(ii_x,ii_y);
            %         end
            %     end
            % end
            %
            % %Plot the space activity maps
            % try
            %     close(figNo)
            % catch
            % end
            %
            %
            % hFig = figure(figNo);
            % set(hFig, 'units','normalized','position',[.1 .1 .75 .75])



            % %Plot dFF timecourses
            % subplot(2,3,1:3)
            % hold on
            % ii_plot=0;

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
                    pffft=1;
                end
            end

            delta_below_zero_ii=max(all_ii_turns(include_trial==1));
            delta_above_zero_ii=max(all_ii_ends(include_trial==1)-all_ii_turns(include_trial==1));
            time_bins=handles_XY.dt*([1:delta_below_zero_ii+delta_above_zero_ii]-delta_below_zero_ii);

            hit1_dFF=[];
            ii_hit1=0;
            miss1_dFF=[];
            ii_miss1=0;
            hit4_dFF=[];
            ii_hit4=0;
            miss4_dFF=[];
            ii_miss4=0;

            %Let's do the glm stats from -3 to 3 sec in 1 sec bins
            trimmed_time_range=[-3 3];
            trimmed_dt=1;
            trimmed_time_bins=[trimmed_time_range(1)+(trimmed_dt/2):trimmed_dt:trimmed_time_range(2)-(trimmed_dt/2)];

            glm_div_ii=0;
            glm_div=[];

            for trNo=1:trials.odor_trNo
                if include_trial(trNo)==1
                    these_dFF=trials.trial(trNo).XdFFtest(:,ii_ROI);
                    these_time_bins=time_bins(delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)));
                    %Okabe_Ito colors
                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            ii_hit1=ii_hit1+1;
                            hit1_dFF(ii_hit1,1:length(time_bins))=min(these_dFF);
                            hit1_dFF(ii_hit1,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
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
                end
            end

           


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

            y_trials_end=max([y_trials1(end) y_trials4(end)])

            % %Plot dFF for lane 1
            % subplot(2, 6, [1 2 3]);
            % hold on
            % 
            % drg_pcolor(repmat(time_bins,length(y_trials1),1),repmat(y_trials1,length(time_bins),1)',hit_miss_1)
            % colormap(this_cmap)
            % if max(hit_miss_1(:))>0
            %     clim([0 max(hit_miss_1(:))]);
            % end
            % 
            % shading flat
            % plot([0 0],[y_trials1(1) y_trials1(end)+1],'-w','LineWidth',3)
            % 
            % %y_trials1=1 to 2 is the first miss
            % %size(miss1_dFF,1) to size(miss1_dFF,1) + 1 is the last miss,
            % %size(miss1_dFF,1)+1 to size(miss1_dFF,1)+1+y_gap is the in between miss
            % %and hit
            % %size(miss1_dFF,1)+1+y_gap to size(miss1_dFF,1)+1+y_gap+size(hit1_dFF,1) is
            % %hit
            % 
            % rectangle('Position', [time_bins(1), size(miss1_dFF,1)+1, 1.03*(time_bins(end)-time_bins(1)), y_gap], ... % [x, y, width, height]
            %     'FaceColor', 'white', ...    % Fill color (white)
            %     'EdgeColor', 'none');        % No border for the rectangle
            % 
            % ylim([0 y_trials_end+2])
            % xlim([time_bins(1)-0.05*(time_bins(end)-time_bins(1)) time_bins(end)+0.05*(time_bins(end)-time_bins(1)) ])
            % xlabel('Time (sec)')
            % yticks([1+(size(miss1_dFF,1)/2) size(miss1_dFF,1)+1+y_gap+(size(hit1_dFF,1)/2)])
            % yticklabels({'Misses','Hits'})
            % title(['Lane 1 dFF'])
            % 
            % %Plot dFF for lane 4
            % subplot(2, 6, [4 5 6]);
            % hold on
            % 
            % drg_pcolor(repmat(time_bins,length(y_trials4),1),repmat(y_trials4,length(time_bins),1)',hit_miss_4)
            % colormap(this_cmap)
            % if max(hit_miss_4(:))>0
            %     clim([0 max(hit_miss_4(:))])
            % end
            % shading flat
            % plot([0 0],[y_trials4(1) y_trials4(end)+1],'-w','LineWidth',3)
            % 
            % rectangle('Position', [time_bins(1), size(miss4_dFF,1)+1, 1.03*(time_bins(end)-time_bins(1)), y_gap], ... % [x, y, width, height]
            %     'FaceColor', 'white', ...    % Fill color (white)
            %     'EdgeColor', 'none');        % No border for the rectangle
            % ylim([0 y_trials_end+2])
            % xlim([time_bins(1)-0.05*(time_bins(end)-time_bins(1)) time_bins(end)+0.05*(time_bins(end)-time_bins(1)) ])
            % xlabel('Time (sec)')
            % yticks([1+(size(miss4_dFF,1)/2) size(miss4_dFF,1)+1+y_gap+(size(hit4_dFF,1)/2)])
            % yticklabels({'Misses','Hits'})
            % title(['Lane 4 dFF'])
            % 
            % %Plot space activity maps as in Moser https://doi.org/10.1126/science.1114037
            % max_this_dFF_activity=max(this_dFF_activity(:));
            % max_this_dFFl4_activity=max(this_dFFl4_activity(:));
            % max_this_dFFl1_activity=max(this_dFFl1_activity(:));
            % 
            % colormap fire
            % this_cmap=colormap;
            % this_cmap(1,:)=[0.3 0.3 0.3];
            % 
            % if max_this_dFFl1_activity>max_this_dFFl4_activity
            %     %Lane 1
            %     subplot(2, 6, [7 8]);
            %     max_activity=max_this_dFFl1_activity;
            %     delta_ac=max_activity/255;
            %     if delta_ac==0
            %         delta_ac=0.000001;
            %     end
            %     this_masked_dFFl1_activity=this_dFFl1_activity;
            %     this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
            %     drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
            %     colormap(this_cmap)
            %     clim([-1.5*delta_ac max_activity])
            %     shading interp
            %     set(gca, 'YDir', 'reverse');
            % 
            %     yticks(0:48:480)
            %     xticks(0:50:500)
            %     xlabel('x (mm)')
            %     ylabel('y (mm)')
            %     title(['Lane 1, SSI= ' num2str(all_info_lane1(all_info_ii))])
            % 
            %     %Lane 4 normalized to lane 1
            %     subplot(2, 6, [9 10]);
            %     max_activity=max_this_dFFl1_activity;
            %     this_masked_dFFl4_activity=this_dFFl4_activity;
            %     this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
            %     drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
            %     colormap(this_cmap)
            %     clim([-1.5*delta_ac max_activity])
            %     shading interp
            %     set(gca, 'YDir', 'reverse');
            % 
            %     yticks(0:48:480)
            %     xticks(0:50:500)
            %     xlabel('x (mm)')
            %     ylabel('y (mm)')
            %     title(['Lane 4 normalized to Lane 1'])
            % 
            %     %Lane 4 normalized to lane 4
            %     subplot(2, 6, [11 12]);
            %     max_activity=max_this_dFFl4_activity;
            %     delta_ac=max_activity/255;
            %     if delta_ac==0
            %         delta_ac=0.000001;
            %     end
            %     this_masked_dFFl4_activity=this_dFFl4_activity;
            %     this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
            %     drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
            %     colormap(this_cmap)
            %     clim([-1.5*delta_ac max_activity])
            %     shading interp
            %     set(gca, 'YDir', 'reverse');
            % 
            %     yticks(0:48:480)
            %     xticks(0:50:500)
            %     xlabel('x (mm)')
            %     ylabel('y (mm)')
            %     title(['Lane 4 , SSI= ' num2str(all_info_lane4(all_info_ii))])
            % 
            % else
            % 
            % 
            %     %Lane 4 normalized to lane 4
            %     subplot(2, 6, [7 8]);
            %     max_activity=max_this_dFFl4_activity;
            %     delta_ac=max_activity/255;
            %     if delta_ac==0
            %         delta_ac=0.000001;
            %     end
            %     this_masked_dFFl4_activity=this_dFFl4_activity;
            %     this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
            %     drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
            %     colormap(this_cmap)
            %     clim([-1.5*delta_ac max_activity])
            %     shading interp
            %     set(gca, 'YDir', 'reverse');
            % 
            %     yticks(0:48:480)
            %     xticks(0:50:500)
            %     xlabel('x (mm)')
            %     ylabel('y (mm)')
            %     title(['Lane 4, SSI= ' num2str(all_info_lane4(all_info_ii))])
            % 
            %     %Lane 1 normalized to lane 4
            %     subplot(2, 6, [9 10]);
            %     max_activity=max_this_dFFl4_activity;
            %     this_masked_dFFl1_activity=this_dFFl1_activity;
            %     this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
            %     drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
            %     colormap(this_cmap)
            %     clim([-1.5*delta_ac max_activity])
            %     shading interp
            %     set(gca, 'YDir', 'reverse');
            % 
            %     yticks(0:48:480)
            %     xticks(0:50:500)
            %     xlabel('x (mm)')
            %     ylabel('y (mm)')
            %     title(['Lane 1 normalized to Lane 4'])
            % 
            %     %Lane 1 normalized to lane 1
            %     subplot(2, 6, [11 12]);
            %     max_activity=max_this_dFFl1_activity;
            %     delta_ac=max_activity/255;
            %     if delta_ac==0
            %         delta_ac=0.00001;
            %     end
            %     this_masked_dFFl1_activity=this_dFFl1_activity;
            %     this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
            %     drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
            %     colormap(this_cmap)
            %     clim([-1.5*delta_ac max_activity])          
            %     shading interp
            %     set(gca, 'YDir', 'reverse');
            % 
            %     yticks(0:48:480)
            %     xticks(0:50:500)
            %     xlabel('x (mm)')
            %     ylabel('y (mm)')
            %     title(['Lane 1, SSI= ' num2str(all_info_lane1(all_info_ii))])
            % end
            % 
            % 
            % sgt_legend=['dFF map file  No ' num2str(fileNo) ' ROI No ' num2str(ii_ROI) ' MI= ' num2str(all_info_mutual_info14(all_info_ii))];
            % 
            % % if sum(this_ROI==imps.file(fileNo).ROIs_conc)>0
            % %     sgt_legend=[sgt_legend ' odor '];
            % % end
            % % 
            % % if sum(this_ROI==imps.file(fileNo).ROIs_x)>0
            % %     sgt_legend=[sgt_legend ' x '];
            % % end
            % % 
            % % if sum(this_ROI==imps.file(fileNo).ROIs_y)>0
            % %     sgt_legend=[sgt_legend ' y '];
            % % end
            % 
            % sgtitle(sgt_legend)

            %Now do glm
            tbl = table(glm_div.data',glm_div.trial_type',glm_div.time',...
                'VariableNames',{'dFF','trial_type','time'});
            mdl = fitglm(tbl,'dFF~trial_type+time'...
                ,'CategoricalVars',[2])
            pffft=1;

        end
    end
end

%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .4 .4])
hold on

[f_aic,x_aic] = drg_ecdf(all_info_lane1);
plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_info_lane4);
plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_info_mutual_info14);
plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)

these_ylim=ylim;
these_xlim=xlim;
text(0.2,0.7,'SSI Lane 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(0.2,0.65,'SSI Lane 4','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(0.2,0.60,'MI shared by lanes 1 and 4','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)


title('Cumulative histograms')
xlabel('Information bits')
ylabel('Cumulative histogram')


data = [all_info_lane1; all_info_lane4; all_info_mutual_info14];

% Step 1 Transpose the data if necessary (N samples x 3 variables)
data = data';

rng('default') % for reproducibility

% Step 2: Run t-SNE
% The default parameters are usually sufficient, but you can adjust them
% For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
% mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
mappedX = tsne(data,'Algorithm','exact','Distance','cosine'); %works better

% mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
% mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK

% Step 3: Plot the results
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE Information Content ' ]);

%Cool, it looks like we have three clear clusters

% Perform k-means clustering
%kmeans misclassifies a small number of points
% k = 3; % Number of clusters
% [idx, centroids] = kmeans(mappedX, k,'Distance','cityblock');

gm = fitgmdist(mappedX, 3); % Assuming 3 clusters
idx = cluster(gm, mappedX);


%Report the mean infos for each cluster
for ii_k=1:3
    mean_lane1=mean(all_info_lane1(idx==ii_k));
    mean_lane4=mean(all_info_lane4(idx==ii_k));
    mean_mutual14=mean(all_info_mutual_info14(idx==ii_k));

    fprintf(1, ['\nInformation content for cluster ' num2str(ii_k) ' lane1, lane 2, mutual: '...
        num2str(mean_lane1) ' ' num2str(mean_lane4) ' ' num2str(mean_mutual14) '\n'])
end

% Step 3: Plot the results
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
for ii_k=1:3
    switch ii_k
        case 1
            plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        case 2
            plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        case 3
            plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
    end
end

text(-60,70,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(-60,60,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(-60,50,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE Information Content ' ]);

%Now calculate 3xSD for the shuffled distributions
sd_all_info_lane1_sh=zeros(all_info_ii,1);
sd_all_info_lane4_sh=zeros(all_info_ii,1);
sd_all_info_mutual_info14_sh=zeros(all_info_ii,1);

mean_all_info_lane1_sh=zeros(all_info_ii,1);
mean_all_info_lane4_sh=zeros(all_info_ii,1);
mean_all_info_mutual_info14_sh=zeros(all_info_ii,1);

all_all_info_lane1_sh=[];
all_all_info_lane4_sh=[];
all_all_info_mutual_info14_sh=[];

for ii_ROI=1:all_info_ii

    this_ii_ROI=ii_ROI

    this_all_info_lane1_sh=zeros(1,n_shuffle_SI);
    this_all_info_lane1_sh(1,:)=all_info_lane1_sh(ii_ROI,:);
    this_all_info_lane4_sh=zeros(1,n_shuffle_SI);
    this_all_info_lane4_sh(1,:)=all_info_lane4_sh(ii_ROI,:);
    this_all_info_mutual_info14_sh=zeros(1,n_shuffle_SI);
    this_all_info_mutual_info14_sh(1,:)=all_info_mutual_info14_sh(ii_ROI,:);

    all_all_info_lane1_sh=[all_all_info_lane1_sh this_all_info_lane1_sh];
    all_all_info_lane4_sh=[all_all_info_lane4_sh this_all_info_lane4_sh];
    all_all_info_mutual_info14_sh=[all_all_info_mutual_info14_sh this_all_info_mutual_info14_sh];

    for jj=1:3
        switch jj
            case 1
                data=this_all_info_lane1_sh;
            case 2
                data=this_all_info_lane4_sh;
            case 3
                data=this_all_info_mutual_info14_sh;
        end
        n_bootstrap = 1000; % Number of bootstrap samples
        n_samples = length(data);

        bootstrap_std = zeros(n_bootstrap, 1);

        for i = 1:n_bootstrap
            % Generate bootstrap sample
            bootstrap_sample = datasample(data, n_samples, 'Replace', true);

            % Calculate standard deviation of the bootstrap sample
            bootstrap_std(i) = std(bootstrap_sample);
        end

        % Estimate of the standard deviation
        estimated_std = mean(bootstrap_std);
        switch jj
            case 1
                sd_all_info_lane1_sh(ii_ROI)=estimated_std;
            case 2
                sd_all_info_lane4_sh(ii_ROI)=estimated_std;
            case 3
                sd_all_info_mutual_info14_sh(ii_ROI)=estimated_std;
        end
    end

    mean_all_info_lane1_sh(ii_ROI)=mean(this_all_info_lane1_sh);
    mean_all_info_lane4_sh(ii_ROI)=mean(this_all_info_lane4_sh);
    mean_all_info_mutual_info14_sh(ii_ROI)=mean(this_all_info_mutual_info14_sh);

end

%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

[f_aic,x_aic] = drg_ecdf(all_info_lane1);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_all_info_lane1_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)

text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)

title('Cumulative histograms info lane 1')
xlabel('Information bits')
ylabel('Cumulative histogram')



%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


[f_aic,x_aic] = drg_ecdf(all_info_lane4);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_all_info_lane4_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)


text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)

title('Cumulative histograms info lane 4')
xlabel('Information bits')
ylabel('Cumulative histogram')

%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on
 

[f_aic,x_aic] = drg_ecdf(all_info_mutual_info14);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_all_info_mutual_info14_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)

text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)


title('Cumulative histograms mutual info lanes 1 and 4')
xlabel('Information bits')
ylabel('Cumulative histogram')

%plot in color those significantly above SD
% Step 3: Plot the results
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
sig_all_info_lane1=zeros(1,all_info_ii);
for ii_ROI=1:all_info_ii
    if all_info_lane1(ii_ROI)>=3*sd_all_info_lane1_sh(ii_ROI)+mean_all_info_lane1_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idx(ii_ROI);
        sig_all_info_lane1(ii_ROI)=1;
        switch ii_k
            case 1
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
            case 2
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
            case 3
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        end
    end
end

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE Information Content,info 1 significant ' ]);

%plot in color those significantly above SD
% Step 3: Plot the results
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
sig_all_info_lane4=zeros(1,all_info_ii);
for ii_ROI=1:all_info_ii
    if all_info_lane4(ii_ROI)>=3*sd_all_info_lane4_sh(ii_ROI)+mean_all_info_lane4_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idx(ii_ROI);
        sig_all_info_lane4(ii_ROI)=1;
        switch ii_k
            case 1
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
            case 2
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
            case 3
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        end
    end
end

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE Information Content,info 4 significant ' ]);


%plot in color those significantly above SD
% Step 3: Plot the results
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
sig_all_info_mutual_info14=zeros(1,all_info_ii);
plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
for ii_ROI=1:all_info_ii
    if all_info_mutual_info14(ii_ROI)>=3*sd_all_info_mutual_info14_sh(ii_ROI)+mean_all_info_mutual_info14_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idx(ii_ROI);
        sig_all_info_mutual_info14(ii_ROI)=1;
        switch ii_k
            case 1
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
            case 2
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
            case 3
                plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        end
    end
end

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE Information Content,mutual info 1,4 significant ' ]);

%Now plot the fraction of significant information for each cluster
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
bar_offset=0;

%lane 1 spatial information 
for ii_k=1:3
    this_fraction=sum(sig_all_info_lane1(idx==ii_k))/sum(idx==ii_k);
       switch ii_k
            case 1
                bar(bar_offset,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+4,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+8,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end
      
end



%lane 4 spatial information 
for ii_k=1:3
    this_fraction=sum(sig_all_info_lane4(idx==ii_k))/sum(idx==ii_k);
       switch ii_k
            case 1
                bar(bar_offset+1,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+5,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+9,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end
  
end



%lanes 1 and 4 mutual spatial information 
for ii_k=1:3
    this_fraction=sum(sig_all_info_mutual_info14(idx==ii_k))/sum(idx==ii_k);
       switch ii_k
            case 1
                bar(bar_offset+2,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+6,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+10,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end

end


text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2 4 5 6 8 9 10])
xticklabels({'Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4'})
xtickangle(45)

title(['Significant spatial information for the different clusters'])
ylabel('Fraction SSI')
ylim([0 1])
xlim([-1 11])

%Plot the spatial information for each cluster

%Now plot the Fusi spatial information for each cluster
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
bar_offset=0;
edges=[-0.1:0.25:42];
rand_offset=0.5;

%lane 1 spatial information
for ii_k=1:3
    these_SSI=all_info_lane1(idx==ii_k);
    switch ii_k
        case 1
            bar(bar_offset,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+4,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+4,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+8,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+8,rand_offset,'k','k',1);
    end
end



%lane 4 spatial information 
for ii_k=1:3
    these_SSI=all_info_lane4(idx==ii_k);
    switch ii_k
        case 1
            bar(bar_offset+1,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+1,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+5,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+5,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+9,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+9,rand_offset,'k','k',1);
    end
end



%both lanes mutual spatial information shared between 1 and 4
for ii_k=1:3
    these_SSI=all_info_mutual_info14(idx==ii_k);
    switch ii_k
        case 1
            bar(bar_offset+2,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+2,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+6,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+6,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+10,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
                ,edges,bar_offset+10,rand_offset,'k','k',1);
    end
end


text(4,14,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(4,13,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(4,12,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2 4 5 6 8 9 10])
xticklabels({'Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4'})
xtickangle(45)

title(['Spatial information for the different clusters'])
ylabel('SI')
ylim([0 0.5])
xlim([-1 11])

%Load information on prediction importance
load([save_PathPredImp save_FilePredImp])

%Now plot the Fusi spatial information for each cluster
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
bar_offset=0;
edges=[-0.1:0.25:42];
rand_offset=0.5;

%lane 1 spatial information
for ii_k=1:3
    these_Fusi_SSI=imps.all_Fusi_SSIl1(idx==ii_k);
    switch ii_k
        case 1
            bar(bar_offset,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+4,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+4,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+8,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+8,rand_offset,'k','k',1);
    end
end



%lane 4 spatial information 
for ii_k=1:3
    these_Fusi_SSI=imps.all_Fusi_SSIl4(idx==ii_k);
    switch ii_k
        case 1
            bar(bar_offset+1,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+1,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+5,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+5,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+9,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+9,rand_offset,'k','k',1);
    end
end



%both lanes mutual spatial information 
for ii_k=1:3
    these_Fusi_SSI=imps.all_Fusi_SSI(idx==ii_k);
    switch ii_k
        case 1
            bar(bar_offset+2,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+2,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+6,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+6,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+10,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
                ,edges,bar_offset+10,rand_offset,'k','k',1);
    end
end


text(4,14,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(4,13,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(4,12,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2 4 5 6 8 9 10])
xticklabels({'Lane 1','Lane 4','Both','Lane 1','Lane 4','Both','Lane 1','Lane 4','Both'})
xtickangle(45)

title(['Fusi significant spatial information for the different clusters'])
ylabel('SSI')
ylim([0 15])
xlim([-1 11])

%Now explore whether there is a relationship between 
%prediction importance and these spatial
%information clusteers




%plot in color prediciton importance categories
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);

for ii=1:length(imps.mean_imp_ROI)
    this_file=imps.mean_imp_fileNo(ii);
    this_ROI=imps.mean_imp_ROI(ii);

    ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
    this_ii_k=idx(ii_ROI);

    %1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s
    switch imps.ii_95_class(ii)
        case 1
            %1 is >=x95
            plot(mappedX(ii_info,1),mappedX(ii_info,2),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0],'MarkerSize',5);
        case 2
            %2 is >=y95
            plot(mappedX(ii_info,1),mappedX(ii_info,2),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9],'MarkerSize',5);
        case 3
            %3 is >=conc95
            plot(mappedX(ii_info,1),mappedX(ii_info,2),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5],'MarkerSize',5);
        case 4
            %4 is >=x95 and conc95
            plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor',[0.95 0.9 0.25],'MarkerSize',5);
        case 5
            %5 is >=y95 and conc95
            plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7],'MarkerSize',5);
        case 6
            %6 is >=x95 and y95
            plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0],'MarkerSize',5);
        case 7
            %7 is >= all 95s
            plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor',[0.8 0.6 0.7],'MarkerSize',5);
    end
end


xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE Information Content,prediction importance ' ]);

%ii_95_class 1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s


%Now plot the fraction of significant predictive importance
sig_pred_imp_x_ii_k=[];
ii_x=0;

sig_pred_imp_y_ii_k=[];
ii_y=0;

sig_pred_imp_conc_ii_k=[];
ii_conc=0;

for ii=1:length(imps.mean_imp_ROI)
    this_file=imps.mean_imp_fileNo(ii);
    this_ROI=imps.mean_imp_ROI(ii);

    ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
    this_ii_k=idx(ii_info);

    if imps.sig_pred_imp_x(ii)==1
        ii_x=ii_x+1;
        sig_pred_imp_x_ii_k(ii_x)=this_ii_k;
    end


    if imps.sig_pred_imp_y(ii)==1
        ii_y=ii_y+1;
        sig_pred_imp_y_ii_k(ii_y)=this_ii_k;
    end


    if imps.sig_pred_imp_conc(ii)==1
        ii_conc=ii_conc+1;
        sig_pred_imp_conc_ii_k(ii_conc)=this_ii_k;
    end

end

figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
bar_offset=0;

%predictive importance for x
for ii_k=1:3
    this_fraction=sum(sig_pred_imp_x_ii_k==ii_k)/length(sig_pred_imp_x_ii_k);

    switch ii_k
        case 1
            bar(bar_offset,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset+4,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset+8,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
    end

end



%predictive importance for y
for ii_k=1:3
    this_fraction=sum(sig_pred_imp_y_ii_k==ii_k)/length(sig_pred_imp_y_ii_k);
       switch ii_k
            case 1
                bar(bar_offset+1,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+5,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+9,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end
  
end



%predictive importance for conc 
for ii_k=1:3
    this_fraction=sum(sig_pred_imp_conc_ii_k==ii_k)/length(sig_pred_imp_conc_ii_k);
       switch ii_k
            case 1
                bar(bar_offset+2,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+6,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+10,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end

end


text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2 4 5 6 8 9 10])
xticklabels({'x','y','odor','x','y','odor','x','y','odor'})
xtickangle(45)

title(['Fraction of prediction important ROIs for the different clusters'])
ylabel('Fraction')
ylim([0 1])
xlim([-1 11])

%Spatial correlations of the spatial pattern
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
bar_offset=0;
edges=[0:0.05:1];
rand_offset=0.5;

%predictive importance for x
for ii_k=1:3
    these_rhos=imps.all_spatial_rhol1l4(idx==ii_k);
    these_rhos=these_rhos(~isnan(these_rhos));
    switch ii_k
        case 1
            bar(bar_offset,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
                ,edges,bar_offset,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+1,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
                ,edges,bar_offset+1,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+2,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
                ,edges,bar_offset+2,rand_offset,'k','k',1);
    end



end




% 
% text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2])
xticklabels({'Cluster 1','Cluster 2','Cluster 3'})
xtickangle(45)

title(['Spatial correlation for the different clusters'])
ylabel('Rho')
ylim([-1 1])
xlim([-1 3])


%Center of mass distances of the spatial pattern
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
bar_offset=0;
edges=[0:25:500];
rand_offset=0.5;

%predictive importance for x
for ii_k=1:3
    these_dcoms=imps.all_delta_center_of_mass(idx==ii_k);
    these_dcoms=these_dcoms(~isnan(these_dcoms));
    switch ii_k
        case 1
            bar(bar_offset,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
                ,edges,bar_offset,rand_offset,'k','k',1);
        case 2
            bar(bar_offset+1,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
                ,edges,bar_offset+1,rand_offset,'k','k',1);
        case 3
            bar(bar_offset+2,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            %Violin plot
            [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
                ,edges,bar_offset+2,rand_offset,'k','k',1);
    end



end

% 
% text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2])
xticklabels({'Cluster 1','Cluster 2','Cluster 3'})
xtickangle(45)

title(['Difference in the center of mass'])
ylabel('delta com')
ylim([0 350])
xlim([-1 3])

%Show a subset of the spatial maps
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Get XY and dFF per trial
        arena_file=handles_XY.arena_file{fileNo};
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        no_neurons=handles_out.no_neurons;


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
        for ii_ROI=1:no_neurons

            ii_ROI_all=find((all_info_fileNo==fileNo)&(all_info_ii_ROI==ii_ROI));

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

            % %Cumulative counts for info
            % cum_xy=zeros(1,10*10);
            % cum_lane=zeros(1,2);
            % cum_bindFF=zeros(1,2);
            % 
            % cum_xy_bindFF_Lane1=zeros(2,10*10);
            % cum_xy_bindFF_Lane4=zeros(2,10*10);
            % cum_xy_bindFF_BothLanes=zeros(2,10*10);
            % 
            % cum_xy_bindFF_lane=zeros(2,10*10,2);
            % 
            % cum_xy_Lane1=zeros(1,10*10);
            % cum_xy_Lane4=zeros(1,10*10);
            % cum_xy_BothLanes=zeros(1,10*10);
            % 
            % cum_bindFF_Lane1=zeros(1,2);
            % cum_bindFF_Lane4=zeros(1,2);
            % cum_bindFF_BothLanes=zeros(1,2);

           


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
            end
   
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
                close(figNo)
            catch
            end


            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .75 .75])



            % %Plot dFF timecourses
            % subplot(2,3,1:3)
            % hold on
            % ii_plot=0;
 
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

            delta_below_zero_ii=max(all_ii_turns(include_trial==1));
            delta_above_zero_ii=max(all_ii_ends(include_trial==1)-all_ii_turns(include_trial==1));
            time_bins=handles_XY.dt*([1:delta_below_zero_ii+delta_above_zero_ii]-delta_below_zero_ii);

            hit1_dFF=[];
            ii_hit1=0;
            miss1_dFF=[];
            ii_miss1=0;
            hit4_dFF=[];
            ii_hit4=0;
            miss4_dFF=[];
            ii_miss4=0;

            %Let's do the glm stats from -3 to 3 sec in 1 sec bins
            trimmed_time_range=[-3 3];
            trimmed_dt=1;
            trimmed_time_bins=[trimmed_time_range(1)+(trimmed_dt/2):trimmed_dt:trimmed_time_range(2)-(trimmed_dt/2)];

            glm_div_ii=0;
            glm_div=[];

            for trNo=1:trials.odor_trNo
                % if include_trial(trNo)==1
                    these_dFF=trials.trial(trNo).XdFFtest(:,ii_ROI);
                    these_time_bins=time_bins(delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)));
                    %Okabe_Ito colors
                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            ii_hit1=ii_hit1+1;
                            hit1_dFF(ii_hit1,1:length(time_bins))=min(these_dFF);
                            hit1_dFF(ii_hit1,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
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
                % end
            end

           


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

            y_trials_end=max([y_trials1(end) y_trials4(end)])

            %Plot dFF for lane 1
            subplot(2, 6, [1 2 3]);
            hold on
            
            drg_pcolor(repmat(time_bins,length(y_trials1),1),repmat(y_trials1,length(time_bins),1)',hit_miss_1)
            colormap(this_cmap)
            if max(max([max(hit_miss_4(:)) max(hit_miss_1(:))]))>0
                clim([0 max([max(hit_miss_4(:)) max(hit_miss_1(:))])]);
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
                clim([0 max([max(hit_miss_4(:)) max(hit_miss_1(:))])])
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

                yticks([50 100 200 300 400 410])
                yticklabels({'Lane 4','100','200','300','400','Lane 1'})

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                

                title_legend=['Lane 1, SI= ' num2str(all_info_lane1(ii_ROI_all))];
                if sig_all_info_lane1(ii_ROI_all)==1
                    title_legend=[title_legend ' S'];
                end

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

                yticks([50 100 200 300 400 410])

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

                yticks([50 100 200 300 400 410])

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
               

                title_legend=['Lane 4, SI= ' num2str(all_info_lane4(ii_ROI_all))];
                if sig_all_info_lane4(ii_ROI_all)==1
                    title_legend=[title_legend ' S'];
                end

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
                
                yticks([50 100 200 300 400 410])
                yticklabels({'Lane 4','100','200','300','400','Lane 1'})

                xlabel('x (mm)')
                ylabel('y (mm)')
                
                title_legend=['Lane 4, SI= ' num2str(all_info_lane4(ii_ROI_all))];
                if sig_all_info_lane4(ii_ROI_all)==1
                    title_legend=[title_legend ' S'];
                end
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

                yticks([50 100 200 300 400 410])
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

                yticks([50 100 200 300 400 410])
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                
                title_legend=['Lane 1, SI= ' num2str(all_info_lane1(ii_ROI_all))];
                if sig_all_info_lane1(ii_ROI_all)==1
                    title_legend=[title_legend ' S'];
                end
                title(title_legend)
            end


            sgt_legend=['dFF map file  No ' num2str(fileNo) ' ROI No ' num2str(ii_ROI) ' MI= ' num2str(all_info_mutual_info14(ii_ROI_all))];

            
          


            if sig_all_info_mutual_info14(ii_ROI_all)==1
                sgt_legend=[sgt_legend ' S']; %MI is significant
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_conc)>0
                sgt_legend=[sgt_legend ' odor '];
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_x)>0
                sgt_legend=[sgt_legend ' x '];
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_y)>0
                sgt_legend=[sgt_legend ' y '];
            end

            switch idx(ii_ROI_all)
                case 1
                    sgt_legend=[sgt_legend ' cl1 '];
                case 2
                    sgt_legend=[sgt_legend ' cl2 '];
                case 3
                    sgt_legend=[sgt_legend ' cl3 '];
            end

            sgtitle(sgt_legend)

            %Now do glm
            tbl = table(glm_div.data',glm_div.trial_type',glm_div.time',...
                'VariableNames',{'dFF','trial_type','time'});
            mdl = fitglm(tbl,'dFF~trial_type+time'...
                ,'CategoricalVars',[2])


            if sig_all_info_lane1(ii_ROI_all)==1
                pffft=1;
            end

            if sig_all_info_lane4(ii_ROI_all)==1
                pffft=1;
            end

            if sig_all_info_mutual_info14(ii_ROI_all)==1
                pfft=1; %MI is significant
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_conc)>0
                pfft=1;
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_x)>0
                pfft=1;
            end

            if sum(ii_ROI==imps.file(fileNo).ROIs_y)>0
                pffft=1;
            end

        end
    end
end

fclose(fileID);

pffft=1;

