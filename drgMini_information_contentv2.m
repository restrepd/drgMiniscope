%drgMini_information_contentv2
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

%Load the odor plumes
load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odor_plume_patternsDR.mat')

op_threshold=-3.2;
for cm_from_floor=1:2
    mean_plume_l1=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1';
    mean_plume_l4=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4';
    x_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume;
    y_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume;

    %Generate binary plumes
    binary_plumel1=zeros(10,10);
    binary_plumel4=zeros(10,10);

    for ii_x=1:10
        for ii_y=1:10
            these_opl1=[];
            these_opl1=mean_plume_l1( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
            this_mean_opl1=mean(these_opl1(:));
            if this_mean_opl1<op_threshold
                binary_plumel1(ii_x,ii_y)=0;
            else
                binary_plumel1(ii_x,ii_y)=1;
            end

             these_opl4=[];
            these_opl4=mean_plume_l4( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
            this_mean_opl4=mean(these_opl4(:));
            if this_mean_opl4<op_threshold
                binary_plumel4(ii_x,ii_y)=0;
            else
                binary_plumel4(ii_x,ii_y)=1;
            end
        end
    end

    odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1=binary_plumel1;
    odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4=binary_plumel4;

end

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
            
            include_xy=zeros(1,10*10);
            include_op=zeros(1,2*10*10);


            cum_xy_Lane1=zeros(1,10*10);
            cum_xy_Lane4=zeros(1,10*10);
            cum_xy_BothLanes=zeros(1,10*10);
            
            cum_lane_bindFF=zeros(2,2);
            

            cum_bindFF_Lane1=zeros(1,2);
            cum_bindFF_Lane4=zeros(1,2);
            cum_bindFF_BothLanes=zeros(1,2);

            % cum_op=zeros(1,2*10*10); %Added op
            % cum_xy_bindFF_op=zeros(2,10*10,2*10*10); %Added op
            % cum_op_bindFF=zeros(2,2*10*10); %Added op

            cum_op_bin=zeros(1,2); %Added op
            cum_xy_bindFF_op_bin=zeros(2,10*10,2); %Added op
            cum_op_bindFF_bin=zeros(2,2); %Added op
            cum_xy_op_bin=zeros(10*10,2);

            % cum_op=zeros(1,2*10*10); %Added op
            % cum_xy_bindFF_op=zeros(2,10*10,10*10); %Added op
            % cum_op_bindFF=zeros(2,10*10); %Added op

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
                    include_xy(this_xy_ii)=1;
                    include_op(this_xy_ii)=1;
                    include_op(this_xy_ii+100)=1;

                    cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
                    cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_BothLanes(this_bin_dFF,this_xy_ii)+1;
                    cum_xy_BothLanes(this_xy_ii)=cum_xy_BothLanes(this_xy_ii)+1;

                    if trials.lane_per_trial(trNo)==1
                        % this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl1_activity=sum_dFFl1_activity+these_dFF(ii_t);
                        cum_lane_bindFF(this_bin_dFF,1)=cum_lane_bindFF(this_bin_dFF,1)+1;
                        cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane1(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,1)+1;
                        cum_xy_Lane1(this_xy_ii)=cum_xy_Lane1(this_xy_ii)+1;
                        cum_bindFF_Lane1(this_bin_dFF)=cum_bindFF_Lane1(this_bin_dFF)+1;
                        cum_lane(1)=cum_lane(1)+1;

                        switch handles_conc.group(fileNo)
                            case 1
                                %2 cm from the floor
                                cm_from_floor=2;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)+1;

                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1;
                                cum_xy_op_bin(this_xy_ii,binary_op+1)=cum_xy_op_bin(this_xy_ii,binary_op+1)+1;
                                cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF,binary_op+1)+1;
                            case 5
                                %1 cm from the floor
                                cm_from_floor=1;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)+1;

                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1; 
                                cum_xy_op_bin(this_xy_ii,binary_op+1)=cum_xy_op_bin(this_xy_ii,binary_op+1)+1;
                                cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF,binary_op+1)+1;
                        end

                    else
                        % this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl4_activity=sum_dFFl4_activity+these_dFF(ii_t);
                        cum_lane_bindFF(this_bin_dFF,2)=cum_lane_bindFF(this_bin_dFF,2)+1;
                        cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_Lane4(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)=cum_xy_bindFF_lane(this_bin_dFF,this_xy_ii,2)+1;
                        cum_xy_Lane4(this_xy_ii)=cum_xy_Lane4(this_xy_ii)+1;
                        cum_bindFF_Lane4(this_bin_dFF)=cum_bindFF_Lane4(this_bin_dFF)+1;
                        cum_lane(2)=cum_lane(2)+1;

                        switch handles_conc.group(fileNo)
                            case 1
                                %2 cm from the floor
                                cm_from_floor=2;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);
                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)+1;

                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1; 
                                cum_xy_op_bin(this_xy_ii,binary_op+1)=cum_xy_op_bin(this_xy_ii,binary_op+1)+1;
                                cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF,binary_op+1)+1;

                            case 5
                                %1 cm from the floor
                                cm_from_floor=1;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);
                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF,binary_op*100+this_xy_ii)+1;

                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1; 
                                cum_xy_op_bin(this_xy_ii,binary_op+1)=cum_xy_op_bin(this_xy_ii,binary_op+1)+1;
                                cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF,this_xy_ii,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF,binary_op+1)+1;
                        end
                    end

                end
            end

            %Calculate information

            p_xy=cum_xy(logical(include_xy))/sum(cum_xy(logical(include_xy)));
            p_lane=cum_lane/sum(cum_lane(:));
            p_bindFF=cum_bindFF/sum(cum_bindFF(:));
            

            p_lane_bindFF=cum_lane_bindFF/sum(cum_lane_bindFF(:));

            included_cum_xy_bindFF_Lane1=cum_xy_bindFF_Lane1(:,logical(include_xy));
            p_xy_bindFF_Lane1=included_cum_xy_bindFF_Lane1/sum(included_cum_xy_bindFF_Lane1(:));

            included_cum_xy_bindFF_Lane4=cum_xy_bindFF_Lane4(:,logical(include_xy));
            p_xy_bindFF_Lane4=included_cum_xy_bindFF_Lane4/sum(included_cum_xy_bindFF_Lane4(:));

            included_cum_xy_bindFF_BothLanes=cum_xy_bindFF_BothLanes(:,logical(include_xy));
            p_xy_bindFF_BothLanes=included_cum_xy_bindFF_BothLanes/sum(included_cum_xy_bindFF_BothLanes(:));

            included_cum_xy_bindFF_lane=cum_xy_bindFF_lane(:,logical(include_xy),:);
            p_xy_bindFF_lane=included_cum_xy_bindFF_lane/sum(included_cum_xy_bindFF_lane(:));
            

            included_cum_xy_Lane1=cum_xy_Lane1(1,logical(include_xy));
            p_xy_Lane1=included_cum_xy_Lane1/sum(included_cum_xy_Lane1(:));

            included_cum_xy_Lane4=cum_xy_Lane4(1,logical(include_xy));
            p_xy_Lane4=included_cum_xy_Lane4/sum(included_cum_xy_Lane4(:));

            included_cum_xy_BothLanes=cum_xy_BothLanes(1,logical(include_xy));
            p_xy_BothLanes=included_cum_xy_BothLanes/sum(included_cum_xy_BothLanes(:));

            p_bindFF_Lane1=cum_bindFF_Lane1/sum(cum_bindFF_Lane1(:));

            p_bindFF_Lane4=cum_bindFF_Lane4/sum(cum_bindFF_Lane4(:));

            p_bindFF_BothLanes=cum_bindFF_BothLanes/sum(cum_bindFF_BothLanes(:));

            % p_op=cum_op(logical(include_op))/sum(cum_op(logical(include_op)));
            % 
            % included_cum_xy_bindFF_op=cum_xy_bindFF_op(:,logical(include_xy),logical(include_op));
            % p_xy_bindFF_op=included_cum_xy_bindFF_op/sum(included_cum_xy_bindFF_op(:));

            included_cum_xy_op_bin=cum_xy_op_bin(logical(include_xy),:);
            p_xy_op_bin=included_cum_xy_op_bin/sum(included_cum_xy_op_bin(:));

            % included_cum_op_bindFF=cum_op_bindFF(:,logical(include_op));
            % p_op_bindFF=included_cum_op_bindFF/sum(included_cum_op_bindFF(:));

            p_op_bin=cum_op_bin/sum(cum_op_bin);

            included_cum_xy_bindFF_op_bin=cum_xy_bindFF_op_bin(:,logical(include_xy),:);
            p_xy_bindFF_op_bin=included_cum_xy_bindFF_op_bin/sum(included_cum_xy_bindFF_op_bin(:));

            p_op_bindFF_bin=cum_op_bindFF_bin/sum(cum_op_bindFF_bin(:));

            all_info_ii=all_info_ii+1;
            all_info_fileNo(all_info_ii)=fileNo;
            all_info_ii_ROI(all_info_ii)=ii_ROI;

            all_info_lane1(all_info_ii)=0;
            all_info_lane4(all_info_ii)=0;
            all_info_mutual_info14(all_info_ii)=0;

            all_info_lane(all_info_ii)=0;

            % all_info_op(all_info_ii)=0;
            % all_info_mutual_info_op(all_info_ii)=0;

            all_info_op_bin(all_info_ii)=0;
            all_info_mutual_info_dFFbin_xy_op_bin(all_info_ii)=0;

            all_info_mutual_xy_op_bin(all_info_ii)=0;

            %Calculate info for lane 1
            for ii_xy=1:size(p_xy_bindFF_Lane1,2)
                % if sum(p_xy_bindFF_Lane1(:,ii_xy))>0
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)~=0
                            all_info_lane1(all_info_ii)=all_info_lane1(all_info_ii)+p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)*...
                                log2(p_xy_bindFF_Lane1(ii_bin_dFF,ii_xy)/(p_xy_Lane1(ii_xy)*p_bindFF_Lane1(ii_bin_dFF)));
                        end
                    end
                % end
            end

            %Calculate info for lane 4
            for ii_xy=1:size(p_xy_bindFF_Lane4,2)
                % if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)~=0
                            all_info_lane4(all_info_ii)=all_info_lane4(all_info_ii)+p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)*...
                                log2(p_xy_bindFF_Lane4(ii_bin_dFF,ii_xy)/(p_xy_Lane4(ii_xy)*p_bindFF_Lane4(ii_bin_dFF)));
                        end
                    end
                % end
            end

            %Calculate mutual info between lanes 1 and 4
            for ii_xy=1:size(p_xy_bindFF_lane,2)
                % if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_lane=1:2
                        for ii_bin_dFF=1:2
                            if p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)~=0
                                all_info_mutual_info14(all_info_ii)=all_info_mutual_info14(all_info_ii)+p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)*...
                                    log2(p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_lane(ii_lane)));
                            end
                        end
                    end
                % end
            end

            %Calculate info for lanes
            for ii_lane=1:2
                for ii_bin_dFF=1:2
                    if p_lane_bindFF(ii_bin_dFF,ii_lane)~=0
                        all_info_lane(all_info_ii)=all_info_lane(all_info_ii)+p_lane_bindFF(ii_bin_dFF,ii_lane)*...
                            log2(p_lane_bindFF(ii_bin_dFF,ii_lane)/(p_lane(ii_lane)*p_bindFF(ii_bin_dFF)));
                    end
                end
            end

            %  %Calculate info for op
            % for ii_op=1:size(p_op_bindFF,2)
            %     for ii_bin_dFF=1:2
            %         if p_op_bindFF(ii_bin_dFF,ii_op)~=0
            %             all_info_op(all_info_ii)=all_info_op(all_info_ii)+p_op_bindFF(ii_bin_dFF,ii_op)*...
            %                 log2(p_op_bindFF(ii_bin_dFF,ii_op)/(p_op(ii_op)*p_bindFF(ii_bin_dFF)));
            %         end
            %     end
            % end
            % 
            % %Calculate mutual info with op
            % for ii_xy=1:size(p_xy_bindFF_op,2)
            %     % if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
            %         for ii_op=1:size(p_xy_bindFF_op,3)
            %             for ii_bin_dFF=1:2
            %                 if p_xy_bindFF_op(ii_bin_dFF,ii_xy,ii_op)~=0
            %                     delta_info=p_xy_bindFF_op(ii_bin_dFF,ii_xy,ii_op)*...
            %                         log2(p_xy_bindFF_op(ii_bin_dFF,ii_xy,ii_op)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_op(ii_op)));
            %                     all_info_mutual_info_op(all_info_ii)=all_info_mutual_info_op(all_info_ii)+delta_info;
            %                 end
            %             end
            %         end
            %     % end
            % end

              %Calculate info for op_bin
            for ii_op=1:2
                for ii_bin_dFF=1:2
                    if p_op_bindFF_bin(ii_bin_dFF,ii_op)~=0
                        all_info_op_bin(all_info_ii)=all_info_op_bin(all_info_ii)+p_op_bindFF_bin(ii_bin_dFF,ii_op)*...
                            log2(p_op_bindFF_bin(ii_bin_dFF,ii_op)/(p_op_bin(ii_op)*p_bindFF(ii_bin_dFF)));
                    end
                end
            end

            %Calculate mutual info with op_bin
            for ii_xy=1:size(p_xy_bindFF_op_bin,2)
                % if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_op=1:2
                        for ii_bin_dFF=1:2
                            if p_xy_bindFF_op_bin(ii_bin_dFF,ii_xy,ii_op)~=0
                                delta_info=p_xy_bindFF_op_bin(ii_bin_dFF,ii_xy,ii_op)*...
                                    log2(p_xy_bindFF_op_bin(ii_bin_dFF,ii_xy,ii_op)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_op_bin(ii_op)));
                                all_info_mutual_info_dFFbin_xy_op_bin(all_info_ii)=all_info_mutual_info_dFFbin_xy_op_bin(all_info_ii)+delta_info;
                            end
                        end
                    end
                % end
            end

              %Calculate mutual info with xy x op_bin
            for ii_xy=1:size(p_xy_op_bin,1)
                % if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_op=1:2
                        for ii_bin_dFF=1:2
                            if p_xy_bindFF_op_bin(ii_bin_dFF,ii_xy,ii_op)~=0
                                delta_info=p_xy_op_bin(ii_xy,ii_op)*...
                                    log2(p_xy_op_bin(ii_xy,ii_op)/(p_xy(ii_xy)*p_op_bin(ii_op)));
                                all_info_mutual_xy_op_bin(all_info_ii)=all_info_mutual_xy_op_bin(all_info_ii)+delta_info;
                            end
                        end
                    end
                % end
            end


            %Now calculate information content with the shuffled dFF
            rng(ii_ROI)
            these_rnd=rand(1,2*n_shuffle_SI);

            for ii_sh=1:n_shuffle_SI

                

                %Flip and roll Stefanini et al 2020 https://doi.org/10.1016/j.neuron.2020.05.022
                

                % We will do a reversal and a circular permutation of dFF
                these_x_reversed=zeros(1,length(all_these_x));
                these_y_reversed=zeros(1,length(all_these_y));

                % With this one you flip and roll x and y
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

                %With this one you flip and roll dFF
                offset_ii=floor(these_rnd(ii_sh+n_shuffle_SI)*length(all_these_x));
                all_these_dFF_reversed=zeros(1,length(all_these_dFF));
                for ii_trl=1:length(all_these_dFF)
                    this_ii_trl=ii_trl+offset_ii;
                    if this_ii_trl>length(all_these_dFF)
                        offset_ii=-ii_trl+1;
                        this_ii_trl=ii_trl+offset_ii;
                    end
                    all_these_dFF_reversed(1,length(all_these_dFF)-ii_trl+1)=all_these_dFF(this_ii_trl);
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

                cum_lane_bindFF=zeros(2,2);

                cum_bindFF_Lane1=zeros(1,2);
                cum_bindFF_Lane4=zeros(1,2);
                cum_bindFF_BothLanes=zeros(1,2);

                % cum_op=zeros(1,2*10*10); %Added op
                % cum_xy_bindFF_op=zeros(2,10*10,2*10*10); %Added op
                % cum_op_bindFF=zeros(2,2*10*10); %Added op


                cum_op_bin=zeros(1,2); %Added op
                cum_xy_bindFF_op_bin=zeros(2,10*10,2); %Added op
                cum_op_bindFF_bin=zeros(2,2); %Added op

                cum_xy_op_bin=zeros(10*10,2); %Added op

                include_xy=zeros(1,10*10);
                include_op=zeros(1,2*10*10);


                for ii_t=1:length(all_these_dFF_reversed)
                    this_x_ii=ceil(all_these_x(ii_t)/50);
                    if this_x_ii==11
                        this_x_ii=10;
                    end

                    this_y_ii=ceil(all_these_y(ii_t)/48);
                    if this_y_ii==11
                        this_y_ii=10;
                    end

                     this_x_ii_rev=ceil(these_x_reversed(ii_t)/50);
                    if this_x_ii_rev==11
                        this_x_ii_rev=10;
                    end

                    this_y_ii_rev=ceil(these_y_reversed(ii_t)/48);
                    if this_y_ii_rev==11
                        this_y_ii_rev=10;
                    end

                    % this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                    % this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
                    % sum_dFF_activity=sum_dFF_activity+these_dFF(ii_t);

                    %Tally info
                    this_bin_dFF_rev=(all_these_dFF_reversed(ii_t)>0)+1;
                    cum_bindFF(this_bin_dFF_rev)=cum_bindFF(this_bin_dFF_rev)+1;
                    cum_bindFF_BothLanes(this_bin_dFF_rev)=cum_bindFF_BothLanes(this_bin_dFF_rev)+1;
                    this_xy_ii=this_x_ii+10*(this_y_ii-1);
                    this_xy_ii_rev=this_x_ii_rev+10*(this_y_ii_rev-1);
                    include_xy(this_xy_ii)=1;
                    include_op(this_xy_ii)=1;
                    include_op(this_xy_ii+100)=1;


                    cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
                    cum_xy_bindFF_BothLanes(this_bin_dFF_rev,this_xy_ii)=cum_xy_bindFF_BothLanes(this_bin_dFF_rev,this_xy_ii)+1;
                    cum_xy_BothLanes(this_xy_ii)=cum_xy_BothLanes(this_xy_ii)+1;

                    if all_lanes(ii_t)==1
                        % this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl1_activity=sum_dFFl1_activity+these_dFF(ii_t);
                        cum_lane_bindFF(this_bin_dFF_rev,1)=cum_lane_bindFF(this_bin_dFF_rev,1)+1;
                        cum_xy_bindFF_Lane1(this_bin_dFF_rev,this_xy_ii)=cum_xy_bindFF_Lane1(this_bin_dFF_rev,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF_rev,this_xy_ii,1)=cum_xy_bindFF_lane(this_bin_dFF_rev,this_xy_ii,1)+1;
                        cum_xy_Lane1(this_xy_ii)=cum_xy_Lane1(this_xy_ii)+1;
                        cum_bindFF_Lane1(this_bin_dFF_rev)=cum_bindFF_Lane1(this_bin_dFF_rev)+1;
                        cum_lane(1)=cum_lane(1)+1;

                        switch handles_conc.group(fileNo)
                            case 1
                                %2 cm from the floor
                                cm_from_floor=2;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % 
                                % cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)+1;
                                % 
                                
                                cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)+1;
                                cum_xy_op_bin(this_xy_ii_rev,binary_op+1)=cum_xy_op_bin(this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1;
                            case 5
                                %1 cm from the floor
                                cm_from_floor=1;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                                % 
                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)+1;
                                
                                cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)+1;
                                cum_xy_op_bin(this_xy_ii_rev,binary_op+1)=cum_xy_op_bin(this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1;
                        end
                    else
                        % this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        % this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
                        % sum_dFFl4_activity=sum_dFFl4_activity+these_dFF(ii_t);
                        cum_lane_bindFF(this_bin_dFF_rev,2)=cum_lane_bindFF(this_bin_dFF_rev,2)+1;
                        cum_xy_bindFF_Lane4(this_bin_dFF_rev,this_xy_ii)=cum_xy_bindFF_Lane4(this_bin_dFF_rev,this_xy_ii)+1;
                        cum_xy_bindFF_lane(this_bin_dFF_rev,this_xy_ii,2)=cum_xy_bindFF_lane(this_bin_dFF_rev,this_xy_ii,2)+1;
                        cum_xy_Lane4(this_xy_ii)=cum_xy_Lane4(this_xy_ii)+1;
                        cum_bindFF_Lane4(this_bin_dFF_rev)=cum_bindFF_Lane4(this_bin_dFF_rev)+1;
                        cum_lane(2)=cum_lane(2)+1;

                        switch handles_conc.group(fileNo)
                            case 1
                                %2 cm from the floor
                                cm_from_floor=2;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);

                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)+1;
                                %
                                cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)+1;
                                cum_xy_op_bin(this_xy_ii_rev,binary_op+1)=cum_xy_op_bin(this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1;
                            case 5
                                %1 cm from the floor
                                cm_from_floor=1;
                                %Find the binary op
                                binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);

                                % cum_op(binary_op*100+this_xy_ii)=cum_op(binary_op*100+this_xy_ii)+1;
                                % cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)=...
                                %     cum_xy_bindFF_op(this_bin_dFF_rev,this_xy_ii,binary_op*100+this_xy_ii)+1;
                                % cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)=cum_op_bindFF(this_bin_dFF_rev,binary_op*100+this_xy_ii)+1;
                                %

                                cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)=...
                                    cum_xy_bindFF_op_bin(this_bin_dFF_rev,this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)=cum_op_bindFF_bin(this_bin_dFF_rev,binary_op+1)+1;
                                cum_xy_op_bin(this_xy_ii_rev,binary_op+1)=cum_xy_op_bin(this_xy_ii_rev,binary_op+1)+1;
                                cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1;
                        end
                    end

                end

                %Calculate information

                p_xy=cum_xy(logical(include_xy))/sum(cum_xy(logical(include_xy)));
                p_lane=cum_lane/sum(cum_lane(:));
                p_bindFF=cum_bindFF/sum(cum_bindFF(:));


                p_lane_bindFF=cum_lane_bindFF/sum(cum_lane_bindFF(:));

                included_cum_xy_bindFF_Lane1=cum_xy_bindFF_Lane1(:,logical(include_xy));
                p_xy_bindFF_Lane1=included_cum_xy_bindFF_Lane1/sum(included_cum_xy_bindFF_Lane1(:));

                included_cum_xy_bindFF_Lane4=cum_xy_bindFF_Lane4(:,logical(include_xy));
                p_xy_bindFF_Lane4=included_cum_xy_bindFF_Lane4/sum(included_cum_xy_bindFF_Lane4(:));

                included_cum_xy_bindFF_BothLanes=cum_xy_bindFF_BothLanes(:,logical(include_xy));
                p_xy_bindFF_BothLanes=included_cum_xy_bindFF_BothLanes/sum(included_cum_xy_bindFF_BothLanes(:));

                included_cum_xy_bindFF_lane=cum_xy_bindFF_lane(:,logical(include_xy),:);
                p_xy_bindFF_lane=included_cum_xy_bindFF_lane/sum(included_cum_xy_bindFF_lane(:));


                included_cum_xy_Lane1=cum_xy_Lane1(1,logical(include_xy));
                p_xy_Lane1=included_cum_xy_Lane1/sum(included_cum_xy_Lane1(:));

                included_cum_xy_Lane4=cum_xy_Lane4(1,logical(include_xy));
                p_xy_Lane4=included_cum_xy_Lane4/sum(included_cum_xy_Lane4(:));

                included_cum_xy_BothLanes=cum_xy_BothLanes(1,logical(include_xy));
                p_xy_BothLanes=included_cum_xy_BothLanes/sum(included_cum_xy_BothLanes(:));

                p_bindFF_Lane1=cum_bindFF_Lane1/sum(cum_bindFF_Lane1(:));

                p_bindFF_Lane4=cum_bindFF_Lane4/sum(cum_bindFF_Lane4(:));

                p_bindFF_BothLanes=cum_bindFF_BothLanes/sum(cum_bindFF_BothLanes(:));

                % p_op=cum_op(logical(include_op))/sum(cum_op(logical(include_op)));

                p_op_bin=cum_op_bin/sum(cum_op_bin);

                % included_cum_xy_bindFF_op=cum_xy_bindFF_op(:,logical(include_xy),logical(include_op));
                % p_xy_bindFF_op=included_cum_xy_bindFF_op/sum(included_cum_xy_bindFF_op(:));
                %
                % included_cum_op_bindFF=cum_op_bindFF(:,logical(include_op));
                % p_op_bindFF=included_cum_op_bindFF/sum(included_cum_op_bindFF(:));

                included_cum_xy_op_bin=cum_xy_op_bin(logical(include_xy),:);
                p_xy_op_bin=included_cum_xy_op_bin/sum(included_cum_xy_op_bin(:));

                included_cum_xy_bindFF_op_bin=cum_xy_bindFF_op_bin(:,logical(include_xy),:);
                p_xy_bindFF_op_bin=included_cum_xy_bindFF_op_bin/sum(included_cum_xy_bindFF_op_bin(:));

                all_info_lane1_sh(all_info_ii,ii_sh)=0;
                all_info_lane4_sh(all_info_ii,ii_sh)=0;
                all_info_mutual_info14_sh(all_info_ii,ii_sh)=0;

                all_info_lane_sh(all_info_ii,ii_sh)=0;

                all_info_op_bin_sh(all_info_ii,ii_sh)=0;
                % all_info_mutual_info_dFFbin_xy_op_bin_sh(all_info_ii,ii_sh)=0;

                all_info_mutual_xy_op_bin_sh(all_info_ii,ii_sh)=0;

                all_info_mutual_info_dFFbin_xy_op_bin_sh(all_info_ii,ii_sh)=0;

                %Calculate info for lane 1
                for ii_xy=1:size(p_xy_bindFF_Lane1,2)
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
                for ii_xy=1:size(p_xy_bindFF_Lane4,2)
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
                for ii_xy=1:size(p_xy_bindFF_lane,2)
                    % if sum(p_xy_bindFF_lane(:,ii_xy))>0
                        for ii_lane=1:2
                            for ii_bin_dFF=1:2
                                if p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)~=0
                                    all_info_mutual_info14_sh(all_info_ii,ii_sh)=all_info_mutual_info14_sh(all_info_ii,ii_sh)+p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)*...
                                        log2(p_xy_bindFF_lane(ii_bin_dFF,ii_xy,ii_lane)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_lane(ii_lane)));
                                end
                            end
                        end
                    % end
                end

                %Calculate info for lanes
                for ii_lane=1:2
                    for ii_bin_dFF=1:2
                        if p_lane_bindFF(ii_bin_dFF,ii_lane)~=0
                            all_info_lane_sh(all_info_ii,ii_sh)=all_info_lane_sh(all_info_ii,ii_sh)+p_lane_bindFF(ii_bin_dFF,ii_lane)*...
                                log2(p_lane_bindFF(ii_bin_dFF,ii_lane)/(p_lane(ii_lane)*p_bindFF(ii_bin_dFF)));
                        end
                    end
                end


                %Calculate info for op
                for ii_op=1:2
                    for ii_bin_dFF=1:2
                        if p_op_bindFF_bin(ii_bin_dFF,ii_op)~=0
                            all_info_op_bin_sh(all_info_ii,ii_sh)=all_info_op_bin_sh(all_info_ii,ii_sh)+p_op_bindFF_bin(ii_bin_dFF,ii_op)*...
                                log2(p_op_bindFF_bin(ii_bin_dFF,ii_op)/(p_op_bin(ii_op)*p_bindFF(ii_bin_dFF)));
                        end
                    end
                end

                %Calculate mutual info with xy, dFFbin, op bin
                for ii_xy=1:size(p_xy_bindFF_op_bin,2)
                    % if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_op=1:2
                        for ii_bin_dFF=1:2
                            if p_xy_bindFF_op_bin(ii_bin_dFF,ii_xy,ii_op)~=0
                                all_info_mutual_info_dFFbin_xy_op_bin_sh(all_info_ii,ii_sh)=all_info_mutual_info_dFFbin_xy_op_bin_sh(all_info_ii)+p_xy_bindFF_op_bin(ii_bin_dFF,ii_xy,ii_op)*...
                                    log2(p_xy_bindFF_op_bin(ii_bin_dFF,ii_xy,ii_op)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)*p_op_bin(ii_op)));
                            end
                        end
                    end
                    % end
                end

                     %Calculate mutual info with op
                for ii_xy=1:size(p_xy_op_bin,1)
                    % if sum(p_xy_bindFF_Lane4(:,ii_xy))>0
                    for ii_op=1:2
                        % for ii_bin_dFF=1:2
                            if p_xy_op_bin(ii_xy,ii_op)~=0
                                all_info_mutual_xy_op_bin_sh(all_info_ii,ii_sh)=all_info_mutual_xy_op_bin_sh(all_info_ii)+p_xy_op_bin(ii_xy,ii_op)*...
                                    log2(p_xy_op_bin(ii_xy,ii_op)/(p_xy(ii_xy)*p_op_bin(ii_op)));
                            end
                        % end
                    end
                    % end
                end

            end

            

            % %Find turn points and time vectors
            % all_ii_turns=zeros(1,trials.odor_trNo);
            % all_ii_ends=zeros(1,trials.odor_trNo);
            % include_trial=zeros(1,trials.odor_trNo);
            % for trNo=1:trials.odor_trNo
            %     this_ii_last_turn=find(angles.trial(trNo).delta_x>=100,1,'last');
            %     if ~isempty(this_ii_last_turn)
            %         all_ii_turns(trNo)=angles.trial(trNo).ii_turns(this_ii_last_turn);
            %         all_ii_ends(trNo)= size(trials.trial(trNo).XYtest,1);
            %         include_trial(trNo)=1;
            %     else
            %         pffft=1;
            %     end
            % end

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

            % delta_below_zero_ii=max(all_ii_turns(include_trial==1));
            % delta_above_zero_ii=max(all_ii_ends(include_trial==1)-all_ii_turns(include_trial==1));

            delta_below_zero_ii=max(all_ii_turns);
            delta_above_zero_ii=max(all_ii_ends-all_ii_turns);

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

            y_trials_end=max([y_trials1(end) y_trials4(end)]);


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

[f_aic,x_aic] = drg_ecdf(all_info_lane);
plot(x_aic,f_aic,'Color',[240/255 228/255 66/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_info_op_bin);
plot(x_aic,f_aic,'Color',[0/255 114/255 178/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_info_mutual_info_dFFbin_xy_op_bin);
plot(x_aic,f_aic,'Color',[213/255 94/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_info_mutual_xy_op_bin);
plot(x_aic,f_aic,'Color',[204/255 121/255 167/255],'LineWidth',3)




these_ylim=ylim;
these_xlim=xlim;
text(0.3,0.4,'MI xy x bindFF Lane 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(0.3,0.35,'MI xy x bindFF Lane 4','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(0.3,0.30,'MI xy x bindFF x lane','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
text(0.3,0.25,'MI lane x bindFF','Color',[240/255 228/255 66/255],'FontWeight','bold','FontSize',16)
text(0.3,0.20,'MI op x bindFF','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
text(0.3,0.15,'MI op x bindFF x xy','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
text(0.3,0.10,'MI op x xy','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)

title('Cumulative histograms for mutual information')
xlabel('Information bits')
ylabel('Cumulative fraction')


%Now calculate 3xSD for the shuffled distributions
sd_all_info_lane1_sh=zeros(all_info_ii,1);
sd_all_info_lane4_sh=zeros(all_info_ii,1);
sd_all_info_mutual_info14_sh=zeros(all_info_ii,1);
sd_all_info_lane_sh=zeros(all_info_ii,1);
sd_all_info_op_bin=zeros(all_info_ii,1);
sd_all_info_mutual_info_dFFbin_xy_op_bin=zeros(all_info_ii,1);

mean_all_info_lane1_sh=zeros(all_info_ii,1);
mean_all_info_lane4_sh=zeros(all_info_ii,1);
mean_all_info_mutual_info14_sh=zeros(all_info_ii,1);
mean_all_info_lane_sh=zeros(all_info_ii,1);
mean_all_info_op_bin=zeros(all_info_ii,1);
mean_all_info_mutual_info_dFFbin_xy_op_bin=zeros(all_info_ii,1);

all_all_info_lane1_sh=[];
all_all_info_lane4_sh=[];
all_all_info_mutual_info14_sh=[];
all_all_info_lane_sh=[];
all_all_info_op_bin_sh=[];
all_all_info_mutual_info_dFFbin_xy_op_bin_sh=[];


for ii_ROI=1:all_info_ii

    this_ii_ROI=ii_ROI

    this_all_info_lane1_sh=zeros(1,n_shuffle_SI);
    this_all_info_lane1_sh(1,:)=all_info_lane1_sh(ii_ROI,:);
    this_all_info_lane4_sh=zeros(1,n_shuffle_SI);
    this_all_info_lane4_sh(1,:)=all_info_lane4_sh(ii_ROI,:);
    this_all_info_mutual_info14_sh=zeros(1,n_shuffle_SI);
    this_all_info_mutual_info14_sh(1,:)=all_info_mutual_info14_sh(ii_ROI,:);
    this_all_info_lane_sh=zeros(1,n_shuffle_SI);
    this_all_info_lane_sh(1,:)=all_info_lane_sh(ii_ROI,:);
    this_all_info_op_bin_sh=zeros(1,n_shuffle_SI);
    this_all_info_op_bin_sh(1,:)=all_info_op_bin_sh(ii_ROI,:);
    this_all_info_mutual_info_dFFbin_xy_op_bin_sh=zeros(1,n_shuffle_SI);
    this_all_info_mutual_info_dFFbin_xy_op_bin_sh(1,:)=all_info_mutual_info_dFFbin_xy_op_bin_sh(ii_ROI,:);

    all_all_info_lane1_sh=[all_all_info_lane1_sh this_all_info_lane1_sh];
    all_all_info_lane4_sh=[all_all_info_lane4_sh this_all_info_lane4_sh];
    all_all_info_mutual_info14_sh=[all_all_info_mutual_info14_sh this_all_info_mutual_info14_sh];
    all_all_info_lane_sh=[all_all_info_lane_sh this_all_info_lane_sh];

    all_all_info_op_bin_sh=[all_all_info_op_bin_sh this_all_info_op_bin_sh];
    all_all_info_mutual_info_dFFbin_xy_op_bin_sh=[all_all_info_mutual_info_dFFbin_xy_op_bin_sh this_all_info_mutual_info_dFFbin_xy_op_bin_sh];

    for jj=1:6
        switch jj
            case 1
                data=this_all_info_lane1_sh;
            case 2
                data=this_all_info_lane4_sh;
            case 3
                data=this_all_info_mutual_info14_sh;
            case 4
                data=this_all_info_lane_sh;
            case 5
                data=this_all_info_op_bin_sh;
            case 6
                data=this_all_info_mutual_info_dFFbin_xy_op_bin_sh;
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
            case 4
                sd_all_info_lane_sh(ii_ROI)=estimated_std;
            case 5
                sd_all_info_op_bin_sh(ii_ROI)=estimated_std;
            case 6
                sd_all_info_mutual_info_dFFbin_xy_op_bin_sh(ii_ROI)=estimated_std;
        end
    end

    mean_all_info_lane1_sh(ii_ROI)=mean(this_all_info_lane1_sh);
    mean_all_info_lane4_sh(ii_ROI)=mean(this_all_info_lane4_sh);
    mean_all_info_mutual_info14_sh(ii_ROI)=mean(this_all_info_mutual_info14_sh);
    mean_all_info_lane_sh(ii_ROI)=mean(this_all_info_lane_sh);
    mean_all_info_op_bin_sh(ii_ROI)=mean(this_all_info_lane_sh);
    mean_all_info_mutual_info_dFFbin_xy_op_bin_sh(ii_ROI)=mean(this_all_info_lane_sh);

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

title('MI xy x bindFF Lane 1, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')


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

title('MI xy x bindFF Lane 4, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')


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


title('MI xy x bindFF x lane, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')



%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on
 

[f_aic,x_aic] = drg_ecdf(all_info_lane);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_all_info_lane_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)

text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)


title('MI lane x bindFF, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')


%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on
 

[f_aic,x_aic] = drg_ecdf(all_info_op_bin);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_all_info_op_bin_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)

text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)


title('MI op x bindFF, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')


%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on
 

[f_aic,x_aic] = drg_ecdf(all_info_mutual_info_dFFbin_xy_op_bin);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(all_all_info_mutual_info_dFFbin_xy_op_bin_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)

text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)


title('MI op x bindFF x xy, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')

%Load information on prediction importance to do the Fusi analysis
load([save_PathPredImp save_FilePredImp])

%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

[f_aic,x_aic] = drg_ecdf(imps.all_Fusi_SSIl1);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(imps.all_Fusi_SSIl1_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)

text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)

title('Fusi MI xy x bindFF Lane 1, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')


%Plot histograms
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


[f_aic,x_aic] = drg_ecdf(imps.all_Fusi_SSIl4);
plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(imps.all_Fusi_SSIl4_sh);
plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)

text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)

title('Fusi MI xy x bindFF Lane 4, grey=shuffled')
xlabel('Information bits')
ylabel('Fraction')

%Now do tsne analysis with all_info_mutual_info14
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
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content ' ]);

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
            % case 4
            % plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255]);
    end
end

text(-60,70,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(-60,60,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(-60,50,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE MI for lanes, xy and dFF vs. MI for xy in each lane' ]);

%Now do tsne analysis with all_info_mutual_info_dFFbin_xy_op_bin
data2 = [all_info_lane1; all_info_lane4; all_info_mutual_info_dFFbin_xy_op_bin]; % all_info_mutual_info_dFFbin_xy_op_bin

% Step 1 Transpose the data if necessary (N samples x 3 variables)
data2 = data2';
 
rng('default') % for reproducibility
 
% Step 2: Run t-SNE
% The default parameters are usually sufficient, but you can adjust them
% For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
% mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
mappedXl = tsne(data2,'Algorithm','exact','Distance','cosine'); %works better

% mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
% mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK

% Step 3: Plot the results
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedXl(:,1),mappedXl(:,2),'.','Color',[0.7 0.7 0.7]);
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content for odor plume, space' ]);

%Cool, it looks like we have three clear clusters

% Perform k-means clustering
%kmeans misclassifies a small number of points
% k = 3; % Number of clusters
% [idx, centroids] = kmeans(mappedX, k,'Distance','cityblock');

gml = fitgmdist(mappedXl, 2); % Assuming 3 clusters
idxl = cluster(gml, mappedXl);


% %Report the mean infos for each cluster
% for ii_k=1:2
%     mean_lane1l=mean(all_info_lane1(idxl==ii_k));
%     mean_lane4l=mean(all_info_lane4(idxl==ii_k));
%     mean_lane=mean(all_info_lane(idxl==ii_k));
% 
%     fprintf(1, ['\nInformation content for cluster ' num2str(ii_k) ' lane1, lane 2, lane: '...
%         num2str(mean_lane1) ' ' num2str(mean_lane4) ' ' num2str(mean_mutual14) '\n'])
% end

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
for ii_k=1:2
    switch ii_k
        case 1
            plot(mappedXl(idxl==ii_k,1),mappedXl(idxl==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        case 2
            plot(mappedXl(idxl==ii_k,1),mappedXl(idxl==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        case 3
            plot(mappedXl(idxl==ii_k,1),mappedXl(idxl==ii_k,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
    end
end

text(-60,70,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(-60,60,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(-60,50,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE MI for odor plume xy and dFF vs. MI for xy in each lane' ]);



%Now do tsne analysis with all_info_mutual_info_dFFbin_xy_op_bin
data3 = [imps.all_Fusi_SSIl1'; imps.all_Fusi_SSIl4'; all_info_mutual_info_dFFbin_xy_op_bin]; % all_info_mutual_info_dFFbin_xy_op_bin

% Step 1 Transpose the data if necessary (N samples x 3 variables)
data3 = data3';
 
rng('default') % for reproducibility
 
% Step 2: Run t-SNE
% The default parameters are usually sufficient, but you can adjust them
% For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
% mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
mappedXF = tsne(data3,'Algorithm','exact','Distance','cosine'); %works better

% mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
% mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK

% Step 3: Plot the results
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedXF(:,1),mappedXF(:,2),'.','Color',[0.7 0.7 0.7]);
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content for odor plume, Fusi space' ]);

%Cool, it looks like we have three clear clusters

% Perform k-means clustering
%kmeans misclassifies a small number of points
% k = 3; % Number of clusters
% [idx, centroids] = kmeans(mappedX, k,'Distance','cityblock');

gmF = fitgmdist(mappedXF, 2); % Assuming 3 clusters
idxF = cluster(gmF, mappedXF);


% %Report the mean infos for each cluster
% for ii_k=1:2
%     mean_lane1l=mean(all_info_lane1(idxl==ii_k));
%     mean_lane4l=mean(all_info_lane4(idxl==ii_k));
%     mean_lane=mean(all_info_lane(idxl==ii_k));
% 
%     fprintf(1, ['\nInformation content for cluster ' num2str(ii_k) ' lane1, lane 2, lane: '...
%         num2str(mean_lane1) ' ' num2str(mean_lane4) ' ' num2str(mean_mutual14) '\n'])
% end

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
for ii_k=1:2
    switch ii_k
        case 1
            plot(mappedXF(idxF==ii_k,1),mappedXF(idxF==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        case 2
            plot(mappedXF(idxF==ii_k,1),mappedXF(idxF==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
    end
end

text(-60,70,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(-60,60,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(-60,50,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE MI for odor plume, xy and dFF vs. MI for Fusi space' ]);

%Estimate those significantly above SD

% Step 3: Plot the results
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
sig_all_info_lane1=zeros(1,all_info_ii);
for ii_ROI=1:all_info_ii
    if all_info_lane1(ii_ROI)>=3*sd_all_info_lane1_sh(ii_ROI)+mean_all_info_lane1_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idx(ii_ROI);
        sig_all_info_lane1(ii_ROI)=1;
        % switch ii_k
        %     case 1
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        %     case 2
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        %     case 3
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        % end
    end
end

% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content,info 1 significant ' ]);

%plot in color those significantly above SD
% Step 3: Plot the results
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
sig_all_info_lane4=zeros(1,all_info_ii);
for ii_ROI=1:all_info_ii
    if all_info_lane4(ii_ROI)>=3*sd_all_info_lane4_sh(ii_ROI)+mean_all_info_lane4_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idx(ii_ROI);
        sig_all_info_lane4(ii_ROI)=1;
        % switch ii_k
        %     case 1
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        %     case 2
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        %     case 3
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        % end
    end
end
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content,info 4 significant ' ]);


%plot in color those significantly above SD
% Step 3: Plot the results
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
sig_all_info_mutual_info14=zeros(1,all_info_ii);
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
for ii_ROI=1:all_info_ii
    if all_info_mutual_info14(ii_ROI)>=3*sd_all_info_mutual_info14_sh(ii_ROI)+mean_all_info_mutual_info14_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idx(ii_ROI);
        sig_all_info_mutual_info14(ii_ROI)=1;
        % switch ii_k
        %     case 1
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        %     case 2
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        %     case 3
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        % end
    end
end

%Now keep track of significnat for dFFbin_xy_op_bin

%lane 1
ii_ROI_positive=0;
sig_all_info_lane1_dFFbin_xy_op_bin=zeros(1,all_info_ii);
for ii_ROI=1:all_info_ii
    if all_info_lane1(ii_ROI)>=3*sd_all_info_lane1_sh(ii_ROI)+mean_all_info_lane1_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idxl(ii_ROI);
        sig_all_info_lane1_dFFbin_xy_op_bin(ii_ROI)=1;
        % switch ii_k
        %     case 1
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        %     case 2
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        %     case 3
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        % end
    end
end

% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content,info 1 significant ' ]);

%plot in color those significantly above SD
% Step 3: Plot the results
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
sig_all_info_lane4_dFFbin_xy_op_bin=zeros(1,all_info_ii);
for ii_ROI=1:all_info_ii
    if all_info_lane4(ii_ROI)>=3*sd_all_info_lane4_sh(ii_ROI)+mean_all_info_lane4_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idxl(ii_ROI);
        sig_all_info_lane4_dFFbin_xy_op_bin(ii_ROI)=1;
        % switch ii_k
        %     case 1
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        %     case 2
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        %     case 3
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        % end
    end
end


%all_info_mutual_info_dFFbin_xy_op_bin
sig_all_info_mutual_dFFbin_xy_op_bin=zeros(1,all_info_ii);
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
ii_ROI_positive=0;
for ii_ROI=1:all_info_ii
    if all_info_mutual_info_dFFbin_xy_op_bin(ii_ROI)>=3*sd_all_info_mutual_info_dFFbin_xy_op_bin_sh(ii_ROI)+mean_all_info_mutual_info_dFFbin_xy_op_bin_sh(ii_ROI)
        ii_ROI_positive=ii_ROI_positive+1;
        ii_k=idxl(ii_ROI);
        sig_all_info_mutual_dFFbin_xy_op_bin(ii_ROI)=1;
        % switch ii_k
        %     case 1
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        %     case 2
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        %     case 3
        %         plot(mappedX(ii_ROI,1),mappedX(ii_ROI,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
        % end
    end
end

% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content,mutual info 1,4 significant ' ]);

%Now plot the fraction of significant information for each cluster for
%info14
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

title(['Significant spatial information for the different clusters, lanes space'])
ylabel('Fraction SSI')
ylim([0 1])
xlim([-1 11])


%Now plot the fraction of significant information for each cluster for
%dFFbin_xy_op_bin
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
for ii_k=1:2
    this_fraction=sum(sig_all_info_lane1_dFFbin_xy_op_bin(idxl==ii_k))/sum(idxl==ii_k);
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
for ii_k=1:2
    this_fraction=sum(sig_all_info_lane4_dFFbin_xy_op_bin(idxl==ii_k))/sum(idxl==ii_k);
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
for ii_k=1:2
    this_fraction=sum(sig_all_info_mutual_dFFbin_xy_op_bin(idxl==ii_k))/sum(idxl==ii_k);
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
% text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2 4 5 6])
xticklabels({'Lane 1','Lane 4','dFF xy op','Lane 1','Lane 4','dFF xy_op'})
xtickangle(45)

title(['Significant spatial information for the different clusters, space odor plume'])
ylabel('Fraction SSI')
ylim([0 1])
xlim([-1 7])

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


%Now plot the spatial information for each cluster
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
for ii_k=1:2
    these_SSI=all_info_lane1(idxl==ii_k);
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
for ii_k=1:2
    these_SSI=all_info_lane4(idxl==ii_k);
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
for ii_k=1:2
    these_SSI=all_info_mutual_info14(idxl==ii_k);
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

xticks([0 1 2 4 5 6 ])
xticklabels({'Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4'})
xtickangle(45)

title(['Spatial information for the different clusters, xy dFF op'])
ylabel('SI')
ylim([0 0.5])
xlim([-1 7])




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

title(['Fusi significant spatial information for the different clusters, space lane'])
ylabel('SSI')
ylim([0 15])
xlim([-1 11])


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
for ii_k=1:2
    these_Fusi_SSI=imps.all_Fusi_SSIl1(idxl==ii_k);
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
for ii_k=1:2
    these_Fusi_SSI=imps.all_Fusi_SSIl4(idxl==ii_k);
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
for ii_k=1:2
    these_Fusi_SSI=imps.all_Fusi_SSI(idxl==ii_k);
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

xticks([0 1 2 4 5 6])
xticklabels({'Lane 1','Lane 4','Both','Lane 1','Lane 4','Both'})
xtickangle(45)

title(['Fusi significant spatial information for the different clusters, xy dFF op'])
ylabel('SSI')
ylim([0 15])
xlim([-1 7])

%Now explore whether there is a relationship between 
%prediction importance and these spatial
%information clusteers



% 
% %plot in color prediciton importance categories
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
% 
% for ii=1:length(imps.mean_imp_ROI)
%     this_file=imps.mean_imp_fileNo(ii);
%     this_ROI=imps.mean_imp_ROI(ii);
% 
%     ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
%     this_ii_k=idx(ii_ROI);
% 
%     %1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s
%     switch imps.ii_95_class(ii)
%         case 1
%             %1 is >=x95
%             plot(mappedX(ii_info,1),mappedX(ii_info,2),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0],'MarkerSize',5);
%         case 2
%             %2 is >=y95
%             plot(mappedX(ii_info,1),mappedX(ii_info,2),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9],'MarkerSize',5);
%         case 3
%             %3 is >=conc95
%             plot(mappedX(ii_info,1),mappedX(ii_info,2),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5],'MarkerSize',5);
%         case 4
%             %4 is >=x95 and conc95
%             plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor',[0.95 0.9 0.25],'MarkerSize',5);
%         case 5
%             %5 is >=y95 and conc95
%             plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7],'MarkerSize',5);
%         case 6
%             %6 is >=x95 and y95
%             plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0],'MarkerSize',5);
%         case 7
%             %7 is >= all 95s
%             plot(mappedX(ii_info,1),mappedX(ii_info,2),'s','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor',[0.8 0.6 0.7],'MarkerSize',5);
%     end
% end
% 
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content,prediction importance ' ]);

%ii_95_class 1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s

% 
% %plot in color prediciton importance categories
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% plot(mappedXl(:,1),mappedXl(:,2),'.','Color',[0.7 0.7 0.7]);
% 
% for ii=1:length(imps.mean_imp_ROI)
%     this_file=imps.mean_imp_fileNo(ii);
%     this_ROI=imps.mean_imp_ROI(ii);
% 
%     ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
%     this_ii_k=idx(ii_ROI);
% 
%     %1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s
%     switch imps.ii_95_class(ii)
%         case 1
%             %1 is >=x95
%             plot(mappedXl(ii_info,1),mappedXl(ii_info,2),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0],'MarkerSize',5);
%         case 2
%             %2 is >=y95
%             plot(mappedXl(ii_info,1),mappedXl(ii_info,2),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9],'MarkerSize',5);
%         case 3
%             %3 is >=conc95
%             plot(mappedXl(ii_info,1),mappedXl(ii_info,2),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5],'MarkerSize',5);
%         case 4
%             %4 is >=x95 and conc95
%             plot(mappedXl(ii_info,1),mappedXl(ii_info,2),'s','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor',[0.95 0.9 0.25],'MarkerSize',5);
%         case 5
%             %5 is >=y95 and conc95
%             plot(mappedXl(ii_info,1),mappedXl(ii_info,2),'s','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7],'MarkerSize',5);
%         case 6
%             %6 is >=x95 and y95
%             plot(mappedXl(ii_info,1),mappedXl(ii_info,2),'s','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0],'MarkerSize',5);
%         case 7
%             %7 is >= all 95s
%             plot(mappedXl(ii_info,1),mappedXl(ii_info,2),'s','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor',[0.8 0.6 0.7],'MarkerSize',5);
%     end
% end
% 
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE Information Content including odor plume,prediction importance ' ]);

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

title(['Fraction of prediction important ROIs for the different clusters space, lane'])
ylabel('Fraction')
ylim([0 1])
xlim([-1 11])


%Now plot the fraction of significant predictive importance
sig_pred_imp_x_ii_k=[];
ii_x=0;

sig_pred_imp_y_ii_k=[];
ii_y=0;

sig_pred_imp_conc_ii_k=[];
ii_conc=0;

all_info_lane1_impx=[];
all_info_lane4_impx=[];
all_info_mutual_info_dFFbin_xy_op_bin_impx=[];

all_info_lane1_impy=[];
all_info_lane4_impy=[];
all_info_mutual_info_dFFbin_xy_op_bin_impy=[];

all_info_lane1_impconc=[];
all_info_lane4_impconc=[];
all_info_mutual_info_dFFbin_xy_op_bin_impconc=[];

for ii=1:length(imps.mean_imp_ROI)
    this_file=imps.mean_imp_fileNo(ii);
    this_ROI=imps.mean_imp_ROI(ii);

    ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
    this_ii_k=idxl(ii_info);

    if imps.sig_pred_imp_x(ii)==1
        ii_x=ii_x+1;
        sig_pred_imp_x_ii_k(ii_x)=this_ii_k;
        all_info_lane1_impx=[all_info_lane1_impx all_info_lane1(ii_info)];
        all_info_lane4_impx=[all_info_lane4_impx all_info_lane4(ii_info)];
        all_info_mutual_info_dFFbin_xy_op_bin_impx=[all_info_mutual_info_dFFbin_xy_op_bin_impx all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];

    end


    if imps.sig_pred_imp_y(ii)==1
        ii_y=ii_y+1;
        sig_pred_imp_y_ii_k(ii_y)=this_ii_k;
        all_info_lane1_impy=[all_info_lane1_impy all_info_lane1(ii_info)];
        all_info_lane4_impy=[all_info_lane4_impy all_info_lane4(ii_info)];
        all_info_mutual_info_dFFbin_xy_op_bin_impy=[all_info_mutual_info_dFFbin_xy_op_bin_impy all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
    end


    if imps.sig_pred_imp_conc(ii)==1
        ii_conc=ii_conc+1;
        sig_pred_imp_conc_ii_k(ii_conc)=this_ii_k;
        all_info_lane1_impconc=[all_info_lane1_impconc all_info_lane1(ii_info)];
        all_info_lane4_impconc=[all_info_lane4_impconc all_info_lane4(ii_info)];
        all_info_mutual_info_dFFbin_xy_op_bin_impconc=[all_info_mutual_info_dFFbin_xy_op_bin_impconc all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
    end

end

%Plot the MIs for the imps

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

edges=[0:0.025:1];
rand_offset=0.5;

%MI for predictive importance for x
bar(bar_offset,mean(all_info_lane1_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_lane1_impx...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar(bar_offset,mean(all_info_lane4_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_lane4_impx...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar(bar_offset,mean(all_info_mutual_info_dFFbin_xy_op_bin_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_mutual_info_dFFbin_xy_op_bin_impx...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar_offset=bar_offset+1;

%MI for predictive importance for y
bar(bar_offset,mean(all_info_lane1_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_lane1_impy...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar(bar_offset,mean(all_info_lane4_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_lane4_impy...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar(bar_offset,mean(all_info_mutual_info_dFFbin_xy_op_bin_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_mutual_info_dFFbin_xy_op_bin_impy...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar_offset=bar_offset+1;


%MI for predictive importance for conc
bar(bar_offset,mean(all_info_lane1_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_lane1_impconc...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar(bar_offset,mean(all_info_lane4_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_lane4_impconc...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;

bar(bar_offset,mean(all_info_mutual_info_dFFbin_xy_op_bin_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

[mean_out, CIout, violin_x]=drgViolinPoint(all_info_mutual_info_dFFbin_xy_op_bin_impconc...
    ,edges,bar_offset,rand_offset,'k','k',2);
bar_offset=bar_offset+1;



text(4,0.7,'MI lane 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(4,0.62,'MI lane 4','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(4,0.54,'MI xy dFF op 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2 4 5 6 8 9 10])
xticklabels({'x','y','odor','x','y','odor','x','y','odor'})
xtickangle(45)

title(['MI for prediction important ROIs including godor plume'])
ylabel('MI')
ylim([0 1])
xlim([-1 11])


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
for ii_k=1:2
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
for ii_k=1:2
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
for ii_k=1:2
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

title(['Fraction of prediction important ROIs including odor plume'])
ylabel('Fraction')
ylim([0 1])
xlim([-1 11])

% 
% %Now plot the fraction of significant predictive importance for space Fusi,
% %odor plume
% sig_predF_imp_x_ii_k=[];
% ii_x=0;
% 
% sig_predF_imp_y_ii_k=[];
% ii_y=0;
% 
% sig_predF_imp_conc_ii_k=[];
% ii_conc=0;
% 
% all_infoF_lane1_impx=[];
% all_infoF_lane4_impx=[];
% all_infoF_mutual_info_dFFbin_xy_op_bin_impx=[];
% 
% all_infoF_lane1_impy=[];
% all_infoF_lane4_impy=[];
% all_infoF_mutual_info_dFFbin_xy_op_bin_impy=[];
% 
% all_infoF_lane1_impconc=[];
% all_infoF_lane4_impconc=[];
% all_infoF_mutual_info_dFFbin_xy_op_bin_impconc=[];
% 
% for ii=1:length(imps.mean_imp_ROI)
%     this_file=imps.mean_imp_fileNo(ii);
%     this_ROI=imps.mean_imp_ROI(ii);
% 
%     ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
%     this_ii_k=idxF(ii_info);
% 
%     if imps.sig_pred_imp_x(ii)==1
%         ii_x=ii_x+1;
%         sig_predF_imp_x_ii_k(ii_x)=this_ii_k;
%         all_infoF_lane1_impx=[all_infoF_lane1_impx all_info_lane1(ii_info)];
%         all_infoF_lane4_impx=[all_infoF_lane4_impx all_info_lane4(ii_info)];
%         all_infoF_mutual_info_dFFbin_xy_op_bin_impx=[all_infoF_mutual_info_dFFbin_xy_op_bin_impx all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
% 
%     end
% 
% 
%     if imps.sig_pred_imp_y(ii)==1
%         ii_y=ii_y+1;
%         sig_predF_imp_y_ii_k(ii_y)=this_ii_k;
%         all_infoF_lane1_impy=[all_infoF_lane1_impy all_info_lane1(ii_info)];
%         all_infoF_lane4_impy=[all_infoF_lane4_impy all_info_lane4(ii_info)];
%         all_infoF_mutual_info_dFFbin_xy_op_bin_impy=[all_infoF_mutual_info_dFFbin_xy_op_bin_impy all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
%     end
% 
% 
%     if imps.sig_pred_imp_conc(ii)==1
%         ii_conc=ii_conc+1;
%         sig_predF_imp_conc_ii_k(ii_conc)=this_ii_k;
%         all_infoF_lane1_impconc=[all_infoF_lane1_impconc all_info_lane1(ii_info)];
%         all_infoF_lane4_impconc=[all_infoF_lane4_impconc all_info_lane4(ii_info)];
%         all_infoF_mutual_info_dFFbin_xy_op_bin_impconc=[all_infoF_mutual_info_dFFbin_xy_op_bin_impconc all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
%     end
% 
% end
% 
% %Plot the MIs for the imps
% 
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% bar_offset=0;
% 
% edges=[0:0.025:1];
% rand_offset=0.5;
% 
% %MI for predictive importance for x
% bar(bar_offset,mean(all_infoF_lane1_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_lane1_impx...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_lane4_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_lane4_impx...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_mutual_info_dFFbin_xy_op_bin_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_mutual_info_dFFbin_xy_op_bin_impx...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar_offset=bar_offset+1;
% 
% %MI for predictive importance for y
% bar(bar_offset,mean(all_infoF_lane1_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_lane1_impy...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_lane4_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_lane4_impy...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_mutual_info_dFFbin_xy_op_bin_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_mutual_info_dFFbin_xy_op_bin_impy...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar_offset=bar_offset+1;
% 
% 
% %MI for predictive importance for conc
% bar(bar_offset,mean(all_infoF_lane1_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_lane1_impconc...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_lane4_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_lane4_impconc...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_mutual_info_dFFbin_xy_op_bin_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_mutual_info_dFFbin_xy_op_bin_impconc...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% 
% 
% text(4,0.7,'MI lane 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,0.62,'MI lane 4','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,0.54,'MI xy dFF op 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2 4 5 6 8 9 10])
% xticklabels({'x','y','odor','x','y','odor','x','y','odor'})
% xtickangle(45)
% 
% title(['MI for prediction important ROIs including godor plume'])
% ylabel('MI')
% ylim([0 1])
% xlim([-1 11])
% 
% 
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% bar_offset=0;
% 
% %predictive importance for x
% for ii_k=1:2
%     this_fraction=sum(sig_predF_imp_x_ii_k==ii_k)/length(sig_predF_imp_x_ii_k);
% 
%     switch ii_k
%         case 1
%             bar(bar_offset,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%         case 2
%             bar(bar_offset+4,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%         case 3
%             bar(bar_offset+8,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%     end
% 
% end
% 
% 
% 
% %predictive importance for y
% for ii_k=1:2
%     this_fraction=sum(sig_predF_imp_y_ii_k==ii_k)/length(sig_predF_imp_y_ii_k);
%        switch ii_k
%             case 1
%                 bar(bar_offset+1,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             case 2
%                 bar(bar_offset+5,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             case 3
%                 bar(bar_offset+9,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%        end
% 
% end
% 
% 
% 
% %predictive importance for conc 
% for ii_k=1:2
%     this_fraction=sum(sig_predF_imp_conc_ii_k==ii_k)/length(sig_predF_imp_conc_ii_k);
%        switch ii_k
%             case 1
%                 bar(bar_offset+2,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             case 2
%                 bar(bar_offset+6,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             case 3
%                 bar(bar_offset+10,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%        end
% 
% end
% 
% 
% text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2 4 5 6 8 9 10])
% xticklabels({'x','y','odor','x','y','odor','x','y','odor'})
% xtickangle(45)
% 
% title(['Fraction of prediction important ROIs including odor plume'])
% ylabel('Fraction')
% ylim([0 1])
% xlim([-1 11])

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


%Spatial correlations of the spatial pattern
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .2 .3])
bar_offset=0;
edges=[0:0.05:1];
rand_offset=0.5;

%predictive importance for x
for ii_k=1:2
    these_rhos=imps.all_spatial_rhol1l4(idxl==ii_k);
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

xticks([0 1])
xticklabels({'Cluster 1','Cluster 2'})
xtickangle(45)

title(['Spatial correlation for the different clusters, xy dFF op'])
ylabel('Rho')
ylim([-1 1])
xlim([-1 2])


%Center of mass distances of the spatial pattern
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .2 .3])
bar_offset=0;
edges=[0:25:500];
rand_offset=0.5;

%predictive importance for x
for ii_k=1:2
    these_dcoms=imps.all_delta_center_of_mass(idxl==ii_k);
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

xticks([0 1])
xticklabels({'Cluster 1','Cluster 2'})
xtickangle(45)

title(['Difference in the center of mass, xy dFF op'])
ylabel('delta com')
ylim([0 350])
xlim([-1 2])



%Show a subset of the spatial maps
figureNo=figureNo+1;
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



            % delta_below_zero_ii=max(all_ii_turns(include_trial==1));
            % delta_above_zero_ii=max(all_ii_ends(include_trial==1)-all_ii_turns(include_trial==1));

            

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

                yticks([70 100 200 300 400 430])
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
    end
end

fclose(fileID);

pffft=1;

