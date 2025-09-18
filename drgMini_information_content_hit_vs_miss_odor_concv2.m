%drgMini_information_contentv3
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
        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc04192024/';
        choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_04192024.m';

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        %Angle file
        % save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        % choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        %This one has the odor encounter
        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle05152025/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_05102025.m';

        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        %This is not used here, all the Moser infor is input through PredImp
        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';


        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        save_PathPredImp='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FilePredImp='outputPredictionImportanceHitMiss.mat';

        %The output of drgMini_information_contentv2 is saved here
        save_PathIC='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FileIC='outputPerROIInformationContentHitMiss.mat';

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_hitmiss_stats.txt','w');
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

lowest_conc=-200;
op_threshold=-3.5;
odor_c_bounds=[-10:8/10:-2];
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

    minC=prctile([mean_plume_l4(:); mean_plume_l1(:)],0.5);
    if minC<lowest_conc
        minC=lowest_conc;
    end
    maxC=max([max(mean_plume_l4(:)) max(mean_plume_l1(:))]);
    if minC==maxC
        maxC=minC+0.1;
    end
    odor_plume_patterns.cm_from_floor(cm_from_floor).minC=minC;
end

%Now calculate mutual information
figNo=figureNo+1;
x=25:50:475;
y=24:48:456;

all_info_ii_max=5000;

all_info_ii=0;
all_info_fileNo=[];
all_info_ii_ROI=[];


all_info_xy_hit=[];
all_info_xy_miss=[];
all_info_xy=[];

all_op_bin_info=[];
all_op_bin_info_hit=[];
all_op_bin_info_miss=[];

all_info_op=[];
all_info_op_hit=[];
all_info_op_miss=[];

all_info_xy_hit_sh=[];
all_info_xy_miss_sh=[];
all_info_xy_sh=[];

all_op_bin_info_sh=[];
all_op_bin_info_hit_sh=[];
all_op_bin_info_miss_sh=[];

all_info_op_sh=[];
all_info_op_hit_sh=[];
all_info_op_miss_sh=[];

handles_outic=[];

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

            this_Hit_activity=zeros(10,10);
            this_Hit_activity_n=zeros(10,10);
            sum_Hit_activity=0;

            this_Miss_activity=zeros(10,10);
            this_Miss_activity_n=zeros(10,10);
            sum_Miss_activity=0;

            %Cumulative counts for info

            %Both hit and miss
            cum_xy=zeros(1,10*10);
            cum_bindFF=zeros(1,2);
            cum_xy_bindFF=zeros(2,10*10);
            cum_op_bin=zeros(1,2);
            cum_op_bin_bindFF=zeros(2,2);
            cum_op=zeros(1,10);
            cum_op_bindFF=zeros(2,10);
            include_xy=zeros(1,10*10);
            include_op=zeros(1,10);

            %Hit
            cum_xy_hit=zeros(1,10*10);
            cum_bindFF_hit=zeros(1,2);
            cum_xy_bindFF_hit=zeros(2,10*10);
            cum_op_bin_hit=zeros(1,2);
            cum_op_bin_bindFF_hit=zeros(2,2);
            cum_op_hit=zeros(1,10);
            cum_op_bindFF_hit=zeros(2,10);
            include_xy_hit=zeros(1,10*10);
            include_op_hit=zeros(1,10);

            %Miss
            cum_xy_miss=zeros(1,10*10);
            cum_bindFF_miss=zeros(1,2);
            cum_xy_bindFF_miss=zeros(2,10*10);
            cum_op_bin_miss=zeros(1,2);
            cum_op_bin_bindFF_miss=zeros(2,2);
            cum_op_miss=zeros(1,10);
            cum_op_bindFF_miss=zeros(2,10);
            include_xy_miss=zeros(1,10*10);
            include_op_miss=zeros(1,10);




            all_these_x=[];
            all_these_y=[];
            all_these_dFF=[];
            all_hit_vs_miss=[];
            all_x_ii=[];
            all_trNo=[];
            all_ii_t=[];


            for trNo=1:trials.odor_trNo
                these_x=trials.trial(trNo).XYtest(:,1);
                these_y=trials.trial(trNo).XYtest(:,2);
                these_dFF=trials.trial(trNo).XdFFtest(:,ii_ROI);

                all_these_x=[all_these_x; these_x];
                all_these_y=[all_these_y; these_y];
                all_these_dFF=[all_these_dFF; these_dFF];
                if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                    is_this_hit=1;
                else
                    is_this_hit=0;
                end
                all_hit_vs_miss=[all_hit_vs_miss; is_this_hit*ones(length(these_x),1)];

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

                    this_xy_ii=this_x_ii+10*(this_y_ii-1);
                    include_xy(this_xy_ii)=1;
                    cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
                    cum_xy_bindFF(this_bin_dFF,this_xy_ii)=cum_xy_bindFF(this_bin_dFF,this_xy_ii)+1;

                    %Where are odor start and end?
                    op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                    op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                    x_ii=op_predictedstart+ii_t;
                    all_x_ii=[all_x_ii x_ii];
                    all_trNo=[all_trNo trNo];
                    all_ii_t=[all_ii_t ii_t];

                    if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                        % Hit
                        cum_xy_bindFF_hit(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_hit(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_hit(this_xy_ii)=cum_xy_hit(this_xy_ii)+1;
                        cum_bindFF_hit(this_bin_dFF)=cum_bindFF_hit(this_bin_dFF)+1;
                        include_xy_hit(this_xy_ii)=1;
                    else
                        %Miss
                        cum_xy_bindFF_miss(this_bin_dFF,this_xy_ii)=cum_xy_bindFF_miss(this_bin_dFF,this_xy_ii)+1;
                        cum_xy_miss(this_xy_ii)=cum_xy_miss(this_xy_ii)+1;
                        cum_bindFF_miss(this_bin_dFF)=cum_bindFF_miss(this_bin_dFF)+1;
                        include_xy_miss(this_xy_ii)=1;
                    end


                    %Now figure out cums for op
                    switch handles_conc.group(fileNo)
                        case 1
                            %2 cm from the floor
                            cm_from_floor=2;
                        case 5
                            %1 cm from the floor
                            cm_from_floor=1;
                    end

                    %Find the binary op
                    if trials.lane_per_trial(trNo)==1
                        binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                    else
                        binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);
                    end
                    %Make odor zero if the time is early or late
                    if x_ii<trials.odor_ii_start(trNo)
                        binary_op=0;
                    else
                        if x_ii>trials.odor_ii_end(trNo)
                            binary_op=0;
                        else
                            this_x=trials.trial(trNo).XYtest(ii_t,1);
                            x_on=(x_ii-trials.odor_ii_start(trNo))*handles_conc.dt*handles_conc.air_flow_speed;
                            if this_x>x_on
                                binary_op=0;
                            end
                        end
                    end

                    %Now do cums for op_bin
                    cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1;
                    cum_op_bin_bindFF(this_bin_dFF,binary_op+1)=cum_op_bin_bindFF(this_bin_dFF,binary_op+1)+1;

                    if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                        %Hit
                        cum_op_bin_hit(binary_op+1)=cum_op_bin_hit(binary_op+1)+1;
                        cum_op_bin_bindFF_hit(this_bin_dFF,binary_op+1)=...
                            cum_op_bin_bindFF_hit(this_bin_dFF,binary_op+1)+1;
                    else
                        %Miss
                        cum_op_bin_miss(binary_op+1)=cum_op_bin_miss(binary_op+1)+1;
                        cum_op_bin_bindFF_miss(this_bin_dFF,binary_op+1)=...
                            cum_op_bin_bindFF_miss(this_bin_dFF,binary_op+1)+1;
                    end

                    %Find the analog op
                    minC=odor_plume_patterns.cm_from_floor(cm_from_floor).minC;
                    this_x=trials.trial(trNo).XYtest(ii_t,1);
                    this_y=trials.trial(trNo).XYtest(ii_t,2);

                    [M,this_x_op_ii]=min(abs(odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume-this_x));
                    [M,this_y_op_ii]=min(abs(odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume-this_y));

                    if trials.lane_per_trial(trNo)==1
                        op_conc=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1(this_y_op_ii,this_x_op_ii);
                    else
                        op_conc=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4(this_y_op_ii,this_x_op_ii);
                    end

                    %Make odor minC if the time is early or late
                    if x_ii<trials.odor_ii_start(trNo)
                        op_conc=minC;
                    else
                        if x_ii>trials.odor_ii_end(trNo)
                            op_conc=minC;
                        else
                            this_x=trials.trial(trNo).XYtest(ii_t,1);
                            x_on=(x_ii-trials.odor_ii_start(trNo))*handles_conc.dt*handles_conc.air_flow_speed;
                            if this_x>x_on
                                op_conc=minC;
                            end
                        end
                    end
                    for ii_b=1:10
                        if (op_conc>=odor_c_bounds(ii_b))&(op_conc<odor_c_bounds(ii_b+1))
                            this_ii_op=ii_b;
                        end
                    end

                    %Now do cums for op analogue
                    cum_op(this_ii_op)=cum_op(this_ii_op)+1;
                    cum_op_bindFF(this_bin_dFF,this_ii_op)=cum_op_bindFF(this_bin_dFF,this_ii_op)+1;
                    include_op(this_ii_op)=1;

                    if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                        %Hit
                        cum_op_hit(binary_op+1)=cum_op_hit(binary_op+1)+1;
                        cum_op_bindFF_hit(this_bin_dFF,binary_op+1)=...
                            cum_op_bindFF_hit(this_bin_dFF,binary_op+1)+1;
                        include_op_hit(this_ii_op)=1;
                    else
                        %Miss
                        cum_op_miss(binary_op+1)=cum_op_miss(binary_op+1)+1;
                        cum_op_bindFF_miss(this_bin_dFF,binary_op+1)=...
                            cum_op_bindFF_miss(this_bin_dFF,binary_op+1)+1;
                        include_op_miss(this_ii_op)=1;
                    end

                end
            end

            %Calculate information

            p_xy=cum_xy(logical(include_xy))/sum(cum_xy(logical(include_xy)));
            p_bindFF=cum_bindFF/sum(cum_bindFF(:));
            included_cum_xy_bindFF=cum_xy_bindFF(:,logical(include_xy));
            p_xy_bindFF=included_cum_xy_bindFF/sum(included_cum_xy_bindFF(:));

            p_xy_hit=cum_xy_hit(logical(include_xy_hit))/sum(cum_xy_hit(logical(include_xy_hit)));
            p_bindFF_hit=cum_bindFF_hit/sum(cum_bindFF_hit(:));
            included_cum_xy_bindFF_hit=cum_xy_bindFF_hit(:,logical(include_xy_hit));
            p_xy_bindFF_hit=included_cum_xy_bindFF_hit/sum(included_cum_xy_bindFF_hit(:));

            p_xy_miss=cum_xy_miss(logical(include_xy_miss))/sum(cum_xy_miss(logical(include_xy_miss)));
            p_bindFF_miss=cum_bindFF_miss/sum(cum_bindFF_miss(:));
            included_cum_xy_bindFF_miss=cum_xy_bindFF_miss(:,logical(include_xy_miss));
            p_xy_bindFF_miss=included_cum_xy_bindFF_miss/sum(included_cum_xy_bindFF_miss(:));


            %Binary op
            p_op_bin=cum_op_bin/sum(cum_op_bin);
            p_op_bin_bindFF=cum_op_bin_bindFF/sum(cum_op_bin_bindFF(:));

            p_op_bin_hit=cum_op_bin/sum(cum_op_bin_hit);
            p_op_bin_bindFF_hit=cum_op_bin_bindFF_hit/sum(cum_op_bin_bindFF_hit(:));

            p_op_bin_miss=cum_op_bin/sum(cum_op_bin_miss);
            p_op_bin_bindFF_miss=cum_op_bin_bindFF_miss/sum(cum_op_bin_bindFF_miss(:));

            %Analog op
            included_cum_op=cum_op(logical(include_op));
            p_op=included_cum_op/sum(included_cum_op);
            included_cum_op_bindFF=cum_op_bindFF(:,logical(include_op));
            p_op_bindFF=included_cum_op_bindFF/sum(included_cum_op_bindFF(:));

            included_cum_op_hit=cum_op_hit(logical(include_op_hit));
            p_op_hit=included_cum_op_hit/sum(included_cum_op_hit);
            included_cum_op_bindFF_hit=cum_op_bindFF_hit(:,logical(include_op_hit));
            p_op_bindFF_hit=included_cum_op_bindFF_hit/sum(included_cum_op_bindFF_hit(:));

            included_cum_op_miss=cum_op_miss(logical(include_op_miss));
            p_op_miss=included_cum_op_miss/sum(included_cum_op_miss);
            included_cum_op_bindFF_miss=cum_op_bindFF_miss(:,logical(include_op_miss));
            p_op_bindFF_miss=included_cum_op_bindFF_miss/sum(included_cum_op_bindFF_miss(:));


            all_info_ii=all_info_ii+1;
            all_info_fileNo(all_info_ii)=fileNo;
            all_info_ii_ROI(all_info_ii)=ii_ROI;

            handles_outic.all_info_ii=all_info_ii;
            handles_outic.fileNo(all_info_ii)=fileNo;
            handles_outic.ii_ROI(all_info_ii)=ii_ROI;

            all_info_xy_hit(all_info_ii)=0;
            all_info_xy_miss(all_info_ii)=0;
            all_info_xy(all_info_ii)=0;

            %Binary op
            all_op_bin_info(all_info_ii)=0;
            all_op_bin_info_hit(all_info_ii)=0;
            all_op_bin_info_miss(all_info_ii)=0;

            %Analogue op
            all_info_op(all_info_ii)=0;
            all_info_op_hit(all_info_ii)=0;
            all_info_op_miss(all_info_ii)=0;

            %Calculate info for hit
            for ii_xy=1:size(p_xy_bindFF_hit,2)
                % if sum(p_xy_bindFF_hit(:,ii_xy))>0
                for ii_bin_dFF=1:2
                    if p_xy_bindFF_hit(ii_bin_dFF,ii_xy)~=0
                        all_info_xy_hit(all_info_ii)=all_info_xy_hit(all_info_ii)+p_xy_bindFF_hit(ii_bin_dFF,ii_xy)*...
                            log2(p_xy_bindFF_hit(ii_bin_dFF,ii_xy)/(p_xy_hit(ii_xy)*p_bindFF_hit(ii_bin_dFF)));
                    end
                end
                % end
            end

            %Calculate info for miss
            for ii_xy=1:size(p_xy_bindFF_miss,2)
                % if sum(p_xy_bindFF_miss(:,ii_xy))>0
                for ii_bin_dFF=1:2
                    if p_xy_bindFF_miss(ii_bin_dFF,ii_xy)~=0
                        all_info_xy_miss(all_info_ii)=all_info_xy_miss(all_info_ii)+p_xy_bindFF_miss(ii_bin_dFF,ii_xy)*...
                            log2(p_xy_bindFF_miss(ii_bin_dFF,ii_xy)/(p_xy_miss(ii_xy)*p_bindFF_miss(ii_bin_dFF)));
                    end
                end
                % end
            end

            %Calculate info for both hit and miss
            for ii_xy=1:size(p_xy_bindFF,2)
                % if sum(p_xy_bindFF_miss(:,ii_xy))>0
                for ii_bin_dFF=1:2
                    if p_xy_bindFF(ii_bin_dFF,ii_xy)~=0
                        all_info_xy(all_info_ii)=all_info_xy(all_info_ii)+p_xy_bindFF(ii_bin_dFF,ii_xy)*...
                            log2(p_xy_bindFF(ii_bin_dFF,ii_xy)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)));
                    end
                end
                % end
            end
 
            %Binary op

            %Calculate info for op_bin
            for ii_op=1:2
                for ii_bin_dFF=1:2
                    if p_op_bin_bindFF_hit(ii_bin_dFF,ii_op)~=0
                        all_op_bin_info_hit(all_info_ii)=all_op_bin_info_hit(all_info_ii)+p_op_bin_bindFF_hit(ii_bin_dFF,ii_op)*...
                            log2(p_op_bin_bindFF_hit(ii_bin_dFF,ii_op)/(p_op_bin_hit(ii_op)*p_bindFF_hit(ii_bin_dFF)));
                    end
                end
            end

            for ii_op=1:2
                for ii_bin_dFF=1:2
                    if p_op_bin_bindFF_miss(ii_bin_dFF,ii_op)~=0
                        all_op_bin_info_miss(all_info_ii)=all_op_bin_info_miss(all_info_ii)+p_op_bin_bindFF_miss(ii_bin_dFF,ii_op)*...
                            log2(p_op_bin_bindFF_miss(ii_bin_dFF,ii_op)/(p_op_bin_miss(ii_op)*p_bindFF_miss(ii_bin_dFF)));
                    end
                end
            end

            for ii_op=1:2
                for ii_bin_dFF=1:2
                    if p_op_bin_bindFF(ii_bin_dFF,ii_op)~=0
                        all_op_bin_info(all_info_ii)=all_op_bin_info(all_info_ii)+p_op_bin_bindFF(ii_bin_dFF,ii_op)*...
                            log2(p_op_bin_bindFF(ii_bin_dFF,ii_op)/(p_op_bin(ii_op)*p_bindFF(ii_bin_dFF)));
                    end
                end
            end


            %Calculate info for op analog
            for ii_op=1:size(p_op_bindFF,2)
                for ii_bin_dFF=1:2
                    if p_op_bindFF(ii_bin_dFF,ii_op)~=0
                        all_info_op(all_info_ii)=all_info_op(all_info_ii)+p_op_bindFF(ii_bin_dFF,ii_op)*...
                            log2(p_op_bindFF(ii_bin_dFF,ii_op)/(p_op(ii_op)*p_bindFF(ii_bin_dFF)));
                    end
                end
            end

            %Calculate info for op analog hit
            for ii_op=1:size(p_op_bindFF_hit,2)
                for ii_bin_dFF=1:2
                    if p_op_bindFF_hit(ii_bin_dFF,ii_op)~=0
                        all_info_op_hit(all_info_ii)=all_info_op_hit(all_info_ii)+p_op_bindFF_hit(ii_bin_dFF,ii_op)*...
                            log2(p_op_bindFF_hit(ii_bin_dFF,ii_op)/(p_op_hit(ii_op)*p_bindFF_hit(ii_bin_dFF)));
                    end
                end
            end

            %Calculate info for op analog miss
            for ii_op=1:size(p_op_bindFF_miss,2)
                for ii_bin_dFF=1:2
                    if p_op_bindFF_miss(ii_bin_dFF,ii_op)~=0
                        all_info_op_miss(all_info_ii)=all_info_op_miss(all_info_ii)+p_op_bindFF_miss(ii_bin_dFF,ii_op)*...
                            log2(p_op_bindFF_miss(ii_bin_dFF,ii_op)/(p_op_miss(ii_op)*p_bindFF_miss(ii_bin_dFF)));
                    end
                end
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

                %Both hit and miss
                cum_xy=zeros(1,10*10);
                cum_bindFF=zeros(1,2);
                cum_xy_bindFF=zeros(2,10*10);
                cum_op_bin=zeros(1,2);
                cum_op_bin_bindFF=zeros(2,2);
                cum_op=zeros(1,10);
                cum_op_bindFF=zeros(2,10);
                include_xy=zeros(1,10*10);
                include_op=zeros(1,10);

                %Hit
                cum_xy_hit=zeros(1,10*10);
                cum_bindFF_hit=zeros(1,2);
                cum_xy_bindFF_hit=zeros(2,10*10);
                cum_op_bin_hit=zeros(1,2);
                cum_op_bin_bindFF_hit=zeros(2,2);
                cum_op_hit=zeros(1,10);
                cum_op_bindFF_hit=zeros(2,10);
                include_xy_hit=zeros(1,10*10);
                include_op_hit=zeros(1,10);

                %Miss
                cum_xy_miss=zeros(1,10*10);
                cum_bindFF_miss=zeros(1,2);
                cum_xy_bindFF_miss=zeros(2,10*10);
                cum_op_bin_miss=zeros(1,2);
                cum_op_bin_bindFF_miss=zeros(2,2);
                cum_op_miss=zeros(1,10);
                cum_op_bindFF_miss=zeros(2,10);
                include_xy_miss=zeros(1,10*10);
                include_op_miss=zeros(1,10);


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

                    this_xy_ii=this_x_ii+10*(this_y_ii-1);
                    this_xy_ii_rev=this_x_ii_rev+10*(this_y_ii_rev-1);
                    include_xy(this_xy_ii)=1;

                    cum_xy(this_xy_ii)=cum_xy(this_xy_ii)+1;
                    cum_xy_bindFF(this_bin_dFF_rev,this_xy_ii)=cum_xy_bindFF(this_bin_dFF_rev,this_xy_ii)+1;


                    %Where are odor start and end?
                    x_ii=all_x_ii(ii_t);
                    trNo=all_trNo(ii_t);
                    this_ii_t=all_ii_t(ii_t);

                    if all_hit_vs_miss(ii_t)==1
                        % Hit
                        cum_xy_bindFF_hit(this_bin_dFF_rev,this_xy_ii)=cum_xy_bindFF_hit(this_bin_dFF_rev,this_xy_ii)+1;
                        cum_xy_hit(this_xy_ii)=cum_xy_hit(this_xy_ii)+1;
                        cum_bindFF_hit(this_bin_dFF_rev)=cum_bindFF_hit(this_bin_dFF_rev)+1;
                        include_xy_hit(this_xy_ii)=1;
                    else
                        %Miss
                        cum_xy_bindFF_miss(this_bin_dFF_rev,this_xy_ii)=cum_xy_bindFF_miss(this_bin_dFF_rev,this_xy_ii)+1;
                        cum_xy_miss(this_xy_ii)=cum_xy_miss(this_xy_ii)+1;
                        cum_bindFF_miss(this_bin_dFF_rev)=cum_bindFF_miss(this_bin_dFF_rev)+1;
                        include_xy_miss(this_xy_ii)=1;
                    end


                    %Now figure out cums for op
                    switch handles_conc.group(fileNo)
                        case 1
                            %2 cm from the floor
                            cm_from_floor=2;
                        case 5
                            %1 cm from the floor
                            cm_from_floor=1;
                    end

                    %Find the binary op
                    if trials.lane_per_trial(trNo)==1
                        binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                    else
                        binary_op=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);
                    end
                    %Make odor zero if the time is early or late
                    if x_ii<trials.odor_ii_start(trNo)
                        binary_op=0;
                    else
                        if x_ii>trials.odor_ii_end(trNo)
                            binary_op=0;
                        else
                            this_x=trials.trial(trNo).XYtest(this_ii_t,1);
                            x_on=(x_ii-trials.odor_ii_start(trNo))*handles_conc.dt*handles_conc.air_flow_speed;
                            if this_x>x_on
                                binary_op=0;
                            end
                        end
                    end

                    %Now do cums for op_bin
                    cum_op_bin(binary_op+1)=cum_op_bin(binary_op+1)+1;
                    cum_op_bin_bindFF(this_bin_dFF_rev,binary_op+1)=cum_op_bin_bindFF(this_bin_dFF_rev,binary_op+1)+1;

                    if all_hit_vs_miss(ii_t)==1
                        %Hit
                        cum_op_bin_hit(binary_op+1)=cum_op_bin_hit(binary_op+1)+1;
                        cum_op_bin_bindFF_hit(this_bin_dFF_rev,binary_op+1)=...
                            cum_op_bin_bindFF_hit(this_bin_dFF_rev,binary_op+1)+1;
                    else
                        %Miss
                        cum_op_bin_miss(binary_op+1)=cum_op_bin_miss(binary_op+1)+1;
                        cum_op_bin_bindFF_miss(this_bin_dFF_rev,binary_op+1)=...
                            cum_op_bin_bindFF_miss(this_bin_dFF_rev,binary_op+1)+1;
                    end

                    %Find the analog op
                    minC=odor_plume_patterns.cm_from_floor(cm_from_floor).minC;
                    this_x=trials.trial(trNo).XYtest(this_ii_t,1);
                    this_y=trials.trial(trNo).XYtest(this_ii_t,2);

                    [M,this_x_op_ii]=min(abs(odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume-this_x));
                    [M,this_y_op_ii]=min(abs(odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume-this_y));

                    if trials.lane_per_trial(trNo)==1
                        op_conc=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1(this_y_op_ii,this_x_op_ii);
                    else
                        op_conc=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4(this_y_op_ii,this_x_op_ii);
                    end

                    %Make odor minC if the time is early or late
                    if x_ii<trials.odor_ii_start(trNo)
                        op_conc=minC;
                    else
                        if x_ii>trials.odor_ii_end(trNo)
                            op_conc=minC;
                        else
                            this_x=trials.trial(trNo).XYtest(this_ii_t,1);
                            x_on=(x_ii-trials.odor_ii_start(trNo))*handles_conc.dt*handles_conc.air_flow_speed;
                            if this_x>x_on
                                op_conc=minC;
                            end
                        end
                    end
                    for ii_b=1:10
                        if (op_conc>=odor_c_bounds(ii_b))&(op_conc<odor_c_bounds(ii_b+1))
                            this_ii_op=ii_b;
                        end
                    end

                    %Now do cums for op analogue
                    cum_op(this_ii_op)=cum_op(this_ii_op)+1;
                    cum_op_bindFF(this_bin_dFF_rev,this_ii_op)=cum_op_bindFF(this_bin_dFF_rev,this_ii_op)+1;
                    include_op(this_ii_op)=1;

                    if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                        %Hit
                        cum_op_hit(binary_op+1)=cum_op_hit(binary_op+1)+1;
                        cum_op_bindFF_hit(this_bin_dFF_rev,binary_op+1)=...
                            cum_op_bindFF_hit(this_bin_dFF_rev,binary_op+1)+1;
                        include_op_hit(this_ii_op)=1;
                    else
                        %Miss
                        cum_op_miss(binary_op+1)=cum_op_miss(binary_op+1)+1;
                        cum_op_bindFF_miss(this_bin_dFF_rev,binary_op+1)=...
                            cum_op_bindFF_miss(this_bin_dFF_rev,binary_op+1)+1;
                        include_op_miss(this_ii_op)=1;
                    end

                end


                %Calculate information

                p_xy=cum_xy(logical(include_xy))/sum(cum_xy(logical(include_xy)));
                p_bindFF=cum_bindFF/sum(cum_bindFF(:));
                included_cum_xy_bindFF=cum_xy_bindFF(:,logical(include_xy));
                p_xy_bindFF=included_cum_xy_bindFF/sum(included_cum_xy_bindFF(:));

                p_xy_hit=cum_xy_hit(logical(include_xy_hit))/sum(cum_xy_hit(logical(include_xy_hit)));
                p_bindFF_hit=cum_bindFF_hit/sum(cum_bindFF_hit(:));
                included_cum_xy_bindFF_hit=cum_xy_bindFF_hit(:,logical(include_xy_hit));
                p_xy_bindFF_hit=included_cum_xy_bindFF_hit/sum(included_cum_xy_bindFF_hit(:));

                p_xy_miss=cum_xy_miss(logical(include_xy_miss))/sum(cum_xy_miss(logical(include_xy_miss)));
                p_bindFF_miss=cum_bindFF_miss/sum(cum_bindFF_miss(:));
                included_cum_xy_bindFF_miss=cum_xy_bindFF_miss(:,logical(include_xy_miss));
                p_xy_bindFF_miss=included_cum_xy_bindFF_miss/sum(included_cum_xy_bindFF_miss(:));


                %Binary op
                p_op_bin=cum_op_bin/sum(cum_op_bin);
                p_op_bin_bindFF=cum_op_bin_bindFF/sum(cum_op_bin_bindFF(:));

                p_op_bin_hit=cum_op_bin/sum(cum_op_bin_hit);
                p_op_bin_bindFF_hit=cum_op_bin_bindFF_hit/sum(cum_op_bin_bindFF_hit(:));

                p_op_bin_miss=cum_op_bin/sum(cum_op_bin_miss);
                p_op_bin_bindFF_miss=cum_op_bin_bindFF_miss/sum(cum_op_bin_bindFF_miss(:));

                %Analog op
                included_cum_op=cum_op(logical(include_op));
                p_op=included_cum_op/sum(included_cum_op);
                included_cum_op_bindFF=cum_op_bindFF(:,logical(include_op));
                p_op_bindFF=included_cum_op_bindFF/sum(included_cum_op_bindFF(:));

                included_cum_op_hit=cum_op_hit(logical(include_op_hit));
                p_op_hit=included_cum_op_hit/sum(included_cum_op_hit);
                included_cum_op_bindFF_hit=cum_op_bindFF_hit(:,logical(include_op_hit));
                p_op_bindFF_hit=included_cum_op_bindFF_hit/sum(included_cum_op_bindFF_hit(:));

                included_cum_op_miss=cum_op_miss(logical(include_op_miss));
                p_op_miss=included_cum_op_miss/sum(included_cum_op_miss);
                included_cum_op_bindFF_miss=cum_op_bindFF_miss(:,logical(include_op_miss));
                p_op_bindFF_miss=included_cum_op_bindFF_miss/sum(included_cum_op_bindFF_miss(:));


                all_info_xy_hit_sh(all_info_ii,ii_sh)=0;
                all_info_xy_miss_sh(all_info_ii,ii_sh)=0;
                all_info_xy_sh(all_info_ii,ii_sh)=0;

                all_op_bin_info_sh(all_info_ii,ii_sh)=0;
                all_op_bin_info_hit_sh(all_info_ii,ii_sh)=0;
                all_op_bin_info_miss_sh(all_info_ii,ii_sh)=0;

                all_info_op_sh(all_info_ii,ii_sh)=0;
                all_info_op_hit_sh(all_info_ii,ii_sh)=0;
                all_info_op_miss_sh(all_info_ii,ii_sh)=0;

                %Calculate info for hit
                for ii_xy=1:size(p_xy_bindFF_hit,2)
                    % if sum(p_xy_bindFF_hit(:,ii_xy))>0
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF_hit(ii_bin_dFF,ii_xy)~=0
                            all_info_xy_hit_sh(all_info_ii,ii_sh)=all_info_xy_hit_sh(all_info_ii,ii_sh)+p_xy_bindFF_hit(ii_bin_dFF,ii_xy)*...
                                log2(p_xy_bindFF_hit(ii_bin_dFF,ii_xy)/(p_xy_hit(ii_xy)*p_bindFF_hit(ii_bin_dFF)));
                        end
                    end
                    % end
                end

                %Calculate info for miss
                for ii_xy=1:size(p_xy_bindFF_miss,2)
                    % if sum(p_xy_bindFF_miss(:,ii_xy))>0
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF_miss(ii_bin_dFF,ii_xy)~=0
                            all_info_xy_miss_sh(all_info_ii,ii_sh)=all_info_xy_miss_sh(all_info_ii,ii_sh)+p_xy_bindFF_miss(ii_bin_dFF,ii_xy)*...
                                log2(p_xy_bindFF_miss(ii_bin_dFF,ii_xy)/(p_xy_miss(ii_xy)*p_bindFF_miss(ii_bin_dFF)));
                        end
                    end
                    % end
                end

                %Calculate info for both hit and miss
                for ii_xy=1:size(p_xy_bindFF,2)
                    % if sum(p_xy_bindFF_miss(:,ii_xy))>0
                    for ii_bin_dFF=1:2
                        if p_xy_bindFF(ii_bin_dFF,ii_xy)~=0
                            all_info_xy_sh(all_info_ii,ii_sh)=all_info_xy_sh(all_info_ii,ii_sh)+p_xy_bindFF(ii_bin_dFF,ii_xy)*...
                                log2(p_xy_bindFF(ii_bin_dFF,ii_xy)/(p_xy(ii_xy)*p_bindFF(ii_bin_dFF)));
                        end
                    end
                    % end
                end

                %Binary op

                %Calculate info for op_bin
                for ii_op=1:2
                    for ii_bin_dFF=1:2
                        if p_op_bin_bindFF_hit(ii_bin_dFF,ii_op)~=0
                            all_op_bin_info_hit_sh(all_info_ii,ii_sh)=all_op_bin_info_hit_sh(all_info_ii,ii_sh)+p_op_bin_bindFF_hit(ii_bin_dFF,ii_op)*...
                                log2(p_op_bin_bindFF_miss(ii_bin_dFF,ii_op)/(p_op_bin_miss(ii_op)*p_bindFF_miss(ii_bin_dFF)));
                        end
                    end
                end

                for ii_op=1:2
                    for ii_bin_dFF=1:2
                        if p_op_bin_bindFF_miss(ii_bin_dFF,ii_op)~=0
                            all_op_bin_info_miss_sh(all_info_ii,ii_sh)=all_op_bin_info_miss_sh(all_info_ii,ii_sh)+p_op_bin_bindFF_miss(ii_bin_dFF,ii_op)*...
                                log2(p_op_bin_bindFF_miss(ii_bin_dFF,ii_op)/(p_op_bin_miss(ii_op)*p_bindFF_miss(ii_bin_dFF)));
                        end
                    end
                end

                for ii_op=1:2
                    for ii_bin_dFF=1:2
                        if p_op_bin_bindFF(ii_bin_dFF,ii_op)~=0
                            all_op_bin_info_sh(all_info_ii,ii_sh)=all_op_bin_info_sh(all_info_ii,ii_sh)+p_op_bin_bindFF(ii_bin_dFF,ii_op)*...
                                log2(p_op_bin_bindFF(ii_bin_dFF,ii_op)/(p_op_bin(ii_op)*p_bindFF(ii_bin_dFF)));
                        end
                    end
                end


                %Calculate info for op analog
                for ii_op=1:size(p_op_bindFF,2)
                    for ii_bin_dFF=1:2
                        if p_op_bindFF(ii_bin_dFF,ii_op)~=0
                            all_info_op_sh(all_info_ii,ii_sh)=all_info_op_sh(all_info_ii,ii_sh)+p_op_bindFF(ii_bin_dFF,ii_op)*...
                                log2(p_op_bindFF(ii_bin_dFF,ii_op)/(p_op(ii_op)*p_bindFF(ii_bin_dFF)));
                        end
                    end
                end

                %Calculate info for op analog hit
                for ii_op=1:size(p_op_bindFF_hit,2)
                    for ii_bin_dFF=1:2
                        if p_op_bindFF_hit(ii_bin_dFF,ii_op)~=0
                            all_info_op_hit_sh(all_info_ii,ii_sh)=all_info_op_hit_sh(all_info_ii,ii_sh)+p_op_bindFF_hit(ii_bin_dFF,ii_op)*...
                                log2(p_op_bindFF_hit(ii_bin_dFF,ii_op)/(p_op_hit(ii_op)*p_bindFF_hit(ii_bin_dFF)));
                        end
                    end
                end

                %Calculate info for op analog miss
                for ii_op=1:size(p_op_bindFF_miss,2)
                    for ii_bin_dFF=1:2
                        if p_op_bindFF_miss(ii_bin_dFF,ii_op)~=0
                            all_info_op_miss_sh(all_info_ii,ii_sh)=all_info_op_miss_sh(all_info_ii,ii_sh)+p_op_bindFF_miss(ii_bin_dFF,ii_op)*...
                                log2(p_op_bindFF_miss(ii_bin_dFF,ii_op)/(p_op_miss(ii_op)*p_bindFF_miss(ii_bin_dFF)));
                        end
                    end
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



%Now calculate 3xSD for the shuffled distributions
sd_all_spatial_info_hit_sh=zeros(all_info_ii,1);
sd_all_spatial_info_miss_sh=zeros(all_info_ii,1);
sd_all_spatial_info=zeros(all_info_ii,1);

sd_all_op_bin_info_hit_sh=zeros(all_info_ii,1);
sd_all_op_bin_info_miss_sh=zeros(all_info_ii,1);
sd_all_op_bin_info_sh=zeros(all_info_ii,1);

sd_all_op_info_hit_sh=zeros(all_info_ii,1);
sd_all_op_info_miss_sh=zeros(all_info_ii,1);
sd_all_op_info_sh=zeros(all_info_ii,1);

mean_all_spatial_info_hit_sh=zeros(all_info_ii,1);
mean_all_spatial_info_miss_sh=zeros(all_info_ii,1);
mean_all_spatial_info_sh=zeros(all_info_ii,1);

mean_all_op_bin_info_sh=zeros(all_info_ii,1);
mean_all_op_bin_info_hit_sh=zeros(all_info_ii,1);
mean_all_op_bin_info_miss_sh=zeros(all_info_ii,1);

mean_all_op_info_sh=zeros(all_info_ii,1);
mean_all_op_info_hit_sh=zeros(all_info_ii,1);
mean_all_op_info_miss_sh=zeros(all_info_ii,1);

all_all_spatial_info_hit_sh=[];
all_all_spatial_info_miss_sh=[];
all_all_spatial_info_sh=[];

all_all_op_bin_info_sh=[];
all_all_op_bin_info_hit_sh=[];
all_all_op_bin_info_miss_sh=[];

all_all_op_info_sh=[];
all_all_op_info_hit_sh=[];
all_all_op_info_miss_sh=[];

ssi_all_spatial_info_hit=zeros(all_info_ii,1);
ssi_all_spatial_info_miss=zeros(all_info_ii,1);
ssi_all_spatial_info=zeros(all_info_ii,1);

ssi_all_op_bin_info=zeros(all_info_ii,1);
ssi_all_op_bin_info_hit=zeros(all_info_ii,1);
ssi_all_op_bin_info_miss=zeros(all_info_ii,1);



for ii_ROI=1:all_info_ii

    this_ii_ROI=ii_ROI

    this_all_spatial_info_hit_sh=zeros(1,n_shuffle_SI);
    this_all_spatial_info_hit_sh(1,:)=all_info_xy_hit_sh(ii_ROI,:);
    this_all_spatial_info_miss_sh=zeros(1,n_shuffle_SI);
    this_all_spatial_info_miss_sh(1,:)=all_info_xy_miss_sh(ii_ROI,:);
    this_all_spatial_info_sh=zeros(1,n_shuffle_SI);
    this_all_spatial_info_sh(1,:)=all_info_xy_sh(ii_ROI,:);

    this_op_bin_info_hit_sh=zeros(1,n_shuffle_SI);
    this_op_bin_info_hit_sh(1,:)=all_op_bin_info_hit_sh(ii_ROI,:);
    this_op_bin_info_miss_sh=zeros(1,n_shuffle_SI);
    this_op_bin_info_miss_sh(1,:)=all_op_bin_info_miss_sh(ii_ROI,:);
    this_op_bin_info_sh=zeros(1,n_shuffle_SI);
    this_op_bin_info_sh(1,:)=all_op_bin_info_sh(ii_ROI,:);

    this_op_info_hit_sh=zeros(1,n_shuffle_SI);
    this_op_info_hit_sh(1,:)=all_info_op_hit_sh(ii_ROI,:);
    this_op_info_miss_sh=zeros(1,n_shuffle_SI);
    this_op_info_miss_sh(1,:)=all_info_op_miss_sh(ii_ROI,:);
    this_op_info_sh=zeros(1,n_shuffle_SI);
    this_op_info_sh(1,:)=all_info_op_sh(ii_ROI,:);

    all_all_spatial_info_hit_sh=[all_all_spatial_info_hit_sh this_all_spatial_info_hit_sh];
    all_all_spatial_info_miss_sh=[all_all_spatial_info_miss_sh this_all_spatial_info_miss_sh];
    all_all_spatial_info_sh=[all_all_spatial_info_sh this_all_spatial_info_sh];

    all_all_op_bin_info_hit_sh=[all_all_op_bin_info_hit_sh this_op_info_hit_sh];
    all_all_op_bin_info_miss_sh=[all_all_op_bin_info_miss_sh this_op_info_miss_sh];
    all_all_op_bin_info_sh=[all_all_op_bin_info_sh this_op_info_sh];

    all_all_op_info_hit_sh=[all_all_op_info_hit_sh this_op_info_hit_sh];
    all_all_op_info_miss_sh=[all_all_op_info_miss_sh this_op_info_miss_sh];
    all_all_op_info_sh=[all_all_op_info_sh this_op_info_sh];


    for jj=1:9
        switch jj
            case 1
                data=this_all_spatial_info_hit_sh;
            case 2
                data=this_all_spatial_info_miss_sh;
            case 3
                data=this_all_spatial_info_sh;
            case 4
                data=this_op_bin_info_hit_sh;
            case 5
                data=this_op_bin_info_miss_sh;
            case 6
                data=this_op_bin_info_sh;
            case 7
                data=this_op_info_hit_sh;
            case 8
                data=this_op_info_miss_sh;
            case 9
                data=this_op_info_sh;
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
                sd_all_spatial_info_hit_sh(ii_ROI)=estimated_std;
            case 2
                sd_all_spatial_info_miss_sh(ii_ROI)=estimated_std;
            case 3
                sd_all_spatial_info(ii_ROI)=estimated_std;
            case 4
                sd_all_op_bin_info_hit_sh(ii_ROI)=estimated_std;
            case 5
                sd_all_op_bin_info_miss_sh(ii_ROI)=estimated_std;
            case 6
                sd_all_op_bin_info_sh(ii_ROI)=estimated_std;
            case 7
                sd_all_op_info_hit_sh(ii_ROI)=estimated_std;
            case 8
                sd_all_op_info_miss_sh(ii_ROI)=estimated_std;
            case 9
                sd_all_op_info_sh(ii_ROI)=estimated_std;
        end
    end

    mean_all_spatial_info_hit_sh(ii_ROI)=mean(this_all_spatial_info_hit_sh);
    mean_all_spatial_info_miss_sh(ii_ROI)=mean(this_all_spatial_info_miss_sh);
    mean_all_spatial_info_sh(ii_ROI)=mean(this_all_spatial_info_sh);

    mean_all_op_bin_info_sh(ii_ROI)=mean(this_op_bin_info_sh);
    mean_all_op_bin_info_hit_sh(ii_ROI)=mean(this_op_bin_info_hit_sh);
    mean_all_op_bin_info_miss_sh(ii_ROI)=mean(this_op_bin_info_miss_sh);

    mean_all_op_info_sh(ii_ROI)=mean(this_op_info_sh);
    mean_all_op_info_hit_sh(ii_ROI)=mean(this_op_info_hit_sh);
    mean_all_op_info_miss_sh(ii_ROI)=mean(this_op_info_miss_sh);

    ssi_all_spatial_info_hit(ii_ROI)=(all_info_xy_hit(ii_ROI)-mean_all_spatial_info_hit_sh(ii_ROI))/sd_all_spatial_info_hit_sh(ii_ROI);
    ssi_all_spatial_info_miss(ii_ROI)=(all_info_xy_miss(ii_ROI)-mean_all_spatial_info_miss_sh(ii_ROI))/sd_all_spatial_info_miss_sh(ii_ROI);
    ssi_all_spatial_info(ii_ROI)=(all_info_xy(ii_ROI)-mean_all_spatial_info_sh(ii_ROI))/sd_all_spatial_info(ii_ROI);

    ssi_all_op_bin_info(ii_ROI)=(all_op_bin_info(ii_ROI)-mean_all_op_bin_info_sh(ii_ROI))/sd_all_op_bin_info_sh(ii_ROI);
    ssi_all_op_bin_info_hit(ii_ROI)=(all_op_bin_info_hit(ii_ROI)-mean_all_op_bin_info_hit_sh(ii_ROI))/sd_all_op_bin_info_hit_sh(ii_ROI);
    ssi_all_op_bin_info_miss(ii_ROI)=(all_op_bin_info_miss(ii_ROI)-mean_all_op_bin_info_miss_sh(ii_ROI))/sd_all_op_bin_info_miss_sh(ii_ROI);

    ssi_all_op_info(ii_ROI)=(all_info_op(ii_ROI)-mean_all_op_info_sh(ii_ROI))/sd_all_op_info_sh(ii_ROI);
    ssi_all_op_info_hit(ii_ROI)=(all_info_op_hit(ii_ROI)-mean_all_op_info_hit_sh(ii_ROI))/sd_all_op_info_hit_sh(ii_ROI);
    ssi_all_op_info_miss(ii_ROI)=(all_info_op_miss(ii_ROI)-mean_all_op_info_miss_sh(ii_ROI))/sd_all_op_info_miss_sh(ii_ROI);

    handles_outic.ssi_all_spatial_info_hit(ii_ROI)=ssi_all_spatial_info_hit(ii_ROI);
    handles_outic.ssi_all_spatial_info_miss(ii_ROI)=ssi_all_spatial_info_miss(ii_ROI);
    handles_outic.ssi_all_spatial_info(ii_ROI)=ssi_all_spatial_info(ii_ROI);

    handles_outic.ssi_all_op_bin_info(ii_ROI)=ssi_all_op_bin_info(ii_ROI);
    handles_outic.ssi_all_op_bin_info_hit(ii_ROI)=ssi_all_op_bin_info_hit(ii_ROI);
    handles_outic.ssi_all_op_bin_info_hit(ii_ROI)=ssi_all_op_bin_info_hit(ii_ROI);


    handles_outic.ssi_all_op_info(ii_ROI)=ssi_all_op_info(ii_ROI);
    handles_outic.ssi_all_op_info_hit(ii_ROI)=ssi_all_op_info_hit(ii_ROI);
    handles_outic.ssi_all_op_info_hit(ii_ROI)=ssi_all_op_info_hit(ii_ROI);

end


%Plot histogram for SSI mutual info for both lanes for xy x bindFF
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

edges=[-4:2:34];
histogram(ssi_all_spatial_info)

plot([3 3],[0 450],'-k','LineWidth',2)

xlim([-5 35])
title('Histogram for SSI MI position x bindFF')
xlabel('Sigma')
ylabel('Count')


%Plot histogram for SSI mutual info for op x bindFF
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

edges=[-4:1:20];
histogram(ssi_all_op_bin_info,edges)

plot([3 3],[0 900],'-k','LineWidth',2)

xlim([-5 20])
title('Histogram for SSI MI binary odor concentration x bindFF')
xlabel('Sigma')
ylabel('Count')


%Plot histogram for SSI mutual info for op x bindFF
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

edges=[-4:1:20];
histogram(ssi_all_op_info,edges)

plot([3 3],[0 900],'-k','LineWidth',2)

xlim([-5 20])
title('Histogram for SSI MI odor concentration x bindFF')
xlabel('Sigma')
ylabel('Count')

%Plot histograms for SSI mutual info for xy x bindFF
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

[f_aic,x_aic] = drg_ecdf(ssi_all_spatial_info_hit);
plot(x_aic,f_aic,'Color',[204/255 121/255 167/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(ssi_all_spatial_info_miss); %xy x bindFF for lane 4
plot(x_aic,f_aic,'Color',[0/255 114/255 178/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(ssi_all_spatial_info);   %xy x bindFF x lane
plot(x_aic,f_aic,'Color',[213/255 94/255 0/255],'LineWidth',3)


plot([3 3],[0 1],'-k','LineWidth',2)

these_ylim=ylim;
these_xlim=xlim;
text(4,0.4,'Hit','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
text(4,0.32,'Miss','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
text(4,0.24,'All trials','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.25,'MI lane x bindFF','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.20,'MI op x bindFF','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.15,'MI op x bindFF x xy','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.10,'MI op x xy','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)

xlim([-2 15])
title('z-scored mutual information (position and binary dF/F)')
xlabel('Sigma')
ylabel('Cumulative fraction')



%Now do tsne analysis with SSI for place cells
% ssi_all_spatial_info_hit_no_nan=ssi_all_spatial_info_hit((~isnan(ssi_all_spatial_info_hit))&(~isnan(ssi_all_spatial_info_miss))&(~isnan(ssi_all_spatial_info)));
% ssi_all_spatial_info_miss_no_nan=ssi_all_spatial_info_miss((~isnan(ssi_all_spatial_info_hit))&(~isnan(ssi_all_spatial_info_miss))&(~isnan(ssi_all_spatial_info)));
% ssi_all_spatial_info_no_nan=ssi_all_spatial_info((~isnan(ssi_all_spatial_info_hit))&(~isnan(ssi_all_spatial_info_miss))&(~isnan(ssi_all_spatial_info)));

datassi = [ssi_all_spatial_info_hit'; ssi_all_spatial_info_miss'; ssi_all_spatial_info'];

%Load correlations, etc
load([save_PathPredImp save_FilePredImp])

ii_all_ROIs=1:length(ssi_all_spatial_info_hit);
not_nan_ii_allROIs=ii_all_ROIs((~isnan(ssi_all_spatial_info_hit))&(~isnan(ssi_all_spatial_info_miss))&(~isnan(ssi_all_spatial_info)));
handles_outic.not_nan_ii_allROIs=not_nan_ii_allROIs;

cropped_ssi_all_spatial_info_hit=zeros(length(not_nan_ii_allROIs),1);
cropped_ssi_all_spatial_info_miss=zeros(length(not_nan_ii_allROIs),1);
cropped_ssi_all_spatial_info=zeros(length(not_nan_ii_allROIs),1);
cropped_all_spatial_rho_hit_vs_miss=zeros(length(not_nan_ii_allROIs),1);
cropped_all_delta_center_of_mass=zeros(length(not_nan_ii_allROIs),1);

cropped_ssi_all_op_info=zeros(length(not_nan_ii_allROIs),1);
cropped_ssi_all_op_info_hit=zeros(length(not_nan_ii_allROIs),1);
cropped_ssi_all_op_info_miss=zeros(length(not_nan_ii_allROIs),1);


for jj=1:length(not_nan_ii_allROIs)
    cropped_ssi_all_spatial_info_hit(jj)=ssi_all_spatial_info_hit(not_nan_ii_allROIs(jj));
    cropped_ssi_all_spatial_info_miss(jj)=ssi_all_spatial_info_miss(not_nan_ii_allROIs(jj));
    cropped_ssi_all_spatial_info(jj)=ssi_all_spatial_info(not_nan_ii_allROIs(jj));
    cropped_all_spatial_rho_hit_vs_miss(jj)=imps.all_spatial_rho_hit_vs_miss(not_nan_ii_allROIs(jj));
    cropped_all_delta_center_of_mass(jj)=imps.all_delta_center_of_mass(not_nan_ii_allROIs(jj));
    cropped_ssi_all_op_info(jj)=ssi_all_op_info(not_nan_ii_allROIs(jj));
    cropped_ssi_all_op_info_hit(jj)=ssi_all_op_info_hit(not_nan_ii_allROIs(jj));
    cropped_ssi_all_op_info_miss(jj)=ssi_all_op_info_miss(not_nan_ii_allROIs(jj));
end

handles_outic.cropped_ssi_all_spatial_info_hit=cropped_ssi_all_spatial_info_hit;
handles_outic.cropped_ssi_all_spatial_info_miss=cropped_ssi_all_spatial_info_miss;
handles_outic.cropped_ssi_all_spatial_info=cropped_ssi_all_spatial_info;
handles_outic.cropped_all_spatial_rho_hit_vs_miss=cropped_all_spatial_rho_hit_vs_miss;
handles_outic.cropped_all_delta_center_of_mass=cropped_all_delta_center_of_mass;



% Step 1 Transpose the data if necessary (N samples x 3 variables)
datassi = datassi';

rng('default') % for reproducibility

% Step 2: Run t-SNE
% The default parameters are usually sufficient, but you can adjust them
% For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% mappedX = tsne(datassi, 'NumPCAComponents', 2, 'Perplexity', 30);
% mappedX = tsne(datassi,'Algorithm','exact','Distance','mahalanobis'); %Works well
% mappedX = tsne(datassi,'Algorithm','exact','Distance','cosine'); %works better
mappedX = tsne(datassi,'Distance','cosine'); %works best

% mappedX = tsne(datassi,'Algorithm','exact','Distance','chebychev'); %Ok, yields two clusters
% mappedX = tsne(datassi,'Algorithm','exact','Distance','euclidean'); %OK, three clusters

% %Step 3: Plot the results
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

no_clusters=3;
gm = fitgmdist(mappedX, no_clusters); % Assuming 2 clusters
idxssi = cluster(gm, mappedX);

handles_outic.idxssi=idxssi;

not_nan_ii_ROI=zeros(1,length(ssi_all_op_info));
ii_not_nan=0;
for ii=1:length(ssi_all_op_info)
    if ~isnan(ssi_all_op_info(ii))
        ii_not_nan=ii_not_nan+1;
        not_nan_ii_ROI(ii)=ii_not_nan;
    else
        not_nan_ii_ROI(ii)=NaN;
    end
end

cropped_ssi_all_op_info=zeros(1,sum(~isnan(ssi_all_op_info)));
cropped_idxssi_all_op_info=zeros(1,sum(~isnan(ssi_all_op_info)));
jj_all=0;
cropped_ssi_all_op_info_hit=zeros(1,sum(~isnan(ssi_all_op_info_hit)));
cropped_idxssi_all_op_info_hit=zeros(1,sum(~isnan(ssi_all_op_info_hit)));
jj_hit=0;
cropped_ssi_all_op_info_miss=zeros(1,sum(~isnan(ssi_all_op_info_miss)));
cropped_idxssi_all_op_info_miss=zeros(1,sum(~isnan(ssi_all_op_info_miss)));
jj_miss=0;
for jj=1:length(ssi_all_op_info)
    if ~isnan(ssi_all_op_info(jj))
        jj_all=jj_all+1;
        cropped_ssi_all_op_info(jj_all)=ssi_all_op_info(jj);
        kk=find(jj==not_nan_ii_allROIs);
        cropped_idxssi_all_op_info(jj_all)=idxssi(kk);
    end
    if ~isnan(ssi_all_op_info_hit(jj))
        jj_hit=jj_hit+1;
        cropped_ssi_all_op_info_hit(jj_hit)=ssi_all_op_info_hit(jj);
        kk=find(jj==not_nan_ii_allROIs);
        cropped_idxssi_all_op_info_hit(jj_hit)=idxssi(kk);
    end
    if ~isnan(ssi_all_op_info_miss(jj))
        jj_miss=jj_miss+1;
        cropped_ssi_all_op_info_miss(jj_miss)=ssi_all_op_info_miss(jj);
        kk=find(jj==not_nan_ii_allROIs);
        cropped_idxssi_all_op_info_miss(jj_miss)=idxssi(kk);
    end
end

handles_outic.cropped_ssi_all_op_info=cropped_ssi_all_op_info;
handles_outic.cropped_ssi_all_op_info_hit=cropped_ssi_all_op_info_hit;
handles_outic.cropped_ssi_all_op_info_miss=cropped_ssi_all_op_info_miss;

%Report the mean infos for each cluster
per_cluster=[];
for ii_k=1:no_clusters
    these_ssi_all_spatial_info_hit=cropped_ssi_all_spatial_info_hit(idxssi==ii_k);
    these_ssi_all_spatial_info_miss=cropped_ssi_all_spatial_info_miss(idxssi==ii_k);
    these_ssi_all_spatial_info=cropped_ssi_all_spatial_info(idxssi==ii_k);
    mean_hi=mean(these_ssi_all_spatial_info_hit(~isnan(these_ssi_all_spatial_info_hit)));
    per_cluster.cluster(ii_k).mean_hi=mean_hi;
    per_cluster.cluster(ii_k).these_ssi_all_spatial_info_hit=these_ssi_all_spatial_info_hit(~isnan(these_ssi_all_spatial_info_hit));
    mean_miss=mean(these_ssi_all_spatial_info_miss(~isnan(these_ssi_all_spatial_info_miss)));
    per_cluster.cluster(ii_k).mean_miss=mean_miss;
    per_cluster.cluster(ii_k).these_ssi_all_spatial_info_miss=these_ssi_all_spatial_info_miss(~isnan(these_ssi_all_spatial_info_miss));
    mean_both=mean(these_ssi_all_spatial_info(~isnan(these_ssi_all_spatial_info)));
    per_cluster.cluster(ii_k).mean_both=mean_both;
    per_cluster.cluster(ii_k).these_ssi_all_spatial_info=these_ssi_all_spatial_info(~isnan(these_ssi_all_spatial_info));

    fprintf(1, ['\nMean SSI for cluster ' num2str(ii_k) ' lane1, lane 4, both: '...
        num2str(mean_hi) ' ' num2str(mean_miss) ' ' num2str(mean_both) '\n'])
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
for ii_k=1:no_clusters
    switch ii_k
        case 1
            plot(mappedX(idxssi==ii_k,1),mappedX(idxssi==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        case 2
            plot(mappedX(idxssi==ii_k,1),mappedX(idxssi==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        case 3
            plot(mappedX(idxssi==ii_k,1),mappedX(idxssi==ii_k,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
            % case 4
            % plot(mappedX(idxssi==ii_k,1),mappedX(idxssi==ii_k,2),'.','MarkerFaceColor',[204/255 121/255 167/255],'MarkerEdgeColor',[204/255 121/255 167/255]);
    end
end

text(-37,30,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(-37,20,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(-37,10,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE for SSI for position and dFF for each lane vs both lanes' ]);

%Plot cumulative histograms for SSIs
for ii_k=1:no_clusters
    %SSIs
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    these_ssis=per_cluster.cluster(ii_k).these_ssi_all_spatial_info_hit;
    [f_aic,x_aic] = drg_ecdf(these_ssis);
    plot(x_aic,f_aic,'Color',[204/255 121/255 167/255],'LineWidth',3)

    these_ssis=per_cluster.cluster(ii_k).these_ssi_all_spatial_info_miss;
    [f_aic,x_aic] = drg_ecdf(these_ssis); %xy x bindFF for lane 4
    plot(x_aic,f_aic,'Color',[0/255 114/255 178/255],'LineWidth',3)

    % these_ssis=per_cluster.cluster(ii_k).these_ssi_all_spatial_info;
    % [f_aic,x_aic] = drg_ecdf(these_ssis);   %xy x bindFF x lane
    % plot(x_aic,f_aic,'Color',[213/255 94/255 0/255],'LineWidth',3)

    plot([3 3],[0 1],'-k','LineWidth',2)

    title(['SSIs for Cluster ' num2str(ii_k)])
    ylabel('Cumulative probability')
    xlabel('SSI')

    text(5,0.4,'Hit','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
    text(5,0.3,'Miss','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
    % text(10,0.2,'SSI MI xy x bindFF both lanes','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
    xlim([-5 15])
end



%Plot cumulative histograms for SSI both lanes
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
for ii_k=1:no_clusters
    %SSIs

    switch ii_k
        case 1
            these_ssis=per_cluster.cluster(ii_k).these_ssi_all_spatial_info;
            [f_aic,x_aic] = drg_ecdf(these_ssis);
            plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)
        case 2

            these_ssis=per_cluster.cluster(ii_k).these_ssi_all_spatial_info;
            [f_aic,x_aic] = drg_ecdf(these_ssis);
            plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)
        case 3
            these_ssis=per_cluster.cluster(ii_k).these_ssi_all_spatial_info;
            [f_aic,x_aic] = drg_ecdf(these_ssis);
            plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)
    end


end


text(10,0.4,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(10,0.3,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(10,0.2,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

title(['SSI for all trials'])
ylabel('Cumulative probability')
xlabel('SSI')
xlim([-5 35])

%Difference in SSIs
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
for ii_k=1:no_clusters

    these_ssis_hit=per_cluster.cluster(ii_k).these_ssi_all_spatial_info_hit;
    these_ssis_miss=per_cluster.cluster(ii_k).these_ssi_all_spatial_info_miss;

    switch ii_k
        case 1
            [f_aic,x_aic] = drg_ecdf(these_ssis_hit-these_ssis_miss);
            plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)
        case 2
            [f_aic,x_aic] = drg_ecdf(these_ssis_hit-these_ssis_miss); %xy x bindFF for lane 4
            plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)
        case 3
            [f_aic,x_aic] = drg_ecdf(these_ssis_hit-these_ssis_miss);   %xy x bindFF x lane
            plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)
    end


end
plot([0 0],[0 1],'-k')
title(['SSI hit -SSI miss'])
ylabel('Cumulative probability')
xlabel('Delta SSI')
xlim([-15 15])

text(5,0.4,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(5,0.3,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(5,0.2,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)





%Correlations
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
for ii_k=1:no_clusters



    these_rhos=cropped_all_spatial_rho_hit_vs_miss(idxssi==ii_k);
    these_rhos=these_rhos(~isnan(these_rhos));

    switch ii_k
        case 1
            [f_aic,x_aic] = drg_ecdf(these_rhos);
            plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)
        case 2
            [f_aic,x_aic] = drg_ecdf(these_rhos); %xy x bindFF for lane 4
            plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)
        case 3
            [f_aic,x_aic] = drg_ecdf(these_rhos);   %xy x bindFF x lane
            plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)
    end


end
plot([0 0],[0 1],'-k')
title(['rho hit vs miss'])
ylabel('Cumulative probability')
xlabel('rho')
xlim([-0.5 1])

text(0.1,0.4,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(0.1,0.3,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(0.1,0.2,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

%Center of mass
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
for ii_k=1:no_clusters

    these_dcoms=cropped_all_delta_center_of_mass(idxssi==ii_k);
    these_dcoms=these_dcoms(~isnan(these_dcoms));
    switch ii_k
        case 1
            [f_aic,x_aic] = drg_ecdf(these_dcoms);
            plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)
        case 2
            [f_aic,x_aic] = drg_ecdf(these_dcoms); %xy x bindFF for lane 4
            plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)
        case 3
            [f_aic,x_aic] = drg_ecdf(these_dcoms);   %xy x bindFF x lane
            plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)
    end


end
plot([0 0],[0 1],'-k')
title(['Distance of center of mass'])
ylabel('Cumulative probability')
xlabel('Distance (mm)')
xlim([0 400])

text(100,0.4,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(100,0.3,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(100,0.2,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

%Plot histograms for spatial correlation and distance for the center of
%mass

%Spatial correlation
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
these_rhos=imps.all_spatial_rho_hit_vs_miss;
these_rhos=these_rhos(~isnan(these_rhos));

edges=[-0.5:0.1:1];
histogram(these_rhos,edges)
title('Histogram for spatial correlation')
xlabel('Spatial correlation')
ylabel('Count')


%Distance center of mass
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
these_dcoms=imps.all_delta_center_of_mass;
these_dcoms=these_dcoms(~isnan(these_dcoms));

edges=[0:20:500];
histogram(these_dcoms,edges)
title('Histogram for distance between center of mass')
xlabel('Distance (mm)')
ylabel('Count')


%Plot histograms for SSI mutual info for op x bindFF
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

[f_aic,x_aic] = drg_ecdf(ssi_all_op_info_hit);
plot(x_aic,f_aic,'Color',[204/255 121/255 167/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(ssi_all_op_info_miss); 
plot(x_aic,f_aic,'Color',[0/255 114/255 178/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(ssi_all_op_info);   
plot(x_aic,f_aic,'Color',[213/255 94/255 0/255],'LineWidth',3)


plot([3 3],[0 1],'-k','LineWidth',2)

these_ylim=ylim;
these_xlim=xlim;
text(4,0.4,'Hit','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
text(4,0.32,'Miss','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
text(4,0.24,'All trials','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.25,'MI lane x bindFF','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.20,'MI op x bindFF','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.15,'MI op x bindFF x xy','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
% text(0.3,0.10,'MI op x xy','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)

xlim([-2 15])
title('z-scored mutual information (odor concentration and binary dF/F)')
xlabel('Sigma')
ylabel('Cumulative fraction')

%Now find out the clusters for the ROIs significant for odor concentration x binary dF/F mutual infomation

%ssi_all_op_info
sig_ssi_all_op_info=zeros(1,all_info_ii);
cluster_ssi_all_op_info=zeros(1,all_info_ii);
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
for ii_ROI=1:length(cropped_ssi_all_op_info)
        if cropped_ssi_all_op_info(ii_ROI)>=3
            cluster_ssi_all_op_info(ii_ROI)=idxssi(ii_ROI);
            sig_ssi_all_op_info(ii_ROI)=1;
        end
end

%ssi_all_op_info_hit
sig_ssi_all_op_info_hit=zeros(1,all_info_ii);
cluster_ssi_all_op_info_hit=zeros(1,all_info_ii);
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
for ii_ROI=1:length(cropped_ssi_all_op_info_hit)
    if cropped_ssi_all_op_info_hit(ii_ROI)>=3
        cluster_ssi_all_op_info_hit(ii_ROI)=idxssi(ii_ROI);
        sig_ssi_all_op_info_hit(ii_ROI)=1;
    end
end

%ssi_all_op_info_miss
sig_ssi_all_op_info_miss=zeros(1,all_info_ii);
cluster_ssi_all_op_info_miss=zeros(1,all_info_ii);
% plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
for ii_ROI=1:length(cropped_ssi_all_op_info_miss)
        if cropped_ssi_all_op_info_miss(ii_ROI)>=3
            cluster_ssi_all_op_info_miss(ii_ROI)=idxssi(ii_ROI);
            sig_ssi_all_op_info_miss(ii_ROI)=1;
        end
end

%Now plot the cluster distribution for ssi_op
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

%Hit spatial information 
for ii_k=1:3
    this_fraction=sum(sig_ssi_all_op_info_hit(cluster_ssi_all_op_info_hit==ii_k))/sum(sig_ssi_all_op_info_hit);
       switch ii_k
            case 1
                bar(bar_offset,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+1,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+2,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end
      
end



%Miss spatial information 
for ii_k=1:3
    this_fraction=sum(sig_ssi_all_op_info_miss(cluster_ssi_all_op_info_miss==ii_k))/sum(sig_ssi_all_op_info_miss);
       switch ii_k
            case 1
                bar(bar_offset+4,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+5,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+6,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end
  
end



%lBoth mutual spatial information 
for ii_k=1:3
    this_fraction=sum(sig_ssi_all_op_info(cluster_ssi_all_op_info==ii_k))/sum(sig_ssi_all_op_info);
       switch ii_k
            case 1
                bar(bar_offset+8,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+9,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+10,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end

end


text(4,0.8,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(4,0.72,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(4,0.64,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([0 1 2 4 5 6 8 9 10])
xticklabels({'','Hit','','','Miss','','','All',''})
xtickangle(45)

title(['Fraction of odor concentration x bin dFF SSI'])
ylabel('Fraction SSI')
ylim([0 1])
xlim([-1 11])


%Now plot the cluster distribution for ssi_op
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


%lBoth mutual spatial information 
for ii_k=1:3
    this_fraction=sum(sig_ssi_all_op_info(cluster_ssi_all_op_info==ii_k))/sum(sig_ssi_all_op_info);
       switch ii_k
            case 1
                bar(bar_offset+1,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset+2,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            case 3
                bar(bar_offset+3,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
       end

end


text(1,0.8,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(1,0.72,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(1,0.64,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xticks([1 2 3])
xticklabels({'Cluster 1','Cluster 2','Cluster 3'})
xtickangle(45)

title(['Cluster distribution for zMI of odor concentration x bin dFF'])
ylabel('Fraction of cells')
ylim([0 1])
xlim([0 4])

%Now do tsne analysis with SSI for odorant concentration
% ssi_all_spatial_info_hit_no_nan=ssi_all_spatial_info_hit((~isnan(ssi_all_spatial_info_hit))&(~isnan(ssi_all_spatial_info_miss))&(~isnan(ssi_all_spatial_info)));
% ssi_all_spatial_info_miss_no_nan=ssi_all_spatial_info_miss((~isnan(ssi_all_spatial_info_hit))&(~isnan(ssi_all_spatial_info_miss))&(~isnan(ssi_all_spatial_info)));
% ssi_all_spatial_info_no_nan=ssi_all_spatial_info((~isnan(ssi_all_spatial_info_hit))&(~isnan(ssi_all_spatial_info_miss))&(~isnan(ssi_all_spatial_info)));

datassi_op = [ssi_all_op_info' ssi_all_spatial_info_hit ssi_all_spatial_info_miss ssi_all_spatial_info];

ii_all_ROIs_op=1:length(ssi_all_op_info_hit);
not_nan_ii_allROIs_op=ii_all_ROIs_op((~isnan(ssi_all_op_info_hit))&(~isnan(ssi_all_op_info_miss))&(~isnan(ssi_all_op_info)));
handles_outic.not_nan_ii_allROIs_op=not_nan_ii_allROIs;

cropped_ssi_all_op_info_hit=zeros(length(not_nan_ii_allROIs_op),1);
cropped_ssi_all_spatial_info_miss=zeros(length(not_nan_ii_allROIs_op),1);
cropped_ssi_all_spatial_info=zeros(length(not_nan_ii_allROIs_op),1);
% cropped_all_spatial_rho_hit_vs_miss=zeros(length(not_nan_ii_allROIs),1);
% cropped_all_delta_center_of_mass=zeros(length(not_nan_ii_allROIs),1);


for jj=1:length(not_nan_ii_allROIs_op)
    cropped_ssi_all_op_info_hit(jj)=ssi_all_op_info_hit(not_nan_ii_allROIs_op(jj));
    cropped_ssi_all_op_info_miss(jj)=ssi_all_op_info_miss(not_nan_ii_allROIs_op(jj));
    cropped_ssi_all_op_info(jj)=ssi_all_op_info(not_nan_ii_allROIs_op(jj));
    % cropped_all_spatial_rho_hit_vs_miss(jj)=imps.all_spatial_rho_hit_vs_miss(not_nan_ii_allROIs(jj));
    % cropped_all_delta_center_of_mass(jj)=imps.all_delta_center_of_mass(not_nan_ii_allROIs(jj));
end

handles_outic.cropped_ssi_all_op_info_hit=cropped_ssi_all_op_info_hit;
handles_outic.cropped_ssi_all_op_info_miss=cropped_ssi_all_op_info_miss;
handles_outic.cropped_ssi_all_op_info=cropped_ssi_all_op_info;
% handles_outic.cropped_all_spatial_rho_hit_vs_miss=cropped_all_spatial_rho_hit_vs_miss;
% handles_outic.cropped_all_delta_center_of_mass=cropped_all_delta_center_of_mass;

% Step 1 Transpose the data if necessary (N samples x 3 variables)
% datassi_op = datassi_op';


% Step 2: Run t-SNE
% The default parameters are usually sufficient, but you can adjust them
% For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% mappedX = tsne(datassi, 'NumPCAComponents', 2, 'Perplexity', 30);
% mappedX = tsne(datassi,'Algorithm','exact','Distance','mahalanobis'); %Works well
% mappedX = tsne(datassi,'Algorithm','exact','Distance','cosine'); %works better
mappedX_op = tsne(datassi_op,'Distance','cosine'); %works best

% mappedX = tsne(datassi,'Algorithm','exact','Distance','chebychev'); %Ok, yields two clusters
% mappedX_op = tsne(datassi,'Algorithm','exact','Distance','euclidean'); %OK, three clusters

% %Step 3: Plot the results
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
 
no_clusters=3;
gm_op = fitgmdist(mappedX_op, no_clusters); % Assuming 2 clusters
idxssi_op = cluster(gm_op, mappedX_op);

handles_outic.idxssi_op=idxssi_op;

% %Report the mean infos for each cluster
% per_cluster_op=[];
% for ii_k=1:no_clusters
%     these_ssi_all_op_info_hit=cropped_ssi_all_op_info_hit(idxssi==ii_k);
%     these_ssi_all_op_info_miss=cropped_ssi_all_op_info_miss(idxssi==ii_k);
%     these_ssi_all_op_info=cropped_ssi_all_op_info(idxssi==ii_k);
%     mean_hi=mean(these_ssi_all_op_info_hit(~isnan(these_ssi_all_op_info_hit)));
%     per_cluster_op.cluster(ii_k).mean_hi=mean_hi;
%     per_cluster_op.cluster(ii_k).these_ssi_all_op_info_hit=these_ssi_all_op_info_hit(~isnan(these_ssi_all_op_info_hit));
%     mean_miss=mean(these_ssi_all_op_info_miss(~isnan(these_ssi_all_op_info_miss)));
%     per_cluster_op.cluster(ii_k).mean_miss=mean_miss;
%     per_cluster_op.cluster(ii_k).these_ssi_all_op_info_miss=these_ssi_all_op_info_miss(~isnan(these_ssi_all_op_info_miss));
%     mean_both=mean(these_ssi_all_op_info(~isnan(these_ssi_all_op_info)));
%     per_cluster_op.cluster(ii_k).mean_both=mean_both;
%     per_cluster_op.cluster(ii_k).these_ssi_all_op_info=these_ssi_all_op_info(~isnan(these_ssi_all_op_info));
% 
%     fprintf(1, ['\nMean SSI for op cluster ' num2str(ii_k) ' lane1, lane 4, both: '...
%         num2str(mean_hi) ' ' num2str(mean_miss) ' ' num2str(mean_both) '\n'])
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
for ii_k=1:no_clusters
    switch ii_k
        case 1
            plot(mappedX_op(idxssi_op==ii_k,1),mappedX_op(idxssi_op==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
        case 2
            plot(mappedX_op(idxssi_op==ii_k,1),mappedX_op(idxssi_op==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
        case 3
            plot(mappedX_op(idxssi_op==ii_k,1),mappedX_op(idxssi_op==ii_k,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
            % case 4
            % plot(mappedX_op(idxssi==ii_k,1),mappedX_op(idxssi==ii_k,2),'.','MarkerFaceColor',[204/255 121/255 167/255],'MarkerEdgeColor',[204/255 121/255 167/255]);
    end
end

text(-37,30,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
text(-37,20,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
text(-37,10,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)

xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE for SSI for odor concentration and dFF and space and dFF (useful?)' ]);

pfff=1;
% 
% %Plot cumulative histograms for SSIs
% for ii_k=1:no_clusters
%     %SSIs
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
%     hFig=figure(figureNo);
%     hold on
% 
%     ax=gca;ax.LineWidth=3;
%     set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
%     these_ssis=per_cluster_op.cluster(ii_k).these_ssi_all_op_info_hit;
%     [f_aic,x_aic] = drg_ecdf(these_ssis);
%     plot(x_aic,f_aic,'Color',[204/255 121/255 167/255],'LineWidth',3)
% 
%     these_ssis=per_cluster_op.cluster(ii_k).these_ssi_all_op_info_miss;
%     [f_aic,x_aic] = drg_ecdf(these_ssis); %xy x bindFF for lane 4
%     plot(x_aic,f_aic,'Color',[0/255 114/255 178/255],'LineWidth',3)
% 
%     % these_ssis=per_cluster.cluster(ii_k).these_ssi_all_op_info;
%     % [f_aic,x_aic] = drg_ecdf(these_ssis);   %xy x bindFF x lane
%     % plot(x_aic,f_aic,'Color',[213/255 94/255 0/255],'LineWidth',3)
% 
%     plot([3 3],[0 1],'-k','LineWidth',2)
% 
%     title(['SSIs for op Cluster ' num2str(ii_k)])
%     ylabel('Cumulative probability')
%     xlabel('SSI')
% 
%     text(5,0.4,'Hit','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
%     text(5,0.3,'Miss','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
%     % text(10,0.2,'SSI MI xy x bindFF both lanes','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
%     xlim([-5 15])
% end
% 
% 
% 
% %Plot cumulative histograms for SSI both lanes
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
% for ii_k=1:no_clusters
%     %SSIs
% 
%     switch ii_k
%         case 1
%             these_ssis=per_cluster_op.cluster(ii_k).these_ssi_all_op_info;
%             [f_aic,x_aic] = drg_ecdf(these_ssis);
%             plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)
%         case 2
% 
%             these_ssis=per_cluster_op.cluster(ii_k).these_ssi_all_op_info;
%             [f_aic,x_aic] = drg_ecdf(these_ssis);
%             plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)
%         case 3
%             these_ssis=per_cluster_op.cluster(ii_k).these_ssi_all_op_info;
%             [f_aic,x_aic] = drg_ecdf(these_ssis);
%             plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)
%     end
% 
% 
% end
% 
% 
% text(10,0.4,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(10,0.3,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(10,0.2,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% title(['SSI odor concentration for all trials'])
% ylabel('Cumulative probability')
% xlabel('SSI')
% xlim([-5 35])
% 
% %Difference in SSIs
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
% for ii_k=1:no_clusters
% 
%     these_ssis_hit=per_cluster.cluster(ii_k).these_ssi_all_op_info_hit;
%     these_ssis_miss=per_cluster.cluster(ii_k).these_ssi_all_op_info_miss;
% 
%     switch ii_k
%         case 1
%             [f_aic,x_aic] = drg_ecdf(these_ssis_hit-these_ssis_miss);
%             plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)
%         case 2
%             [f_aic,x_aic] = drg_ecdf(these_ssis_hit-these_ssis_miss); %xy x bindFF for lane 4
%             plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)
%         case 3
%             [f_aic,x_aic] = drg_ecdf(these_ssis_hit-these_ssis_miss);   %xy x bindFF x lane
%             plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)
%     end
% 
% 
% end
% plot([0 0],[0 1],'-k')
% title(['Odor concentration SSI hit -SSI miss'])
% ylabel('Cumulative probability')
% xlabel('Delta SSI')
% xlim([-15 15])
% 
% text(5,0.4,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(5,0.3,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(5,0.2,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 

% 
% %Plot histograms
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% hold on
% 
% [f_aic,x_aic] = drg_ecdf(all_info_both_hit_miss);
% plot(x_aic,f_aic,'Color',[230/255 159/255 0/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_spatial_info_sh);
% plot(x_aic,f_aic,'-.','Color',[230/255 159/255 0/255],'LineWidth',3)
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_op_bin_info); %xy x bindFF for lane 4
% plot(x_aic,f_aic,'Color',[86/255 180/255 233/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_op_bin_info_sh);
% plot(x_aic,f_aic,'-.','Color',[86/255 180/255 233/255],'LineWidth',3)
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_info_mutual_xy_op_bin);
% plot(x_aic,f_aic,'Color',[0/255 158/255 115/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_info_mutual_xy_op_bin_sh);
% plot(x_aic,f_aic,'-.','Color',[0/255 158/255 115/255],'LineWidth',3)
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_info_mutual_info_dFFbin_xy_op_bin);
% plot(x_aic,f_aic,'Color',[0/255 114/255 178/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_info_mutual_info_dFFbin_xy_op_bin_sh);
% plot(x_aic,f_aic,'-.','Color',[0/255 114/255 178/255],'LineWidth',3)
% 
% 
% these_ylim=ylim;
% these_xlim=xlim;
% text(0.45,0.4,'MI xy x bindFF both lanes','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(0.45,0.35,'MI op x bindFF','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(0.45,0.30,'MI op x xy','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% text(0.45,0.25,'MI op x bindFF x xy','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
% % text(0.3,0.20,'MI op x bindFF','Color',[0/255 114/255 178/255],'FontWeight','bold','FontSize',16)
% % text(0.3,0.15,'MI op x bindFF x xy','Color',[213/255 94/255 0/255],'FontWeight','bold','FontSize',16)
% % text(0.3,0.10,'MI op x xy','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
% 
% title('Cumulative histograms for mutual information')
% xlabel('Information bits')
% ylabel('Cumulative fraction')


% %Plot histograms
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% hold on
% 
% [f_aic,x_aic] = drg_ecdf(all_info_hit);
% plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_spatial_info_hit_sh);
% plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)
% 
% text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
% text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)
% 
% title('MI xy x bindFF Lane 1, grey=shuffled')
% xlabel('Information bits')
% ylabel('Fraction')
% 
% 
% %Plot histograms
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% hold on
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_info_miss);
% plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_spatial_info_miss_sh);
% plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)
% 
% 
% text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
% text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)
% 
% title('MI xy x bindFF Lane 4, grey=shuffled')
% xlabel('Information bits')
% ylabel('Fraction')
% 
% 
% %Plot histograms
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% hold on
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_info_mutual_info14);
% plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_info_mutual_info14_sh);
% plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)
% 
% text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
% text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)
% 
% 
% title('MI xy x bindFF x lane, grey=shuffled')
% xlabel('Information bits')
% ylabel('Fraction')
% 
% 
% 
% %Plot histograms
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% hold on
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_info_lane);
% plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_info_lane_sh);
% plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)
% 
% text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
% text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)
% 
% 
% title('MI lane x bindFF, grey=shuffled')
% xlabel('Information bits')
% ylabel('Fraction')
% 
% 
% %Plot histograms
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% hold on
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_op_bin_info);
% plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_op_bin_info_sh);
% plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)
% 
% text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
% text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)
% 
% 
% title('MI op x bindFF, grey=shuffled')
% xlabel('Information bits')
% ylabel('Fraction')
% 
% 
% %Plot histograms
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
% hold on
% 
% 
% [f_aic,x_aic] = drg_ecdf(all_info_mutual_info_dFFbin_xy_op_bin);
% plot(x_aic,f_aic,'Color',[0/255 0/255 0/255],'LineWidth',3)
% 
% [f_aic,x_aic] = drg_ecdf(all_all_info_mutual_info_dFFbin_xy_op_bin_sh);
% plot(x_aic,f_aic,'Color',[0.7 0.7 0.7],'LineWidth',3)
% 
% text(0.2,0.7,'Original','Color','k','FontWeight','bold','FontSize',16)
% text(0.2,0.62,'Shuffled','Color',[0.7 0.7 0.7],'FontWeight','bold','FontSize',16)
% 
% 
% title('MI op x bindFF x xy, grey=shuffled')
% xlabel('Information bits')
% ylabel('Fraction')

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
% 
% %Now do tsne analysis with all_info_mutual_info14
% data = [all_info_hit; all_info_miss; all_info_mutual_info14];
% 
% % Step 1 Transpose the data if necessary (N samples x 3 variables)
% data = data';
% 
% rng('default') % for reproducibility
% 
% % Step 2: Run t-SNE
% % The default parameters are usually sufficient, but you can adjust them
% % For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% % mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
% % mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
% mappedX = tsne(data,'Algorithm','exact','Distance','cosine'); %works better
% 
% % mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
% % mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK
% 
% % Step 3: Plot the results
% % figureNo=figureNo+1;
% % try
% %     close(figureNo)
% % catch
% % end
% % hFig=figure(figureNo);
% % hold on
% %
% % ax=gca;ax.LineWidth=3;
% % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% %
% % plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
% %
% % xlabel('t-SNE Component 1');
% % ylabel('t-SNE Component 2');
% % title(['t-SNE Information Content ' ]);
% 
% %Cool, it looks like we have three clear clusters
% 
% % Perform k-means clustering
% %kmeans misclassifies a small number of points
% % k = 3; % Number of clusters
% % [idx, centroids] = kmeans(mappedX, k,'Distance','cityblock');
% 
% gm = fitgmdist(mappedX, 3); % Assuming 3 clusters
% idx = cluster(gm, mappedX);
% 
% 
% %Report the mean infos for each cluster
% for ii_k=1:3
%     mean_hi=mean(all_info_hit(idx==ii_k));
%     mean_miss=mean(all_info_miss(idx==ii_k));
%     mean_mutual14=mean(all_info_mutual_info14(idx==ii_k));
% 
%     fprintf(1, ['\nInformation content for cluster ' num2str(ii_k) ' lane1, lane 2, mutual: '...
%         num2str(mean_hi) ' ' num2str(mean_miss) ' ' num2str(mean_mutual14) '\n'])
% end
% 
% % Step 3: Plot the results
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
% for ii_k=1:3
%     switch ii_k
%         case 1
%             plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
%         case 2
%             plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
%         case 3
%             plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
%             % case 4
%             % plot(mappedX(idx==ii_k,1),mappedX(idx==ii_k,2),'.','MarkerFaceColor',[204/255 121/255 167/255],'MarkerEdgeColor',[204/255 121/255 167/255]);
%     end
% end
% 
% text(-60,70,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(-60,60,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(-60,50,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE MI for lanes, xy and dFF vs. MI for xy in each lane' ]);
% 
% %Now do tsne analysis with all_info_mutual_info_dFFbin_xy_op_bin
% data2 = [all_info_hit; all_info_miss; all_info_mutual_info_dFFbin_xy_op_bin]; % all_info_mutual_info_dFFbin_xy_op_bin
% 
% % Step 1 Transpose the data if necessary (N samples x 3 variables)
% data2 = data2';
% 
% rng('default') % for reproducibility
% 
% % Step 2: Run t-SNE
% % The default parameters are usually sufficient, but you can adjust them
% % For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% % mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
% % mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
% mappedXl = tsne(data2,'Algorithm','exact','Distance','cosine'); %works better
% 
% % mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
% % mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK
% 
% % Step 3: Plot the results
% % figureNo=figureNo+1;
% % try
% %     close(figureNo)
% % catch
% % end
% % hFig=figure(figureNo);
% % hold on
% %
% % ax=gca;ax.LineWidth=3;
% % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% %
% % plot(mappedXl(:,1),mappedXl(:,2),'.','Color',[0.7 0.7 0.7]);
% %
% % xlabel('t-SNE Component 1');
% % ylabel('t-SNE Component 2');
% % title(['t-SNE Information Content for odor plume, space' ]);
% 
% %Cool, it looks like we have three clear clusters
% 
% % Perform k-means clustering
% %kmeans misclassifies a small number of points
% % k = 3; % Number of clusters
% % [idx, centroids] = kmeans(mappedX, k,'Distance','cityblock');
% 
% nclusXl=3;
% gml = fitgmdist(mappedXl, nclusXl); % Assuming 3 clusters
% idxl = cluster(gml, mappedXl);
% 
% 
% % %Report the mean infos for each cluster
% % for ii_k=1:2
% %     mean_hil=mean(all_info_hit(idxl==ii_k));
% %     mean_missl=mean(all_info_miss(idxl==ii_k));
% %     mean_lane=mean(all_info_lane(idxl==ii_k));
% %
% %     fprintf(1, ['\nInformation content for cluster ' num2str(ii_k) ' lane1, lane 2, lane: '...
% %         num2str(mean_hi) ' ' num2str(mean_miss) ' ' num2str(mean_mutual14) '\n'])
% % end
% 
% % Step 3: Plot the results
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
% for ii_k=1:nclusXl
%     switch ii_k
%         case 1
%             plot(mappedXl(idxl==ii_k,1),mappedXl(idxl==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
%         case 2
%             plot(mappedXl(idxl==ii_k,1),mappedXl(idxl==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
%         case 3
%             plot(mappedXl(idxl==ii_k,1),mappedXl(idxl==ii_k,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
%     end
% end
% 
% text(-60,70,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(-60,60,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% if nclusXl>2
%     text(-60,50,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% end
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE MI for odor plume xy and dFF vs. MI for xy in each lane' ]);
% 

% 
% %Now do tsne analysis with Fusi mutual info
% data3 = [imps.all_Fusi_SSIl1'; imps.all_Fusi_SSIl4'; all_info_mutual_info_dFFbin_xy_op_bin]; % all_info_mutual_info_dFFbin_xy_op_bin
% 
% % Step 1 Transpose the data if necessary (N samples x 3 variables)
% data3 = data3';
% 
% rng('default') % for reproducibility
% 
% % Step 2: Run t-SNE
% % The default parameters are usually sufficient, but you can adjust them
% % For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% % mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
% % mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
% mappedXF = tsne(data3,'Algorithm','exact','Distance','cosine'); %works better
% 
% % mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
% % mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK
% 
% % Step 3: Plot the results
% % figureNo=figureNo+1;
% % try
% %     close(figureNo)
% % catch
% % end
% % hFig=figure(figureNo);
% % hold on
% %
% % ax=gca;ax.LineWidth=3;
% % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% %
% % plot(mappedXF(:,1),mappedXF(:,2),'.','Color',[0.7 0.7 0.7]);
% %
% % xlabel('t-SNE Component 1');
% % ylabel('t-SNE Component 2');
% % title(['t-SNE Information Content for odor plume, Fusi space' ]);
% 
% %Cool, it looks like we have three clear clusters
% 
% % Perform k-means clustering
% %kmeans misclassifies a small number of points
% % k = 3; % Number of clusters
% % [idx, centroids] = kmeans(mappedX, k,'Distance','cityblock');
% 
% nclusXF=3;
% gmF = fitgmdist(mappedXF, nclusXF); % Assuming 3 clusters
% idxF = cluster(gmF, mappedXF);
% 
% 
% % %Report the mean infos for each cluster
% % for ii_k=1:2
% %     mean_hil=mean(all_info_hit(idxl==ii_k));
% %     mean_missl=mean(all_info_miss(idxl==ii_k));
% %     mean_lane=mean(all_info_lane(idxl==ii_k));
% %
% %     fprintf(1, ['\nInformation content for cluster ' num2str(ii_k) ' lane1, lane 2, lane: '...
% %         num2str(mean_hi) ' ' num2str(mean_miss) ' ' num2str(mean_mutual14) '\n'])
% % end
% 
% % Step 3: Plot the results
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
% for ii_k=1:nclusXF
%     switch ii_k
%         case 1
%             plot(mappedXF(idxF==ii_k,1),mappedXF(idxF==ii_k,2),'.','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255]);
%         case 2
%             plot(mappedXF(idxF==ii_k,1),mappedXF(idxF==ii_k,2),'.','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255]);
%         case 3
%             plot(mappedXF(idxF==ii_k,1),mappedXF(idxF==ii_k,2),'.','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255]);
%     end
% end
% 
% text(-60,70,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(-60,60,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% if nclusXF>2
%     text(-60,50,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% end
% 
% xlabel('t-SNE Component 1');
% ylabel('t-SNE Component 2');
% title(['t-SNE MI for odor plume, xy and dFF vs. MI for Fusi space' ]);

%Estimate those significantly above SD

% 
% ii_ROI_positive=0;
% sig_all_info_hit=zeros(1,all_info_ii);
% for ii_ROI=1:all_info_ii
%     if all_info_hit(ii_ROI)>=3*sd_all_spatial_info_hit_sh(ii_ROI)+mean_all_spatial_info_hit_sh(ii_ROI)
%         ii_ROI_positive=ii_ROI_positive+1;
%         ii_k=idx(ii_ROI);
%         sig_all_info_hit(ii_ROI)=1;
%     end
% end
% 
% 
% ii_ROI_positive=0;
% sig_all_info_miss=zeros(1,all_info_ii);
% for ii_ROI=1:all_info_ii
%     if all_info_miss(ii_ROI)>=3*sd_all_spatial_info_miss_sh(ii_ROI)+mean_all_spatial_info_miss_sh(ii_ROI)
%         ii_ROI_positive=ii_ROI_positive+1;
%         ii_k=idx(ii_ROI);
%         sig_all_info_miss(ii_ROI)=1;
% 
%     end
% end
% 
% sig_all_info_mutual_info14=zeros(1,all_info_ii);
% % plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
% ii_ROI_positive=0;
% for ii_ROI=1:all_info_ii
%     if all_info_mutual_info14(ii_ROI)>=3*sd_all_info_mutual_info14_sh(ii_ROI)+mean_all_info_mutual_info14_sh(ii_ROI)
%         ii_ROI_positive=ii_ROI_positive+1;
%         ii_k=idx(ii_ROI);
%         sig_all_info_mutual_info14(ii_ROI)=1;
% 
%     end
% end
% 
% %Now keep track of significnat for dFFbin_xy_op_bin
% 
% %lane 1
% ii_ROI_positive=0;
% sig_all_info_hit_dFFbin_xy_op_bin=zeros(1,all_info_ii);
% for ii_ROI=1:all_info_ii
%     if all_info_hit(ii_ROI)>=3*sd_all_spatial_info_hit_sh(ii_ROI)+mean_all_spatial_info_hit_sh(ii_ROI)
%         ii_ROI_positive=ii_ROI_positive+1;
%         ii_k=idxl(ii_ROI);
%         sig_all_info_hit_dFFbin_xy_op_bin(ii_ROI)=1;
% 
%     end
% end
% 
% 
% ii_ROI_positive=0;
% sig_all_info_miss_dFFbin_xy_op_bin=zeros(1,all_info_ii);
% for ii_ROI=1:all_info_ii
%     if all_info_miss(ii_ROI)>=3*sd_all_spatial_info_miss_sh(ii_ROI)+mean_all_spatial_info_miss_sh(ii_ROI)
%         ii_ROI_positive=ii_ROI_positive+1;
%         ii_k=idxl(ii_ROI);
%         sig_all_info_miss_dFFbin_xy_op_bin(ii_ROI)=1;
% 
%     end
% end
% 
% 
% %all_info_mutual_info_dFFbin_xy_op_bin
% sig_all_info_mutual_dFFbin_xy_op_bin=zeros(1,all_info_ii);
% % plot(mappedX(:,1),mappedX(:,2),'.','Color',[0.7 0.7 0.7]);
% ii_ROI_positive=0;
% for ii_ROI=1:all_info_ii
%     if all_info_mutual_info_dFFbin_xy_op_bin(ii_ROI)>=3*sd_all_info_mutual_info_dFFbin_xy_op_bin_sh(ii_ROI)+mean_all_info_mutual_info_dFFbin_xy_op_bin_sh(ii_ROI)
%         ii_ROI_positive=ii_ROI_positive+1;
%         ii_k=idxl(ii_ROI);
%         sig_all_info_mutual_dFFbin_xy_op_bin(ii_ROI)=1;
% 
%     end
% end

% 
% 
% %Now plot the fraction of significant information for each cluster for
% %info14
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
% %lane 1 spatial information 
% for ii_k=1:3
%     this_fraction=sum(sig_all_info_hit(idx==ii_k))/sum(idx==ii_k);
%        switch ii_k
%             case 1
%                 bar(bar_offset,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             case 2
%                 bar(bar_offset+4,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             case 3
%                 bar(bar_offset+8,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%        end
% 
% end
% 
% 
% 
% %lane 4 spatial information 
% for ii_k=1:3
%     this_fraction=sum(sig_all_info_miss(idx==ii_k))/sum(idx==ii_k);
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
% %lanes 1 and 4 mutual spatial information 
% for ii_k=1:3
%     this_fraction=sum(sig_all_info_mutual_info14(idx==ii_k))/sum(idx==ii_k);
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
% xticklabels({'Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4'})
% xtickangle(45)
% 
% title(['Significant spatial information for the different clusters, lanes space'])
% ylabel('Fraction SSI')
% ylim([0 1])
% xlim([-1 11])
% 
% 
% %Now plot the fraction of significant information for each cluster for
% %dFFbin_xy_op_bin
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
% %lane 1 spatial information 
% for ii_k=1:nclusXl
%     this_fraction=sum(sig_all_info_hit_dFFbin_xy_op_bin(idxl==ii_k))/sum(idxl==ii_k);
%        switch ii_k
%             case 1
%                 bar(bar_offset,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             case 2
%                 bar(bar_offset+4,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             case 3
%                 bar(bar_offset+8,this_fraction,'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%        end
% 
% end
% 
% 
% 
% %lane 4 spatial information 
% for ii_k=1:nclusXl
%     this_fraction=sum(sig_all_info_miss_dFFbin_xy_op_bin(idxl==ii_k))/sum(idxl==ii_k);
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
% %lanes 1 and 4 mutual spatial information 
% for ii_k=1:nclusXl
%     this_fraction=sum(sig_all_info_mutual_dFFbin_xy_op_bin(idxl==ii_k))/sum(idxl==ii_k);
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
% if nclusXl>2
%     text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% end
% 
% xticks([0 1 2 4 5 6 8 9 10])
% xticklabels({'Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4'})
% xtickangle(45)
% 
% title(['Significant spatial information for the different clusters, space odor plume'])
% ylabel('Fraction SSI')
% ylim([0 1])
% xlim([-1 11])
% 
% %Plot the spatial information for each cluster
% 
% %Now plot the Fusi spatial information for each cluster
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
% edges=[-0.1:0.25:42];
% rand_offset=0.5;
% 
% %lane 1 spatial information
% for ii_k=1:3
%     these_SSI=all_info_hit(idx==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+4,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+4,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+8,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+8,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %lane 4 spatial information 
% for ii_k=1:3
%     these_SSI=all_info_miss(idx==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+1,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+5,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+5,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+9,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+9,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %both lanes mutual spatial information shared between 1 and 4
% for ii_k=1:3
%     these_SSI=all_info_mutual_info14(idx==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+2,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+6,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+6,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+10,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+10,rand_offset,'k','k',1);
%     end
% end
% 
% 
% text(4,14,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,13,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,12,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2 4 5 6 8 9 10])
% xticklabels({'Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4'})
% xtickangle(45)
% 
% title(['Spatial information for the different clusters'])
% ylabel('SI')
% ylim([0 0.5])
% xlim([-1 11])

% 
% %Now plot the spatial information for each cluster
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
% edges=[-0.1:0.25:42];
% rand_offset=0.5;
% 
% %lane 1 spatial information
% for ii_k=1:2
%     these_SSI=all_info_hit(idxl==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+4,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+4,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+8,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+8,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %lane 4 spatial information 
% for ii_k=1:2
%     these_SSI=all_info_miss(idxl==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+1,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+5,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+5,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+9,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+9,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %both lanes mutual spatial information shared between 1 and 4
% for ii_k=1:2
%     these_SSI=all_info_mutual_info14(idxl==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+2,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+6,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+6,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+10,mean(these_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_SSI...
%                 ,edges,bar_offset+10,rand_offset,'k','k',1);
%     end
% end
% 
% 
% text(4,14,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,13,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,12,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2 4 5 6 ])
% xticklabels({'Lane 1','Lane 4','Mutual 1 4','Lane 1','Lane 4','Mutual 1 4'})
% xtickangle(45)
% 
% title(['Spatial information for the different clusters, xy dFF op'])
% ylabel('SI')
% ylim([0 0.5])
% xlim([-1 7])


% 
% 
% %Now plot the Fusi spatial information for each cluster
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
% edges=[-0.1:0.25:42];
% rand_offset=0.5;
% 
% %lane 1 spatial information
% for ii_k=1:3
%     these_Fusi_SSI=imps.all_Fusi_SSIl1(idx==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+4,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+4,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+8,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+8,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %lane 4 spatial information 
% for ii_k=1:3
%     these_Fusi_SSI=imps.all_Fusi_SSIl4(idx==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+1,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+5,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+5,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+9,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+9,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %both lanes mutual spatial information 
% for ii_k=1:3
%     these_Fusi_SSI=imps.all_Fusi_SSI(idx==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+2,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+6,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+6,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+10,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+10,rand_offset,'k','k',1);
%     end
% end
% 
% 
% text(4,14,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,13,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,12,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2 4 5 6 8 9 10])
% xticklabels({'Lane 1','Lane 4','Both','Lane 1','Lane 4','Both','Lane 1','Lane 4','Both'})
% xtickangle(45)
% 
% title(['Fusi significant spatial information for the different clusters, space lane'])
% ylabel('SSI')
% ylim([0 15])
% xlim([-1 11])
% 
% 
% %Now plot the Fusi spatial information for each cluster
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
% edges=[-0.1:0.25:42];
% rand_offset=0.5;
% 
% %lane 1 spatial information
% for ii_k=1:2
%     these_Fusi_SSI=imps.all_Fusi_SSIl1(idxl==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+4,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+4,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+8,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+8,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %lane 4 spatial information 
% for ii_k=1:2
%     these_Fusi_SSI=imps.all_Fusi_SSIl4(idxl==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+1,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+5,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+5,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+9,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+9,rand_offset,'k','k',1);
%     end
% end
% 
% 
% 
% %both lanes mutual spatial information 
% for ii_k=1:2
%     these_Fusi_SSI=imps.all_Fusi_SSI(idxl==ii_k);
%     switch ii_k
%         case 1
%             bar(bar_offset+2,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+6,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+6,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+10,mean(these_Fusi_SSI),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_Fusi_SSI...
%                 ,edges,bar_offset+10,rand_offset,'k','k',1);
%     end
% end
% 
% 
% text(4,14,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% text(4,13,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% text(4,12,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2 4 5 6])
% xticklabels({'Lane 1','Lane 4','Both','Lane 1','Lane 4','Both'})
% xtickangle(45)
% 
% title(['Fusi significant spatial information for the different clusters, xy dFF op'])
% ylabel('SSI')
% ylim([0 15])
% xlim([-1 7])

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
% 
% %Now plot the fraction of significant predictive importance
% sig_pred_imp_x_ii_k=[];
% ii_x=0;
% 
% sig_pred_imp_y_ii_k=[];
% ii_y=0;
% 
% sig_pred_imp_conc_ii_k=[];
% ii_conc=0;
% 
% for ii=1:length(imps.mean_imp_ROI)
%     this_file=imps.mean_imp_fileNo(ii);
%     this_ROI=imps.mean_imp_ROI(ii);
% 
%     ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
%     this_ii_k=idx(ii_info);
% 
%     if imps.sig_pred_imp_x(ii)==1
%         ii_x=ii_x+1;
%         sig_pred_imp_x_ii_k(ii_x)=this_ii_k;
%     end
% 
% 
%     if imps.sig_pred_imp_y(ii)==1
%         ii_y=ii_y+1;
%         sig_pred_imp_y_ii_k(ii_y)=this_ii_k;
%     end
% 
% 
%     if imps.sig_pred_imp_conc(ii)==1
%         ii_conc=ii_conc+1;
%         sig_pred_imp_conc_ii_k(ii_conc)=this_ii_k;
%     end
% 
% end
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
% for ii_k=1:3
%     this_fraction=sum(sig_pred_imp_x_ii_k==ii_k)/length(sig_pred_imp_x_ii_k);
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
% for ii_k=1:3
%     this_fraction=sum(sig_pred_imp_y_ii_k==ii_k)/length(sig_pred_imp_y_ii_k);
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
% for ii_k=1:3
%     this_fraction=sum(sig_pred_imp_conc_ii_k==ii_k)/length(sig_pred_imp_conc_ii_k);
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
% title(['Fraction of prediction important ROIs for the different clusters space, lane'])
% ylabel('Fraction')
% ylim([0 1])
% xlim([-1 11])

% 
% %Now plot the fraction of significant predictive importance
% sig_pred_imp_x_ii_k=[];
% ii_x=0;
% 
% sig_pred_imp_y_ii_k=[];
% ii_y=0;
% 
% sig_pred_imp_conc_ii_k=[];
% ii_conc=0;
% 
% all_info_hit_impx=[];
% all_info_miss_impx=[];
% all_info_mutual_info_dFFbin_xy_op_bin_impx=[];
% 
% all_info_hit_impy=[];
% all_info_miss_impy=[];
% all_info_mutual_info_dFFbin_xy_op_bin_impy=[];
% 
% all_info_hit_impconc=[];
% all_info_miss_impconc=[];
% all_info_mutual_info_dFFbin_xy_op_bin_impconc=[];
% 
% for ii=1:length(imps.mean_imp_ROI)
%     this_file=imps.mean_imp_fileNo(ii);
%     this_ROI=imps.mean_imp_ROI(ii);
% 
%     ii_info=find((all_info_fileNo==this_file)&(all_info_ii_ROI==this_ROI));
%     this_ii_k=idxl(ii_info);
% 
%     if imps.sig_pred_imp_x(ii)==1
%         ii_x=ii_x+1;
%         sig_pred_imp_x_ii_k(ii_x)=this_ii_k;
%         all_info_hit_impx=[all_info_hit_impx all_info_hit(ii_info)];
%         all_info_miss_impx=[all_info_miss_impx all_info_miss(ii_info)];
%         all_info_mutual_info_dFFbin_xy_op_bin_impx=[all_info_mutual_info_dFFbin_xy_op_bin_impx all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
% 
%     end
% 
% 
%     if imps.sig_pred_imp_y(ii)==1
%         ii_y=ii_y+1;
%         sig_pred_imp_y_ii_k(ii_y)=this_ii_k;
%         all_info_hit_impy=[all_info_hit_impy all_info_hit(ii_info)];
%         all_info_miss_impy=[all_info_miss_impy all_info_miss(ii_info)];
%         all_info_mutual_info_dFFbin_xy_op_bin_impy=[all_info_mutual_info_dFFbin_xy_op_bin_impy all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
%     end
% 
% 
%     if imps.sig_pred_imp_conc(ii)==1
%         ii_conc=ii_conc+1;
%         sig_pred_imp_conc_ii_k(ii_conc)=this_ii_k;
%         all_info_hit_impconc=[all_info_hit_impconc all_info_hit(ii_info)];
%         all_info_miss_impconc=[all_info_miss_impconc all_info_miss(ii_info)];
%         all_info_mutual_info_dFFbin_xy_op_bin_impconc=[all_info_mutual_info_dFFbin_xy_op_bin_impconc all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
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
% bar(bar_offset,mean(all_info_hit_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_hit_impx...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_info_miss_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_miss_impx...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_info_mutual_info_dFFbin_xy_op_bin_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_mutual_info_dFFbin_xy_op_bin_impx...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar_offset=bar_offset+1;
% 
% %MI for predictive importance for y
% bar(bar_offset,mean(all_info_hit_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_hit_impy...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_info_miss_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_miss_impy...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_info_mutual_info_dFFbin_xy_op_bin_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_mutual_info_dFFbin_xy_op_bin_impy...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar_offset=bar_offset+1;
% 
% 
% %MI for predictive importance for conc
% bar(bar_offset,mean(all_info_hit_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_hit_impconc...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_info_miss_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_miss_impconc...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_info_mutual_info_dFFbin_xy_op_bin_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_info_mutual_info_dFFbin_xy_op_bin_impconc...
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
%     this_fraction=sum(sig_pred_imp_x_ii_k==ii_k)/length(sig_pred_imp_x_ii_k);
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
%     this_fraction=sum(sig_pred_imp_y_ii_k==ii_k)/length(sig_pred_imp_y_ii_k);
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
%     this_fraction=sum(sig_pred_imp_conc_ii_k==ii_k)/length(sig_pred_imp_conc_ii_k);
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
% all_infoF_hi_impx=[];
% all_infoF_miss_impx=[];
% all_infoF_mutual_info_dFFbin_xy_op_bin_impx=[];
% 
% all_infoF_hi_impy=[];
% all_infoF_miss_impy=[];
% all_infoF_mutual_info_dFFbin_xy_op_bin_impy=[];
% 
% all_infoF_hi_impconc=[];
% all_infoF_miss_impconc=[];
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
%         all_infoF_hi_impx=[all_infoF_hi_impx all_info_hit(ii_info)];
%         all_infoF_miss_impx=[all_infoF_miss_impx all_info_miss(ii_info)];
%         all_infoF_mutual_info_dFFbin_xy_op_bin_impx=[all_infoF_mutual_info_dFFbin_xy_op_bin_impx all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
% 
%     end
% 
% 
%     if imps.sig_pred_imp_y(ii)==1
%         ii_y=ii_y+1;
%         sig_predF_imp_y_ii_k(ii_y)=this_ii_k;
%         all_infoF_hi_impy=[all_infoF_hi_impy all_info_hit(ii_info)];
%         all_infoF_miss_impy=[all_infoF_miss_impy all_info_miss(ii_info)];
%         all_infoF_mutual_info_dFFbin_xy_op_bin_impy=[all_infoF_mutual_info_dFFbin_xy_op_bin_impy all_info_mutual_info_dFFbin_xy_op_bin(ii_info)];
%     end
% 
% 
%     if imps.sig_pred_imp_conc(ii)==1
%         ii_conc=ii_conc+1;
%         sig_predF_imp_conc_ii_k(ii_conc)=this_ii_k;
%         all_infoF_hi_impconc=[all_infoF_hi_impconc all_info_hit(ii_info)];
%         all_infoF_miss_impconc=[all_infoF_miss_impconc all_info_miss(ii_info)];
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
% bar(bar_offset,mean(all_infoF_hi_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_hi_impx...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_miss_impx),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_miss_impx...
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
% bar(bar_offset,mean(all_infoF_hi_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_hi_impy...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_miss_impy),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_miss_impy...
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
% bar(bar_offset,mean(all_infoF_hi_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_hi_impconc...
%     ,edges,bar_offset,rand_offset,'k','k',2);
% bar_offset=bar_offset+1;
% 
% bar(bar_offset,mean(all_infoF_miss_impconc),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
% 
% [mean_out, CIout, violin_x]=drgViolinPoint(all_infoF_miss_impconc...
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
% 
% %Spatial correlations of the spatial pattern
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
% edges=[0:0.05:1];
% rand_offset=0.5;
% 
% %predictive importance for x
% for ii_k=1:3
%     these_rhos=imps.all_spatial_rho_hit_vs_miss(idx==ii_k);
%     these_rhos=these_rhos(~isnan(these_rhos));
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+1,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+2,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%     end
% 
% 
% 
% end
% 
% 
% 
% 
% % 
% % text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% % text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% % text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2])
% xticklabels({'Cluster 1','Cluster 2','Cluster 3'})
% xtickangle(45)
% 
% title(['Spatial correlation for the different clusters'])
% ylabel('Rho')
% ylim([-1 1])
% xlim([-1 3])
% 
% 
% %Center of mass distances of the spatial pattern
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
% edges=[0:25:500];
% rand_offset=0.5;
% 
% %predictive importance for x
% for ii_k=1:3
%     these_dcoms=imps.all_delta_center_of_mass(idx==ii_k);
%     these_dcoms=these_dcoms(~isnan(these_dcoms));
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+1,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+2,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%     end
% 
% 
% 
% end
% 
% % 
% % text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% % text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% % text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1 2])
% xticklabels({'Cluster 1','Cluster 2','Cluster 3'})
% xtickangle(45)
% 
% title(['Difference in the center of mass'])
% ylabel('delta com')
% ylim([0 350])
% xlim([-1 3])
% 
% 
% %Spatial correlations of the spatial pattern
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .2 .3])
% bar_offset=0;
% edges=[0:0.05:1];
% rand_offset=0.5;
% 
% %predictive importance for x
% for ii_k=1:2
%     these_rhos=imps.all_spatial_rho_hit_vs_miss(idxl==ii_k);
%     these_rhos=these_rhos(~isnan(these_rhos));
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+1,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+2,mean(these_rhos),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_rhos...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%     end
% 
% 
% 
% end
% 
% 
% 
% 
% % 
% % text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% % text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% % text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1])
% xticklabels({'Cluster 1','Cluster 2'})
% xtickangle(45)
% 
% title(['Spatial correlation for the different clusters, xy dFF op'])
% ylabel('Rho')
% ylim([-1 1])
% xlim([-1 2])
% 
% 
% %Center of mass distances of the spatial pattern
% figureNo=figureNo+1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .2 .3])
% bar_offset=0;
% edges=[0:25:500];
% rand_offset=0.5;
% 
% %predictive importance for x
% for ii_k=1:2
%     these_dcoms=imps.all_delta_center_of_mass(idxl==ii_k);
%     these_dcoms=these_dcoms(~isnan(these_dcoms));
%     switch ii_k
%         case 1
%             bar(bar_offset,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
%         case 2
%             bar(bar_offset+1,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
%                 ,edges,bar_offset+1,rand_offset,'k','k',1);
%         case 3
%             bar(bar_offset+2,mean(these_dcoms),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             %Violin plot
%             [mean_out, CIout,violin_x]=drgViolinPoint(these_dcoms...
%                 ,edges,bar_offset+2,rand_offset,'k','k',1);
%     end
% 
% 
% 
% end
% 
% % 
% % text(4,0.7,'Cluster 1','Color',[230/255 159/255 0/255],'FontWeight','bold','FontSize',16)
% % text(4,0.62,'Cluster 2','Color',[86/255 180/255 233/255],'FontWeight','bold','FontSize',16)
% % text(4,0.54,'Cluster 3','Color',[0/255 158/255 115/255],'FontWeight','bold','FontSize',16)
% 
% xticks([0 1])
% xticklabels({'Cluster 1','Cluster 2'})
% xtickangle(45)
% 
% title(['Difference in the center of mass, xy dFF op'])
% ylabel('delta com')
% ylim([0 350])
% xlim([-1 2])
% 
% 
% save([save_PathIC save_FileIC],'handles_outic','-v7.3')

%Odor cell placeholders
oc_dFF_norm_map_all=zeros(sum(ssi_all_op_info_hit>=3),11);
oc_dFF_norm_map_all_n=zeros(sum(ssi_all_op_info_hit>=3),11);

oc_dFF_norm_map_fileNo_all=zeros(1,sum(ssi_all_op_info_hit>=3));
oc_dFF_norm_map_iiROI_all=zeros(1,sum(ssi_all_op_info_hit>=3));
oc_ii_all=0;

oc_dFF_norm_map_hit=zeros(sum(ssi_all_op_info_hit>=3),11);
oc_dFF_norm_map_hit_n=zeros(sum(ssi_all_op_info_hit>=3),11);

oc_dFF_norm_map_miss=zeros(sum(ssi_all_op_info_hit>=3),11);
oc_dFF_norm_map_miss_n=zeros(sum(ssi_all_op_info_hit>=3),11);

oc_dFF_norm_map_fileNo=zeros(1,sum(ssi_all_op_info_hit>=3));
oc_dFF_norm_map_iiROI=zeros(1,sum(ssi_all_op_info_hit>=3));
oc_ii=0;

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
            % spatial_rhol_hit_vs_miss=imps.file(fileNo).spatial_rhol_hit_vs_miss;
            % delta_center_of_mass=imps.file(fileNo).delta_center_of_mass;
            % SSI=imps.file(fileNo).information_content;
            % SSIl1=imps.file(fileNo).information_contentl1;
            % SSIl4=imps.file(fileNo).information_contentl4;
            % sparsity=imps.file(fileNo).sparsity;

            %Initialize variables
            this_dFF_activity=zeros(10,10);
            this_dFF_activity_n=zeros(10,10);
            sum_dFF_activity=0;

            this_Hit_activity=zeros(10,10);
            this_Hit_activity_n=zeros(10,10);
            sum_Hit_activity=0;

            this_Miss_activity=zeros(10,10);
            this_Miss_activity_n=zeros(10,10);
            sum_Miss_activity=0;

            this_dFF_op_activity=zeros(1,length(odor_c_bounds));
            this_dFF_op_activity_n=zeros(1,length(odor_c_bounds));

            these_dFF_op_activity=zeros(trials.odor_trNo,length(odor_c_bounds));
            these_dFF_op_activity_n=zeros(trials.odor_trNo,length(odor_c_bounds));
            

            this_Hit_op_activity=zeros(1,length(odor_c_bounds));
            this_Hit_op_activity_n=zeros(1,length(odor_c_bounds));
            

            this_Miss_op_activity=zeros(1,length(odor_c_bounds));
            this_Miss_op_activity_n=zeros(1,length(odor_c_bounds));
           



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

                    %Figure out odor concentration bin
                    switch handles_conc.group(fileNo)
                        case 1
                            %2 cm from the floor
                            cm_from_floor=2;
                        case 5
                            %1 cm from the floor
                            cm_from_floor=1;
                    end

                    %Find the analog op
                    minC=odor_plume_patterns.cm_from_floor(cm_from_floor).minC;
                    this_x=trials.trial(trNo).XYtest(ii_t,1);
                    this_y=trials.trial(trNo).XYtest(ii_t,2);

                    [M,this_x_op_ii]=min(abs(odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume-this_x));
                    [M,this_y_op_ii]=min(abs(odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume-this_y));

                    if trials.lane_per_trial(trNo)==1
                        op_conc=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1(this_y_op_ii,this_x_op_ii);
                    else
                        op_conc=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4(this_y_op_ii,this_x_op_ii);
                    end

                    %Where are odor start and end?
                    op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                    op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                    x_ii=op_predictedstart+ii_t;

                    %Make odor minC if the time is early or late
                    if x_ii<trials.odor_ii_start(trNo)
                        op_conc=minC;
                    else
                        if x_ii>trials.odor_ii_end(trNo)
                            op_conc=minC;
                        else
                            this_x=trials.trial(trNo).XYtest(ii_t,1);
                            x_on=(x_ii-trials.odor_ii_start(trNo))*handles_conc.dt*handles_conc.air_flow_speed;
                            if this_x>x_on
                                op_conc=minC;
                            end
                        end
                    end
                    delta_oc_bound=(odor_c_bounds(2)-odor_c_bounds(1))/2;
                    for ii_b=1:11
                        if (op_conc>=odor_c_bounds(ii_b)-delta_oc_bound)&(op_conc<odor_c_bounds(ii_b)+delta_oc_bound)
                            this_ii_op=ii_b;
                        end
                    end

                    %Get dFF
                    this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                    this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
                    sum_dFF_activity=sum_dFF_activity+these_dFF(ii_t);

                    this_dFF_op_activity(this_ii_op)=this_dFF_op_activity(this_ii_op)+these_dFF(ii_t);
                    this_dFF_op_activity_n(this_ii_op)=this_dFF_op_activity_n(this_ii_op)+1;

                    these_dFF_op_activity(trNo,this_ii_op)=these_dFF_op_activity(trNo,this_ii_op)+these_dFF(ii_t);
                    these_dFF_op_activity_n(trNo,this_ii_op)=these_dFF_op_activity_n(trNo,this_ii_op)+1;
                    

                    %Tally info
                    this_bin_dFF=(these_dFF(ii_t)>0)+1;
                    this_xy_ii=this_x_ii+10*(this_y_ii-1);
            
                    if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                        this_Hit_activity(this_x_ii,this_y_ii)=this_Hit_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        this_Hit_activity_n(this_x_ii,this_y_ii)=this_Hit_activity_n(this_x_ii,this_y_ii)+1;
                        sum_Hit_activity=sum_Hit_activity+these_dFF(ii_t);

                        this_Hit_op_activity(this_ii_op)=this_Hit_op_activity(this_ii_op)+these_dFF(ii_t);
                        this_Hit_op_activity_n(this_ii_op)=this_Hit_op_activity_n(this_ii_op)+1;
                      
                    else
                        this_Miss_activity(this_x_ii,this_y_ii)=this_Miss_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        this_Miss_activity_n(this_x_ii,this_y_ii)=this_Miss_activity_n(this_x_ii,this_y_ii)+1;
                        sum_Miss_activity=sum_Miss_activity+these_dFF(ii_t);

                        this_Miss_op_activity(this_ii_op)=this_Miss_op_activity(this_ii_op)+these_dFF(ii_t);
                        this_Miss_op_activity_n(this_ii_op)=this_Miss_op_activity_n(this_ii_op)+1;
                     
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
                    if this_Hit_activity_n(ii_x,ii_y)~=0
                        this_Hit_activity(ii_x,ii_y)=this_Hit_activity(ii_x,ii_y)/this_Hit_activity_n(ii_x,ii_y);
                    end
                end
            end

            for ii_x=1:10
                for ii_y=1:10
                    if this_Miss_activity_n(ii_x,ii_y)~=0
                        this_Miss_activity(ii_x,ii_y)=this_Miss_activity(ii_x,ii_y)/this_Miss_activity_n(ii_x,ii_y);
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

            %Plot dFFs for lanes 1 and 4
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

            %Plot dFF for lane1
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
            title(['Lane 1 odor'])

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
            title(['Lane 4 odor'])

            %Plot space activity maps as in Moser https://doi.org/10.1126/science.1114037
            max_this_dFF_activity=max(this_dFF_activity(:));
            max_this_Miss_activity=max(this_Miss_activity(:));
            max_this_Hit_activity=max(this_Hit_activity(:));

            colormap fire
            this_cmap=colormap;
            this_cmap(1,:)=[0.3 0.3 0.3];

            if max_this_Hit_activity>max_this_Miss_activity
                %Lane 1
                subplot(2, 6, [7 8]);
                max_activity=max_this_Hit_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_Hit_activity=this_Hit_activity;
                this_masked_Hit_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_Hit_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([70 100 200 300 400 430])
                yticklabels({'Lane 4','100','200','300','400','Lane 1'})

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')


                title_legend=['Hits, SSI= ' num2str(ssi_all_spatial_info_hit(ii_ROI_all))];
               
                title(title_legend)

                %Lane 4 normalized to lane 1
                subplot(2, 6, [9 10]);
                max_activity=max_this_Hit_activity;
                this_masked_Miss_activity=this_Miss_activity;
                this_masked_Miss_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_Miss_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Miss normalized to Hits'])

                %Lane 4 normalized to lane 4
                subplot(2, 6, [11 12]);
                max_activity=max_this_Miss_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_Miss_activity=this_Miss_activity;
                this_masked_Miss_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_Miss_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])

                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')


                title_legend=['Miss, SSI= ' num2str(ssi_all_spatial_info_miss(ii_ROI_all))];
              

                title(title_legend)
            else


                %Lane 4 normalized to lane 4
                subplot(2, 6, [7 8]);
                max_activity=max_this_Miss_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_Miss_activity=this_Miss_activity;
                this_masked_Miss_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_Miss_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                xticks(0:50:500)

                yticks([70 100 200 300 400 430])
                yticklabels({'Lane 4','100','200','300','400','Lane 1'})

                xlabel('x (mm)')
                ylabel('y (mm)')

                title_legend=['Miss, SSI= ' num2str(ssi_all_spatial_info_miss(ii_ROI_all))];
          
                title(title_legend)

                %Lane 1 normalized to lane 4
                subplot(2, 6, [9 10]);
                max_activity=max_this_Miss_activity;
                this_masked_Hit_activity=this_Hit_activity;
                this_masked_Hit_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_Hit_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Hits normalized to Miss'])

                %Lane 1 normalized to lane 1
                subplot(2, 6, [11 12]);
                max_activity=max_this_Hit_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.00001;
                end
                this_masked_Hit_activity=this_Hit_activity;
                this_masked_Hit_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_Hit_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks([100 200 300 400])
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')

                title_legend=['Hits, SSI= ' num2str(ssi_all_spatial_info_hit(ii_ROI_all))];
               
                title(title_legend)
            end


            sgt_legend=['dFF map file  No ' num2str(fileNo) ' ROI No ' num2str(ii_ROI) ' SSI both= ' num2str(ssi_all_spatial_info(ii_ROI_all))...
                ' rho= ' num2str(imps.all_spatial_rho_hit_vs_miss(ii_ROI_all)) ' dCOM= ' num2str(imps.all_delta_center_of_mass(ii_ROI_all)) ' mm '...
                'SSIop= ' num2str(ssi_all_op_bin_info(ii_ROI_all))];


       
            if ~isnan(not_nan_ii_ROI(ii_ROI_all))
                switch idxssi(not_nan_ii_ROI(ii_ROI_all))
                    case 1
                        sgt_legend=[sgt_legend ' cl1 '];
                    case 2
                        sgt_legend=[sgt_legend ' cl2 '];
                    case 3
                        sgt_legend=[sgt_legend ' cl3 '];
                end
            else
                sgt_legend=[sgt_legend ' nan '];
            end

            sgtitle(sgt_legend)



            %Now do glm
            tbl = table(glm_div.data',glm_div.trial_type',glm_div.time',...
                'VariableNames',{'dFF','trial_type','time'});
            mdl = fitglm(tbl,'dFF~trial_type+time'...
                ,'CategoricalVars',[2])

       



            %Normalize the dFF
            for ii_op=1:11
                 if this_dFF_op_activity_n(ii_op)~=0
                    this_dFF_op_activity(ii_op)=this_dFF_op_activity(ii_op)/this_dFF_op_activity_n(ii_op);
                end
            end

            for ii_op=1:11
                if this_Hit_op_activity_n(ii_op)~=0
                    this_Hit_op_activity(ii_op)=this_Hit_op_activity(ii_op)/this_Hit_op_activity_n(ii_op);
                end
            end

            for ii_op=1:11
                if this_Miss_op_activity_n(ii_op)~=0
                    this_Miss_op_activity(ii_op)=this_Miss_op_activity(ii_op)/this_Miss_op_activity_n(ii_op);
                end
            end

            % max_op_activity=max([this_Hit_op_activity(this_Hit_op_activity_n>0) this_Miss_op_activity(this_Miss_op_activity_n>0)]);
            % min_op_hit_activity=min([this_Hit_op_activity(this_Hit_op_activity_n>0) this_Miss_op_activity(this_Miss_op_activity_n>0)]);

            max_op_hit_activity=max(this_Hit_op_activity(this_Hit_op_activity_n>0));
            min_op_hit_activity=min(this_Hit_op_activity(this_Hit_op_activity_n>0));

            max_op_dFF_activity=max(this_dFF_op_activity(this_dFF_op_activity_n>0));
            min_op_dFF_activity=min(this_dFF_op_activity(this_dFF_op_activity_n>0));
            
            %Plot the odor activity maps
            try
                close(figureNo+1)
            catch
            end


            hFig = figure(figureNo+1);
            set(hFig, 'units','normalized','position',[.02 .3 .3 .4])

            subplot(3,1,1)
            hold on
            for ii_op=1:11
                if (this_Hit_op_activity_n(ii_op)>0)
                    bar(ii_op,(this_Hit_op_activity(ii_op)-min_op_hit_activity)/(max_op_hit_activity-min_op_hit_activity),'EdgeColor','none','FaceColor',[204/255 121/255 167/255])
                else
                    bar(ii_op,0.05,'EdgeColor','none','FaceColor',[150/255 150/255 150/255])
                end
            end
            xticks([1 6 11])
            ylim([0 1.2])
            xticklabels({'-10','-6','-2'})
            title('Hit')
            xlabel('Odor concentration')
            ylabel('dFF activity (AU)')

            subplot(3,1,2)
            hold on
            for ii_op=1:11
                if (this_Miss_op_activity_n(ii_op)>0)
                    bar(ii_op,(this_Miss_op_activity(ii_op)-min_op_hit_activity)/(max_op_hit_activity-min_op_hit_activity),'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
                else
                    bar(ii_op,0.05,'EdgeColor','none','FaceColor',[150/255 150/255 150/255])
                end
            end
            xticks([1 6 11])
            ylim([0 1.2])
            xticklabels({'-10','-6','-2'})
            title('Miss')
            xlabel('Odor concentration')
            ylabel('dFF activity (AU)')

              subplot(3,1,3)
            hold on
            for ii_op=1:11
                if (this_dFF_op_activity_n(ii_op)>0)
                    bar(ii_op,(this_dFF_op_activity(ii_op)-min_op_dFF_activity)/(max_op_dFF_activity-min_op_dFF_activity),'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
                else
                    bar(ii_op,0.05,'EdgeColor','none','FaceColor',[150/255 150/255 150/255])
                end
            end
            xticks([1 6 11])
            xticklabels({'-10','-6','-2'})
            ylim([0 1.2])
            title('All')
            xlabel('Odor concentration')
            ylabel('dFF activity (AU)')

            sgt_legend=['dFF map file  No ' num2str(fileNo) ' ROI No ' num2str(ii_ROI) ' SSI op= ' num2str(ssi_all_op_info(ii_ROI_all))...
                ' SSI op hit= ' num2str(ssi_all_op_info_hit(ii_ROI_all))...
                ' ssi op miss= ' num2str(ssi_all_op_info_miss(ii_ROI_all))];

            sgtitle(sgt_legend)

            if ssi_all_op_info_hit(ii_ROI_all)>=3
                oc_ii=oc_ii+1;

                oc_dFF_norm_map_hit(oc_ii,:)=(this_Hit_op_activity-min_op_hit_activity)/(max_op_hit_activity-min_op_hit_activity);
                oc_dFF_norm_map_hit(oc_ii,this_Hit_op_activity_n==0)=0;
                oc_dFF_norm_map_hit_n(oc_ii,:)=this_Hit_op_activity_n;

                oc_dFF_norm_map_miss(oc_ii,:)=(this_Miss_op_activity-min_op_hit_activity)/(max_op_hit_activity-min_op_hit_activity);
                oc_dFF_norm_map_miss(oc_ii,this_Miss_op_activity_n==0)=0;
                oc_dFF_norm_map_miss_n(oc_ii,:)=this_Miss_op_activity_n;

                oc_dFF_norm_map_fileNo(oc_ii)=fileNo;
                oc_dFF_norm_map_iiROI(oc_ii)=ii_ROI;
            end

            if ssi_all_op_info(ii_ROI_all)>=3
                oc_ii_all=oc_ii_all+1;

                oc_dFF_norm_map_all(oc_ii_all,:)=(this_dFF_op_activity-min_op_dFF_activity)/(max_op_dFF_activity-min_op_dFF_activity);
                oc_dFF_norm_map_all(oc_ii_all,this_dFF_op_activity_n==0)=0;
                oc_dFF_norm_map_all_n(oc_ii_all,:)=this_dFF_op_activity_n;

                oc_dFF_norm_map_fileNo_all(oc_ii_all)=fileNo;
                oc_dFF_norm_map_iiROI_all(oc_ii_all)=ii_ROI;
            end

            %Now plot the per trial pseudocolor plot
            all_these_dFF_op_activity=[];
            for trNo=1:trials.odor_trNo
                for ii_occ=1:11
                    if these_dFF_op_activity_n(trNo,ii_occ)~=0
                        these_dFF_op_activity(trNo,ii_occ)=these_dFF_op_activity(trNo,ii_occ)/these_dFF_op_activity_n(trNo,ii_occ);
                        all_these_dFF_op_activity=[all_these_dFF_op_activity these_dFF_op_activity(trNo,ii_occ)];
                    end
                end
            end
            max_these_dFF_op_activity=max(all_these_dFF_op_activity);
            min_these_dFF_op_activity=min(all_these_dFF_op_activity);
            these_dFF_op_activity=(these_dFF_op_activity-min_these_dFF_op_activity)/(max_these_dFF_op_activity-min_these_dFF_op_activity);
            for trNo=1:trials.odor_trNo
                for ii_occ=1:11
                    if these_dFF_op_activity_n(trNo,ii_occ)==0
                        these_dFF_op_activity(trNo,ii_occ)=-0.1;
                    end
                end
            end

            try
                close(figureNo+2)
            catch
            end


            hFig = figure(figureNo+2);
            set(hFig, 'units','normalized','position',[.35 .3 .3 .4])
            hold on

            drg_pcolor(repmat(odor_c_bounds',1,trials.odor_trNo),repmat([1:trials.odor_trNo],length(odor_c_bounds),1),these_dFF_op_activity')
            colormap(this_cmap)
            clim([0 1])


            shading flat


            ylim([1 trials.odor_trNo])
            xlim([-10 -2 ])
            xlabel('log(odor concentration) AU')
            ylabel('Trial #')
            title(['Neural activity vs. odor conc, all trials'])

            try
                close(figureNo+3)
            catch
            end


            hFig = figure(figureNo+3);
            set(hFig, 'units','normalized','position',[.35 .3 .3 .4])
            hold on

            hit_trials=logical((trials.hit1==1)|(trials.hit4==1));
            no_hits=sum( (trials.hit1==1)|(trials.hit4==1));
            drg_pcolor(repmat(odor_c_bounds',1,no_hits),repmat([1:no_hits],length(odor_c_bounds),1),these_dFF_op_activity(hit_trials,:)')
            colormap(this_cmap)
            clim([0 1])



            shading flat


            ylim([1 sum(hit_trials)])
            xlim([-10 -2 ])
            xlabel('log(odor concentration) AU')
            ylabel('Trial #')
            title(['Neural activity vs. odor conc, hit trials'])

            try
                close(figureNo+4)
            catch
            end


            hFig = figure(figureNo+4);
            set(hFig, 'units','normalized','position',[.67 .3 .3 .4])
            hold on

            miss_trials=logical((trials.miss1==1)|(trials.miss4==1));
            no_miss=sum( (trials.miss1==1)|(trials.miss4==1));
            drg_pcolor(repmat(odor_c_bounds',1,no_miss),repmat([1:no_miss],length(odor_c_bounds),1),these_dFF_op_activity(miss_trials,:)')
            colormap(this_cmap)
            clim([0 1])



            shading flat


            ylim([1 sum(miss_trials)])
            xlim([-10 -2 ])
            xlabel('log(odor concentration) AU')
            ylabel('Trial #')
            title(['Neural activity vs. odor conc, miss trials'])

             if (fileNo==7)
                pffft=1;
            end

            if (fileNo==7)&(ii_ROI==23)
                pffft=1;
            end

            if (fileNo==7)&(ii_ROI==99)
                pffft=1;
            end

            if (fileNo==7)&(ii_ROI==154)
                pffft=1;
            end

            if (fileNo==13)&(ii_ROI==92)
                pffft=1;
            end

            %For Figure 3
            if (fileNo==7)&(ii_ROI==121)
                pffft=1;
            end

            if (fileNo==7)&(ii_ROI==176)
                pffft=1;
            end

            if (fileNo==7)&(ii_ROI==180)
                pffft=1;
            end

            if (fileNo==13)&(ii_ROI==143)
                pffft=1;
            end
        end
    end
end

%Now plot the response of the odor concentration cells


%Report the results for hits and misses
%Find the peak for oc_dFF_norm_map_hit
ii_oc_peaks=zeros(1,oc_ii);
oc_dFF_norm_peaks=zeros(1,oc_ii);
odor_c_peak=zeros(1,oc_ii);
for ii_oc=1:oc_ii
    this_oc_dFF=zeros(1,11);
    this_oc_dFF(1,:)=oc_dFF_norm_map_hit(ii_oc,:);
    this_oc_dFF_n=zeros(1,11);
    this_oc_dFF_n(1,:)=oc_dFF_norm_map_hit_n(ii_oc,:);
    [oc_dFF_norm_peaks(ii_oc), ii_oc_peaks(ii_oc)]=max(this_oc_dFF);
    odor_c_peak(ii_oc)=odor_c_bounds(ii_oc_peaks(ii_oc));
end

%Plot the histogram for odor_c peaks
figureNo=figureNo+4;
try
    close(figureNo)
catch
end


hFig = figure(figureNo);
set(hFig, 'units','normalized','position',[.3 .3 .3 .4])

hold on
histogram(odor_c_peak,[-11:1:-1])
xlim([odor_c_bounds(1)-0.5*(odor_c_bounds(2)-odor_c_bounds(1)) odor_c_bounds(end)+0.5*(odor_c_bounds(2)-odor_c_bounds(1))])
xlabel('Log(odor concentraiton) AU')

%Plot cumulative histogram for spatial SSI vs oc SSI
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

[f_aic,x_aic] = drg_ecdf(ssi_all_spatial_info_hit);
plot(x_aic,f_aic,'Color',[204/255 121/255 167/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(ssi_all_spatial_info_hit(ssi_all_op_info_hit>=3));
plot(x_aic,f_aic,'Color',[150/255 150/255 150/255],'LineWidth',3)

text(6,0.4,'All ROIs','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
text(6,0.35,'Odor concentraion SSI>=3','Color',[150/255 150/255 150/255],'FontWeight','bold','FontSize',16)
xlabel('SSI spatial for hits')
ylabel('Cumulative prbability')

%Sort by peak ii
oc_dFF_norm_map_hit_sorted=zeros(oc_ii,11);
oc_dFF_norm_map_miss_sorted=zeros(oc_ii,11);

to_sort_oc=[ii_oc_peaks' [1:oc_ii]'];
sorted_oc=sortrows(to_sort_oc);

for ii_oc=1:oc_ii
    this_oc_dFF_norm_map_hit=zeros(1,11);
    this_oc_dFF_norm_map_hit(1,:)=oc_dFF_norm_map_hit(sorted_oc(ii_oc,2),:);
    this_oc_dFF_norm_map_hit_n=zeros(1,11);
    this_oc_dFF_norm_map_hit_n(1,:)=oc_dFF_norm_map_hit_n(sorted_oc(ii_oc,2),:);
    this_oc_dFF_norm_map_hit(1,oc_dFF_norm_map_hit_n(sorted_oc(ii_oc,2),:)==0)=-0.1;
    oc_dFF_norm_map_hit_sorted(ii_oc,:)=this_oc_dFF_norm_map_hit;

    this_oc_dFF_norm_map_miss=zeros(1,11);
    this_oc_dFF_norm_map_miss(1,:)=oc_dFF_norm_map_miss(sorted_oc(ii_oc,2),:);
    this_oc_dFF_norm_map_miss_n=zeros(1,11);
    this_oc_dFF_norm_map_miss_n(1,:)=oc_dFF_norm_map_miss_n(sorted_oc(ii_oc,2),:);
    this_oc_dFF_norm_map_miss(1,oc_dFF_norm_map_miss_n(sorted_oc(ii_oc,2),:)==0)=-0.1;
    oc_dFF_norm_map_miss_sorted(ii_oc,:)=this_oc_dFF_norm_map_miss;
end


%Odor concentration activity map for hits
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

colormap fire
this_cmap=colormap;
this_cmap(1,:)=[0.8 0.8 0.8];

drg_pcolor(repmat(odor_c_bounds,oc_ii,1),repmat([1:oc_ii]',1,11),oc_dFF_norm_map_hit_sorted)

colormap(this_cmap)
shading flat

%caxis([prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),1) prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),99.9)])

caxis([-0.01,max(oc_dFF_norm_map_hit_sorted(:))])



xlim([-10 -1.2])
ylim([1 oc_ii])
title('Odor concentration activity map, hits')
xlabel('log(odor concentrtion) AU')
ylabel('ROI number')

%Odor concentration activity map for miss
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

colormap fire
this_cmap=colormap;
this_cmap(1,:)=[0.8 0.8 0.8];

drg_pcolor(repmat(odor_c_bounds,oc_ii,1),repmat([1:oc_ii]',1,11),oc_dFF_norm_map_miss_sorted)

colormap(this_cmap)
shading flat

%caxis([prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),1) prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),99.9)])

caxis([-0.01,max(oc_dFF_norm_map_hit_sorted(:))])



xlim([-10 -1.2])
ylim([1 oc_ii])
title('Odor concentration activity map, miss')
xlabel('log(odor concentrtion) AU')
ylabel('ROI number')

%Now report the results for all trials

%Find the peak for oc_dFF_norm_map_all
ii_oc_peaks_all=zeros(1,oc_ii);
oc_dFF_norm_peaks_all=zeros(1,oc_ii);
odor_c_peak_all=zeros(1,oc_ii);
for ii_oc=1:oc_ii_all
    this_oc_dFF=zeros(1,11);
    this_oc_dFF(1,:)=oc_dFF_norm_map_all(ii_oc,:);
    this_oc_dFF_n=zeros(1,11);
    this_oc_dFF_n(1,:)=oc_dFF_norm_map_all_n(ii_oc,:);
    [oc_dFF_norm_peaks_all(ii_oc), ii_oc_peaks_all(ii_oc)]=max(this_oc_dFF);
    odor_c_peak_all(ii_oc)=odor_c_bounds(ii_oc_peaks_all(ii_oc));
end

%Plot the histogram for odor_c peaks for all trials
figureNo=figureNo+2;
try
    close(figureNo)
catch
end


hFig = figure(figureNo);
set(hFig, 'units','normalized','position',[.3 .3 .3 .4])

hold on
histogram(odor_c_peak_all,[-11:1:-1])
xlim([-10 -1])
xlabel('Log(odor concentraiton) AU')
title('Peak neural activity odor concentration for all trials')

%Plot cumulative histogram for spatial SSI vs oc SSI for all trials
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

[f_aic,x_aic] = drg_ecdf(ssi_all_spatial_info);
plot(x_aic,f_aic,'Color',[204/255 121/255 167/255],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(ssi_all_spatial_info(ssi_all_op_info>=3));
plot(x_aic,f_aic,'Color',[150/255 150/255 150/255],'LineWidth',3)

plot([3 3],[0 1],'-k','LineWidth',2)

text(6,0.4,'All ROIs','Color',[204/255 121/255 167/255],'FontWeight','bold','FontSize',16)
text(6,0.35,'Odor concentraion SSI>=3','Color',[150/255 150/255 150/255],'FontWeight','bold','FontSize',16)
xlabel('SSI spatial for all trials')
ylabel('Cumulative prbability')

%Sort by peak ii
oc_dFF_norm_map_all_sorted=zeros(oc_ii_all,11);

to_sort_oc_all=[ii_oc_peaks_all' [1:oc_ii_all]'];
sorted_oc_all=sortrows(to_sort_oc_all);

for ii_oc=1:oc_ii_all
    this_oc_dFF_norm_map_all=zeros(1,11);
    this_oc_dFF_norm_map_all(1,:)=oc_dFF_norm_map_all(sorted_oc_all(ii_oc,2),:);
    this_oc_dFF_norm_map_all_n=zeros(1,11);
    this_oc_dFF_norm_map_all_n(1,:)=oc_dFF_norm_map_all_n(sorted_oc_all(ii_oc,2),:);
    this_oc_dFF_norm_map_all(1,oc_dFF_norm_map_all_n(sorted_oc_all(ii_oc,2),:)==0)=-0.1;
    oc_dFF_norm_map_all_sorted(ii_oc,:)=this_oc_dFF_norm_map_all;
end


%Odor concentration activity map for all trials
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

colormap fire
this_cmap=colormap;
this_cmap(1,:)=[0.8 0.8 0.8];

drg_pcolor(repmat(odor_c_bounds,oc_ii_all,1),repmat([1:oc_ii_all]',1,11),oc_dFF_norm_map_all_sorted)

colormap(this_cmap)
shading flat

%caxis([prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),1) prctile(sorted_handles_out.ii_comp(ii_comp).all_div_hit1(:),99.9)])

caxis([-0.01,max(oc_dFF_norm_map_all_sorted(:))])



xlim([-10 -1.2])
ylim([1 oc_ii_all])
title('Odor concentration activity map, all trials')
xlabel('log(odor concentrtion) AU')
ylabel('ROI number')

%Plot rainbow
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


prain=[0:1/99:1];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
%             colormap jet
colormap(this_cmap)
shading interp
ax=gca;
set(ax,'XTickLabel','')

fclose(fileID);

pffft=1;

