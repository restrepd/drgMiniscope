%drgMini_analyze_batch_DecodeCorrXYConcSubROIs
close all
clear all

is_sphgpu=0;

%Use R2ER hat?
%Popsil and Bair, PLOS Computational Biology, 2021
%https://doi.org/10.1371/journal.pcbi.1009212
is_hat=0;

switch is_sphgpu
    case 0

        %There was a bug and the shuffled runs were the same for all shuffled runs
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';
        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        %For the files below the shuffled runs should be different
        %
        % %Trained with all trials
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01062025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01062025.m';

        %Trained with hits only using old code v3 that was lost
        %  save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

        %Trained with miss only using new drgMini_DecodeOdorConcv4.m code
        %  save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc02232025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_02232025.m'

        % %Trained with hit only using new drgMini_DecodeOdorConcv4.m code
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc02242025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_02242025.m'

        % %Trained with hits only predicting the odor cone for the entire path
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

        %Odor concentraiton decoding for paper
        %Trained with hits only taking on account when mouse detects the odor

        %Files for the different decoding algorithms
        
        %BINARY TREE
        
        %binary tree for paper before_bins=0

        %Odor decoing
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc04192024/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_04192024.m';

        %Trained with all trials and mult=0.5 taking on account when mouse detects the odor
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc11112025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_11112025.m';

        %Trained with all trials and mult=0.1 taking on account when mouse detects the odor
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc0p1_11212025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_0p1_11212025.m';

        %Trained with all trials and mult=1 taking on account when mouse detects the odor
        %Please note that .imp is too long because this has multiple time points
        %for the decoding
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc_1_11222025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_1_11222025.m';

          %Trained with all trials and mult=1 taking on account when mouse detects the odor
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc_1_01012026/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_1_01012026.m';

        %Trained with all trials and mult=1 taking on account when mouse detects the odor
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcSub_1_01082026/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_subROIs_1_01082026.m';

        %Trained with all trials and mult=1 taking on account when mouse detects the odor
        save_PathConc{1}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcSub_1_01122026/';
        choiceOdorConcFileName{1}='drgDynamicOdorConcChoices_Fabio_subROIs_1_01122026.m';

        save_PathConc{2}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcSub_1_01182026/';
        choiceOdorConcFileName{2}='DecodeDynOdorConcSub_1_01182026.m';

        %Trained with all trials and mult=0.3 taking on account when mouse detects the odor
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc_0p3_11223025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_0p3_11223025.m';

        %Trained with all trials and mult=0.6 taking on account when mouse detects the odor
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc_0p6_11224025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_0p6_11224025.m';

        % %Position decoding 
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        %Position decoding 
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01102926/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_SubROIs_01102026.m';

        %Position decoding 
        save_PathXY{1}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01132926/';
        choiceXYFileName{1}='drgOdorArenaChoices_Fabio_Good_SubROIs_01132026.m';

        save_PathXY{2}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01192026/';
        choiceXYFileName{2}='drgOdorArenaChoices_Fabio_Good_SubROIs_01192026.m';
        
        %binary tree before_bins=5

        %Odor decoding
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc06042025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_06042025.m'
        % 
        % %Position decoding 
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput_bin_b5_06052025/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_bin_b5_06052025.m';

        %SMV

        %before_bins=0
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcSVM06042025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_SVMb0_06052025.m';

        %before_bins=5
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcSVM0b06052025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_SVMb0_06052025.m'


        %ANN
        
        %ann for paper before_bins=0

        % %Odor decoing
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcANN0b06202025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_ANNb0_06202025.m';
        % 
        % %Position decoding, note: this is currenlty binary tree 
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';


        %GLM
        
        %glm for paper before_bins=0

        % %Odor decoing
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcAGLM0b06232025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_glmb0_06232025.m';
        % 
        % %Position decoding, note: this is currenlty binary tree 
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeBayesOdorConc05072025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Bayes_05072025.m';

        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';
        %
        % save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        % choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        

        %Angle processing
        % save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        % choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        %This one has the odor encounter
        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle05152025/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_05102025.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        saveFile='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/SubROI_analysis.mat';

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end


%6/27/2024 Exclude files with p_R1_hits(fileNo) > 0.05
thr_rho=1; %Updated for a final setting 10/25/25
% 
% R_thr=0.1;
% R2ER_thr=0.05;

%These are used to classify angle approach
low_angle=-130;
high_angle=-50;


%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4, 2 cm from floor
%Group 2 is rewarded, with odor lane 4, no odor in lane 1choiceAngleFileName
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
this_choiceOdorConcFileName=choiceOdorConcFileName{1};
eval(['handles_conc=' this_choiceOdorConcFileName(1:end-2) ';'])
this_choiceXYFileName=choiceXYFileName{1};
eval(['handles_XY=' this_choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])

figureNo=0;
%Exclude one file with high fraction_other_angle
%We will exclude fraction_other_angle>thr_froa
thr_froa=0.2;

%Calculate fraction of horizontal angle approach for each file
fraction_other_angle=[];
all_meanAngles=[];
for fileNo=1:length(handles_conc.arena_file)
    % if sum(handles_conc.group(fileNo)==these_groups)>0
    angle_file=handles_Angle.arena_file{fileNo};
    %load the ouptut file
    load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
    trials=handles_out.trials;

    meanAngles=[];
    for trNo=1:length(handles_out.angles.trial)
        if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
            meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
        else
            meanAngles(trNo)=0; %I need to fix this
        end
        all_meanAngles=[all_meanAngles meanAngles(trNo)];
    end

    no_hit90=0;
    no_hito=0;
    for trNo=1:trials.odor_trNo


        %Okabe_Ito colors
        switch trials.odor_trial_type(trNo)
            case 1
                %Hit 90 degrees
                if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                    no_hit90=no_hit90+1;
                end

                %Horizontal hit
                if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                    no_hito=no_hito+1;
                end

            case 3

                %Hit 90 degrees
                if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                    no_hit90=no_hit90+1;
                end

                %Horizontal hit
                if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                    no_hito=no_hito+1;
                end
        end

    end

    fraction_other_angle=[fraction_other_angle no_hito/(no_hito+no_hit90)];

    % end
end


%Now calcualte R1, R2 and R2ER and p values 
%Note that R2ER is an unbiased estimation of the fraction of variance explained
%calculated as in Pospisil and Blair PLOS Comp. Biol. 2021
ii_run=1;


r2ER_all_trials=[];



these_groups=[1:5];
handles_anal_conc=[];

for fileNo=1:length(handles_conc.arena_file)
    if sum(handles_conc.group(fileNo)==these_groups)>0
        %Load conc data
        arena_file=handles_conc.arena_file{fileNo};
        ii_sub_cum=0;
        ii_sub_cum_sh=0;
        load([save_PathConc{1} arena_file(1:end-4) handles_conc.save_tag{ii_run} '_subROI.mat'])
        %Please note that no_sub_ROIs must be the same between files 
        for ROIsubset=1:length(handles_choices.no_sub_ROIs)
            handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum=0;
            handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum_sh=0;
        end

        for conc_fileNo=1:length(save_PathConc)
            %load the ouptut file
            load([save_PathConc{conc_fileNo} arena_file(1:end-4) handles_conc.save_tag{ii_run} '_subROI.mat'])
            trials=handles_out.trials;
            odor_plume_template=handles_out.odor_plume_template;
            for ROIsubset=1:length(handles_choices.no_sub_ROIs)
                for ii_sub=1:handles_choices.no_runs__per_sub_ROI(ROIsubset)
                    op_all_trials=[];
                    op_decod_all_trials=[];



                    % %Load angle file
                    % angle_file=handles_Angle.arena_file{fileNo};
                    % %load the ouptut file
                    % load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
                    % meanAngles=[];
                    % for trNo=1:length(handles_out.angles.trial)
                    %     if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
                    %         meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
                    %     else
                    %         meanAngles(trNo)=0; %I need to fix this
                    %     end
                    % end


                    op_predicted=handles_out.ROIsubset(ROIsubset).ii_sub(ii_sub).op_predicted;

                    no_conv_points=11;
                    % conv_win=ones(1,no_conv_points)/no_conv_points;
                    conv_win_gauss = gausswin(no_conv_points);
                    conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

                    op_predicted_conv=conv(op_predicted,conv_win_gauss,'same');

                    %Now limit the x and y to max and min
                    minop=min(odor_plume_template);
                    op_predicted_conv(op_predicted_conv<minop)=minop;
                    maxop=max(odor_plume_template);
                    op_predicted_conv(op_predicted_conv>maxop)=maxop;


                    last_op_predictedend=1;
                    ii_start=0;

                    % for trNo=1:trials.odor_trNo
                    % 
                    %     op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                    %     op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                    %     if op_predictedend>length(odor_plume_template)
                    %         op_predictedend=length(odor_plume_template);
                    %     end
                    % 
                    %     op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % 
                    % 
                    %     last_op_predictedend=op_predictedend;
                    % 
                    %     % ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;
                    %     % ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
                    % end

                    nTrials = trials.odor_trNo;

                    % Preallocate cell outputs (variable segment lengths)[web:16][web:19]
                    op_all_trials_cell       = cell(nTrials,1);
                    op_decod_all_trials_cell = cell(nTrials,1);
                    last_op_predictedend_vec = zeros(nTrials,1);  % if you still need this per trial

                    parfor trNo = 1:nTrials
                        op_predictedstart = trials.odor_ii_start(trNo) + handles_choices.trial_start_offset;
                        op_predictedend   = trials.odor_ii_end(trNo)   + handles_choices.trial_end_offset;
                        if op_predictedend > length(odor_plume_template)
                            op_predictedend = length(odor_plume_template);
                        end

                        seg_op  = odor_plume_template(op_predictedstart:op_predictedend).';
                        seg_dec = op_predicted_conv(op_predictedstart:op_predictedend);

                        % Store into sliced cell variables[web:23][web:28]
                        op_all_trials_cell{trNo}       = seg_op;
                        op_decod_all_trials_cell{trNo} = seg_dec;

                        last_op_predictedend_vec(trNo) = op_predictedend;
                    end

                    % Concatenate after parfor[web:16][web:29]
                    op_all_trials       = vertcat(op_all_trials_cell{:});
                    op_decod_all_trials = vertcat(op_decod_all_trials_cell{:});

                    %Calculate r2, the fraction of variance explained
                    handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum=handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum+1;
                    ii_sub_cum=handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum;
                    if ~isempty(op_all_trials)
                        % r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
                        [this_r2ER,this_P_r2ER]=drgMini_r2ER_bootstrap(op_all_trials,op_decod_all_trials,10000,is_hat);
                        handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).r2ER_all_trials=this_r2ER;
                        handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).p_r2ER_all=this_P_r2ER;
                    else
                        handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).r2ER_all_trials=NaN;
                        handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).p_r2ER_all=NaN;
                    end

                    %Calculate r2, the fraction of variance explained for shifted traces
                    if isfield(handles_out.ROIsubset(ROIsubset).ii_sub(ii_sub),'op_predicted_sh')==1
                        op_all_trials=[];
                        op_decod_all_trials=[];


                        op_predicted=handles_out.ROIsubset(ROIsubset).ii_sub(ii_sub).op_predicted_sh;

                        no_conv_points=11;
                        % conv_win=ones(1,no_conv_points)/no_conv_points;
                        conv_win_gauss = gausswin(no_conv_points);
                        conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

                        op_predicted_conv=conv(op_predicted,conv_win_gauss,'same');

                        %Now limit the x and y to max and min
                        minop=min(odor_plume_template);
                        op_predicted_conv(op_predicted_conv<minop)=minop;
                        maxop=max(odor_plume_template);
                        op_predicted_conv(op_predicted_conv>maxop)=maxop;


                        last_op_predictedend=1;
                        ii_start=0;

                        % for trNo=1:trials.odor_trNo
                        %
                        %     op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                        %     op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                        %     if op_predictedend>length(odor_plume_template)
                        %         op_predictedend=length(odor_plume_template);
                        %     end
                        %
                        %     op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
                        %     op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];
                        %
                        %
                        %     last_op_predictedend=op_predictedend;
                        %
                        %     % ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;
                        %     % ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
                        % end

                        nTrials = trials.odor_trNo;

                        % Preallocate cell outputs (variable segment lengths)[web:16][web:19]
                        op_all_trials_cell       = cell(nTrials,1);
                        op_decod_all_trials_cell = cell(nTrials,1);
                        last_op_predictedend_vec = zeros(nTrials,1);  % if you still need this per trial

                        parfor trNo = 1:nTrials
                            op_predictedstart = trials.odor_ii_start(trNo) + handles_choices.trial_start_offset;
                            op_predictedend   = trials.odor_ii_end(trNo)   + handles_choices.trial_end_offset;
                            if op_predictedend > length(odor_plume_template)
                                op_predictedend = length(odor_plume_template);
                            end

                            seg_op  = odor_plume_template(op_predictedstart:op_predictedend).';
                            seg_dec = op_predicted_conv(op_predictedstart:op_predictedend);

                            % Store into sliced cell variables[web:23][web:28]
                            op_all_trials_cell{trNo}       = seg_op;
                            op_decod_all_trials_cell{trNo} = seg_dec;

                            last_op_predictedend_vec(trNo) = op_predictedend;
                        end

                        % Concatenate after parfor[web:16][web:29]
                        op_all_trials       = vertcat(op_all_trials_cell{:});
                        op_decod_all_trials = vertcat(op_decod_all_trials_cell{:});

                        %Calculate R2, the fraction of variance explained
                        handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum_sh=handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum_sh+1;
                        ii_sub_cum_sh=handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum_sh;
                        if ~isempty(op_all_trials)
                            % r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
                            [this_r2ER,this_P_r2ER]=drgMini_r2ER_bootstrap(op_all_trials,op_decod_all_trials,10000,is_hat);
                            handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).r2ER_all_trials_sh=this_r2ER;
                            handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).p_r2ER_all_sh=this_P_r2ER;
                        else
                            handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).r2ER_all_trials_sh=NaN;
                            handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).p_r2ER_all_sh=NaN;
                        end
                    end

                end
            end
        end
        fprintf(1, ['\nr2ER odor concentration processing done file No ' num2str(fileNo) '\n']);
    end
end

%Now plot r2ER vs the number of ROIs for odor runs
figureNo = figureNo + 1;
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

glm_r2ER=[];
glm_r2ER_ii=0;

id_r2ER_ii=0;
input_r2ER_data=[];


these_groups=[1 5]; %2 cm and 1 cm all trials
these_r2ER_means=[];
these_r2ERs_conc=[];
these_r2ER_sh_means=[];
these_r2ER_sh_conc=[];

for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    %First get the r2ERs for each number of ROIs
    these_r2ERs=[];
    for ii_sub=1:handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ERs=[these_r2ERs handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_all_trials];
            end
        end
    end

    these_r2ERs_conc=[these_r2ERs_conc  these_r2ERs];
    these_r2ER_means=[these_r2ER_means mean(these_r2ERs)];
    

    %Then get the r2ERs_sh for each number of ROIs
    these_r2ER_sh=[];
    for ii_sub=1:handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum_sh
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ER_sh=[these_r2ER_sh handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_all_trials_sh];
            end
        end
    end
    
    
    these_r2ER_sh_conc=[these_r2ER_sh_conc  these_r2ER_sh];
    these_r2ER_sh_means=[these_r2ER_sh_means mean(these_r2ER_sh)];
    

    %plot bar
    bar(bar_offset,mean(these_r2ERs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_x]=drgViolinPoint(these_r2ERs...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;



    %plot bar
    bar(bar_offset,mean(these_r2ER_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_x]=drgViolinPoint(these_r2ER_sh...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=these_r2ERs;
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=ROIsubset*ones(1,length(these_r2ERs));
    glm_r2ER.shuffled(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=0*ones(1,length(these_r2ERs));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ERs);

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ERs;
    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=num2str(handles_choices.no_sub_ROIs(ROIsubset));
        ROIsubset_labels{ROIsubset}=num2str(handles_choices.no_sub_ROIs(ROIsubset));
    else
        input_r2ER_data(id_r2ER_ii).description='All';
        ROIsubset_labels{ROIsubset}='All';
    end

      glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=these_r2ER_sh;
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=ROIsubset*ones(1,length(these_r2ER_sh));
    glm_r2ER.shuffled(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=1*ones(1,length(these_r2ER_sh));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ER_sh);

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ERs;
    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=[num2str(handles_choices.no_sub_ROIs(ROIsubset)) ' sh'];
        ROIsubset_labels{ROIsubset}=num2str(handles_choices.no_sub_ROIs(ROIsubset));
    else
        input_r2ER_data(id_r2ER_ii).description='sh';
        ROIsubset_labels{ROIsubset}='sh';
    end

end

bar_offsets=[0:2:bar_offset-2];
plot(bar_offsets,these_r2ER_means,'-o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerFaceColor',[0.7 0.7 0.7])

xlim([-0.5 bar_offset-0.5])
xticklabels(ROIsubset_labels);   % tick labels from your cell array
xlabel('No ROIs')
ylabel('r2ER for odor concentration')
title('Decoding with different subsets of ROIs')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for r2ER vs ROI number\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for r2ER vs ROI number\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_r2ER.data',glm_r2ER.ROIsubset',glm_r2ER.shuffled',...
    'VariableNames',{'r2ES','ROI_number','shuffled'});
mdl = fitglm(tbl,'r2ES~ROI_number+shuffled'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for r2ER vs ROI number\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for r2ER vs ROI number\n']);

 
[output_data] = drgMutiRanksumorTtest(input_r2ER_data, fileID,0);

%Now exclude all the r2ERs below 3xSD of shuffled mean r2ER
figureNo = figureNo + 1;
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

glm_r2ER=[];
glm_r2ER_ii=0;

id_r2ER_ii=0;
input_r2ER_data=[];


these_groups=[1 5]; %2 cm and 1 cm all trials
these_r2ER_means=[];
these_r2ERs_conc=[];
these_r2ER_sh_means=[];
these_r2ER_sh_conc=[];
these_SD_significant_r2ERs_conc=[];
for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    %First get the r2ERs for each number of ROIs
    these_r2ERs=[];
    for ii_sub=1:handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ERs=[these_r2ERs handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_all_trials];
            end
        end
    end

    these_r2ERs_conc=[these_r2ERs_conc  these_r2ERs];
    these_r2ER_means=[these_r2ER_means mean(these_r2ERs)];
    

    %Then get the r2ERs_sh for each number of ROIs
    these_r2ER_sh=[];
    for ii_sub=1:handles_anal_conc.ROIsubset(ROIsubset).ii_sub_cum_sh
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ER_sh=[these_r2ER_sh handles_anal_conc.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_all_trials_sh];
            end
        end
    end
    
    
    these_r2ER_sh_conc=[these_r2ER_sh_conc  these_r2ER_sh];
    these_r2ER_sh_means=[these_r2ER_sh_means mean(these_r2ER_sh)];
    
    %Estimate significant r2ERs
    SD_these_r2ERs_sh=std(these_r2ER_sh);
    mean_these_r2ERs_sh=mean(these_r2ER_sh);
    these_SD_significant_r2ERs_conc=[these_SD_significant_r2ERs_conc these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh];

    %plot bar
    bar(bar_offset,mean(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_x]=drgViolinPoint(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;


    glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=these_r2ERs;
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=ROIsubset*ones(1,length(these_r2ERs));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ERs);

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ERs;
    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=num2str(handles_choices.no_sub_ROIs(ROIsubset));
        ROIsubset_labels{ROIsubset}=num2str(handles_choices.no_sub_ROIs(ROIsubset));
    else
        input_r2ER_data(id_r2ER_ii).description='All';
        ROIsubset_labels{ROIsubset}='All';
    end

end

bar_offsets=[0:bar_offset-1];
plot(bar_offsets,these_r2ER_means,'-o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerFaceColor',[0.7 0.7 0.7])

xlim([-0.5 bar_offset-0.5])
xticklabels(ROIsubset_labels);   % tick labels from your cell array
xlabel('No ROIs')
ylabel('Significant r2ER for odor concentration')
title('Decoding with different subsets of ROIs')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for r2ER vs ROI number, significant r2ERs\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for r2ER vs ROI number, significant r2ERs\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_r2ER.data',glm_r2ER.ROIsubset',...
    'VariableNames',{'r2ES','ROI_number'});
mdl = fitglm(tbl,'r2ES~ROI_number'...
    ,'CategoricalVars',[2])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for r2ER vs ROI number\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for r2ER vs ROI number\n']);

 
[output_data] = drgMutiRanksumorTtest(input_r2ER_data, fileID,0);

pffft=1;


%Now calculate r2ER for XY 
these_groups=[1:5];

ii_run=1;

r2ER_x_all_trials=[];
r2ER_y_all_trials=[];
handles_anal_XY=[];

for fileNo=1:length(handles_XY.arena_file)
    if sum(handles_XY.group(fileNo)==these_groups)>0

        %Load XY navigation data
        arena_file=handles_XY.arena_file{fileNo};

        ii_sub_XY=0;
        ii_sub_XY_sh=0;
        load([save_PathXY{1} arena_file(1:end-4) handles_XY.save_tag{ii_run} '_subROI.mat'])
        %Please note that no_sub_ROIs must be the same between files
        for ROIsubset=1:length(handles_choices.no_sub_ROIs)
            handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum=0;
            handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh=0;
        end

        for XY_fileNo=1:length(save_PathXY)
            %load the ouptut file
            load([save_PathXY{XY_fileNo} arena_file(1:end-4) handles_XY.save_tag{ii_run} '_subROI.mat'])
            trials=handles_outXY.trials;
            XYtest=handles_outXY.XYtest;

            for ROIsubset=1:length(handles_choices.no_sub_ROIs)
                for ii_sub=1:handles_choices.no_runs__per_sub_ROI(ROIsubset)
                    x_all_trials=[];
                    y_all_trials=[];

                    x_decod_all_trials=[];
                    y_decod_all_trials=[];

                    % %Load angle file
                    % angle_file=handles_Angle.arena_file{fileNo};
                    % %load the ouptut file
                    % load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
                    % meanAngles=[];
                    % for trNo=1:length(handles_out.angles.trial)
                    %     if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
                    %         meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
                    %     else
                    %         meanAngles(trNo)=0; %I need to fix this
                    %     end
                    % end

                    % %Load conc data
                    % arena_file=handles_XY.arena_file{fileNo};
                    % %load the ouptut file
                    % load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
                    % trials=handles_out.trials;
                    x_predicted=handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).x_predicted;
                    y_predicted=handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).y_predicted;

                    no_conv_points=11;
                    % conv_win=ones(1,no_conv_points)/no_conv_points;
                    conv_win_gauss = gausswin(no_conv_points);
                    conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

                    x_predicted_conv=conv(x_predicted,conv_win_gauss,'same');
                    y_predicted_conv=conv(y_predicted,conv_win_gauss,'same');

                    %Now limit the x and y to max and min
                    minY1=min(XYtest(:,1));
                    x_predicted_conv(x_predicted_conv<minY1)=minY1;
                    maxY1=max(XYtest(:,1));
                    x_predicted_conv(x_predicted_conv>maxY1)=maxY1;

                    minY2=min(XYtest(:,2));
                    y_predicted_conv(y_predicted_conv<minY2)=minY2;
                    maxY2=max(XYtest(:,2));
                    y_predicted_conv(y_predicted_conv>maxY2)=maxY2;




                    % for trNo=1:trials.odor_trNo
                    % 
                    %     ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                    %     ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                    %     if ii_predictedend>size(XYtest,1)
                    %         ii_predictedend=size(XYtest,1);
                    %     end
                    % 
                    %     x_all_trials=[x_all_trials; XYtest(ii_predictedstart:ii_predictedend,1)];
                    %     x_decod_all_trials=[x_decod_all_trials; x_predicted_conv(ii_predictedstart:ii_predictedend)];
                    % 
                    %     y_all_trials=[y_all_trials; XYtest(ii_predictedstart:ii_predictedend,2)];
                    %     y_decod_all_trials=[y_decod_all_trials; y_predicted_conv(ii_predictedstart:ii_predictedend)];
                    % 
                    % 
                    % end


                    nTrials=trials.odor_trNo;

                    % Per-trial storage (cells because lengths may differ)[web:4]
                    x_all_trials_cell       = cell(nTrials,1);
                    x_decod_all_trials_cell = cell(nTrials,1);
                    y_all_trials_cell       = cell(nTrials,1);
                    y_decod_all_trials_cell = cell(nTrials,1);

                    parfor trNo = 1:nTrials
                        ii_predictedstart = trials.odor_ii_start(trNo) + handles_choices.trial_start_offset;
                        ii_predictedend   = trials.odor_ii_end(trNo)   + handles_choices.trial_end_offset;
                        ii_predictedend   = min(ii_predictedend, size(XYtest,1));

                        x_all_trials_cell{trNo}       = XYtest(ii_predictedstart:ii_predictedend, 1);
                        x_decod_all_trials_cell{trNo} = x_predicted_conv(ii_predictedstart:ii_predictedend);
                        y_all_trials_cell{trNo}       = XYtest(ii_predictedstart:ii_predictedend, 2);
                        y_decod_all_trials_cell{trNo} = y_predicted_conv(ii_predictedstart:ii_predictedend);
                    end

                    % Concatenate after parfor[web:4]
                    x_all_trials       = vertcat(x_all_trials_cell{:});
                    x_decod_all_trials = vertcat(x_decod_all_trials_cell{:});
                    y_all_trials       = vertcat(y_all_trials_cell{:});
                    y_decod_all_trials = vertcat(y_decod_all_trials_cell{:});



                    %Now calculate R2ER for position
                    %Calculate R2, the fraction of variance explained

                    handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum=handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum+1;
                    ii_sub_cum=handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum;

                    %x
                    if ~isempty(x_all_trials)
                        % r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
                        [this_r2ER,this_P_r2ER]=drgMini_r2ER_bootstrap(x_all_trials,x_decod_all_trials,10000,is_hat);
                        handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).r2ER_x_all_trials=this_r2ER;
                        % p_r2ER_all(fileNo)=this_P_r2ER;
                    else
                        handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).r2ER_x_all_trials=NaN;
                    end


                    %y
                    if ~isempty(y_all_trials)
                        % r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
                        [this_r2ER,this_P_r2ER]=drgMini_r2ER_bootstrap(y_all_trials,y_decod_all_trials,10000,is_hat);
                        handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).r2ER_y_all_trials=this_r2ER;
                        % p_r2ER_all(fileNo)=this_P_r2ER;
                    else
                        handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum).file(fileNo).r2ER_y_all_trials=NaN;
                    end

                    %Calculate r2, the fraction of variance explained for shifted traces
                    if isfield(handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub),'x_predicted_sh')==1
                        x_all_trials=[];
                        y_all_trials=[];

                        x_decod_all_trials=[];
                        y_decod_all_trials=[];

                        % %Load angle file
                        % angle_file=handles_Angle.arena_file{fileNo};
                        % %load the ouptut file
                        % load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
                        % meanAngles=[];
                        % for trNo=1:length(handles_out.angles.trial)
                        %     if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
                        %         meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
                        %     else
                        %         meanAngles(trNo)=0; %I need to fix this
                        %     end
                        % end

                        % %Load conc data
                        % arena_file=handles_XY.arena_file{fileNo};
                        % %load the ouptut file
                        % load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
                        % trials=handles_out.trials;
                        x_predicted=handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).x_predicted_sh;
                        y_predicted=handles_outXY.ROIsubset(ROIsubset).ii_sub(ii_sub).y_predicted_sh;

                        no_conv_points=11;
                        % conv_win=ones(1,no_conv_points)/no_conv_points;
                        conv_win_gauss = gausswin(no_conv_points);
                        conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

                        x_predicted_conv=conv(x_predicted,conv_win_gauss,'same');
                        y_predicted_conv=conv(y_predicted,conv_win_gauss,'same');

                        %Now limit the x and y to max and min
                        minY1=min(XYtest(:,1));
                        x_predicted_conv(x_predicted_conv<minY1)=minY1;
                        maxY1=max(XYtest(:,1));
                        x_predicted_conv(x_predicted_conv>maxY1)=maxY1;

                        minY2=min(XYtest(:,2));
                        y_predicted_conv(y_predicted_conv<minY2)=minY2;
                        maxY2=max(XYtest(:,2));
                        y_predicted_conv(y_predicted_conv>maxY2)=maxY2;



                        % 
                        % for trNo=1:trials.odor_trNo
                        % 
                        %     ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                        %     ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                        %     if ii_predictedend>size(XYtest,1)
                        %         ii_predictedend=size(XYtest,1);
                        %     end
                        % 
                        %     x_all_trials=[x_all_trials; XYtest(ii_predictedstart:ii_predictedend,1)];
                        %     x_decod_all_trials=[x_decod_all_trials; x_predicted_conv(ii_predictedstart:ii_predictedend)];
                        % 
                        %     y_all_trials=[y_all_trials; XYtest(ii_predictedstart:ii_predictedend,2)];
                        %     y_decod_all_trials=[y_decod_all_trials; y_predicted_conv(ii_predictedstart:ii_predictedend)];
                        % 
                        % 
                        % end

                        nTrials=trials.odor_trNo;

                        % Per-trial storage (cells because lengths may differ)[web:4]
                        x_all_trials_cell       = cell(nTrials,1);
                        x_decod_all_trials_cell = cell(nTrials,1);
                        y_all_trials_cell       = cell(nTrials,1);
                        y_decod_all_trials_cell = cell(nTrials,1);

                        parfor trNo = 1:nTrials
                            ii_predictedstart = trials.odor_ii_start(trNo) + handles_choices.trial_start_offset;
                            ii_predictedend   = trials.odor_ii_end(trNo)   + handles_choices.trial_end_offset;
                            ii_predictedend   = min(ii_predictedend, size(XYtest,1));

                            x_all_trials_cell{trNo}       = XYtest(ii_predictedstart:ii_predictedend, 1);
                            x_decod_all_trials_cell{trNo} = x_predicted_conv(ii_predictedstart:ii_predictedend);
                            y_all_trials_cell{trNo}       = XYtest(ii_predictedstart:ii_predictedend, 2);
                            y_decod_all_trials_cell{trNo} = y_predicted_conv(ii_predictedstart:ii_predictedend);
                        end

                        % Concatenate after parfor[web:4]
                        x_all_trials       = vertcat(x_all_trials_cell{:});
                        x_decod_all_trials = vertcat(x_decod_all_trials_cell{:});
                        y_all_trials       = vertcat(y_all_trials_cell{:});
                        y_decod_all_trials = vertcat(y_decod_all_trials_cell{:});



                        %Now calculate R2ER for position
                        %Calculate R2, the fraction of variance explained

                        handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh=handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh+1;
                        ii_sub_cum_sh=handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh;

                        %x
                        if ~isempty(x_all_trials)
                            % r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
                            [this_r2ER,this_P_r2ER]=drgMini_r2ER_bootstrap(x_all_trials,x_decod_all_trials,10000,is_hat);
                            handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).r2ER_x_all_trials_sh=this_r2ER;
                            % p_r2ER_all(fileNo)=this_P_r2ER;
                        else
                            handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).r2ER_x_all_trials_sh=NaN;
                        end


                        %y
                        if ~isempty(y_all_trials)
                            % r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
                            [this_r2ER,this_P_r2ER]=drgMini_r2ER_bootstrap(y_all_trials,y_decod_all_trials,10000,is_hat);
                            handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).r2ER_y_all_trials_sh=this_r2ER;
                            % p_r2ER_all(fileNo)=this_P_r2ER;
                        else
                            handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub_cum_sh).file(fileNo).r2ER_y_all_trials_sh=NaN;
                        end

                
                    end
                end
            end
        end
        fprintf(1, ['\nr2ER odor XY processing done file No ' num2str(fileNo) '\n']);
    end
end



%Now plot r2ER vs the number of ROIs for x, including shuffled
figureNo = figureNo + 1;
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

glm_r2ER=[];
glm_r2ER_ii=0;

id_r2ER_ii=0;
input_r2ER_data=[];


these_groups=[1 5]; %2 cm and 1 cm all trials
these_r2ER_x_means=[];
these_r2ERs_x=[];
these_r2ER_x_sh_means=[];
these_r2ER_x_sh=[];
for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    %First get the r2ERs for each number of ROIs
    these_r2ERs=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ERs=[these_r2ERs handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_x_all_trials];
            end
        end
    end
    
    these_r2ERs_x=[these_r2ERs_x  these_r2ERs];
    these_r2ER_x_means=[these_r2ER_x_means mean(these_r2ERs)];

    

    %Then get the r2ERs_sh for each number of ROIs
    these_r2ER_sh=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ER_sh=[these_r2ER_sh handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_x_all_trials_sh];
            end
        end
    end
    
    
    these_r2ER_x_sh=[these_r2ER_x_sh  these_r2ER_sh];
    these_r2ER_x_sh_means=[these_r2ER_x_sh_means mean(these_r2ER_sh)];
    

    %plot bar
    bar(bar_offset,mean(these_r2ERs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_x]=drgViolinPoint(these_r2ERs...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    %plot bar
    bar(bar_offset,mean(these_r2ER_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_x]=drgViolinPoint(these_r2ER_sh...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    %r2ERs_x
    glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=these_r2ERs;
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=ROIsubset*ones(1,length(these_r2ERs));
    glm_r2ER.shuffled(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=0*ones(1,length(these_r2ERs));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ERs);

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ERs;

    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=num2str(handles_choices.no_sub_ROIs(ROIsubset));
        ROIsubset_labels{ROIsubset}=num2str(handles_choices.no_sub_ROIs(ROIsubset));
    else
        input_r2ER_data(id_r2ER_ii).description='All';
        ROIsubset_labels{ROIsubset}='All';
    end

     %r2ER_x_sh
     glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=these_r2ER_sh;
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=ROIsubset*ones(1,length(these_r2ER_sh));
    glm_r2ER.shuffled(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=1*ones(1,length(these_r2ER_sh));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ER_sh);

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ER_sh;

    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=[num2str(handles_choices.no_sub_ROIs(ROIsubset)) ' sh'];
    else
        input_r2ER_data(id_r2ER_ii).description='All sh';
    end

end

bar_offsets=[0:2:bar_offset-2];
plot(bar_offsets,these_r2ER_means,'-o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerFaceColor',[0.7 0.7 0.7])

xlim([-0.5 bar_offset-0.5])
xticklabels(ROIsubset_labels);   % tick labels from your cell array
xlabel('No ROIs')
ylabel('r2ER for x')
title('Decoding x with different subsets of ROIs')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for r2ER x vs ROI number\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for r2ER x vs ROI number\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_r2ER.data',glm_r2ER.ROIsubset',glm_r2ER.shuffled',...
    'VariableNames',{'r2ES','ROI_number','shuffled'});
mdl = fitglm(tbl,'r2ES~ROI_number+shuffled'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for r2ER x vs ROI number\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for r2ER x vs ROI number\n']);

 
[output_data] = drgMutiRanksumorTtest(input_r2ER_data, fileID,0);


%Now plot r2ER vs the number of ROIs for y, including shuffled
figureNo = figureNo + 1;
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

glm_r2ER=[];
glm_r2ER_ii=0;

id_r2ER_ii=0;
input_r2ER_data=[];


these_groups=[1 5]; %2 cm and 1 cm all trials
these_r2ER_y_means=[];
these_r2ERs_y=[];
these_r2ER_y_sh_means=[];
these_r2ER_y_sh=[];
for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    %First get the r2ERs for each number of ROIs
    these_r2ERs=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ERs=[these_r2ERs handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_y_all_trials];
            end
        end
    end
    
    these_r2ERs_y=[these_r2ERs_y  these_r2ERs];
    these_r2ER_y_means=[these_r2ER_y_means mean(these_r2ERs)];

    

    %Then get the r2ERs_sh for each number of ROIs
    these_r2ER_sh=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ER_sh=[these_r2ER_sh handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_y_all_trials_sh];
            end
        end
    end
    
    
    these_r2ER_y_sh=[these_r2ER_y_sh  these_r2ER_sh];
    these_r2ER_y_sh_means=[these_r2ER_y_sh_means mean(these_r2ER_sh)];
    

    %plot bar
    bar(bar_offset,mean(these_r2ERs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_y]=drgViolinPoint(these_r2ERs...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    %plot bar
    bar(bar_offset,mean(these_r2ER_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_y]=drgViolinPoint(these_r2ER_sh...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    %r2ERs_y
    glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=these_r2ERs;
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=ROIsubset*ones(1,length(these_r2ERs));
    glm_r2ER.shuffled(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs))=0*ones(1,length(these_r2ERs));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ERs);

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ERs;

    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=num2str(handles_choices.no_sub_ROIs(ROIsubset));
        ROIsubset_labels{ROIsubset}=num2str(handles_choices.no_sub_ROIs(ROIsubset));
    else
        input_r2ER_data(id_r2ER_ii).description='All';
        ROIsubset_labels{ROIsubset}='All';
    end

     %r2ER_y_sh
     glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=these_r2ER_sh;
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=ROIsubset*ones(1,length(these_r2ER_sh));
    glm_r2ER.shuffled(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ER_sh))=1*ones(1,length(these_r2ER_sh));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ER_sh);

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ER_sh;

    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=[num2str(handles_choices.no_sub_ROIs(ROIsubset)) ' sh'];
    else
        input_r2ER_data(id_r2ER_ii).description='All sh';
    end

end

bar_offsets=[0:2:bar_offset-2];
plot(bar_offsets,these_r2ER_means,'-o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerFaceColor',[0.7 0.7 0.7])

xlim([-0.5 bar_offset-0.5])
xticklabels(ROIsubset_labels);   % tick labels from your cell array
xlabel('No ROIs')
ylabel('r2ER for y')
title('Decoding y with different subsets of ROIs')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for r2ER y vs ROI number\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for r2ER y vs ROI number\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_r2ER.data',glm_r2ER.ROIsubset',glm_r2ER.shuffled',...
    'VariableNames',{'r2ES','ROI_number','shuffled'});
mdl = fitglm(tbl,'r2ES~ROI_number+shuffled'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for r2ER x vs ROI number\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for r2ER x vs ROI number\n']);

 
[output_data] = drgMutiRanksumorTtest(input_r2ER_data, fileID,0);


%Now plot r2ER vs the number of ROIs for x, including shuffled
figureNo = figureNo + 1;
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

glm_r2ER=[];
glm_r2ER_ii=0;

id_r2ER_ii=0;
input_r2ER_data=[];


these_groups=[1 5]; %2 cm and 1 cm all trials
these_r2ER_x_means=[];
these_r2ERs_x=[];
these_r2ER_x_sh_means=[];
these_r2ER_x_sh=[];
these_SD_significant_r2ERs_x=[];
for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    %First get the r2ERs for each number of ROIs
    these_r2ERs=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ERs=[these_r2ERs handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_x_all_trials];
            end
        end
    end
    
    these_r2ERs_x=[these_r2ERs_x  these_r2ERs];
    these_r2ER_x_means=[these_r2ER_x_means mean(these_r2ERs)];

    

    %Then get the r2ERs_sh for each number of ROIs
    these_r2ER_sh=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ER_sh=[these_r2ER_sh handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_x_all_trials_sh];
            end
        end
    end
    
    
    these_r2ER_x_sh=[these_r2ER_x_sh  these_r2ER_sh];
    these_r2ER_x_sh_means=[these_r2ER_x_sh_means mean(these_r2ER_sh)];
    

    %Estimate significant r2ERs
    SD_these_r2ERs_sh=std(these_r2ER_sh);
    mean_these_r2ERs_sh=mean(these_r2ER_sh);
    these_SD_significant_r2ERs_x=[these_SD_significant_r2ERs_x these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh];

    %plot bar
    bar(bar_offset,mean(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_x]=drgViolinPoint(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;


    %r2ERs_x
    glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)))=these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh);
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)))=ROIsubset*ones(1,length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh));

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh);

    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=num2str(handles_choices.no_sub_ROIs(ROIsubset));
        ROIsubset_labels{ROIsubset}=num2str(handles_choices.no_sub_ROIs(ROIsubset));
    else
        input_r2ER_data(id_r2ER_ii).description='All';
        ROIsubset_labels{ROIsubset}='All';
    end

end

bar_offsets=[0:bar_offset-1];
plot(bar_offsets,these_r2ER_x_means,'-o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerFaceColor',[0.7 0.7 0.7])

xlim([-0.5 bar_offset-0.5])
xticklabels(ROIsubset_labels);   % tick labels from your cell array
xlabel('No ROIs')
ylabel('Significant r2ER for x')
title('Decoding x with different subsets of ROIs')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for r2ER x significant vs ROI number\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for r2ER x significant vs ROI number\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_r2ER.data',glm_r2ER.ROIsubset',...
    'VariableNames',{'r2ES','ROI_number'});
mdl = fitglm(tbl,'r2ES~ROI_number'...
    ,'CategoricalVars',[2])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for r2ER x significant vs ROI number\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for r2ER x significant vs ROI number\n']);

 
[output_data] = drgMutiRanksumorTtest(input_r2ER_data, fileID,0);


%Now plot r2ER vs the number of ROIs for y, including shuffled
figureNo = figureNo + 1;
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

glm_r2ER=[];
glm_r2ER_ii=0;

id_r2ER_ii=0;
input_r2ER_data=[];


these_groups=[1 5]; %2 cm and 1 cm all trials
these_r2ER_y_means=[];
these_r2ERs_y=[];
these_r2ER_y_sh_means=[];
these_r2ER_y_sh=[];
these_SD_significant_r2ERs_y=[];
for ROIsubset=1:length(handles_choices.no_sub_ROIs)
    %First get the r2ERs for each number of ROIs
    these_r2ERs=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ERs=[these_r2ERs handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_y_all_trials];
            end
        end
    end
    
    these_r2ERs_y=[these_r2ERs_y  these_r2ERs];
    these_r2ER_y_means=[these_r2ER_y_means mean(these_r2ERs)];

    

    %Then get the r2ERs_sh for each number of ROIs
    these_r2ER_sh=[];
    for ii_sub=1:handles_anal_XY.ROIsubset(ROIsubset).ii_sub_cum_sh
        for fileNo=1:length(handles_conc.arena_file)
            if sum(handles_conc.group(fileNo)==these_groups)>0
                these_r2ER_sh=[these_r2ER_sh handles_anal_XY.ROIsubset(ROIsubset).ii_sub(ii_sub).file(fileNo).r2ER_y_all_trials_sh];
            end
        end
    end
    
    
    these_r2ER_y_sh=[these_r2ER_y_sh  these_r2ER_sh];
    these_r2ER_y_sh_means=[these_r2ER_y_sh_means mean(these_r2ER_sh)];
    

    %Estimate significant r2ERs
    SD_these_r2ERs_sh=std(these_r2ER_sh);
    mean_these_r2ERs_sh=mean(these_r2ER_sh);
    these_SD_significant_r2ERs_y=[these_SD_significant_r2ERs_y these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh];

    %plot bar
    bar(bar_offset,mean(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
           
    %Violin plot
    try
    [mean_out, CIout,violin_y]=drgViolinPoint(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)...
        ,edges,bar_offset,rand_offset,'k','k',1);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    

    %r2ERs_y
    glm_r2ER.data(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)))=these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh);
    glm_r2ER.ROIsubset(glm_r2ER_ii+1:glm_r2ER_ii+length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)))=ROIsubset*ones(1,length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh)));
    glm_r2ER_ii=glm_r2ER_ii+length(these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh));

    id_r2ER_ii=id_r2ER_ii+1;
    input_r2ER_data(id_r2ER_ii).data=these_r2ERs(these_r2ERs>mean_these_r2ERs_sh+3*SD_these_r2ERs_sh);

    if ROIsubset<length(handles_choices.no_sub_ROIs)
        input_r2ER_data(id_r2ER_ii).description=num2str(handles_choices.no_sub_ROIs(ROIsubset));
        ROIsubset_labels{ROIsubset}=num2str(handles_choices.no_sub_ROIs(ROIsubset));
    else
        input_r2ER_data(id_r2ER_ii).description='All';
        ROIsubset_labels{ROIsubset}='All';
    end


end

bar_offsets=[0:1:bar_offset-1];
plot(bar_offsets,these_r2ER_y_means,'-o','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerFaceColor',[0.7 0.7 0.7])

xlim([-0.5 bar_offset-0.5])
xticklabels(ROIsubset_labels);   % tick labels from your cell array
xlabel('No ROIs')
ylabel('Significant r2ER for y')
title('Decoding y with different subsets of ROIs')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for significant r2ER y vs ROI number\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for significant r2ER y vs ROI number\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_r2ER.data',glm_r2ER.ROIsubset',...
    'VariableNames',{'r2ES','ROI_number'});
mdl = fitglm(tbl,'r2ES~ROI_number'...
    ,'CategoricalVars',[2])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for r2ER x vs ROI number\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for r2ER x vs ROI number\n']);

 
[output_data] = drgMutiRanksumorTtest(input_r2ER_data, fileID,0);


%Now plot r2ER x vs r2ER_conc
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

%Calculate mean
r2ER_edges=[0 0.025 0.05 0.1 0.15 0.25 0.35 0.45 0.6 0.9];
mean_these_r2ERs_x=zeros(1,length(r2ER_edges));
SD_these_r2ERs_x=zeros(1,length(r2ER_edges));
r2ER_edge_markers=zeros(1,length(r2ER_edges));
mean_these_r2ERs_x(1)=0;
SD_these_r2ERs_x(1)=0;
r2ER_edge_markers(1)=0;
for ii_edge=1:length(r2ER_edges)-1
    mean_these_r2ERs_x(ii_edge+1)=mean(these_r2ERs_x((these_r2ERs_conc>=r2ER_edges(ii_edge))&(these_r2ERs_conc<r2ER_edges(ii_edge+1))));
    SD_these_r2ERs_x(ii_edge+1)=std(these_r2ERs_x((these_r2ERs_conc>=r2ER_edges(ii_edge))&(these_r2ERs_conc<r2ER_edges(ii_edge+1))));
    r2ER_edge_markers(ii_edge+1)=mean(r2ER_edges(ii_edge)+r2ER_edges(ii_edge+1))/2;
end
[hlpvl, hppvl] = boundedline(r2ER_edge_markers,mean_these_r2ERs_x, SD_these_r2ERs_x,'cmap',[1 0 0]);


plot(these_r2ERs_conc,these_r2ERs_x,'.k')
plot(r2ER_edge_markers,mean_these_r2ERs_x,'-r','LineWidth',2)

xlim([0 0.6])
ylim([0 0.9])


xlabel('r2ER odor concentration')
ylabel('r2ER for x')
title('r2ERx vs r2ER odor concentration')


%Now plot r2ER y vs r2ER_conc
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

%Calculate mean
r2ER_edges=[0 0.025 0.05 0.1 0.15 0.25 0.35 0.45 0.6 0.9];
mean_these_r2ERs_y=zeros(1,length(r2ER_edges));
SD_these_r2ERs_y=zeros(1,length(r2ER_edges));
r2ER_edge_markers=zeros(1,length(r2ER_edges));
mean_these_r2ERs_y(1)=0;
SD_these_r2ERs_y(1)=0;
r2ER_edge_markers(1)=0;
for ii_edge=1:length(r2ER_edges)-1
    mean_these_r2ERs_y(ii_edge+1)=mean(these_r2ERs_y((these_r2ERs_conc>=r2ER_edges(ii_edge))&(these_r2ERs_conc<r2ER_edges(ii_edge+1))));
    SD_these_r2ERs_y(ii_edge+1)=std(these_r2ERs_y((these_r2ERs_conc>=r2ER_edges(ii_edge))&(these_r2ERs_conc<r2ER_edges(ii_edge+1))));
    r2ER_edge_markers(ii_edge+1)=mean(r2ER_edges(ii_edge)+r2ER_edges(ii_edge+1))/2;
end
[hlpvl, hppvl] = boundedline(r2ER_edge_markers,mean_these_r2ERs_y, SD_these_r2ERs_y,'cmap',[1 0 0]);


plot(these_r2ERs_conc,these_r2ERs_y,'.k')
plot(r2ER_edge_markers,mean_these_r2ERs_y,'-r','LineWidth',2)

xlim([0 0.6])
ylim([0 0.9])


xlabel('r2ER odor concentration')
ylabel('r2ER for y')
title('r2ERy vs r2ER odor concentration')

%Show histograms for r2ERs
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);


ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

%r2ER odor concentration
subplot(3,1,1)
hold on

edges=[0:0.05:0.9];

histogram(these_r2ERs_conc,edges)
ax = gca;                               % get current axes
ax.YScale = 'log';                      % logarithmic counts axis
ylabel('Count (log scale)');

title('r2ER odor concentration')

%r2ER odor x
subplot(3,1,2)
hold on

histogram(these_r2ERs_x,edges)

ax = gca;                               % get current axes
ax.YScale = 'log';                      % logarithmic counts axis
ylabel('Count (log scale)');

title('log10(r2ERx)')

%r2ER odor y
subplot(3,1,3)
hold on

histogram(these_r2ERs_y,edges)

ax = gca;                               % get current axes
ax.YScale = 'log';                      % logarithmic counts axis
ylabel('Count (log scale)');

title('log10(r2ERy)')



%Now plot r2ER x vs r2ER_conc for r2ERs that are significant for either
%conc or x
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

%Calculate mean
these_r2ERs_x_sig=these_r2ERs_x(these_SD_significant_r2ERs_x|these_SD_significant_r2ERs_conc);
these_r2ERs_conc_sig=these_r2ERs_conc(these_SD_significant_r2ERs_x|these_SD_significant_r2ERs_conc);
r2ER_edges=[0 0.025 0.05 0.1 0.15 0.25 0.35 0.45 0.6 0.9];
mean_these_r2ERs_x=zeros(1,length(r2ER_edges));
SD_these_r2ERs_x=zeros(1,length(r2ER_edges));
r2ER_edge_markers=zeros(1,length(r2ER_edges));
mean_these_r2ERs_x(1)=0;
SD_these_r2ERs_x(1)=0;
r2ER_edge_markers(1)=0;
for ii_edge=1:length(r2ER_edges)-1
    mean_these_r2ERs_x(ii_edge+1)=mean(these_r2ERs_x((these_r2ERs_conc>=r2ER_edges(ii_edge))&(these_r2ERs_conc<r2ER_edges(ii_edge+1))));
    SD_these_r2ERs_x(ii_edge+1)=std(these_r2ERs_x((these_r2ERs_conc>=r2ER_edges(ii_edge))&(these_r2ERs_conc<r2ER_edges(ii_edge+1))));
    r2ER_edge_markers(ii_edge+1)=mean(r2ER_edges(ii_edge)+r2ER_edges(ii_edge+1))/2;
end
[hlpvl, hppvl] = boundedline(r2ER_edge_markers,mean_these_r2ERs_x, SD_these_r2ERs_x,'cmap',[1 0 0]);


plot(these_r2ERs_conc,these_r2ERs_x,'.k')
plot(r2ER_edge_markers,mean_these_r2ERs_x,'-r','LineWidth',2)

xlim([0 0.6])
ylim([0 0.9])


xlabel('sig r2ER odor concentration')
ylabel('sog r2ER for x')
title('sig r2ERx vs r2ER odor concentration')


%Determine whether this is best fit by a linear relationship
x=these_r2ERs_conc_sig';
y=these_r2ERs_x_sig';

tbl = table(x, y);

% x, y are column vectors of the same length
% Example dummy data (replace with your own):
% x = linspace(0.1, 10, 20)';
% y = (2 * x ./ (1 + x)) + 0.1*randn(size(x));   % Vmax=2, Km=1 + noise

tbl = table(x, y);

% 1) Linear model: y ~ 1 + x
mdl_lin = fitlm(tbl, 'y ~ x');      % intercept + linear term

% 2) Quadratic model: y ~ 1 + x + x^2
% Adding x^2 in the formula automatically includes x as well.[web:114]
mdl_quad = fitlm(tbl, 'y ~ x + x^2');
% (Equivalently: mdl_quad = fitlm(tbl, 'y ~ x^2'); in recent MATLAB.[web:114])

% 3) Extra–sum-of-squares F-test: linear (reduced) vs quadratic (full)

% Residual sum of squares
RSS_lin  = sum(mdl_lin.Residuals.Raw.^2);
RSS_quad = sum(mdl_quad.Residuals.Raw.^2);

n      = numel(y);
p_lin  = numel(mdl_lin.Coefficients.Estimate);   % typically 2: intercept, x
p_quad = numel(mdl_quad.Coefficients.Estimate);  % typically 3: intercept, x, x^2

% Degrees of freedom
df1 = p_quad - p_lin;          % number of added parameters (should be 1)
df2 = n - p_quad;              % df for full model residuals

F = ((RSS_lin - RSS_quad) / df1) / (RSS_quad / df2);
p_F = 1 - fcdf(F, df1, df2);

fprintf('Extra-SS F-test (linear vs quadratic) for r2ER x vs conc: F = %.3f, p = %.3g\n', F, p_F);

if p_F < 0.05
    disp('Quadratic term significantly improves fit over linear (F-test).');
else
    disp('No significant improvement of quadratic over linear (F-test).');
end

% 4) Model comparison via AIC
%
%Akaike Information Criterion (AIC) is an information-theoretic criterion
%Smaller AIC is better; it trades off goodness of fit (via RSS) against model 
%complexity (via the number of parameters)

AIC_lin  = n*log(RSS_lin/n)  + 2*p_lin;
AIC_quad = n*log(RSS_quad/n) + 2*p_quad;

fprintf('RSS linear = %.4g, RSS quadratic = %.4g\n', RSS_lin, RSS_quad);
fprintf('AIC linear = %.4g, AIC quadratic = %.4g\n', AIC_lin, AIC_quad);

if AIC_quad < AIC_lin
    disp('Quadratic model preferred by AIC for r2ER y vs r2ER OC');
else
    disp('Linear model preferred or similar by AIC for r2ER y vs r2ER OC');
end


%Now plot r2ER y vs r2ER_conc
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

%Calculate mean
these_r2ERs_y_sig=these_r2ERs_y(these_SD_significant_r2ERs_y|these_SD_significant_r2ERs_conc);
these_r2ERs_conc_sig=these_r2ERs_conc(these_SD_significant_r2ERs_y|these_SD_significant_r2ERs_conc);
r2ER_edges=[0 0.025 0.05 0.1 0.15 0.25 0.35 0.45 0.6 0.9];
mean_these_r2ERs_y=zeros(1,length(r2ER_edges));
SD_these_r2ERs_y=zeros(1,length(r2ER_edges));
r2ER_edge_markers=zeros(1,length(r2ER_edges));
mean_these_r2ERs_y(1)=0;
SD_these_r2ERs_y(1)=0;
r2ER_edge_markers(1)=0;
for ii_edge=1:length(r2ER_edges)-1
    mean_these_r2ERs_y(ii_edge+1)=mean(these_r2ERs_y_sig((these_r2ERs_conc_sig>=r2ER_edges(ii_edge))&(these_r2ERs_conc_sig<r2ER_edges(ii_edge+1))));
    SD_these_r2ERs_y(ii_edge+1)=std(these_r2ERs_y_sig((these_r2ERs_conc_sig>=r2ER_edges(ii_edge))&(these_r2ERs_conc_sig<r2ER_edges(ii_edge+1))));
    r2ER_edge_markers(ii_edge+1)=mean(r2ER_edges(ii_edge)+r2ER_edges(ii_edge+1))/2;
end
[hlpvl, hppvl] = boundedline(r2ER_edge_markers,mean_these_r2ERs_y, SD_these_r2ERs_y,'cmap',[1 0 0]);


plot(these_r2ERs_conc_sig,these_r2ERs_y_sig,'.k')
plot(r2ER_edge_markers,mean_these_r2ERs_y,'-r','LineWidth',2)

xlim([0 0.6])
ylim([0 0.9])


xlabel('sig r2ER odor concentration')
ylabel('sig r2ER for y')
title('sig r2ERy vs r2ER odor concentration')



%Determine whether this is best fit by a linear relationship
x=these_r2ERs_conc_sig';
y=these_r2ERs_y_sig';

tbl = table(x, y);

% x, y are column vectors of the same length
% Example dummy data (replace with your own):
% x = linspace(0.1, 10, 20)';
% y = (2 * x ./ (1 + x)) + 0.1*randn(size(x));   % Vmax=2, Km=1 + noise

tbl = table(x, y);

% 1) Linear model: y ~ 1 + x
mdl_lin = fitlm(tbl, 'y ~ x');      % intercept + linear term

% 2) Quadratic model: y ~ 1 + x + x^2
% Adding x^2 in the formula automatically includes x as well.[web:114]
mdl_quad = fitlm(tbl, 'y ~ x + x^2');
% (Equivalently: mdl_quad = fitlm(tbl, 'y ~ x^2'); in recent MATLAB.[web:114])

% 3) Extra–sum-of-squares F-test: linear (reduced) vs quadratic (full)

% Residual sum of squares
RSS_lin  = sum(mdl_lin.Residuals.Raw.^2);
RSS_quad = sum(mdl_quad.Residuals.Raw.^2);

n      = numel(y);
p_lin  = numel(mdl_lin.Coefficients.Estimate);   % typically 2: intercept, x
p_quad = numel(mdl_quad.Coefficients.Estimate);  % typically 3: intercept, x, x^2

% Degrees of freedom
df1 = p_quad - p_lin;          % number of added parameters (should be 1)
df2 = n - p_quad;              % df for full model residuals

F = ((RSS_lin - RSS_quad) / df1) / (RSS_quad / df2);
p_F = 1 - fcdf(F, df1, df2);

fprintf('Extra-SS F-test (linear vs quadratic) for r2ER y vs conc: F = %.3f, p = %.3g\n', F, p_F);

if p_F < 0.05
    disp('Quadratic term significantly improves fit over linear (F-test).');
else
    disp('No significant improvement of quadratic over linear (F-test).');
end

% 4) Model comparison via AIC
%
%Akaike Information Criterion (AIC) is an information-theoretic criterion
%Smaller AIC is better; it trades off goodness of fit (via RSS) against model 
%complexity (via the number of parameters)

AIC_lin  = n*log(RSS_lin/n)  + 2*p_lin;
AIC_quad = n*log(RSS_quad/n) + 2*p_quad;

fprintf('RSS linear = %.4g, RSS quadratic = %.4g\n', RSS_lin, RSS_quad);
fprintf('AIC linear = %.4g, AIC quadratic = %.4g\n', AIC_lin, AIC_quad);

if AIC_quad < AIC_lin
    disp('Quadratic model preferred by AIC for r2ER y vs r2ER OC');
else
    disp('Linear model preferred or similar by AIC for r2ER y vs r2ER OC');
end


%Show histograms for sig r2ERs
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);


ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

%r2ER odor concentration
subplot(3,1,1)
hold on

edges=[0:0.05:0.9];

histogram(these_r2ERs_conc_sig,edges)
ax = gca;                               % get current axes
ax.YScale = 'log';                      % logarithmic counts axis
ylabel('Count (log scale)');

title('r2ER odor concentration')

%r2ER odor x
subplot(3,1,2)
hold on

histogram(these_r2ERs_x_sig,edges)

ax = gca;                               % get current axes
ax.YScale = 'log';                      % logarithmic counts axis
ylabel('Count (log scale)');

title('log10(r2ERx)')

%r2ER odor y
subplot(3,1,3)
hold on

histogram(these_r2ERs_y_sig,edges)

ax = gca;                               % get current axes
ax.YScale = 'log';                      % logarithmic counts axis
ylabel('Count (log scale)');

title('log10(r2ERy)')

save(saveFile,'handles_anal_XY','handles_anal_conc','-v7.3')

fclose(fileID);

pffft=1;