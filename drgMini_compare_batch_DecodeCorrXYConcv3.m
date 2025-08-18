%drgMini_compare_batch_DecodeCorrXYConcv3
close all
clear all

is_sphgpu=0;
% is_pearson=1; %If this is 1 Pearson correlation is calculated, otherwise Spearman

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
          %1=fitrnet
    %2=fitrgp
    %3=fitrtree
    %4=fitglm
    %5=fitrsvm

        algo_legend{3}='ANN';
        algo_legend{5}='GPR';
        algo_legend{1}='BT';
        algo_legend{2}='GLM';
        algo_legend{4}='SVM';

        %GPR
        
        %before_bins=0

        %Odor decoding
        save_PathConc{5}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcGP0b06062025/';
        choiceOdorConcFileName{5}='drgDynamicOdorConcChoices_Fabio_GPb0_06062025.m';

        %Position decoding 
        save_PathXY{5}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName{5}='drgOdorArenaChoices_Fabio_Good_01122025.m';
        
        %BINARY TREE
        
        %binary tree for paper before_bins=0

        %Odor decoding
        save_PathConc{1}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc04192024/';
        choiceOdorConcFileName{1}='drgDynamicOdorConcChoices_Fabio_Good_04192024.m';

        %Position decoding 
        save_PathXY{1}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName{1}='drgOdorArenaChoices_Fabio_Good_01122025.m';

        %binary tree before_bins=5

        %Odor decoding
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc06042025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_06042025.m'
        % 
        % %Position decoding 
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput_bin_b5_06052025/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_bin_b5_06052025.m';

        %SVM

        %before_bins=0
        save_PathConc{4}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcSVM06042025/';
        choiceOdorConcFileName{4}='drgDynamicOdorConcChoices_Fabio_SVMb0_06052025.m';

        %Position decoding, note: this is currenlty binary tree 
        save_PathXY{4}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName{4}='drgOdorArenaChoices_Fabio_Good_01122025.m';

        %before_bins=5
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcSVM0b06052025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_SVMb0_06052025.m'


        %ANN
        
        %ann for paper before_bins=0

        %Odor decoing
        save_PathConc{3}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcANN0b06202025/';
        choiceOdorConcFileName{3}='drgDynamicOdorConcChoices_Fabio_ANNb0_06202025.m';

        %Position decoding, note: this is currenlty binary tree 
        save_PathXY{3}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName{3}='drgOdorArenaChoices_Fabio_Good_01122025.m';


        %GLM
        
        %glm for paper before_bins=0

        %Odor decoding
        save_PathConc{2}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConcAGLM0b06232025/';
        choiceOdorConcFileName{2}='drgDynamicOdorConcChoices_Fabio_glmb0_06232025.m';

        %Position decoding, note: this is currenlty binary tree 
        save_PathXY{2}='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName{2}='drgOdorArenaChoices_Fabio_Good_01122025.m';

        
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeBayesOdorConc05072025/';
        % choiceOdorConcFileName='drgDynamicOdorConcChoices_Bayes_05072025.m';

        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';
        %
        % save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        % choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        

        %Angle processing
        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end

%Exclude files with p value for R1_all_trials > 0.05

exclusion_criterion=3;
%1 p<=thr_rho
%2 p<pDFR
%3 p<Bonferroni
%4 R>=0.1

thr_rho=0.05;

R_thr=0.1;
R2ER_thr=0.05;

%These are used to classify angle approach
low_angle=-130;
high_angle=-50;


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

%Calcualte R1, R2 and R2ER and p values for each decoding algo
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

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

for ii_algo=1:length(save_PathConc)

    this_choiceOdorConcFileName=choiceOdorConcFileName{ii_algo};
    eval(['handles_conc=' this_choiceOdorConcFileName(1:end-2) ';'])
    this_choiceXYFileName=choiceXYFileName{ii_algo};
    eval(['handles_XY=' this_choiceXYFileName(1:end-2) ';'])
    %Now calcualte R1, R2 and R2ER and p values

    ii_run=1;


    R1_all_trials=[];
    R1_hits=[];
    R1_miss=[];
    R1_all_trials_sh=[];

    p_R1_hits=[];

    R2_rect_all_trials=[];
    R2_rect_hits=[];
    R2_rect_miss=[];
    R2_rect_all_trials_sh=[];

    R2sh_all_trials=[];
    R2sh_hits=[];
    R2sh_miss=[];
    R2sh_all_trials_sh=[];

    r2ER_all_trials=[];
    r2ER_hits=[];
    r2ER_miss=[];
    r2ER_all_trials_sh=[];

    p_r2ER_hits=[];

    Var_per_point_all=[];
    Var_per_point_hits=[];
    Var_per_point_miss=[];
    Var_per_point_all_sh=[];

    these_groups=[1:5];

    for fileNo=1:length(handles_conc.arena_file)
        if sum(handles_conc.group(fileNo)==these_groups)>0
            op_all_trials=[];
            op_decod_all_trials=[];



            op_all_hits=[];
            op_decod_all_hits=[];
            op_all_hits90=[];
            op_decod_all_hits90=[];
            op_all_hitso=[];
            op_decod_all_hitso=[];
            % op_all_hits_sh=[];
            % op_decod_all_hits_sh=[];
            op_all_miss=[];
            op_decod_all_miss=[];
            % op_all_miss_sh=[];
            % op_decod_all_miss_sh=[];
            % op_decod_all_trials_sh=[];
            op_between_trials=[];
            op_decod_between_trials=[];
            % op_between_trials_sh=[];
            % op_decod_between_trials_sh=[];

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

            %Load conc data
            arena_file=handles_conc.arena_file{fileNo};
            %load the ouptut file
            this_save_PathConc=save_PathConc{ii_algo};
            load([this_save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            trials=handles_out.trials;
            odor_plume_template=handles_out.odor_plume_template;
            op_predicted=handles_out.op_predicted;

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

            op_predicted_sh=handles_out.op_predicted_sh;
            op_predicted_sh_conv=zeros(size(op_predicted_sh,1),size(op_predicted_sh,2));
            for ii_sh=1:size(op_predicted_sh,2)
                this_op_predicted_sh=zeros(size(op_predicted_sh,1),1);
                this_op_predicted_sh(:,1)=op_predicted_sh(:,ii_sh);
                this_op_predicted_conv_sh=conv(this_op_predicted_sh,conv_win_gauss,'same');

                %Now limit the x and y to max and min
                minop=min(odor_plume_template);
                this_op_predicted_conv_sh(this_op_predicted_conv_sh<minop)=minop;
                maxop=max(odor_plume_template);
                this_op_predicted_conv_sh(this_op_predicted_conv_sh>maxop)=maxop;
                op_predicted_sh_conv(:,ii_sh)=this_op_predicted_conv_sh;
            end

            op_decod_trials_sh=[];
            for ii_sh=1:size(op_predicted_sh_conv,2)
                op_decod_trials_sh.ii_sh(ii_sh).oppsh=[];
            end

            last_op_predictedend=1;
            ii_start=0;

            for trNo=1:trials.odor_trNo

                op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if op_predictedend>length(odor_plume_template)
                    op_predictedend=length(odor_plume_template);
                end

                op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
                op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];

                % op_between_trials=[op_between_trials; odor_plume_template(last_op_predictedend:op_predictedstart)'];
                % op_decod_between_trials=[op_decod_between_trials; op_predicted_conv(last_op_predictedend:op_predictedstart)];
                last_op_predictedend=op_predictedend;

                ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;


                for ii_sh=1:size(op_predicted_sh_conv,2)
                    this_op_decod_trials_sh=zeros(size(op_predicted_sh_conv,1),1);
                    this_op_decod_trials_sh(:,1)=op_predicted_sh_conv(:,ii_sh);
                    op_decod_trials_sh.ii_sh(ii_sh).oppsh=[op_decod_trials_sh.ii_sh(ii_sh).oppsh; this_op_decod_trials_sh(op_predictedstart:op_predictedend)];
                end


                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
                        op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];

                    case 2
                        %Lane 1 miss orange
                        op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];

                    case 3
                        %Lane 4 hit blue
                        op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];

                    case 4
                        %Lane 4 miss sky blue
                        op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
                end
                ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
            end

            these_fileNos(fileNo)=fileNo;

            if ~isempty(op_all_trials)

                [R1,this_p_R1_hits]=corrcoef(op_all_trials,op_decod_all_trials);
                R1_all_trials(fileNo)=R1(1,2);
                p_R1_hits(fileNo)=this_p_R1_hits(1,2);

            else
                R1_all_trials(fileNo)=NaN;
            end

            if ~isempty(op_all_hits)

                R1=corrcoef(op_all_hits,op_decod_all_hits);
                R1_hits(fileNo)=R1(1,2);

            else
                R1_hits(fileNo)=NaN;
            end



            if ~isempty(op_all_miss)

                R1=corrcoef(op_all_miss,op_decod_all_miss);
                R1_miss(fileNo)=R1(1,2);

            else
                R1_miss(fileNo)=NaN;
            end

            %Shuffled
            if ~isempty(op_predicted_sh_conv)

                these_R2shs=[];
                these_R1shs=[];

                for ii_sh=1:size(op_predicted_sh_conv,2)

                    this_op_all_trials_sh=op_decod_trials_sh.ii_sh(ii_sh).oppsh;
                    this_R1=corrcoef(op_all_trials,this_op_all_trials_sh);
                    these_R1shs=[these_R1shs this_R1(1,2)];

                    if sum((op_all_trials-mean(op_all_trials)).^2)>=sum((this_op_all_trials_sh-mean(this_op_all_trials_sh)).^2)
                        this_R2sh=1-sum((this_op_all_trials_sh-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
                    else
                        this_R2sh=1-sum((this_op_all_trials_sh-op_all_trials).^2)/sum((this_op_all_trials_sh-mean(this_op_all_trials_sh)).^2);
                    end
                    these_R2shs=[these_R2shs this_R2sh];

                end
                R2sh_all_trials_sh(fileNo)=mean(these_R2shs);
                R1_all_trials_sh(fileNo)=mean(these_R1shs);
            else
                R2sh_all_trials_sh(fileNo)=NaN;
                R1_all_trials_sh(fileNo)==NaN;
            end



            %Calculate R2, the fraction of variance explained
            if ~isempty(op_all_trials)
                this_R2=1-sum((op_decod_all_trials-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
                R2_rect_all_trials(fileNo)=drgMini_rectifyR2(this_R2);
                % r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
                Var_per_point_all(fileNo)=sum((op_all_trials-mean(op_all_trials)).^2)/length(op_all_trials);
            else
                R2_rect_all_trials(fileNo)=NaN;
            end

            if ~isempty(op_all_hits)
                this_R2=1-sum((op_decod_all_hits-op_all_hits).^2)/(sum((op_all_hits-mean(op_all_hits)).^2));
                R2_rect_hits(fileNo)=drgMini_rectifyR2(this_R2);
                % r2ER_hits(fileNo)=drgMini_r2ER(op_all_hits,op_decod_all_hits);
                Var_per_point_hits(fileNo)=(sum((op_all_hits-mean(op_all_hits)).^2))/length(op_all_hits);
            else
                R2_rect_hits(fileNo)=NaN;
            end

            if ~isempty(op_all_miss)
                this_R2=1-sum((op_decod_all_miss-op_all_miss).^2)/(sum((op_all_miss-mean(op_all_miss)).^2));
                R2_rect_miss(fileNo)=drgMini_rectifyR2(this_R2);
                % r2ER_miss(fileNo)=drgMini_r2ER(op_all_miss,op_decod_all_miss);
                Var_per_point_miss(fileNo)=(sum((op_all_miss-mean(op_all_miss)).^2))/length(op_all_miss);
            else
                R2_rect_miss(fileNo)=NaN;
            end

            %Shuffled
            if ~isempty(op_predicted_sh_conv)
                these_R2s=[];
                these_Vars=[];
                these_r2ERs=[];
                for ii_sh=1:size(op_predicted_sh_conv,2)
                    this_op_all_trials_sh=op_decod_trials_sh.ii_sh(ii_sh).oppsh;
                    this_R2=1-sum((this_op_all_trials_sh-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
                    these_R2s=[these_R2s drgMini_rectifyR2(this_R2)];
                    % these_r2ERs=[these_r2ERs drgMini_r2ER(op_all_trials,this_op_all_trials_sh)];
                    these_Vars=[these_Vars sum((op_all_trials-mean(op_all_trials)).^2)/length(op_all_trials)];
                end
                R2_rect_all_trials_sh(fileNo)=mean(these_R2s);
                % r2ER_all_trials_sh(fileNo)=mean(these_r2ERs);
                Var_per_point_all_sh(fileNo)=mean(these_Vars);
            else
                R2_rect_all_trials_sh(fileNo)=NaN;
            end

            %Calculate R2sh, the fraction of variance shared
            if ~isempty(op_all_trials)
                if sum((op_all_trials-mean(op_all_trials)).^2)>=sum((op_decod_all_trials-mean(op_decod_all_trials)).^2)
                    this_R2sh=1-sum((op_decod_all_trials-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
                else
                    this_R2sh=1-sum((op_decod_all_trials-op_all_trials).^2)/sum((op_decod_all_trials-mean(op_decod_all_trials)).^2);
                end
                R2sh_all_trials(fileNo)=this_R2sh;
            else
                R2sh_all_trials(fileNo)=NaN;
            end

            if ~isempty(op_all_hits)
                if (sum((op_all_hits-mean(op_all_hits)).^2))>=(sum((op_decod_all_hits-mean(op_decod_all_hits)).^2))
                    this_R2sh=1-sum((op_decod_all_hits-op_all_hits).^2)/(sum((op_all_hits-mean(op_all_hits)).^2));
                else
                    this_R2sh=1-sum((op_decod_all_hits-op_all_hits).^2)/(sum((op_decod_all_hits-mean(op_decod_all_hits)).^2));
                end
                R2sh_hits(fileNo)=this_R2sh;
            else
                R2_rect_hits(fileNo)=NaN;
            end

            if ~isempty(op_all_miss)
                if (sum((op_all_miss-mean(op_all_miss)).^2))>=(sum((op_decod_all_miss-mean(op_decod_all_miss)).^2))
                    this_R2sh=1-sum((op_decod_all_miss-op_all_miss).^2)/(sum((op_all_miss-mean(op_all_miss)).^2));
                else
                    this_R2sh=1-sum((op_decod_all_miss-op_all_miss).^2)/(sum((op_decod_all_miss-mean(op_decod_all_miss)).^2));
                end
                R2sh_miss(fileNo)=this_R2sh;
            else
                R2sh_miss(fileNo)=NaN;
            end

            %Shuffled
            if ~isempty(op_predicted_sh_conv)
                these_R2shs=[];

                for ii_sh=1:size(op_predicted_sh_conv,2)
                    this_op_all_trials_sh=op_decod_trials_sh.ii_sh(ii_sh).oppsh;
                    if sum((op_all_trials-mean(op_all_trials)).^2)>=sum((this_op_all_trials_sh-mean(this_op_all_trials_sh)).^2)
                        this_R2sh=1-sum((this_op_all_trials_sh-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
                    else
                        this_R2sh=1-sum((this_op_all_trials_sh-op_all_trials).^2)/sum((this_op_all_trials_sh-mean(this_op_all_trials_sh)).^2);
                    end
                    these_R2shs=[these_R2shs this_R2sh];

                end
                R2sh_all_trials_sh(fileNo)=mean(these_R2shs);
            else
                R2sh_all_trials_sh(fileNo)=NaN;
            end

        end
    end

    for fileNo=1:length(handles_conc.arena_file)
        if sum(handles_conc.group(fileNo)==these_groups)>0
            op_all_trials=[];
            op_decod_all_trials=[];



            op_all_hits=[];
            op_decod_all_hits=[];
            op_all_hits90=[];
            op_decod_all_hits90=[];
            op_all_hitso=[];
            op_decod_all_hitso=[];
            % op_all_hits_sh=[];
            % op_decod_all_hits_sh=[];
            op_all_miss=[];
            op_decod_all_miss=[];
            % op_all_miss_sh=[];
            % op_decod_all_miss_sh=[];
            % op_decod_all_trials_sh=[];
            op_between_trials=[];
            op_decod_between_trials=[];
            % op_between_trials_sh=[];
            % op_decod_between_trials_sh=[];

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

            %Load conc data
            arena_file=handles_conc.arena_file{fileNo};
            %load the ouptut file
            this_save_PathConc=save_PathConc{ii_algo};
            load([this_save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            trials=handles_out.trials;
            odor_plume_template=handles_out.odor_plume_template;
            op_predicted=handles_out.op_predicted;

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

            op_predicted_sh=handles_out.op_predicted_sh;
            op_predicted_sh_conv=zeros(size(op_predicted_sh,1),size(op_predicted_sh,2));
            for ii_sh=1:size(op_predicted_sh,2)
                this_op_predicted_sh=zeros(size(op_predicted_sh,1),1);
                this_op_predicted_sh(:,1)=op_predicted_sh(:,ii_sh);
                this_op_predicted_conv_sh=conv(this_op_predicted_sh,conv_win_gauss,'same');

                %Now limit the x and y to max and min
                minop=min(odor_plume_template);
                this_op_predicted_conv_sh(this_op_predicted_conv_sh<minop)=minop;
                maxop=max(odor_plume_template);
                this_op_predicted_conv_sh(this_op_predicted_conv_sh>maxop)=maxop;
                op_predicted_sh_conv(:,ii_sh)=this_op_predicted_conv_sh;
            end

            op_decod_trials_sh=[];
            for ii_sh=1:size(op_predicted_sh_conv,2)
                op_decod_trials_sh.ii_sh(ii_sh).oppsh=[];
            end

            last_op_predictedend=1;
            ii_start=0;

            for trNo=1:trials.odor_trNo

                op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if op_predictedend>length(odor_plume_template)
                    op_predictedend=length(odor_plume_template);
                end

                op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
                op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];

                op_between_trials=[op_between_trials; odor_plume_template(last_op_predictedend:op_predictedstart)'];
                op_decod_between_trials=[op_decod_between_trials; op_predicted_conv(last_op_predictedend:op_predictedstart)];
                last_op_predictedend=op_predictedend;

                ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;


                for ii_sh=1:size(op_predicted_sh_conv,2)
                    this_op_decod_trials_sh=zeros(size(op_predicted_sh_conv,1),1);
                    this_op_decod_trials_sh(:,1)=op_predicted_sh_conv(:,ii_sh);
                    op_decod_trials_sh.ii_sh(ii_sh).oppsh=[op_decod_trials_sh.ii_sh(ii_sh).oppsh; this_op_decod_trials_sh(op_predictedstart:op_predictedend)];
                end


                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
                        op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];

                        %Hit 90 degrees
                        if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                            op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
                            op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
                        end

                        %Horizontal hit
                        if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                            op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
                            op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
                        end
                    case 2
                        %Lane 1 miss orange
                        op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];

                    case 3
                        %Lane 4 hit blue
                        op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];

                        %Hit 90 degrees
                        if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                            op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
                            op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
                        end

                        %Horizontal hit
                        if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                            op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
                            op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
                        end
                    case 4
                        %Lane 4 miss sky blue
                        op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                        op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
                end
                ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
            end


            %Calculate R2, the fraction of variance explained
            if ~isempty(op_all_trials)
                r2ER_all_trials(fileNo)=drgMini_r2ER(op_all_trials,op_decod_all_trials);
            else
                r2ER_all_trials(fileNo)=NaN;
            end

            if ~isempty(op_all_hits)
                [this_r2ER,this_P_r2ER]=drgMini_r2ER_bootstrap(op_all_hits,op_decod_all_hits,10000);
                r2ER_hits(fileNo)=this_r2ER;
                p_r2ER_hits(fileNo)=this_P_r2ER;
            else
                r2ER_hits(fileNo)=NaN;
            end

            if ~isempty(op_all_miss)
                r2ER_miss(fileNo)=drgMini_r2ER(op_all_miss,op_decod_all_miss);
            else
                r2ER_miss(fileNo)=NaN;
            end

            %Shuffled
            if ~isempty(op_predicted_sh_conv)
                these_R2s=[];
                these_Vars=[];
                these_r2ERs=[];
                for ii_sh=1:size(op_predicted_sh_conv,2)
                    this_op_all_trials_sh=op_decod_trials_sh.ii_sh(ii_sh).oppsh;
                    these_r2ERs=[these_r2ERs drgMini_r2ER(op_all_trials,this_op_all_trials_sh)];
                end
                r2ER_all_trials_sh(fileNo)=mean(these_r2ERs);
            else
                r2ER_all_trials_sh(fileNo)=NaN;
            end


        end
    end

    pFDRR1 = drsFDRpval(p_R1_hits);
    pFDRr2ER = drsFDRpval(p_r2ER_hits);

    p_bonf_R1 = 0.05/length(p_R1_hits);
    p_bonf_r2ER = 0.05/length(p_r2ER_hits);

    %Now cacluate file exclusions
    file_exclusions=zeros(1,length(handles_conc.arena_file));

    for fileNo=1:length(handles_conc.arena_file)
        switch exclusion_criterion
            case 1
                %1 p<=thr_rho
                if p_r2ER_hits(fileNo)<=thr_rho
                    exclude=0;
                else
                    exclude=1;
                end
            case 2
                %2 p<pDFR
                if p_r2ER_hits(fileNo)<=pFDRr2ER
                    exclude=0;
                else
                    exclude=1;
                end
            case 3
                %3 p<Bonferroni
                if p_r2ER_hits(fileNo)<=p_bonf_r2ER
                    exclude=0;
                else
                    exclude=1;
                end
            case 4
                %4 R>=0.1
                if r2ER_hits(fileNo)>=R2ER_thr
                    exclude=0;
                else
                    exclude=1;
                end
        end
        file_exclusions(fileNo)=exclude;
    end

    include_no_odor=[];
    these_groups=[4];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(file_exclusions(fileNo)==0)
            include_no_odor(fileNo)=1;
        end
    end

    include_odor=[];
    these_groups=[1 5];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(file_exclusions(fileNo)==0)
            include_odor(fileNo)=1;
        end
    end

    %Now calculate the hits et al for XY with no odor
    these_groups=[1:5];

    ii_run=1;
    ii_for_corr=0;
    R1_XY_all_trials=[];
    R1_XY_hits=[];
    R1_XY_miss=[];
    R1_XY_all_trials_sh=[];


    R1_x_all_trials=[];
    R1_x_hits=[];
    R1_x_miss=[];
    R1_x_all_trials_sh=[];

    R1_y_all_trials=[];
    R1_y_hits=[];
    R1_y_miss=[];
    R1_y_all_trials_sh=[];

    these_XYfileNos=[];


    for fileNo=1:length(handles_XY.arena_file)
        if sum(handles_XY.group(fileNo)==these_groups)>0
            x_all_trials=[];
            x_decod_all_trials=[];
            x_all_hits=[];
            x_decod_all_hits=[];
            x_all_miss=[];
            x_decod_all_miss=[];
            P_rho_x_all_trials=[];

            y_all_trials=[];
            y_decod_all_trials=[];
            y_all_hits=[];
            y_decod_all_hits=[];
            y_all_miss=[];
            y_decod_all_miss=[];
            P_rho_x_all_trials=[];

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

            %Load conc data
            arena_file=handles_XY.arena_file{fileNo};
            %load the ouptut file
            this_save_PathXY=save_PathXY{ii_algo};
            load([this_save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
            trials=handles_out.trials;
            x_predicted=handles_out.x_predicted;
            y_predicted=handles_out.y_predicted;
            x_predicted_sh=handles_out.x_predicted_sh;
            y_predicted_sh=handles_out.y_predicted_sh;
            XYtest=handles_out.XYtest;

            for ii_sh=1:size(x_predicted_sh,2)
                x_decod_all_trials_sh.ii_sh(ii_sh).xdec=[];
            end

            for ii_sh=1:size(y_predicted_sh,2)
                y_decod_all_trials_sh.ii_sh(ii_sh).ydec=[];
            end

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

            x_predicted_sh_conv=zeros(size(x_predicted_sh,1),size(x_predicted_sh,2));
            y_predicted_sh_conv=zeros(size(y_predicted_sh,1),size(y_predicted_sh,2));

            for ii_sh=1:size(x_predicted_sh,2)

                this_x_predicted_sh=zeros(size(x_predicted_sh,1),1);
                this_x_predicted_sh(:,1)=x_predicted_sh(:,ii_sh);
                this_x_predicted_sh_conv=zeros(size(x_predicted_sh,1),1);

                this_y_predicted_sh=zeros(size(y_predicted_sh,1),1);
                this_y_predicted_sh(:,1)=y_predicted_sh(:,ii_sh);
                this_y_predicted_sh_conv=zeros(size(y_predicted_sh,1),1);

                this_x_predicted_sh_conv=conv(this_x_predicted_sh,conv_win_gauss,'same');
                this_y_predicted_sh_conv=conv(this_y_predicted_sh,conv_win_gauss,'same');

                minY1=min(XYtest(:,1));
                this_x_predicted_sh_conv(this_x_predicted_sh_conv<minY1)=minY1;
                maxY1=max(XYtest(:,1));
                this_x_predicted_sh_conv(this_x_predicted_sh_conv>maxY1)=maxY1;

                minY2=min(XYtest(:,2));
                this_y_predicted_sh_conv(this_y_predicted_sh_conv<minY2)=minY2;
                maxY2=max(XYtest(:,2));
                this_y_predicted_sh_conv(this_y_predicted_sh_conv>maxY2)=maxY2;

                x_predicted_sh_conv(:,ii_sh)=this_x_predicted_sh_conv;
                y_predicted_sh_conv(:,ii_sh)=this_y_predicted_sh_conv;
            end


            for trNo=1:trials.odor_trNo

                ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if ii_predictedend>size(XYtest,1)
                    ii_predictedend=size(XYtest,1);
                end

                x_all_trials=[x_all_trials; XYtest(ii_predictedstart:ii_predictedend,1)];
                x_decod_all_trials=[x_decod_all_trials; x_predicted_conv(ii_predictedstart:ii_predictedend)];

                y_all_trials=[y_all_trials; XYtest(ii_predictedstart:ii_predictedend,2)];
                y_decod_all_trials=[y_decod_all_trials; y_predicted_conv(ii_predictedstart:ii_predictedend)];


                for ii_sh=1:size(x_predicted_sh,2)
                    this_x_predicted_sh_conv=zeros(size(x_predicted_sh_conv,1),1);
                    this_x_predicted_sh_conv(:,1)=x_predicted_sh_conv(:,ii_sh);
                    x_decod_all_trials_sh.ii_sh(ii_sh).xdec=[x_decod_all_trials_sh.ii_sh(ii_sh).xdec; this_x_predicted_sh_conv(ii_predictedstart:ii_predictedend)];

                    this_y_predicted_sh_conv=zeros(size(y_predicted_sh_conv,1),1);
                    this_y_predicted_sh_conv(:,1)=y_predicted_sh_conv(:,ii_sh);
                    y_decod_all_trials_sh.ii_sh(ii_sh).ydec=[y_decod_all_trials_sh.ii_sh(ii_sh).ydec; this_y_predicted_sh_conv(ii_predictedstart:ii_predictedend)];
                end
                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
                        x_all_hits=[x_all_hits; XYtest(ii_predictedstart:ii_predictedend,1)];
                        x_decod_all_hits=[x_decod_all_hits; x_predicted_conv(ii_predictedstart:ii_predictedend)];

                        y_all_hits=[y_all_hits; XYtest(ii_predictedstart:ii_predictedend,2)];
                        y_decod_all_hits=[y_decod_all_hits; y_predicted_conv(ii_predictedstart:ii_predictedend)];

                    case 2
                        %Lane 1 miss orange
                        x_all_miss=[x_all_miss; XYtest(ii_predictedstart:ii_predictedend,1)];
                        x_decod_all_miss=[x_decod_all_miss; x_predicted_conv(ii_predictedstart:ii_predictedend)];

                        y_all_miss=[y_all_miss; XYtest(ii_predictedstart:ii_predictedend,2)];
                        y_decod_all_miss=[y_decod_all_miss; y_predicted_conv(ii_predictedstart:ii_predictedend)];

                    case 3
                        %Lane 4 hit blue
                        x_all_hits=[x_all_hits; XYtest(ii_predictedstart:ii_predictedend,1)];
                        x_decod_all_hits=[x_decod_all_hits; x_predicted_conv(ii_predictedstart:ii_predictedend)];

                        y_all_hits=[y_all_hits; XYtest(ii_predictedstart:ii_predictedend,2)];
                        y_decod_all_hits=[y_decod_all_hits; y_predicted_conv(ii_predictedstart:ii_predictedend)];
                    case 4
                        %Lane 4 miss sky blue
                        x_all_miss=[x_all_miss; XYtest(ii_predictedstart:ii_predictedend,1)];
                        x_decod_all_miss=[x_decod_all_miss; x_predicted_conv(ii_predictedstart:ii_predictedend)];

                        y_all_miss=[y_all_miss; XYtest(ii_predictedstart:ii_predictedend,2)];
                        y_decod_all_miss=[y_decod_all_miss; y_predicted_conv(ii_predictedstart:ii_predictedend)];
                end

            end




            if ~isempty(x_all_trials)

                [R1,this_p_R1_hits]=corrcoef(x_all_trials,x_decod_all_trials);
                R1_x_all_trials(fileNo)=R1(1,2);
                P_rho_x_all_trials(fileNo)=this_p_R1_hits(1,2);

            else
                R1_x_all_trials(fileNo)=NaN;
            end

            if ~isempty(y_all_trials)

                [R1,this_p_R1_hits]=corrcoef(y_all_trials,y_decod_all_trials);
                R1_y_all_trials(fileNo)=R1(1,2);
                P_rho_y_all_trials(fileNo)=this_p_R1_hits(1,2);

            else
                R1_y_all_trials(fileNo)=NaN;
            end

            if ~isempty(x_all_trials)

                R1XY=corr2([x_all_trials y_all_trials],[x_decod_all_trials y_decod_all_trials]);
                R1_XY_all_trials(fileNo)=R1XY;

            else
                R1_XY_all_trials(fileNo)=NaN;
            end

            if ~isempty(x_all_hits)

                R1=corrcoef(x_all_hits,x_decod_all_hits);
                R1_x_hits(fileNo)=R1(1,2);

            else
                R1_x_hits(fileNo)=NaN;
            end

            if ~isempty(y_all_hits)

                R1=corrcoef(y_all_hits,y_decod_all_hits);
                R1_y_hits(fileNo)=R1(1,2);

            else
                R1_y_hits(fileNo)=NaN;
            end

            if ~isempty(x_all_trials)
                R1XY=corr2([x_all_hits y_all_hits],[x_decod_all_hits y_decod_all_hits]);
                R1_XY_hits(fileNo)=R1XY;
                % P_rho_XY_all_trials(fileNo)=this_p_R1_hits(1,2);
            else
                R1_XY_hits(fileNo)=NaN;
            end

            if ~isempty(x_all_miss)

                R1=corrcoef(x_all_miss,x_decod_all_miss);
                R1_x_miss(fileNo)=R1(1,2);

            else
                R1_x_miss(fileNo)=NaN;
            end

            if ~isempty(y_all_miss)

                R1=corrcoef(y_all_miss,y_decod_all_miss);
                R1_y_miss(fileNo)=R1(1,2);

            else
                R1_y_miss(fileNo)=NaN;
            end

            if ~isempty(x_all_trials)
                R1XY=corr2([x_all_miss y_all_miss],[x_decod_all_miss y_decod_all_miss]);
                R1_XY_miss(fileNo)=R1XY;
                % P_rho_XY_all_trials(fileNo)=this_p_R1_hits(1,2);
            else
                R1_XY_miss(fileNo)=NaN;
            end


            if ~isempty(x_all_trials)
                these_R1s=[];
                for ii_sh=1:size(x_predicted_sh,2)

                    [R1,this_p_R1_hits]=corrcoef(x_all_trials,x_decod_all_trials_sh.ii_sh(ii_sh).xdec);
                    these_R1s=[these_R1s R1(1,2)];

                end
                R1_x_all_trials_sh(fileNo)=mean(these_R1s);
            else
                R1_x_all_trials_sh(fileNo)=NaN;
            end

            if ~isempty(y_all_trials)
                these_R1s=[];
                for ii_sh=1:size(y_predicted_sh,2)

                    [R1,this_p_R1_hits]=corrcoef(y_all_trials,y_decod_all_trials_sh.ii_sh(ii_sh).ydec);
                    these_R1s=[these_R1s R1(1,2)];

                end
                R1_y_all_trials_sh(fileNo)=mean(these_R1s);
            else
                R1_y_all_trials_sh(fileNo)=NaN;
            end



            if ~isempty(x_all_trials)
                these_R1XYs=[];
                for ii_sh=1:size(y_predicted_sh,2)
                    these_R1XYs=[these_R1XYs corr2([x_all_trials y_all_trials],[x_decod_all_trials_sh.ii_sh(ii_sh).xdec y_decod_all_trials_sh.ii_sh(ii_sh).ydec])];
                end
                R1_XY_all_trials_sh(fileNo)=mean(these_R1XYs);
                % P_rho_XY_all_trials(fileNo)=this_p_R1_hits(1,2);
            else
                R1_XY_all_trials_sh(fileNo)=NaN;
            end

        end
    end

    % figureNo = figureNo + 1;
    % try
    %     close(figureNo)
    % catch
    % end
    % hFig=figure(figureNo);
    % hold on
    % 
    % ax=gca;ax.LineWidth=3;
    % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    % edges=[-180:10:180];
    % histogram(all_meanAngles,edges)
    % title('Angle for final approach for all files')
    % xlabel('Angle (degrees)')
    % 
    % figureNo = figureNo + 1;
    % try
    %     close(figureNo)
    % catch
    % end
    % hFig=figure(figureNo);
    % hold on
    % 
    % ax=gca;ax.LineWidth=3;
    % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    % edges=[0:0.05:1];
    % histogram(fraction_other_angle,edges)
    % title('Fraction of horizontal approaches for all files')
    % xlabel('Fraction')

    %Calculate the percent of sessions above the exclusion criterion
    %Odor

    %r2ER
    total_sessions=0;
    sessions_p_less=0;
    these_groups=[1 5]; %2 cm and 1 cm all tirals
    for grNo=these_groups
        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            arena_file=handles_conc.arena_file{fileNo};
            if (handles_conc.group(fileNo)==grNo)&(fraction_other_angle(fileNo)<thr_froa)
                total_sessions=total_sessions+1;
                if file_exclusions(fileNo)==0
                    sessions_p_less=sessions_p_less+1;
                end
            end

        end
    end

    fprintf(1, ['\nPercent sessions with r2ER above exclusion criterion for odor trials ' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n'])
    fprintf(fileID, ['\nPercent sessions with r2ER above exclusion criterion for odor trials' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n']);

    odor_r2ER_percent_above_thr=100*sessions_p_less/total_sessions;

    %R1
    total_sessions=0;
    sessions_p_less=0;
    these_groups=[1 5]; %2 cm and 1 cm all tirals
    for grNo=these_groups
        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            arena_file=handles_conc.arena_file{fileNo};
            if (handles_conc.group(fileNo)==grNo)&(fraction_other_angle(fileNo)<thr_froa)
                total_sessions=total_sessions+1;
                if file_exclusions(fileNo)==0
                    sessions_p_less=sessions_p_less+1;
                end
            end

        end
    end

    fprintf(1, ['\nPercent sessions with R1 above exclusion criterion for odor trials ' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n'])
    fprintf(fileID, ['\nPercent sessions with R1 above exclusion criterion for odor trials' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n']);

    odor_percent_above_thr=100*sessions_p_less/total_sessions;

    %No odor

    %r2ER
    total_sessions=0;
    sessions_p_less=0;
    these_groups=[4];
    for grNo=these_groups
        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            arena_file=handles_conc.arena_file{fileNo};
            if (handles_conc.group(fileNo)==grNo)&(fraction_other_angle(fileNo)<thr_froa)
                total_sessions=total_sessions+1;
                if file_exclusions(fileNo)==0
                    sessions_p_less=sessions_p_less+1;
                end
            end

        end
    end

    fprintf(1, ['\nPercent sessions for r2ER with above exclusion criterion for no odor trials' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n'])
    fprintf(fileID, ['\nPercent sessions for r2ER with above exclusion criterion for no odor trials' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n']);
    no_odor_r2ER_percent_above_thr=100*sessions_p_less/total_sessions;

    %R1
    total_sessions=0;
    sessions_p_less=0;
    these_groups=[4];
    for grNo=these_groups
        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            arena_file=handles_conc.arena_file{fileNo};
            if (handles_conc.group(fileNo)==grNo)&(fraction_other_angle(fileNo)<thr_froa)
                total_sessions=total_sessions+1;
                if file_exclusions(fileNo)==0
                    sessions_p_less=sessions_p_less+1;
                end
            end

        end

    end

    fprintf(1, ['\nPercent sessions for R1 above exclusion criterion for no odor trials' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n'])
    fprintf(fileID, ['\nPercent sessions for R1 above exclusion criterion for no odor trials' num2str(100*sessions_p_less/total_sessions) ' , total sessions ' num2str(total_sessions)  '\n']);
    no_odor_percent_above_thr=100*sessions_p_less/total_sessions;
    % 
    % %Do the bar comparing percent above threshold
    % figureNo = figureNo + 1;
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
    % bar_offset=0;
    % bar(bar_offset,odor_percent_above_thr,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
    % 
    % bar_offset=1;
    % bar(bar_offset,no_odor_percent_above_thr,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
    % 
    % xticks([0 1])
    % xticklabels({'odor','No odor'})
    % 
    % 
    % title(['R1 for prediction of odor concentration'])
    % ylabel('percent above thr')
    % 
    % %Do the bar comparing percent above threshold for r2ER
    % figureNo = figureNo + 1;
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
    % bar_offset=0;
    % bar(bar_offset,odor_r2ER_percent_above_thr,'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
    % 
    % bar_offset=1;
    % bar(bar_offset,no_odor_r2ER_percent_above_thr,'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
    % 
    % xticks([0 1])
    % xticklabels({'odor','no odor'})
    % 
    % 
    % title(['r2ER for prediction of odor concentration'])
    % ylabel('Percent above thr')
    % 
    % %Do the R1 conc bar comparing 1 cm vs 2 cm with 0 or 16 bins before
    % figureNo = figureNo + 1;
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
    % bar_offset=0;
    % 
    % edges=[0:0.05:1];
    % rand_offset=0.5;
    % 
    % glm_r1=[];
    % glm_r1_ii=0;
    % 
    % id_r1_ii=0;
    % input_r1_data=[];
    % 
    % 
    % ii_run=1; %In the past I used to run multiple runs with different numbers of before bins
    % 
    % these_groups=[1 5]; %2 cm and 1 cm all tirals
    % 
    % 
    % for grNo=these_groups
    %     these_R1s=[];
    %     for fileNo=1:length(handles_conc.arena_file)
    %         if (handles_conc.group(fileNo)==grNo)&(fraction_other_angle(fileNo)<thr_froa)&(file_exclusions(fileNo)==0)
    %             these_R1s=[these_R1s R1_hits(fileNo)];
    %         end
    %     end
    % 
    %     %plot bar
    %     switch grNo
    %         case 1
    %             bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
    %         case 5
    %             bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
    %     end
    % 
    %     %Violin plot
    %     try
    %         [mean_out, CIout,violin_x]=drgViolinPoint(these_R1s...
    %             ,edges,bar_offset,rand_offset,'k','k',4);
    %     catch
    %         pffft=1;
    %     end
    %     bar_offset=bar_offset+1;
    % 
    %     glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
    %     glm_r1.group(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=grNo*ones(1,length(these_R1s));
    %     glm_r1.run(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_run*ones(1,length(these_R1s));
    %     glm_r1_ii=glm_r1_ii+length(these_R1s);
    % 
    %     id_r1_ii=id_r1_ii+1;
    %     input_r1_data(id_r1_ii).data=these_R1s;
    %     input_r1_data(id_r1_ii).description=[group_label{grNo} ' ' run_label{ii_run}];
    % end
    % bar_offset=bar_offset+1;
    % 
    % 
    % text(1.5,0.38,'2 cm','Color',[230/255 159/255 0/255],'FontWeight','bold')
    % text(1.5,0.35,'1 cm','Color',[86/255 180/255 233/255],'FontWeight','bold')
    % 
    % 
    % 
    % xticks([0 1])
    % xticklabels({'2 cm','1 cm'})
    % 
    % 
    % title(['R1 for prediction of odor concentration for hits'])
    % ylabel('R1')
    % ylim([-0.2 1])
    % xlim([-1 2])
    % 
    % %Perform the glm  for errors
    % fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1 cm from floor and bins before\n'])
    % fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1 cm from floor and bins before\n']);
    % %
    % % fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
    % % fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);
    % 
    % tbl = table(glm_r1.data',glm_r1.group',glm_r1.run',...
    %     'VariableNames',{'R1','group','run'});
    % mdl = fitglm(tbl,'R1~group+run'...
    %     ,'CategoricalVars',[2,3])
    % 
    % txt = evalc('mdl');
    % txt=regexp(txt,'<strong>','split');
    % txt=cell2mat(txt);
    % txt=regexp(txt,'</strong>','split');
    % txt=cell2mat(txt);
    % 
    % fprintf(fileID,'%s\n', txt);
    % 
    % 
    % %Do the ranksum/t-test
    % fprintf(1, ['\n\nRanksum or t-test p values for R1 cm from floor and bins before\n'])
    % fprintf(fileID, ['\n\nRanksum or t-test p values for R1 cm from floor and bins before\n']);
    % 
    % 
    % [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);
    % 
    % 
    % %R1 for conc vs XY
    % ii_run=1; %0 bins before
    % 
    % figureNo = figureNo + 1;
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
    % these_groups=[1 5];
    % these_new_groups.gr(1).groups=[1 5]; %2 cm and 1 cm
    % these_new_groups.gr(2).groups=[2 3]; %One lane odor
    % edges=[0:0.05:1];
    % rand_offset=1;
    % R1_per_mouse=[];
    % 
    % for grNo=1:1
    %     for ii_mouse=unique(handles_conc.mouse)
    %         R1_per_mouse.group(grNo).mouse(ii_mouse).R1_hits=[];
    %         R1_per_mouse.group(grNo).mouse(ii_mouse).pc=[];
    %         R1_per_mouse.group(grNo).mouse(ii_mouse).R1_XY_hits=[];
    %     end
    % end
    % 
    % grNo=1;
    % 
    % %R1 for conc all hits
    % 
    % these_R1_conc_hits=[];
    % percent_correct=[];
    % these_files_excluded=[];
    % 
    % for fileNo=1:length(handles_conc.arena_file)
    %     if (sum(handles_conc.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
    %         arena_file=handles_conc.arena_file{fileNo};
    %         %load the ouptut file
    %         load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
    %         these_R1_conc_hits=[these_R1_conc_hits R1_hits(fileNo)];
    %         R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).R1_hits=...
    %             [R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).R1_hits R1_hits(fileNo)];
    %         this_pc=100*(sum(handles_out.trials.hit1)+sum(handles_out.trials.hit4))/length(handles_out.trials.hit1);
    %         R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).pc=...
    %             [R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).pc this_pc];
    %         percent_correct=[percent_correct this_pc];
    %         if file_exclusions(fileNo)==0
    %             these_files_excluded=[these_files_excluded 0];
    %         else
    %             these_files_excluded=[these_files_excluded 1];
    %         end
    %     end
    % end
    % these_new_groups.gr(grNo).R1_hits=these_R1_conc_hits;
    % these_new_groups.gr(grNo).percent_correct=percent_correct;
    % 
    % %R1 for XY all trials
    % these_R1_XY_hits=[];
    % 
    % for fileNo=1:length(handles_XY.arena_file)
    %     if (sum(handles_XY.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
    %         these_R1_XY_hits=[these_R1_XY_hits R1_XY_hits(fileNo)];
    %         R1_per_mouse.group(grNo).mouse(handles_XY.mouse(fileNo)).R1_XY_hits=...
    %             [R1_per_mouse.group(grNo).mouse(handles_XY.mouse(fileNo)).R1_XY_hits R1_XY_hits(fileNo)];
    %     end
    % end
    % 
    % these_new_groups.gr(grNo).R1_XY_hits=these_R1_XY_hits;
    % 
    % for ii_excluded=0:1
    %     switch ii_excluded
    %         case 0
    %             plot(these_R1_conc_hits(these_files_excluded==0),these_R1_XY_hits(these_files_excluded==0),'ok','MarkerFaceColor','k')
    %         case 1
    %             plot(these_R1_conc_hits(these_files_excluded==1),these_R1_XY_hits(these_files_excluded==1),'ok','MarkerFaceColor',[0.7 0.7 0.7])
    %     end
    % end
    % 
    % p = polyfit(these_R1_conc_hits, these_R1_XY_hits, 1);        % Fit a first-degree polynomial (a line)
    % R1_XY_fit = polyval(p, these_R1_conc_hits);        % Evaluate the fitted line at your x data
    % 
    % plot(these_R1_conc_hits, R1_XY_fit, '-k','LineWidth',2)              % Plot original data as circles
    % 
    % 
    % % text(0.5,0.5,'Both spouts','Color','k')
    % % text(0.5,0.4,'One spout','Color',[0.7 0.7 0.7])
    % title('R1 for hits conc vs. R1 XY, gray=excluded')
    % xlabel('R1 conc')
    % ylabel('R1 XY')
    % 
    % [this_R1,this_P]=corrcoef(these_R1_conc_hits,these_R1_XY_hits);
    % fprintf(1, ['\nFor Fig ' num2str(figureNo) ' R1 between conc and XY is ' num2str(this_R1(1,2)) ' with p = ' num2str(this_P(1,2)) ', slope = ' num2str(p(1)) '\n'])
    % 
    % %R1 per mouse
    % figureNo = figureNo + 1;
    % try
    %     close(figureNo)
    % catch
    % end
    % hFig=figure(figureNo);
    % hold on
    % 
    % marker_per_group{1}=['' 'o' ''];
    % marker_per_group{2}=['' 's' ''];
    % 
    % color_okabe_ito{1}='[230/255 159/255 0/255]';
    % color_okabe_ito{2}='[86/255 180/255 233/255]';
    % color_okabe_ito{3}='[0/255 158/255 115/255]';
    % color_okabe_ito{4}='[240/255 228/255 66/255]';
    % color_okabe_ito{5}='[0/255 114/255 178/255]';
    % 
    % ax=gca;ax.LineWidth=3;
    % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    % 
    % 
    % grNo=1
    % 
    % for ii_mouse=unique(handles_conc.mouse)
    % 
    %     if ~isempty(R1_per_mouse.group(grNo).mouse(ii_mouse).R1_hits)
    %         these_R1conc=R1_per_mouse.group(grNo).mouse(ii_mouse).R1_hits;
    %         these_R1XY=R1_per_mouse.group(grNo).mouse(ii_mouse).R1_XY_hits;
    %         this_mean_R1conc=mean(these_R1conc);
    %         this_mean_R1XY=mean(these_R1XY);
    %         eval(['plot(this_mean_R1conc, this_mean_R1XY,''' marker_per_group{grNo} ''',''MarkerSize'',12,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_mouse} ')'])
    %         for ii_point=1:length(these_R1conc)
    %             eval(['plot(these_R1conc(ii_point), these_R1XY(ii_point),''' marker_per_group{grNo} ''',''MarkerSize'',6,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_mouse} ')'])
    %         end
    % 
    %     end
    % end
    % 
    % 
    % 
    % % text(0.5,0.5,'Both spouts','Color','k')
    % % text(0.5,0.4,'One spout','Color',[0.7 0.7 0.7])
    % title('R1 per mouse')
    % xlabel('R1 conc for hits')
    % ylabel('R1 XY')
    % 
    % 
    % %Plot R1 vs percent correct behavior
    % % group_labels{1}='odor in both spouts';
    % % group_labels{2}='odor in one spout';
    % for grNo=1:1
    %     figureNo = figureNo + 1;
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
    %     R1_conc=these_new_groups.gr(grNo).R1_hits;
    %     R1_XY=these_new_groups.gr(grNo).R1_XY_hits;
    %     percent_correct=these_new_groups.gr(grNo).percent_correct;
    % 
    %     plot(percent_correct,R1_conc,'o','MarkerEdgeColor','none','MarkerFaceColor',[230/255 159/255 0/255])
    %     plot(percent_correct,R1_XY,'o','MarkerEdgeColor','none','MarkerFaceColor',[86/255 180/255 233/255])
    % 
    %     text(60,0.7,'Odor','Color',[230/255 159/255 0/255])
    %     text(60,0.65,'XY','Color',[86/255 180/255 233/255])
    % 
    %     title(['R1 vs. percent correct behavior '])
    %     ylabel('R1')
    %     xlabel('Percent correct')
    % end
    % 

    %Calculate R1 bar graph for 1 and 2 cm all trials vs hit, miss
    ii_run=1;
    R1type=[];
    % R1type{1}='R1_all_trials';
    R1type{1}='R1_hits';
    % R1type{3}='R1_miss';
    R1type{2}='R1_all_trials_sh';



    R1type_label=[];
    % R1type_label{1}='all';
    R1type_label{1}='hits';
    % R1type_label{3}='miss';
    R1type_label{2}='shuffled';

    iiR1type_reorder=[1 2];
    % 
    % figureNo = figureNo + 1;
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
    % bar_offset=0;
    % 
    % edges=[0:0.05:1];
    % rand_offset=0.5;
    % 
    % glm_r1=[];
    % glm_r1_ii=0;
    % 
    % id_r1_ii=0;
    % input_r1_data=[];

    all_R1s=[];

    %Plot the different R1s
    these_groups=[1 5]; %1 and 2 cm all trials


    fprintf(1, ['\nFiles for R1 for prediction of odor all, hit, miss, shuffled:\n'])
    % include_odor=zeros(1,length(handles_conc.arena_file));
    for ii_R1type=1:length(R1type)

        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(file_exclusions(fileNo)==0)
                if ii_R1type==1
                    fprintf(1, [' ' num2str(fileNo)])
                end
                eval(['this_R1=' R1type{ii_R1type} '(fileNo);'])
                these_R1s=[these_R1s this_R1];
                % if ii_R1type==1
                %     include_odor(fileNo)=1;
                % end
            end

        end
        if ii_R1type==1
            fprintf(1, ['\n\n'])
        end
        % all_R1s.R1type(ii_R1type).R1=these_R1s;
        %plot bar
        switch ii_R1type
            case 1
            %     bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            % case 2
                bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
            % case 3
            %     bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
            case 2
                bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
            % case 5
            %     bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
        end

        %Violin plot
        [mean_out, CIout, all_R1s.R1type(ii_R1type).violin_x]=drgViolinPoint(these_R1s...
            ,edges,bar_offset,rand_offset,'k','k',4);
        bar_offset=bar_offset+1;

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
        glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=iiR1type_reorder(ii_R1type)*ones(1,length(these_R1s));
        glm_r1.algo(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_algo*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[algo_legend{ii_algo} ' ' R1type_label{ii_R1type}];

    end
    bar_offset=bar_offset+1;

    % %Plot lines between points
    % for ii=1:length(these_R1s)
    %     these_R1s=[];
    %     these_x=[];
    %     for ii_R1type=1:length(R1type)
    %         these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)];
    %         these_x=[these_x all_R1s.R1type(ii_R1type).violin_x(ii)];
    %     end
    %     plot(these_x,these_R1s,'-','Color',[0.7 0.7 0.7])
    % end


    % xticks([0 1 2 3])
    % xticklabels({'all','hit','miss','shuffled'})
    % 
    % 
    % title(['R1 for prediction of odor all, hit, miss, shuffled'])
    % ylabel('R1')
    % ylim([-0.25 1])
    % xlim([-1 4])
    % 
    % %Perform the glm  for prediction of odor all, hit, miss, shuffled
    % fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for R1 vs trial type\n']);
    % fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
    % fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for R1 vs trial type\n']);
    % fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
    % 
    % tbl = table(glm_r1.data',glm_r1.trial_type',...
    %     'VariableNames',{'R1','trial_type'});
    % mdl = fitglm(tbl,'R1~trial_type'...
    %     ,'CategoricalVars',[2])
    % 
    % txt = evalc('mdl');
    % txt=regexp(txt,'<strong>','split');
    % txt=cell2mat(txt);
    % txt=regexp(txt,'</strong>','split');
    % txt=cell2mat(txt);
    %
    % fprintf(fileID,'%s\n', txt);
    %
    %
    % %Do the ranksum/t-test
    % fprintf(1, ['\n\nRanksum or t-test p values for R1 vs trial typ\n'])
    % fprintf(fileID, ['\n\nRanksum or t-test p values for R1 vs. trial type\n']);
    %
    %
    % [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

end

xticks([0.5 3.5 6.5 9.5 12.5])
xticklabels({'ANN','GPR','BT','GLM','SVM'})


title(['R1 for prediction of odor (hit and shuffled)'])
ylabel('R1')
ylim([-0.4 1])
xlim([-1 14])

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for R1 vs algorithm\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for R1 vs algorithm\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);

tbl = table(glm_r1.data',glm_r1.trial_type',glm_r1.algo',...
    'VariableNames',{'R1','trial_type','algo'});
mdl = fitglm(tbl,'R1~trial_type+algo+trial_type*algo'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for R1 vs trial typ\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for R1 vs. trial type\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

fclose(fileID);

pffft=1;