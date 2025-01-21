%drgMini_batch_bootstrap_R1_constrained_range

close all
clear all

is_sphgpu=0;
is_pearson=1; %If this is 1 Pearson correlation is calculated, otherwise Spearman

switch is_sphgpu
    case 0

        %There was a bug and the shuffled runs were the same for all shuffled runs
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';
        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        %For the files below the shuffled runs should be different
        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01062025/';
        choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01062025.m';

        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end

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
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])

figureNo=0;

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);


%Now re-calculate the hits
these_groups=[1 5];
ii_run=1;
ii_for_corr=0;
R1_all_trials=[];
R1_hits=[];
% R1_hits90=[];
% R1_hitso=[];
R1_miss=[];
% R1_between=[];
R1_all_trials_sh=[];
% R2_all_trials=[];
% R2_hits=[];
% R2_miss=[];
% R2mod_hits=[];
% R2mod_miss=[];
% R2_all_trials_sh=[];
% Var_per_point_all=[];
% Var_per_point_hits=[];
% Var_per_point_miss=[];
% Var_per_point_all_sh=[];
% P_rho_all_trials=[];
these_fileNos=[];
no_neurons_per_file=[];
no_trials_per_file=[];

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
        op_all_trials=[];
        op_decod_all_trials=[];



        op_all_hits=[];
        op_decod_all_hits=[];
        % op_all_hits90=[];
        % op_decod_all_hits90=[];
        % op_all_hitso=[];
        % op_decod_all_hitso=[];
        % op_all_hits_sh=[];
        % op_decod_all_hits_sh=[];
        op_all_miss=[];
        op_decod_all_miss=[];
        % op_all_miss_sh=[];
        % op_decod_all_miss_sh=[];
        % op_decod_all_trials_sh=[];
        % op_between_trials=[];
        % op_decod_between_trials=[];
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
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        no_trials_per_file=[no_trials_per_file handles_out.trials.odor_trNo];
        odor_plume_template=handles_out.odor_plume_template;
        no_neurons_per_file=[no_neurons_per_file handles_out.no_neurons];
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
            % last_op_predictedend=op_predictedend;
            %
            % ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;


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

                    % %Hit 90 degrees
                    % if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                    %     op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                    %
                    % %Horizontal hit
                    % if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                    %     op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                case 2
                    %Lane 1 miss orange
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];

                case 3
                    %Lane 4 hit blue
                    op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];

                    %    %Hit 90 degrees
                    % if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                    %     op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                    %
                    % %Horizontal hit
                    % if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                    %     op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                case 4
                    %Lane 4 miss sky blue
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
            end
            ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
        end

        ii_for_corr=ii_for_corr+1;
        these_fileNos(ii_for_corr)=fileNo;

        if ~isempty(op_all_trials)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(op_all_trials,op_decod_all_trials);
                R1_all_trials(ii_for_corr)=R1(1,2);
                P_rho_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(op_all_trials,op_decod_all_trials,'Type','Spearman');
                R1_all_trials(ii_for_corr)=R1;
                P_rho_all_trials(ii_for_corr)=this_P_rho_all_trials;
            end
        else
            R1_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits)
            if is_pearson==1
                R1=corrcoef(op_all_hits,op_decod_all_hits);
                R1_hits(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_hits,op_decod_all_hits,'Type','Spearman');
                R1_hits(ii_for_corr)=R1;
            end
        else
            R1_hits(ii_for_corr)=NaN;
        end

        % if ~isempty(op_all_hits90)
        %     if is_pearson==1
        %         R1=corrcoef(op_all_hits90,op_decod_all_hits90);
        %         R1_hits90(ii_for_corr)=R1(1,2);
        %     else
        %         R1=corr(op_all_hits90,op_decod_all_hits90,'Type','Spearman');
        %         R1_hits90(ii_for_corr)=R1;
        %     end
        % else
        %     R1_hits90(ii_for_corr)=NaN;
        % end
        %
        % if ~isempty(op_all_hitso)
        %     if is_pearson==1
        %         R1=corrcoef(op_all_hitso,op_decod_all_hitso);
        %         R1_hitso(ii_for_corr)=R1(1,2);
        %     else
        %         R1=corr(op_all_hitso,op_decod_all_hitso,'Type','Spearman');
        %         R1_hitso(ii_for_corr)=R1;
        %     end
        % else
        %     R1_hitso(ii_for_corr)=NaN;
        % end

        if ~isempty(op_all_miss)
            if is_pearson==1
                R1=corrcoef(op_all_miss,op_decod_all_miss);
                R1_miss(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_miss,op_decod_all_miss,'Type','Spearman');
                R1_miss(ii_for_corr)=R1;
            end
        else
            R1_miss(ii_for_corr)=NaN;
        end

        %Shuffled
        if ~isempty(op_predicted_sh_conv)
            these_R1s=[];

            for ii_sh=1:size(op_predicted_sh_conv,2)
                this_op_decode_all_trials_sh=op_decod_trials_sh.ii_sh(ii_sh).oppsh;
                this_R1=corrcoef(op_all_trials,this_op_decode_all_trials_sh);
                these_R1s=[these_R1s this_R1(1,2)];
            end
            R1_all_trials_sh(ii_for_corr)=mean(these_R1s);
        else
            R1_all_trials_sh(ii_for_corr)=NaN;
        end
        % if ~isempty(op_between_trials)
        %     if is_pearson==1
        %         R1=corrcoef(op_between_trials,op_decod_between_trials);
        %         R1_between(ii_for_corr)=R1(1,2);
        %     else
        %         R1=corr(op_between_trials,op_decod_between_trials,'Type','Spearman');
        %         R1_between(ii_for_corr)=R1;
        %     end
        % else
        %     R1_between(ii_for_corr)=NaN;
        % end
        %
        % %Calculate R2
        % if ~isempty(op_all_trials)
        %     this_R2=1-sum((op_decod_all_trials-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
        %     R2_all_trials(ii_for_corr)=drgMini_rectifyR2(this_R2);
        %     Var_per_point_all(ii_for_corr)=sum((op_all_trials-mean(op_all_trials)).^2)/length(op_all_trials);
        % else
        %     R2_all_trials(ii_for_corr)=NaN;
        % end
        %
        % % if ~isempty(op_all_hits)
        % %     this_R2=1-sum((op_decod_all_hits-op_all_hits).^2)/( ((length(op_all_hits)/length(op_all_trials))...
        % %         *sum((op_all_trials-mean(op_all_trials)).^2)));
        % %     R2mod_hits(ii_for_corr)=drgMini_rectifyR2(this_R2);
        % % else
        % %     R2mod_hits(ii_for_corr)=NaN;
        % % end
        % %
        % % if ~isempty(op_all_miss)
        % %     this_R2=1-sum((op_decod_all_miss-op_all_miss).^2)/( ((length(op_all_miss)/length(op_all_trials))...
        % %         *sum((op_all_trials-mean(op_all_trials)).^2)));
        % %     R2mod_miss(ii_for_corr)=drgMini_rectifyR2(this_R2);
        % % else
        % %     R2mod_miss(ii_for_corr)=NaN;
        % % end
        %
        % if ~isempty(op_all_hits)
        %     this_R2=1-sum((op_decod_all_hits-op_all_hits).^2)/(sum((op_all_hits-mean(op_all_hits)).^2));
        %     R2_hits(ii_for_corr)=drgMini_rectifyR2(this_R2);
        %     Var_per_point_hits(ii_for_corr)=(sum((op_all_hits-mean(op_all_hits)).^2))/length(op_all_hits);
        % else
        %     R1_hits(ii_for_corr)=NaN;
        % end
        %
        % if ~isempty(op_all_miss)
        %     this_R2=1-sum((op_decod_all_miss-op_all_miss).^2)/(sum((op_all_miss-mean(op_all_miss)).^2));
        %     R2_miss(ii_for_corr)=drgMini_rectifyR2(this_R2);
        %     Var_per_point_miss(ii_for_corr)=(sum((op_all_miss-mean(op_all_miss)).^2))/length(op_all_miss);
        % else
        %     R2_miss(ii_for_corr)=NaN;
        % end


    end
end

%Do the R1 bar graph for 1 and 2 cm all trials vs hit, miss
ii_run=1;
R1type=[];
R1type{1}='R1_all_trials';
R1type{2}='R1_hits';
R1type{3}='R1_miss';
R1type{4}='R1_all_trials_sh';

R1type_label=[];
R1type_label{1}='all';
R1type_label{2}='hits';
R1type_label{3}='miss';
R1type_label{4}='shuffled';

iiR1type_reorder=[2 1 3 4];

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

all_R1s=[];
%Plot the different R1s
these_groups=[1 5]; %1 and 2 cm all trials
R1_all_trials_orig=[];
R1_hits_orig=[];
R1_miss_orig=[];
R1_shuffle_orig=[];

for ii_R1type=1:length(R1type)

    these_R1s=[];

    eval(['these_R1s=' R1type{ii_R1type} ';'])

    all_R1s.R1type(ii_R1type).R1=these_R1s;
    %plot bar
    switch ii_R1type
        case 1
            bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, all_R1s.R1type(ii_R1type).violin_x]=drgViolinPoint(these_R1s...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=iiR1type_reorder(ii_R1type)*ones(1,length(these_R1s));
    glm_r1_ii=glm_r1_ii+length(these_R1s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R1s;
    input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

end
bar_offset=bar_offset+1;

%Plot lines between points
for ii=1:length(these_R1s)
    these_R1s=[];
    these_x=[];
    for ii_R1type=1:length(R1type)
        these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)];
        these_x=[these_x all_R1s.R1type(ii_R1type).violin_x(ii)];
    end
    plot(these_x,these_R1s,'-','Color',[0.7 0.7 0.7])
end


xticks([0 1 2 3])
xticklabels({'within','hit','miss','shuffled'})


title(['R1 for prediction of odor hit, miss, shuffled'])
ylabel('R1')
ylim([-0.2 0.65])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig 11 for R1 vs trial type\n'])
fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig 11 for R1 vs trial type\n']);
fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n'])
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.trial_type',...
    'VariableNames',{'R1','trial_type'});
mdl = fitglm(tbl,'R1~trial_type'...
    ,'CategoricalVars',[2])

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

%Now perform range-constrained bootstrapping

%Now re-calculate the hits
these_groups=[1 5];
ii_run=1;
ii_for_corr=0;
R1_all_trials=[];
R1_hits=[];
% R1_hits90=[];
% R1_hitso=[];
R1_miss=[];
% R1_between=[];
R1_all_trials_sh=[];
% R2_all_trials=[];
% R2_hits=[];
% R2_miss=[];
% R2mod_hits=[];
% R2mod_miss=[];
% R2_all_trials_sh=[];
% Var_per_point_all=[];
% Var_per_point_hits=[];
% Var_per_point_miss=[];
% Var_per_point_all_sh=[];
% P_rho_all_trials=[];
these_fileNos=[];


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
        op_all_trials=[];
        op_decod_all_trials=[];



        op_all_hits=[];
        op_decod_all_hits=[];
        % op_all_hits90=[];
        % op_decod_all_hits90=[];
        % op_all_hitso=[];
        % op_decod_all_hitso=[];
        % op_all_hits_sh=[];
        % op_decod_all_hits_sh=[];
        op_all_miss=[];
        op_decod_all_miss=[];
        % op_all_miss_sh=[];
        % op_decod_all_miss_sh=[];
        % op_decod_all_trials_sh=[];
        % op_between_trials=[];
        % op_decod_between_trials=[];
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
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
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
            % last_op_predictedend=op_predictedend;
            %
            % ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;


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

                    % %Hit 90 degrees
                    % if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                    %     op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                    %
                    % %Horizontal hit
                    % if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                    %     op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                case 2
                    %Lane 1 miss orange
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];

                case 3
                    %Lane 4 hit blue
                    op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];

                    %    %Hit 90 degrees
                    % if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                    %     op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                    %
                    % %Horizontal hit
                    % if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                    %     op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
                    %     op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
                    % end
                case 4
                    %Lane 4 miss sky blue
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
            end
            ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
        end

        ii_for_corr=ii_for_corr+1;
        these_fileNos(ii_for_corr)=fileNo;




     

        %calculate the sum of squared residuals per time
        SSR_hits=sum((op_all_hits-op_decod_all_hits).^2)/length(op_all_hits);
        SSR_misss=sum((op_all_miss-op_decod_all_miss).^2)/length(op_all_miss);
        SSR_mean_hits=sum((mean(op_all_hits)-op_all_hits).^2)/length(op_all_hits);
        SSR_mean_miss=sum((mean(op_all_miss)-op_all_miss).^2)/length(op_all_miss);
        SSR_mean_all=sum((mean(op_all_trials)-op_all_trials).^2)/length(op_all_trials);

        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .5 .6 .3])

        plot([1:length(op_all_hits)],op_all_hits,'-b')
        plot([1:length(op_all_hits)],op_decod_all_hits,'-k')

        plot([1+1.2*length(op_all_hits):1.2*length(op_all_hits)+length(op_all_miss)],op_all_miss,'-b')
        plot([1+1.2*length(op_all_hits):1.2*length(op_all_hits)+length(op_all_miss)],op_decod_all_miss,'-k')

        min_all=min(op_all_trials);
        max_all=max(op_all_trials);
        delta_all=(max_all-min_all)/10;
        odor_steps=zeros(1,10);
        ii_odor=min_all;
        prediction_mean_hits=[];
        odor_hit=[]
        prediction_mean_miss=[];
        odor_miss=[];
        ii_included_hit=0;
        ii_included_miss=0;
        for ii=1:10
            if sum((op_all_hits>=ii_odor)&(op_all_hits<ii_odor+delta_all))>0
                ii_included_hit=ii_included_hit+1;
                prediction_mean_hits(ii_included_hit)=mean(op_decod_all_hits((op_all_hits>=ii_odor)&(op_all_hits<ii_odor+delta_all)));
                odor_hit(ii_included_hit)=ii_odor+(delta_all/2);
            end
            if sum((op_all_miss>=ii_odor)&(op_all_miss<ii_odor+delta_all))
                ii_included_miss=ii_included_miss+1;
                prediction_mean_miss(ii_included_miss)=mean(op_decod_all_miss((op_all_miss>=ii_odor)&(op_all_miss<ii_odor+delta_all)));
                odor_miss(ii_included_miss)=ii_odor+(delta_all/2);
            end
            ii_odor=ii_odor+delta_all;
        end

        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .1 .3 .3])

        plot(op_all_hits,op_decod_all_hits,'.b')
        plot(odor_hit,prediction_mean_hits,'ob','MarkerSize',10,'MarkerFaceColor','b')
        title('Hits')

        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.5 .1 .3 .3])

        plot(op_all_miss,op_decod_all_miss,'.b')
        plot(odor_miss,prediction_mean_miss,'ob','MarkerSize',10,'MarkerFaceColor','b')
        title('Misses')

         figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .1 .3 .3])

        edges=[min_all:(max_all-min_all)/20:max_all];
        histogram(op_all_hits,edges)
        this_ylim=ylim;
        text(-6,this_ylim(1)+0.5*(this_ylim(2)-this_ylim(1)),['Var = ' num2str(var(op_all_hits))])
        title('Odor concentraiton hits')

         figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.5 .1 .3 .3])

        edges=[min_all:(max_all-min_all)/20:max_all];
        histogram(op_all_miss,edges)
        this_ylim=ylim;
        text(-6,this_ylim(1)+0.5*(this_ylim(2)-this_ylim(1)),['Var = ' num2str(var(op_all_miss))])
        title('Odor concentraiton miss')

        [p_value(ii_for_corr), ii_attempted_boots, no_boots] = drgMin_range_restricted_bootstrap(op_all_hits,op_decod_all_hits,...
            op_all_miss,op_decod_all_miss,0);

        fprintf(['\n\nFile number ' num2str(fileNo) '\n']);
        fprintf('P-value (lowvar correlation >= bootstrapped highvar correlations): %.4f\n', p_value);
        fprintf(['Attempted bootstraps ' num2str(ii_attempted_boots) '\n']);
        fprintf(['Succesful bootstraps ' num2str(ii_attempted_boots) '\n']);

        pffft=1;

    end
end
pffft=1;