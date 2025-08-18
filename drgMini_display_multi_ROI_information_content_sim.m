%drgMini_display_multi_ROI_information_content_sim
close all
clear all

% %Trained with hits only
% save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
% choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

%Trained with hits only and taking on accoount odor on
save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc04192024/';
choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_04192024.m'

% save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
% choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

%This one has the dFF per trial
save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';


choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
fileID = fopen([choiceBatchPathName 'multiROI_info_stats.txt'],'w');

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])

% outputFile='multi_ROI_info_content_best.mat';
outputPath='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
outputFile='multi_ROI_info_content_dynamic_sim.mat';
load([outputPath outputFile])

figureNo=0;
these_groups=[1 5];
ii_run=1;
R1s=[];

% all_n_bits=[1 3 6 12 24 2000]; %Number of bits (ROIs) for most and least important prediction ROIs
%                             %And for best

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

%Now do analysis for z-scored MI
miss_div_hit=[];
for ii_alpha=1:length(handles_out2.alpha)
    all_hit_z_MIs=[];
    all_miss_z_MIs=[];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)



            %hit MI
            data=handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).sh_hit_mi;
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

            handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).ssi_hit_mi=(handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).hit_mi-mean(data))/estimated_std;

            all_hit_z_MIs=[all_hit_z_MIs handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).ssi_hit_mi];
            %miss MI
            data=handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).sh_miss_mi;
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

            handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).ssi_miss_mi=(handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).miss_mi-mean(data))/estimated_std;
            all_miss_z_MIs=[all_miss_z_MIs handles_out2.hit_miss.file(fileNo).ii_alpha(ii_alpha).ssi_miss_mi];

        end

    end

    miss_div_hit=[miss_div_hit mean(all_miss_z_MIs(~isnan(all_miss_z_MIs)))/mean(all_hit_z_MIs(~isnan(all_hit_z_MIs)))];

    if ii_alpha<11
        %Now plot MIs comparing all trials with hits, miss and shuffled

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

        these_MIs=all_hit_z_MIs;
        edges=[0:0.05*max(these_MIs(~isnan(these_MIs))):max(these_MIs(~isnan(these_MIs)))];
        rand_offset=0.5;

        glm_r1=[];
        glm_r1_ii=0;

        id_r1_ii=0;
        input_r1_data=[];

        MI_labels{1}='All trials';
        MI_labels{2}='Hits';
        MI_labels{3}='Misses';
        MI_labels{4}='Shuffled';





        bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_MIs(~isnan(these_MIs))...
            ,edges,bar_offset,rand_offset,'k','k',4);
        bar_offset=bar_offset+1;

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_MIs))=these_MIs;
        glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_MIs))=0*ones(1,length(these_MIs));
        glm_r1_ii=glm_r1_ii+length(these_MIs);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_MIs;
        input_r1_data(id_r1_ii).description=['Hit'];

        %Miss
        these_MIs=all_miss_z_MIs;

        bar(bar_offset,mean(these_MIs(~isnan(these_MIs))),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_MIs(~isnan(these_MIs))...
            ,edges,bar_offset,rand_offset,'k','k',4);
        bar_offset=bar_offset+1;

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_MIs(~isnan(these_MIs))))=these_MIs(~isnan(these_MIs));
        glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_MIs(~isnan(these_MIs))))=1*ones(1,length(these_MIs(~isnan(these_MIs))));
        glm_r1_ii=glm_r1_ii+length(these_MIs(~isnan(these_MIs)));

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_MIs(~isnan(these_MIs));
        input_r1_data(id_r1_ii).description=['Miss'];


        % %Plot lines between points
        % for ii=1:length(all_miss_z_MIs)
        %     plot([0 1],[all_hit_z_MIs(ii) all_miss_z_MIs(ii)],'-','Color',[0.7 0.7 0.7])
        %     plot(zeros(1,length(all_hit_z_MIs)),all_hit_z_MIs,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4)
        %     plot(ones(1,length(all_miss_z_MIs)),all_miss_z_MIs,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4)
        % end


        xticks([0 1])
        xticklabels({'hit','miss'})


        title(['multi ROI MI betweeen op and binary dFF, alpha= ' num2str(handles_out2.alpha(ii_alpha))])
        ylabel('z-scored MI')
        % ylim([-0.2 2.5])
        % xlim([-1 4])

        %Perform the glm  for MIs
        fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for MI vs trial type\n'])
        fprintf(1, ['\n1 Hit, 2 Miss, 3 All, 4 Shuffled\n'])
        fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for MI vs trial type\n']);
        fprintf(fileID, ['\n1 Hit, 2 Miss, 3 All, 4 Shuffled\n'])


        tbl = table(glm_r1.data',glm_r1.trial_type',...
            'VariableNames',{'MI','trial_type'});
        mdl = fitglm(tbl,'MI~trial_type'...
            ,'CategoricalVars',[2])

        txt = evalc('mdl');
        txt=regexp(txt,'<strong>','split');
        txt=cell2mat(txt);
        txt=regexp(txt,'</strong>','split');
        txt=cell2mat(txt);

        fprintf(fileID,'%s\n', txt);


        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for M1 vs trial typ\n'])
        fprintf(fileID, ['\n\nRanksum or t-test p values for M1 vs. trial type\n']);


        [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);
    end
end

experimental_value = 0.093648; % your experimental value here

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

histogram(miss_div_hit)
% xlim([0 0.5])
% plot([experimental_value experimental_value],[0 50],'-r','LineWidth',2)
title('Ratio of Miss zMI divided by Hit zMI')
xlabel('Miss zMI/Hit z MI')
ylabel('Count')

% Assume miss_div_hit is your 200x1 vector of values
nBoot = 1000; % Number of bootstrap resamples

% Bootstrap resampling: each row in bootstat is a bootstrap sample
bootstat = bootstrp(nBoot, @(x) prctile(x,5), miss_div_hit); % [1]

% The bootstrapped 5th percentiles are in bootstat
% To get the estimated 5th percentile (the percentile of percentiles):
percentile_5 = prctile(bootstat, 5); % [6]

% To compare an experimental value:

is_below_5th = experimental_value < percentile_5;

fprintf('Bootstrapped 5th percentile: %.4f\n', percentile_5);
if is_below_5th
    fprintf('Experimental value is below the bootstrapped 5th percentile.\n');
else
    fprintf('Experimental value is NOT below the bootstrapped 5th percentile.\n');
end

pffft=1;

