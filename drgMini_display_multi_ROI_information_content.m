%drgMini_display_multi_ROI_information_content
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


choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
fileID = fopen([choiceBatchPathName 'multiROI_info_stats.txt'],'w');

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])

% outputFile='multi_ROI_info_content_best.mat';
outputFile='multi_ROI_info_content_dynamic.mat';
load([choiceBatchPathName outputFile])

figureNo=0;
these_groups=[1 5];
ii_run=1;
R1s=[];

all_n_bits=[1 3 6 12 24 2000]; %Number of bits (ROIs) for most and least important prediction ROIs
                            %And for best

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

%Plot the best combined mi

%Plot best mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.best.file(fileNo).ii_neuron(ii).best_mi;
        end

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})
 
xlabel('Number of ROIs')
ylabel('Mutual Information')
title('Best multi ROI MI betweeen op and binary dFF')

%Plot worst mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.best.file(fileNo).ii_neuron(ii).worst_mi;
        end

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})

xlabel('Number of ROIs')
ylabel('Mutual Information')
title('Worst multi ROI MI between op and binary dFF')


%Plot most predicitive importance mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).most_imp_mi;
        end

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})

xlabel('Number of ROIs')
ylabel('Mutual Information')
title('Most pred imp multi ROI MI betweeen op and binary dFF')

%Plot least predicitve importance mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).least_imp_mi;
        end

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})

xlabel('Number of ROIs')
ylabel('Mutual Information')
title('Least pred imp multi ROI MI between op and binary dFF')

%Plot the relationship between most and least predicitve importance mutual information between bindFF and op for most important
%Exclude all ROIs
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_least_info_op_bin=[];
        these_most_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_most_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).most_imp_mi;
            these_least_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).least_imp_mi;
        end

        ii_color=ii_color+1;

        for ii=3:length(all_n_bits)-1
            switch ii_color
                case 1
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
        end
            if ii_color==5
                ii_color=0;
            end
            
    end
end

plot([0 1],[0 1],'-k')

xlim([0 0.3])
ylim([0 0.3])

xlabel('MI least imp')
ylabel('MI most imp')
title('multi ROI MI between op and binary dFF most vs. least important')

%Now plot the relationship between conc R1 and all ROI best mutual
        

%Plot relationship between R1 for decoding and MIs
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

MI_max=1.5;

ii_color=0;
R1s_all=[];
R1s_hit=[];
R1s_miss=[];
MIs=[];
MIs_hit=[];
MIs_miss=[];
MIs_sh=[];
all_no_neurons=[];
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        arena_file=handles_conc.arena_file{fileNo};
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        R1s_all=[R1s_all handles_out.R1.all_trials];
        R1s_hit=[R1s_hit handles_out.R1.all_hits];
        R1s_miss=[R1s_miss handles_out.R1.all_miss];

        ii=length(all_n_bits);
        MIs=[MIs handles_out2.best.file(fileNo).ii_neuron(ii).best_mi/handles_out2.file(fileNo).no_neurons];
        MIs_hit=[MIs_hit handles_out2.hit_miss.file(fileNo).hit_mi/handles_out2.file(fileNo).no_neurons];
        MIs_miss=[MIs_miss handles_out2.hit_miss.file(fileNo).miss_mi/handles_out2.file(fileNo).no_neurons];
        MIs_sh=[MIs_sh mean(handles_out2.best.file(fileNo).ii_neuron(length(all_n_bits)).sh_best_mi)/handles_out2.file(fileNo).no_neurons];
        all_no_neurons=[all_no_neurons handles_out2.file(fileNo).no_neurons];
    end
end

MIs=mean(all_no_neurons)*MIs;
MIs_hit=mean(all_no_neurons)*MIs_hit;
MIs_miss=mean(all_no_neurons)*MIs_miss;
MIs_sh=mean(all_no_neurons)*MIs_sh;


all_MIs.type(1).mi=MIs;
all_MIs.type(2).mi=MIs_hit;
all_MIs.type(3).mi=MIs_miss;
all_MIs.type(4).mi=MIs_sh;

plot(MIs,R1s_all,'o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
p = polyfit(MIs, R1s_all, 1);
R1s_all_fit = polyval(p, [0 MI_max]);
plot([0 MI_max],R1s_all_fit,'-','Color',[230/255 159/255 0/255],'LineWidth',3)
[R1, P]=corrcoef(MIs,R1s_all);
fprintf(1, 'Correlation coefficient for all trials %d, p value %d\n\n',R1(1,2),P(1,2));

plot(MIs,R1s_hit,'o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
p = polyfit(MIs, R1s_hit, 1);
R1s_all_fit = polyval(p, [0 MI_max]);
plot([0 MI_max],R1s_all_fit,'-','Color',[86/255 180/255 233/255],'LineWidth',3)
[R1,P]=corrcoef(MIs,R1s_hit);
fprintf(1, 'Correlation coefficient for hits %d, p value %d\n\n',R1(1,2),P(1,2));

plot(MIs,R1s_miss,'o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
p = polyfit(MIs, R1s_miss, 1);
R1s_all_fit = polyval(p, [0 MI_max]);
plot([0 MI_max],R1s_all_fit,'-','Color',[0/255 158/255 115/255],'LineWidth',3)
[R1,P]=corrcoef(MIs,R1s_miss);
fprintf(1, 'Correlation coefficient for miss %d, p value %d\n\n',R1(1,2),P(1,2));

xlabel('MIs')
ylabel('R1')
title('R1 for decoding vs. multi ROI MI betweeen op and binary dFF')

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

edges=[0:0.05:1];
rand_offset=0.5;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

MI_labels{1}='All trials';
MI_labels{2}='Hits';
MI_labels{3}='Misses';
MI_labels{4}='Shuffled';

type_order=[3 1 2 4];

for ii_type=1:length(all_MIs.type)
    these_MIs=all_MIs.type(ii_type).mi;
    %plot bar
    switch ii_type
        case 1
            bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, all_MIs.type(ii_type).violin_x]=drgViolinPoint(these_MIs...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_MIs))=these_MIs;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_MIs))=type_order(ii_type)*ones(1,length(these_MIs));
    glm_r1_ii=glm_r1_ii+length(these_MIs);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_MIs;
    input_r1_data(id_r1_ii).description=[MI_labels{type_order(ii_type)}];

end

%Plot lines between points
for ii=1:length(MIs)
    these_MIs=[];
    these_x=[];
    for ii_type=1:length(all_MIs.type)
        these_MIs=[these_MIs all_MIs.type(ii_type).mi(ii)];
        these_x=[these_x all_MIs.type(ii_type).violin_x(ii)];
    end
    plot(these_x,these_MIs,'-','Color',[0.7 0.7 0.7])
end


xticks([0 1 2 3])
xticklabels({'within','hit','miss','shuffled'})


title(['multi ROI MI betweeen op and binary dFF'])
ylabel('MI in bits')
ylim([-0.2 2.5])
xlim([-1 4])

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

%Now do the same analysis for z-scored MI
all_hit_z_MIs=[];
all_miss_z_MIs=[];
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)


        for ii=1:length(all_n_bits)

            %Best MI
            data=handles_out2.best.file(fileNo).ii_neuron(ii).sh_best_mi;
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

            handles_out2.best.file(fileNo).ii_neuron(ii).ssi_best_mi=(handles_out2.best.file(fileNo).ii_neuron(ii).best_mi-mean(data))/estimated_std;

             if ii==length(all_n_bits) 
                %Do hit z score
                this_hit_z_mi=(handles_out2.hit_miss.file(fileNo).hit_mi-mean(data))/estimated_std;
                all_hit_z_MIs=[all_hit_z_MIs this_hit_z_mi];

                %Do miss z score
                this_miss_z_mi=(handles_out2.hit_miss.file(fileNo).miss_mi-mean(data))/estimated_std;
                all_miss_z_MIs=[all_miss_z_MIs this_miss_z_mi];

            end
            %Worst MI
             data=handles_out2.best.file(fileNo).ii_neuron(ii).sh_worst_mi;
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

            handles_out2.best.file(fileNo).ii_neuron(ii).ssi_worst_mi=(handles_out2.best.file(fileNo).ii_neuron(ii).worst_mi-mean(data))/estimated_std;

            %Most important MI
            data=handles_out2.imps.file(fileNo).ii_neuron(ii).sh_most_imp_mi;
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

            handles_out2.imps.file(fileNo).ii_neuron(ii).ssi_most_imp_mi=(handles_out2.imps.file(fileNo).ii_neuron(ii).most_imp_mi-mean(data))/estimated_std;

            %Least important MI
             data=handles_out2.imps.file(fileNo).ii_neuron(ii).sh_least_imp_mi;
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

            handles_out2.imps.file(fileNo).ii_neuron(ii).ssi_least_imp_mi=(handles_out2.imps.file(fileNo).ii_neuron(ii).least_imp_mi-mean(data))/estimated_std;

         
        end

    end
end


%Plot best mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;
max_best_mi=[];
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.best.file(fileNo).ii_neuron(ii).ssi_best_mi;
        end

        max_best_mi=[max_best_mi max(these_info_op_bin)];

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

plot([1 length(all_n_bits)],[3 3],'-k','LineWidth',2)
plot([1 length(all_n_bits)],[1 1],'-k')
xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})
ylim([-5 60]) 
xlim([0 7])
xlabel('Number of ROIs')
ylabel('z-scored MI')
title('Best multi ROI MI betweeen op and binary dFF')

%Plot worst mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.best.file(fileNo).ii_neuron(ii).ssi_worst_mi;
        end

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

plot([1 length(all_n_bits)],[3 3],'-k','LineWidth',2)
plot([1 length(all_n_bits)],[1 1],'-k')
xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})
ylim([-5 60]) 
xlim([0 7])
xlabel('Number of ROIs')
ylabel('z-scored MI')
title('Worst multi ROI MI between op and binary dFF')


%Plot most predicitive importance mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).ssi_most_imp_mi;
        end

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

plot([1 length(all_n_bits)],[3 3],'-k','LineWidth',2)
plot([1 length(all_n_bits)],[1 1],'-k')
xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})
ylim([-5 60]) 
xlim([0 7])
xlabel('Number of ROIs')
ylabel('z-scored MI')
title('Most predictive importance multi ROI MI betweeen op and binary dFF')

%Plot least predicitve importance mutual information between bindFF and op for most important
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).ssi_least_imp_mi;
        end

        ii_color=ii_color+1;

        switch ii_color
            case 1
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot([1:length(all_n_bits)],these_info_op_bin,'-o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
            if ii_color==5
                ii_color=0;
            end
    end
end

plot([1 length(all_n_bits)],[3 3],'-k','LineWidth',2)
plot([1 length(all_n_bits)],[1 1],'-k')
xticks([1:length(all_n_bits)])
xticklabels({'1','3','6','12','24','All'})
ylim([-5 60]) 
xlim([0 7])

xlabel('Number of ROIs')
ylabel('z-scored MI')
title('Least predictive importance multi ROI MI between op and binary dFF')

%Plot the relationship between most and least predicitve importance mutual information between bindFF and op for most important
%Do this only for 3:end-1
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on


ii_color=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        these_least_info_op_bin=[];
        these_most_info_op_bin=[];
        for ii=1:length(all_n_bits)
            these_most_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).ssi_most_imp_mi;
            these_least_info_op_bin(ii)=handles_out2.imps.file(fileNo).ii_neuron(ii).ssi_least_imp_mi;
        end

        ii_color=ii_color+1;

        for ii=3:length(all_n_bits)-1
            switch ii_color
                case 1
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
                case 2
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)
                case 3
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)
                case 4
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[240/255 228/255 66/255],'MarkerEdgeColor',[240/255 228/255 66/255],'Color',[240/255 228/255 66/255],'MarkerSize',8)
                case 5
                    plot(these_least_info_op_bin(ii),these_most_info_op_bin(ii),'o','MarkerFaceColor',[0/255 114/255 178/255],'MarkerEdgeColor',[0/255 114/255 178/255],'Color',[0/255 114/255 178/255],'MarkerSize',8)
            end
        end
            if ii_color==5
                ii_color=0;
            end
            
    end
end

plot([0 30],[0 30],'-k')

xlim([0 30])
ylim([0 30])

xlabel('MI least imp')
ylabel('MI most imp')
title('multi ROI MI between op and binary dFF most vs. least important')

%Now plot the relationship between conc R1 and all ROI best mutual
        

%Plot relationship between R1 for decoding and z scored MIs
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

MI_max=1.5;

ii_color=0;
R1s_all=[];
R1s_hit=[];
R1s_miss=[];
MIs=[];
MIs_hit=[];
MIs_miss=[];
MIs_sh=[];
all_no_neurons=[];
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        arena_file=handles_conc.arena_file{fileNo};
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        R1s_all=[R1s_all handles_out.R1.all_trials];
        R1s_hit=[R1s_hit handles_out.R1.all_hits];
        R1s_miss=[R1s_miss handles_out.R1.all_miss];

        % ii=length(all_n_bits);
        % MIs=[MIs handles_out2.best.file(fileNo).ii_neuron(ii).ssi_best_mi/handles_out2.file(fileNo).no_neurons];
        % MIs_hit=[MIs_hit handles_out2.hit_miss.file(fileNo).hit_mi/handles_out2.file(fileNo).no_neurons];
        % MIs_miss=[MIs_miss handles_out2.hit_miss.file(fileNo).miss_mi/handles_out2.file(fileNo).no_neurons];
        % MIs_sh=[MIs_sh mean(handles_out2.best.file(fileNo).ii_neuron(length(all_n_bits)).sh_best_mi)/handles_out2.file(fileNo).no_neurons];
        % all_no_neurons=[all_no_neurons handles_out2.file(fileNo).no_neurons];
    end
end

% MIs=mean(all_no_neurons)*MIs;
% MIs_hit=mean(all_no_neurons)*MIs_hit;
% MIs_miss=mean(all_no_neurons)*MIs_miss;
% MIs_sh=mean(all_no_neurons)*MIs_sh;
% 
% 
% all_MIs.type(1).mi=MIs;
% all_MIs.type(2).mi=MIs_hit;
% all_MIs.type(3).mi=MIs_miss;
% all_MIs.type(4).mi=MIs_sh;

% plot(MIs,R1s_all,'o','MarkerFaceColor',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'Color',[230/255 159/255 0/255],'MarkerSize',8)
% p = polyfit(MIs, R1s_all, 1);
% R1s_all_fit = polyval(p, [0 MI_max]);
% plot([0 MI_max],R1s_all_fit,'-','Color',[230/255 159/255 0/255],'LineWidth',3)
% [R1, P]=corrcoef(MIs,R1s_all);
% fprintf(1, 'Correlation coefficient for all trials %d, p value %d\n\n',R1(1,2),P(1,2));
MI_max=max([max(all_hit_z_MIs) max(all_miss_z_MIs)]);

plot(all_hit_z_MIs, R1s_hit,'o','MarkerFaceColor',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'Color',[86/255 180/255 233/255],'MarkerSize',8)


plot(all_miss_z_MIs, R1s_miss,'o','MarkerFaceColor',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'Color',[0/255 158/255 115/255],'MarkerSize',8)


p = polyfit([all_miss_z_MIs all_hit_z_MIs], [R1s_miss R1s_hit], 1);
R1s_all_fit = polyval(p, [-5 130]);
plot([-5 130],R1s_all_fit,'-','Color',[0/255 0/255 0/255],'LineWidth',2)
[R1,P]=corrcoef(all_miss_z_MIs,R1s_miss);
fprintf(1, 'Correlation coefficient for R1 vs zMIs %d, p value %d\n\n',R1(1,2),P(1,2));

xlabel('zMIs')
ylabel('R1')
title('R1 for decoding vs. multi ROI zMI betweeen op and binary dFF')

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

edges=[0:0.05:1];
rand_offset=0.5;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

MI_labels{1}='All trials';
MI_labels{2}='Hits';
MI_labels{3}='Misses';
MI_labels{4}='Shuffled';

 

these_MIs=all_hit_z_MIs;

bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])

%Violin plot
[mean_out, CIout, all_MIs.type(ii_type).violin_x]=drgViolinPoint(these_MIs...
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

bar(bar_offset,mean(these_MIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

%Violin plot
[mean_out, CIout, all_MIs.type(ii_type).violin_x]=drgViolinPoint(these_MIs...
    ,edges,bar_offset,rand_offset,'k','k',4);
bar_offset=bar_offset+1;

glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_MIs))=these_MIs;
glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_MIs))=1*ones(1,length(these_MIs));
glm_r1_ii=glm_r1_ii+length(these_MIs);

id_r1_ii=id_r1_ii+1;
input_r1_data(id_r1_ii).data=these_MIs;
input_r1_data(id_r1_ii).description=['Miss'];
       

%Plot lines between points
for ii=1:length(all_miss_z_MIs)
    plot([0 1],[all_hit_z_MIs(ii) all_miss_z_MIs(ii)],'-','Color',[0.7 0.7 0.7])
    plot(zeros(1,length(all_hit_z_MIs)),all_hit_z_MIs,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4)
    plot(ones(1,length(all_miss_z_MIs)),all_miss_z_MIs,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4)
end


xticks([0 1])
xticklabels({'hit','miss'})


title(['multi ROI MI betweeen op and binary dFF'])
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

fprintf(1, ['\nMean miss divided by mean hit ' num2str(mean(all_miss_z_MIs)/mean(all_hit_z_MIs)) '\n'])

pfft=1;

