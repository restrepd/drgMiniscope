%drgMini_analyze_batch_DecodeCorrXYConc
close all
clear all

is_sphgpu=0;

switch is_sphgpu
    case 0

        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';

        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

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

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
edges=[-180:10:180];
histogram(all_meanAngles,edges)
title('Angle for final approach for all files')
xlabel('Angle (degrees)')

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
edges=[0:0.05:1];
histogram(fraction_other_angle,edges)
title('Fraction of horizontal approaches for all files')
xlabel('Fraction')


%Do the R1 conc bar comparing 1 cm vs 2 cm with 0 or 16 bins before
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


these_ii_run=[1 5]; %0 or 16 bins
these_groups=[1 5]; %2 cm and 1 cm all tirals

for ii_run=these_ii_run
    for grNo=these_groups
        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            arena_file=handles_conc.arena_file{fileNo};
            if (handles_conc.group(fileNo)==grNo)&(fraction_other_angle(fileNo)<thr_froa)
                %load the ouptut file
                load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
                these_R1s=[these_R1s handles_out.R1.all_trials];
            end

        end
        
        %plot bar
        switch grNo
            case 1
                bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 5
                bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        end

        %Violin plot
        [mean_out, CIout,violin_x]=drgViolinPoint(these_R1s...
            ,edges,bar_offset,rand_offset,'k','k',4);
        bar_offset=bar_offset+1;

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
        glm_r1.group(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=grNo*ones(1,length(these_R1s));
        glm_r1.run(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_run*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[group_label{grNo} ' ' run_label{ii_run}];
    end
    bar_offset=bar_offset+1;
end

text(1.5,0.38,'2 cm','Color',[230/255 159/255 0/255],'FontWeight','bold')
text(1.5,0.35,'1 cm','Color',[86/255 180/255 233/255],'FontWeight','bold')
% %Plot lines between points
% for ii=1:length(these_R1s)
%     these_R1s=[];
%     for ii_R1type=1:length(R1type)
%         these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)]
%     end
%     plot([0:length(R1type)-1],these_R1s,'-','Color',[0.7 0.7 0.7])
% end


xticks([0.5 3.5])
xticklabels({'0 bins','16 bins'})


title(['R1 for predcition of odor concentration'])
ylabel('R1')
ylim([-0.2 0.5])
xlim([-1 5])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig 3 R1 cm from floor and bins before\n'])
fprintf(fileID, ['\nglm for Fig 3 R1 cm from floor and bins before\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.group',glm_r1.run',...
    'VariableNames',{'R1','group','run'});
mdl = fitglm(tbl,'R1~group+run'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for R1 cm from floor and bins before\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for R1 cm from floor and bins before\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

%Let's look at dependence of correlation on bins before for all 1 and 2 cm files
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

these_groups=[1 5];


edges=[0:0.05:1];
rand_offset=1;

%R1 for conc all trials
meanR1_conc=[];
bins_before=handles_conc.bins_before;

for ii_run=1:length(handles_conc.bins_before)
    these_R1s=[];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
            arena_file=handles_conc.arena_file{fileNo};
            %load the ouptut file
            load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            these_R1s=[these_R1s handles_out.R1.all_trials];
        end
    end

    [mean_out, CIout]=drgViolinPoint(these_R1s...
        ,edges,bins_before(ii_run),rand_offset,[0.9 0.6 0],[0.9 0.6 0],4);
    meanR1_conc=[meanR1_conc mean(these_R1s)];
end
 
plot(bins_before,meanR1_conc,'o-','MarkerSize',10,'MarkerEdgeColor',[0.9 0.6 0],'MarkerFaceColor',[0.9 0.6 0],'Color',[0.9 0.6 0])

%R1 for XY all trials
meanR1_XY=[];
bins_before=handles_XY.bins_before;

for ii_run=1:length(handles_XY.bins_before)
    these_R1s=[];
    for fileNo=1:length(handles_XY.arena_file)
        if (sum(handles_XY.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
            arena_file=handles_conc.arena_file{fileNo};
            %load the ouptut file
            try
                load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
                these_R1s=[these_R1s handles_out.R1.within_XY];
            catch
            end
        end
    end

    [mean_out, CIout]=drgViolinPoint(these_R1s...
        ,edges,bins_before(ii_run),rand_offset,[0.35 0.7 0.9],[0.35 0.7 0.9],4);
    meanR1_XY=[meanR1_XY mean(these_R1s)];
end

plot(bins_before,meanR1_XY,'o-','MarkerSize',10,'MarkerEdgeColor',[0.35 0.7 0.9],'MarkerFaceColor',[0.35 0.7 0.9],'Color',[0.35 0.7 0.9])
title('Within trial R1')
xlabel('Bins before')
ylabel('R1')

%For 0 bins before plot R1 for conc vs XY
ii_run=1; %0 bins before 

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

these_groups=[1 5];
these_new_groups.gr(1).groups=[1 5]; %2 cm and 1 cm
these_new_groups.gr(2).groups=[2 3]; %One lane odor
edges=[0:0.05:1];
rand_offset=1;
R1_per_mouse=[];

for grNo=1:2
    for ii_mouse=unique(handles_conc.mouse)
        R1_per_mouse.group(grNo).mouse(ii_mouse).R1_conc=[];
        R1_per_mouse.group(grNo).mouse(ii_mouse).pc=[];
        R1_per_mouse.group(grNo).mouse(ii_mouse).R1_XY=[];
    end
end

for grNo=1:2
    %R1 for conc all trials

    R1_conc=[];
    percent_correct=[];

    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
            arena_file=handles_conc.arena_file{fileNo};
            %load the ouptut file
            load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            R1_conc=[R1_conc handles_out.R1.all_trials];
            R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).R1_conc=...
                [R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).R1_conc handles_out.R1.all_trials];
            this_pc=100*(sum(handles_out.trials.hit1)+sum(handles_out.trials.hit4))/length(handles_out.trials.hit1);
            R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).pc=...
                [R1_per_mouse.group(grNo).mouse(handles_conc.mouse(fileNo)).pc this_pc];
            percent_correct=[percent_correct this_pc];
        end
    end
    these_new_groups.gr(grNo).R1_conc=R1_conc;
    these_new_groups.gr(grNo).percent_correct=percent_correct;

    %R1 for XY all trials
    R1_XY=[];

    for fileNo=1:length(handles_XY.arena_file)
        if (sum(handles_XY.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
            arena_file=handles_conc.arena_file{fileNo};
            %load the ouptut file
            load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
            R1_XY=[R1_XY handles_out.R1.within_XY];
            R1_per_mouse.group(grNo).mouse(handles_XY.mouse(fileNo)).R1_XY=...
                [R1_per_mouse.group(grNo).mouse(handles_XY.mouse(fileNo)).R1_XY handles_out.R1.within_XY];
        end
    end

    these_new_groups.gr(grNo).R1_XY=R1_XY;


    switch grNo
        case 1
            plot(R1_conc,R1_XY,'ok','MarkerFaceColor','k')
        case 2
            plot(R1_conc,R1_XY,'ok','MarkerFaceColor',[0.7 0.7 0.7])
    end

end

text(0.5,0.5,'Both spouts','Color','k')
text(0.5,0.4,'One spout','Color',[0.7 0.7 0.7])
title('Within trial R1')
xlabel('R1 conc')
ylabel('R1 XY')

%R1 0 bins before per mouse

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

marker_per_group{1}=['' 'o' ''];
marker_per_group{2}=['' 's' ''];

color_okabe_ito{1}='[230/255 159/255 0/255]';
color_okabe_ito{2}='[86/255 180/255 233/255]';
color_okabe_ito{3}='[0/255 158/255 115/255]';
color_okabe_ito{4}='[240/255 228/255 66/255]';
color_okabe_ito{5}='[0/255 114/255 178/255]';

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])


for grNo=1:2

    for ii_mouse=unique(handles_conc.mouse)

        if ~isempty(R1_per_mouse.group(grNo).mouse(ii_mouse).R1_conc)
            these_R1conc=R1_per_mouse.group(grNo).mouse(ii_mouse).R1_conc;
            these_R1XY=R1_per_mouse.group(grNo).mouse(ii_mouse).R1_XY;
            this_mean_R1conc=mean(these_R1conc);
            this_mean_R1XY=mean(these_R1XY);
            eval(['plot(this_mean_R1conc, this_mean_R1XY,''' marker_per_group{grNo} ''',''MarkerSize'',12,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_mouse} ')'])
            for ii_point=1:length(these_R1conc)
                eval(['plot(these_R1conc(ii_point), these_R1XY(ii_point),''' marker_per_group{grNo} ''',''MarkerSize'',6,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_mouse} ')'])
            end
            
        end
    end

end

text(0.5,0.5,'Both spouts','Color','k')
text(0.5,0.4,'One spout','Color',[0.7 0.7 0.7])
title('Within trial R1 per mouse')
xlabel('R1 conc')
ylabel('R1 XY')


%Plot R1 vs percent correct behavior
group_labels{1}='odor in both spouts';
group_labels{2}='odor in one spout';
for grNo=1:2
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    R1_conc=these_new_groups.gr(grNo).R1_conc;
    R1_XY=these_new_groups.gr(grNo).R1_XY;
    percent_correct=these_new_groups.gr(grNo).percent_correct;

    plot(percent_correct,R1_conc,'o','MarkerEdgeColor','none','MarkerFaceColor',[230/255 159/255 0/255])
    plot(percent_correct,R1_XY,'o','MarkerEdgeColor','none','MarkerFaceColor',[86/255 180/255 233/255])

    text(60,0.7,'Odor','Color',[230/255 159/255 0/255])
    text(60,0.65,'XY','Color',[86/255 180/255 233/255])

    title(['R1 vs. percent correct behavior ' group_labels{grNo}])
    ylabel('R1')
    xlabel('Percent correct')
end

%Compare groups with fittree, weber frechner and only current data
%point
ii_run=1;
groups=unique(handles_conc.group);

group_order=[1 4 5 2 3];

R1type=[];
R1type{1}='all_trials';
R1type{2}='between_trials';
R1type{3}='all_hits';
R1type{4}='all_miss';
R1type{5}='all_trials_sh';

R1type_label=[];
R1type_label{1}='within';
R1type_label{2}='between';
R1type_label{3}='hits';
R1type_label{4}='miss';
R1type_label{5}='shuffled';

%Do the R1 bar graph
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

bar_offset=0;

edges=[0:0.05:1];
rand_offset=0.5;

%Plot the different R1s
for ii_group=1:3
    for ii_R1type=1:length(R1type)

        this_group=group_order(ii_group);

        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            arena_file=handles_conc.arena_file{fileNo};
            if (handles_conc.group(fileNo)==this_group)&(fraction_other_angle(fileNo)<thr_froa)
                %load the ouptut file
                load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
                eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
                % these_R1s=[these_R1s handles_out.R1.all_trials];
            end

        end

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
        [mean_out, CIout]=drgViolinPoint(these_R1s...
            ,edges,bar_offset,rand_offset,'k','k',4);
        bar_offset=bar_offset+1;
    end
    bar_offset=bar_offset+1;
end

x_pos=15.3;
text(x_pos,0.57,'within','Color',[230/255 159/255 0/255])
text(x_pos,0.53,'between','Color',[86/255 180/255 233/255])
text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([2 8 14])
xticklabels({'2 cm','No odor','1 cm','Pre','Odor','Reinf'})


title(['R1, correlation between decoded odor concentration and odor plume mean'])
ylabel('R1')
ylim([-0.2 0.65])
xlim([-1 17])


%Do the R1 bar graph for odor vs one lane vs shuffled


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

ii_run=1; %0 bins before
these_new_groups.gr(1).groups=[1 5]; %2 cm and 1 cm
these_new_groups.gr(2).groups=[2 3]; %One lane odor

R1type_label=[];
R1type_label{1}='all_trials';
R1type_label{2}='all_trials_sh';

R1type{1}='all_trials';
R1type{2}='all_trials_sh';

group_labels{1}='Both spouts';
group_labels{2}='One spout';

%Plot the different R1s
all_R1s=[];

for grNo=1:2
    for ii_R1type=1:length(R1type_label)

        these_R1s=[];
        for fileNo=1:length(handles_conc.arena_file)
            arena_file=handles_conc.arena_file{fileNo};
            if (sum(handles_conc.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
                %load the ouptut file
                load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
                eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
            end
        end
        all_R1s.group(grNo).R1type(ii_R1type).R1=these_R1s;
        %plot bar
        switch ii_R1type
            case 1
                bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
            case 2
                bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        end

        %Violin plot
        [mean_out, CIout, all_R1s.group(grNo).ii_R1type(ii_R1type).violin_x]=drgViolinPoint(these_R1s...
            ,edges,bar_offset,rand_offset,'k','k',4);
        bar_offset=bar_offset+1;

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
        glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_R1type*ones(1,length(these_R1s));
        glm_r1.group(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=grNo*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type} ' ' group_labels{grNo}];

    end
    bar_offset=bar_offset+1;
end


%Plot lines between points
ii_R1type=1;
for grNo=1:2
    for ii=1:length(all_R1s.group(grNo).R1type(ii_R1type).R1)
        these_R1s=[];
        for ii_R1type=1:2
            these_R1s=[these_R1s all_R1s.group(grNo).R1type(ii_R1type).R1(ii)];   
        end
        plot([all_R1s.group(grNo).ii_R1type(1).violin_x(ii) all_R1s.group(grNo).ii_R1type(2).violin_x(ii)],these_R1s,'-','Color',[0.7 0.7 0.7])
    end
end


% x_pos=3;
text(0,0.9,'within','Color',[230/255 159/255 0/255])
text(0,0.8,'shuffled','Color',[86/255 180/255 233/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 3])
xticklabels({'both spouts','one spout'})


title(['R1 for odor in either both spouts or one spout'])
ylabel('R1')
ylim([-0.2 1])
xlim([-1 5])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig 9 odor in one or both spouts\n'])
fprintf(fileID, ['\nglm for Fig 9 odor in one or both spouts\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.trial_type',glm_r1.group',...
    'VariableNames',{'R1','trial_type','group'});
mdl = fitglm(tbl,'R1~trial_type+group'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for odor in one or both spouts\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for odor in one or both spouts\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

%Do the R1 bar graph for 1 and 2 cm all trials vs hit, miss
ii_run=1;
R1type=[];
R1type{1}='all_trials';
R1type{2}='all_hits';
R1type{3}='all_miss';
R1type{4}='all_trials_sh';

R1type_label=[];
R1type_label{1}='within';
R1type_label{2}='hits';
R1type_label{3}='miss';
R1type_label{4}='shuffled';

iiR1type_reorder=[2 3 1 4];

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
    for fileNo=1:length(handles_conc.arena_file)
        arena_file=handles_conc.arena_file{fileNo};
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
            %load the ouptut file
            load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
            switch ii_R1type
                case 1
                    R1_all_trials_orig=[R1_all_trials_orig handles_out.R1.all_trials];
                case 2
                    R1_hits_orig=[R1_hits_orig handles_out.R1.all_hits];

                case 3
                    R1_miss_orig=[R1_miss_orig handles_out.R1.all_miss];
                    if handles_out.R1.all_miss>0.5
                        pffft=1; 
                    end
                case 4
                    R1_shuffle_orig=[R1_shuffle_orig handles_out.R1.all_trials_sh];
            end
        end

    end
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


title(['R1 for prediction of odor hit,miss, shuffled'])
ylabel('R1')
ylim([-0.2 0.65])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig 10 trial type\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig 10 trial type\n']);
fprintf(fileID, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
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

%Let's get the angle for final approach
%Let's look at dependence of correlation on bins before for all 2 cm files
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

these_groups=[1 5];


%Calculate fraction of horizontal angle approach
meanAngles=[];

these_meanAngles=[];
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        for trNo=1:length(handles_out.angles.trial)
            meanAngles=[meanAngles handles_out.angles.trial(trNo).mean_end_angle];
        end
    end
end


edges=[-180:10:180];
histogram(meanAngles,edges)
title('Angle for final approach for 1 cm and 2 cm odor in both spouts')
xlabel('Angle (degrees)')

% low_angle=-130;
% high_angle=-50;

%Now re-calculate the hits
these_groups=[1 5];
ii_run=1;
ii_for_corr=0;
R1_all_trials=[];
R1_hits=[];
R1_hits90=[];
R1_hitso=[];
R1_miss=[];
R1_between=[];


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
        op_all_trials=[];
        op_decod_all_trials=[];
        % op_all_trials_sh=[];
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
        op_all_miss_sh=[];
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

        ii_for_corr=ii_for_corr+1;
        if ~isempty(op_all_trials)
            R1=corrcoef(op_all_trials,op_decod_all_trials);
            R1_all_trials(ii_for_corr)=R1(1,2);
        else
            R1_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits)
            R1=corrcoef(op_all_hits,op_decod_all_hits);
            R1_hits(ii_for_corr)=R1(1,2);
        else
            R1_hits(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits90)
            R1=corrcoef(op_all_hits90,op_decod_all_hits90);
            R1_hits90(ii_for_corr)=R1(1,2);
        else
            R1_hits90(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hitso)
            R1=corrcoef(op_all_hitso,op_decod_all_hitso);
            R1_hitso(ii_for_corr)=R1(1,2);
        else
            R1_hitso(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_miss)
            R1=corrcoef(op_all_miss,op_decod_all_miss);
            R1_miss(ii_for_corr)=R1(1,2);
        else
            R1_miss(ii_for_corr)=NaN;
        end

        if ~isempty(op_between_trials)
            R1=corrcoef(op_between_trials,op_decod_between_trials);
            R1_between(ii_for_corr)=R1(1,2);
        else
            R1_between(ii_for_corr)=NaN;
        end

    end
end

%Let's plot a bar graph for hits with different approach angles
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


%Plot the different R1s
these_groups=[1 5]; %1,2 cm all tirals

R1type=[];
R1type{1}='R1_hits90';
R1type{2}='R1_hitso';
R1type{3}='R1_miss';


R1type_label=[];
R1type_label{1}='hit90';
R1type_label{2}='hito';
R1type_label{3}='miss';

allR1s=[];
for ii_R1type=1:length(R1type)

    these_R1s=[];
    eval(['these_R1s=' R1type{ii_R1type} ';'])
   
    for ii=1:length(these_R1s)
        all_R1s.R1type(ii_R1type).R1(ii)=these_R1s(ii);
    end

    these_R1s=these_R1s(~isnan(these_R1s));

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
    [mean_out, CIout, violin_x]=drgViolinPoint(these_R1s...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;
 
    ii_next=0;
    for ii=1:length(all_R1s.R1type(ii_R1type).R1)
        if ~isnan(all_R1s.R1type(ii_R1type).R1(ii))
            ii_next=ii_next+1;
            all_R1s.R1type(ii_R1type).violin_x(ii)=violin_x(ii_next);
        end
    end

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_R1type*ones(1,length(these_R1s));
    glm_r1_ii=glm_r1_ii+length(these_R1s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R1s;
    input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

end

%Plot lines between points
for ii=1:length(all_R1s.R1type(ii_R1type).R1)
    these_R1s=[];
    these_violin_x=[];
    no_nans=1;
    for ii_R1type=1:length(R1type)
        if isnan(all_R1s.R1type(ii_R1type).R1(ii))
            no_nans=0;
        end
    end

    if no_nans==1
        for ii_R1type=1:length(R1type)
            these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)];
            these_violin_x=[these_violin_x all_R1s.R1type(ii_R1type).violin_x(ii)];
        end

        plot(these_violin_x,these_R1s,'-','Color',[0.7 0.7 0.7])
    end
end

% x_pos=3;
% text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 2 3])
xticklabels({'hit90','hito','miss'})


title(['R1 for predcition of odor concentration'])
ylabel('R1')
ylim([-0.2 0.65])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig 12 trial type\n'])
fprintf(fileID, ['\nglm for Fig 12 trial type\n']);
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
 
%Now re-calculate R1 within before and after last turn
ii_run=1;
these_groups=[1 5];
ii_for_corr=0;
R1_all_trials_bt=[];
R1_all_trials_at=[];
R1_hits_bt=[];
R1_hits_at=[];
R1_miss_bt=[];
R1_miss_at=[];


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
        op_all_trials_bt=[];
        op_decod_all_trials_bt=[];

        op_all_trials_at=[];
        op_decod_all_trials_at=[];

        op_all_hits_bt=[];
        op_decod_all_hits_bt=[];

        op_all_hits_at=[];
        op_decod_all_hits_at=[];

        op_all_miss_bt=[];
        op_decod_all_miss_bt=[];

        op_all_miss_at=[];
        op_decod_all_miss_at=[];

        op_all_trials_bt_sh=[];
        op_decod_all_trials_bt_sh=[];

        op_all_trials_at_sh=[];
        op_decod_all_trials_at_sh=[];

        for ii_sh=1:handles_conc.n_shuffle
            op_all_trials_bt_sh.ii_sh(ii_sh).allt=[];
            op_decod_all_trials_bt_sh.ii_sh(ii_sh).allt=[];

            op_all_trials_at_sh.ii_sh(ii_sh).allt=[];
            op_decod_all_trials_at_sh.ii_sh(ii_sh).allt=[];
        end


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

        op_predicted_sh=handles_out.op_predicted_sh;
    

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

        op_predicted_conv_sh=zeros(size(op_predicted_sh,1),size(op_predicted_sh,2));
        for ii_sh=1:size(op_predicted_sh,2)
            this_op_predicted_sh=zeros(size(op_predicted_sh,1),1);
            this_op_predicted_sh(:,1)=op_predicted_sh(:,ii_sh);
            this_op_predicted_conv_sh=conv(this_op_predicted_sh,conv_win_gauss,'same');
            this_op_predicted_conv_sh(this_op_predicted_conv_sh<minop)=minop;
            this_op_predicted_conv_sh(this_op_predicted_conv_sh>maxop)=maxop;
            op_predicted_conv_sh(:,ii_sh)=this_op_predicted_conv_sh;
        end


        last_op_predictedend=1;
        ii_start=0;

        for trNo=1:trials.odor_trNo



            %Find the last turn
            ii_turn=find(handles_out_angle.angles.trial(trNo).delta_x>100,1,'last');
            if ~isempty(ii_turn)
                op_predictedstart_bt=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                op_predictedend_bt=op_predictedstart_bt+handles_out_angle.angles.trial(trNo).ii_turns(ii_turn)-1;
                if op_predictedend_bt>length(odor_plume_template)
                    op_predictedend_bt=length(odor_plume_template);
                end

                op_predictedstart_at=op_predictedstart_bt+handles_out_angle.angles.trial(trNo).ii_turns(ii_turn);
                op_predictedend_at=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if op_predictedend_at>length(odor_plume_template)
                    op_predictedend_at=length(odor_plume_template);
                end

                op_all_trials_bt=[op_all_trials_bt; odor_plume_template(op_predictedstart_bt:op_predictedend_bt)'];
                op_decod_all_trials_bt=[op_decod_all_trials_bt; op_predicted_conv(op_predictedstart_bt:op_predictedend_bt)];

                op_all_trials_at=[op_all_trials_at; odor_plume_template(op_predictedstart_at:op_predictedend_at)'];
                op_decod_all_trials_at=[op_decod_all_trials_at; op_predicted_conv(op_predictedstart_at:op_predictedend_at)];

                for ii_sh=1:size(op_predicted_sh,2)
                    op_all_trials_bt_sh.ii_sh(ii_sh).allt=[op_all_trials_bt_sh.ii_sh(ii_sh).allt; odor_plume_template(op_predictedstart_bt:op_predictedend_bt)'];
                    op_decod_all_trials_bt_sh.ii_sh(ii_sh).allt=[op_decod_all_trials_bt_sh.ii_sh(ii_sh).allt; op_predicted_conv_sh(op_predictedstart_bt:op_predictedend_bt,ii_sh)];

                    op_all_trials_at_sh.ii_sh(ii_sh).allt=[op_all_trials_at_sh.ii_sh(ii_sh).allt; odor_plume_template(op_predictedstart_at:op_predictedend_at)'];
                    op_decod_all_trials_at_sh.ii_sh(ii_sh).allt=[op_decod_all_trials_at_sh.ii_sh(ii_sh).allt; op_predicted_conv_sh(op_predictedstart_at:op_predictedend_at,ii_sh)];
                end


                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
                        op_all_hits_bt=[op_all_hits_bt; odor_plume_template(op_predictedstart_bt:op_predictedend_bt)'];
                        op_decod_all_hits_bt=[op_decod_all_hits_bt; op_predicted_conv(op_predictedstart_bt:op_predictedend_bt)];

                        op_all_hits_at=[op_all_hits_at; odor_plume_template(op_predictedstart_at:op_predictedend_at)'];
                        op_decod_all_hits_at=[op_decod_all_hits_at; op_predicted_conv(op_predictedstart_at:op_predictedend_at)];
                
                    case 2
                        %Lane 1 miss orange
                        op_all_miss_bt=[op_all_miss_bt; odor_plume_template(op_predictedstart_bt:op_predictedend_bt)'];
                        op_decod_all_miss_bt=[op_decod_all_miss_bt; op_predicted_conv(op_predictedstart_bt:op_predictedend_bt)];

                        op_all_miss_at=[op_all_miss_at; odor_plume_template(op_predictedstart_at:op_predictedend_at)'];
                        op_decod_all_miss_at=[op_decod_all_miss_at; op_predicted_conv(op_predictedstart_at:op_predictedend_at)];

                    case 3
                        %Lane 4 hit blue
                         op_all_hits_bt=[op_all_hits_bt; odor_plume_template(op_predictedstart_bt:op_predictedend_bt)'];
                        op_decod_all_hits_bt=[op_decod_all_hits_bt; op_predicted_conv(op_predictedstart_bt:op_predictedend_bt)];

                        op_all_hits_at=[op_all_hits_at; odor_plume_template(op_predictedstart_at:op_predictedend_at)'];
                        op_decod_all_hits_at=[op_decod_all_hits_at; op_predicted_conv(op_predictedstart_at:op_predictedend_at)];
                    case 4
                        %Lane 4 miss sky blue
                        op_all_miss_bt=[op_all_miss_bt; odor_plume_template(op_predictedstart_bt:op_predictedend_bt)'];
                        op_decod_all_miss_bt=[op_decod_all_miss_bt; op_predicted_conv(op_predictedstart_bt:op_predictedend_bt)];

                        op_all_miss_at=[op_all_miss_at; odor_plume_template(op_predictedstart_at:op_predictedend_at)'];
                        op_decod_all_miss_at=[op_decod_all_miss_at; op_predicted_conv(op_predictedstart_at:op_predictedend_at)];
                end
            end

            % op_between_trials=[op_between_trials; odor_plume_template(last_op_predictedend:op_predictedstart)'];
            % op_decod_between_trials=[op_decod_between_trials; op_predicted_conv(last_op_predictedend:op_predictedstart)];  
            % last_op_predictedend=op_predictedend;

            % ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;


            % 
            % %Okabe_Ito colors
            % switch trials.odor_trial_type(trNo)
            %     case 1
            %         %Lane 1 hits vermillion
            %         op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
            %         op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];
            % 
            %         %Hit 90 degrees
            %         if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
            %             op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
            %             op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
            %         end
            % 
            %         %Horizontal hit
            %         if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
            %             op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
            %             op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
            %         end
            %     case 2
            %         %Lane 1 miss orange
            %         op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
            %         op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
            % 
            %     case 3
            %         %Lane 4 hit blue
            %         op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
            %         op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];
            % 
            %            %Hit 90 degrees
            %         if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
            %             op_all_hits90=[op_all_hits90; odor_plume_template(op_predictedstart:op_predictedend)'];
            %             op_decod_all_hits90=[op_decod_all_hits90; op_predicted_conv(op_predictedstart:op_predictedend)];
            %         end
            % 
            %         %Horizontal hit
            %         if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
            %             op_all_hitso=[op_all_hitso; odor_plume_template(op_predictedstart:op_predictedend)'];
            %             op_decod_all_hitso=[op_decod_all_hitso; op_predicted_conv(op_predictedstart:op_predictedend)];
            %         end
            %     case 4
            %         %Lane 4 miss sky blue
            %         op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
            %         op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
            % end
            % ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
        end

        ii_for_corr=ii_for_corr+1;
        if ~isempty(op_all_trials_bt)
            R1=corrcoef(op_all_trials_bt,op_decod_all_trials_bt);
            R1_all_trials_bt(ii_for_corr)=R1(1,2);
        else
            R1_all_trials_bt(ii_for_corr)=NaN;
        end
 
        if ~isempty(op_all_trials_at)
            R1=corrcoef(op_all_trials_at,op_decod_all_trials_at);
            R1_all_trials_at(ii_for_corr)=R1(1,2);
        else
            R1_all_trials_at(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits_bt)
            R1=corrcoef(op_all_hits_bt,op_decod_all_hits_bt);
            R1_hits_bt(ii_for_corr)=R1(1,2);
        else
            R1_hits_bt(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits_at)
            R1=corrcoef(op_all_hits_at,op_decod_all_hits_at);
            R1_hits_at(ii_for_corr)=R1(1,2);
        else
            R1_hits_at(ii_for_corr)=NaN;
        end

         if ~isempty(op_all_miss_bt)
            R1=corrcoef(op_all_miss_bt,op_decod_all_miss_bt);
            R1_miss_bt(ii_for_corr)=R1(1,2);
        else
            R1_miss_bt(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_miss_at)
            R1=corrcoef(op_all_miss_at,op_decod_all_miss_at);
            R1_miss_at(ii_for_corr)=R1(1,2);
        else
            R1_miss_at(ii_for_corr)=NaN;
        end

        these_R1s=[];
        for ii_sh=1:handles_conc.n_shuffle
            if ~isempty(op_all_trials_bt_sh.ii_sh(ii_sh).allt)
                R1=corrcoef(op_all_trials_bt_sh.ii_sh(ii_sh).allt,op_decod_all_trials_bt_sh.ii_sh(ii_sh).allt);
                these_R1s=[these_R1s R1(1,2)];
            else
                these_R1s=[these_R1s NaN];
            end
        end
        R1_all_trials_bt_sh(ii_for_corr)=mean(these_R1s);

         these_R1s=[];
        for ii_sh=1:handles_conc.n_shuffle
            if ~isempty(op_all_trials_at_sh.ii_sh(ii_sh).allt)
                R1=corrcoef(op_all_trials_at_sh.ii_sh(ii_sh).allt,op_decod_all_trials_at_sh.ii_sh(ii_sh).allt);
                these_R1s=[these_R1s R1(1,2)];
            else
                these_R1s=[these_R1s NaN];
            end
        end
        R1_all_trials_at_sh(ii_for_corr)=mean(these_R1s);
 
        pffft=1;

    end
end


%Let's plot within vs, within before and after last turn
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


%Plot the different R1s
these_groups=[1 5]; %1,2 cm all tirals

R1type=[];
R1type{1}='R1_all_trials';
R1type{2}='R1_all_trials_bt';
R1type{3}='R1_all_trials_at';


R1type_label=[];
R1type_label{1}='within';
R1type_label{2}='before';
R1type_label{3}='after';

allR1s=[];
for ii_R1type=1:length(R1type)

    these_R1s=[];
    eval(['these_R1s=' R1type{ii_R1type} ';'])
   
    for ii=1:length(these_R1s)
        all_R1s.R1type(ii_R1type).R1(ii)=these_R1s(ii);
    end

    these_R1s=these_R1s(~isnan(these_R1s));

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
    [mean_out, CIout, violin_x]=drgViolinPoint(these_R1s...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;
 
    ii_next=0;
    for ii=1:length(all_R1s.R1type(ii_R1type).R1)
        if ~isnan(all_R1s.R1type(ii_R1type).R1(ii))
            ii_next=ii_next+1;
            all_R1s.R1type(ii_R1type).violin_x(ii)=violin_x(ii_next);
        end
    end

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_R1type*ones(1,length(these_R1s));
    glm_r1_ii=glm_r1_ii+length(these_R1s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R1s;
    input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

end

%Plot lines between points
for ii=1:length(all_R1s.R1type(ii_R1type).R1)
    these_R1s=[];
    these_violin_x=[];
    no_nans=1;
    for ii_R1type=1:length(R1type)
        if isnan(all_R1s.R1type(ii_R1type).R1(ii))
            no_nans=0;
        end
    end

    if no_nans==1
        for ii_R1type=1:length(R1type)
            these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)];
            these_violin_x=[these_violin_x all_R1s.R1type(ii_R1type).violin_x(ii)];
        end

        plot(these_violin_x,these_R1s,'-','Color',[0.7 0.7 0.7])
    end
end

% x_pos=3;
% text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 2 3])
xticklabels({'within','before','after'})


title(['R1 for prediction of odor concentration before and after last turn'])
ylabel('R1')
ylim([-0.2 0.65])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig 14 R1 before after last turn\n'])
fprintf(fileID, ['\nglm for Fig 14 R1 before after last turn\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.trial_type',...
    'VariableNames',{'R1','turn'});
mdl = fitglm(tbl,'R1~turn'...
    ,'CategoricalVars',[2])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig 14 R1 before after last turn\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig 14 R1 before after last turn\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

 
%Do the R1 bar graph for 1 and 2 cm all trials vs hit, miss, after last
%turn
ii_run=1;
R1type=[];
R1type{1}='all_trials';
R1type{2}='all_hits';
R1type{3}='all_miss';
R1type{4}='all_trials_sh';

R1type_label=[];
R1type_label{1}='within';
R1type_label{2}='hits';
R1type_label{3}='miss';
R1type_label{4}='shuffled';

iiR1type_reorder=[2 3 1 4];

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
    for fileNo=1:length(handles_conc.arena_file)
        arena_file=handles_conc.arena_file{fileNo};
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)
            %load the ouptut file
            load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
            switch ii_R1type
                case 1
                    R1_all_trials_orig=[R1_all_trials_orig handles_out.R1.all_trials];
                case 2
                    R1_hits_orig=[R1_hits_orig handles_out.R1.all_hits];

                case 3
                    R1_miss_orig=[R1_miss_orig handles_out.R1.all_miss];
                    if handles_out.R1.all_miss>0.5
                        pffft=1; 
                    end
                case 4
                    R1_shuffle_orig=[R1_shuffle_orig handles_out.R1.all_trials_sh];
            end
        end

    end
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


title(['R1 for prediction of odor hit,miss, shuffled'])
ylabel('R1')
ylim([-0.2 0.65])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig 10 trial type\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig 10 trial type\n']);
fprintf(fileID, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
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
pffft=1;
%
% 
% %Let's look at imp
% for this_group=1:5
%     ii_run=1;
%     all_mean_conc_imps=[];
%     above_thr_conc_imps=[];
%     all_mean_x_imps=[];
%     above_thr_x_imps=[];
%     all_mean_y_imps=[];
%     above_thr_y_imps=[];
%     thr_x_imps=zeros(1,length(handles_conc.arena_file));
%     thr_y_imps=zeros(1,length(handles_conc.arena_file));
%     thr_conc_imps=zeros(1,length(handles_conc.arena_file));
%     file_numbers=[];
% 
% 
%     figureNo=figureNo+1;
%     these_R1s=[];
%     for fileNo=1:length(handles_conc.arena_file)
%         if handles_conc.group(fileNo)==this_group
%             arena_file=handles_conc.arena_file{fileNo};
% 
%             %Get predictor importance for conc decoding
%             load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
%             all_conc_imps=[];
%             for trNo=1:length(handles_out.imp.trial)
%                 these_imps=handles_out.imp.trial(trNo).imp;
%                 all_conc_imps=[all_conc_imps these_imps'];
%             end
%             sum_all_conc_imps=sum(all_conc_imps);
%             mat_sum_all_conc_imps=repmat(sum_all_conc_imps,size(all_conc_imps,1),1);
%             all_conc_imps=all_conc_imps./mat_sum_all_conc_imps; %Normalize to all added imps=1
%             mean_conc_imps=mean(all_conc_imps,2);
%             thr_conc_imps(fileNo)=prctile(mean_conc_imps,95);
%             all_mean_conc_imps=[all_mean_conc_imps; mean_conc_imps];
%             above_thr_conc_imps=[above_thr_conc_imps; mean_conc_imps>=thr_conc_imps(fileNo)];
%             file_numbers=[file_numbers; fileNo*ones(length(mean_conc_imps),1)];
% 
%             %Get predictor importance for x and y decoding
%             load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
%             all_x_imps=[];
%             all_y_imps=[];
%             for trNo=1:length(handles_out.imp.trial)
%                 these_x_imps=handles_out.imp.trial(trNo).impx;
%                 all_x_imps=[all_x_imps these_x_imps'];
%                 these_y_imps=handles_out.imp.trial(trNo).impy;
%                 all_y_imps=[all_y_imps these_y_imps'];
%             end
%             sum_all_x_imps=sum(all_x_imps);
%             mat_sum_all_x_imps=repmat(sum_all_x_imps,size(all_x_imps,1),1);
%             all_x_imps=all_x_imps./mat_sum_all_x_imps; %Normalize to all added imps=1
%             mean_x_imps=mean(all_x_imps,2);
%             thr_x_imps(fileNo)=prctile(mean_x_imps,95);
%             all_mean_x_imps=[all_mean_x_imps; mean_x_imps];
%             above_thr_x_imps=[above_thr_x_imps; mean_x_imps>=thr_x_imps(fileNo)];
% 
% 
%             sum_all_y_imps=sum(all_y_imps);
%             mat_sum_all_y_imps=repmat(sum_all_y_imps,size(all_y_imps,1),1);
%             all_y_imps=all_y_imps./mat_sum_all_y_imps; %Normalize to all added imps=1
%             mean_y_imps=mean(all_y_imps,2);
%             thr_y_imps(fileNo)=prctile(mean_y_imps,95);
%             all_mean_y_imps=[all_mean_y_imps; mean_y_imps];
%             above_thr_y_imps=[above_thr_y_imps; mean_y_imps>=thr_y_imps(fileNo)];
% 
%             %Let's look at importance relationships between the different decoders
% 
%             %x vs conc
%             try
%                 close(figureNo)
%             catch
%             end
%             hFig=figure(figureNo);
%             hold on
% 
%             ax=gca;ax.LineWidth=3;
%             set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
%             max_imp=max([max(mean_x_imps) max(mean_conc_imps)]);
%             min_imp_non_zero=min([min(mean_x_imps(mean_x_imps~=0)) min(mean_conc_imps(mean_conc_imps~=0))]);
%             min_imp=0.9*min_imp_non_zero;
%             mean_x_imps(mean_x_imps==0)=min_imp;
%             mean_conc_imps(mean_conc_imps==0)=min_imp;
% 
%             plot(log10(mean_x_imps),log10(mean_conc_imps),'ob','MarkerFaceColor','b');
%             plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
%             xlim([log10(min_imp) log10(max_imp)])
%             ylim([log10(min_imp) log10(max_imp)])
%             plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
%             plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
%             xlabel('Importance x (log10)')
%             ylabel('Importance conc (log10)')
%             title(['Prediction importance for ' arena_file(1:14)], 'Interpreter', 'none')
% 
%             %y vs conc
%             try
%                 close(figureNo+1)
%             catch
%             end
%             hFig=figure(figureNo+1);
%             hold on
% 
%             ax=gca;ax.LineWidth=3;
%             set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
%             max_imp=max([max(mean_y_imps) max(mean_conc_imps)]);
%             min_imp_non_zero=min([min(mean_y_imps(mean_y_imps~=0)) min(mean_conc_imps(mean_conc_imps~=0))]);
%             min_imp=0.9*min_imp_non_zero;
%             mean_y_imps(mean_y_imps==0)=min_imp;
%             mean_conc_imps(mean_conc_imps==0)=min_imp;
% 
%             plot(log10(mean_y_imps),log10(mean_conc_imps),'ob','MarkerFaceColor','b')
%             plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
%             xlim([log10(min_imp) log10(max_imp)])
%             ylim([log10(min_imp) log10(max_imp)])
%             plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
%             plot([log10(min_imp) log10(max_imp)],[log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],'-k')
% 
%             xlabel('Importance y (log10)')
%             ylabel('Importance conc (log10)')
%             title(['Prediction importance for ' arena_file(1:14)], 'Interpreter', 'none')
% 
%             %x vs y
%             try
%                 close(figureNo+2)
%             catch
%             end
%             hFig=figure(figureNo+2);
%             hold on
% 
%             ax=gca;ax.LineWidth=3;
%             set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
%             max_imp=max([max(mean_y_imps) max(mean_x_imps)]);
%             min_imp_non_zero=min([min(mean_y_imps(mean_y_imps~=0)) min(mean_x_imps(mean_x_imps~=0))]);
%             min_imp=0.9*min_imp_non_zero;
%             mean_y_imps(mean_y_imps==0)=min_imp;
%             mean_x_imps(mean_x_imps==0)=min_imp;
% 
%             plot(log10(mean_x_imps),log10(mean_y_imps),'ob','MarkerFaceColor','b')
%             plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
%             xlim([log10(min_imp) log10(max_imp)])
%             ylim([log10(min_imp) log10(max_imp)])
%             plot([log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
%             plot([log10(min_imp) log10(max_imp)],[log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],'-k')
%             xlabel('Importance x (log10)')
%             ylabel('Importance y (log10)')
%             title(['Prediction importance for ' arena_file(1:14)], 'Interpreter', 'none')
% 
%             pffft=1;
%         end
%     end
%     figureNo=figureNo+2;
%     %
%     % %Now plot x vs conc for all above thr
%     % try
%     %     close(figureNo)
%     % catch
%     % end
%     % hFig=figure(figureNo);
%     % hold on
%     %
%     % ax=gca;ax.LineWidth=3;
%     % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
%     %
%     % mean_x_imps_x95=all_mean_x_imps(logical(above_thr_x_imps));
%     % mean_conc_imps_x95=all_mean_conc_imps(logical(above_thr_x_imps));
%     %
%     % mean_conc_imps_conc95=all_mean_conc_imps(logical(above_thr_conc_imps));
%     % mean_x_imps_conc95=all_mean_x_imps(logical(above_thr_conc_imps));
%     %
%     % mean_conc_imps_xconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));
%     % mean_x_imps_xconc95=all_mean_x_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));
%     %
%     % mean_conc_imps_yconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));
%     % mean_y_imps_yconc95=all_mean_y_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));
%     %
%     % plot(log10(mean_x_imps_x95),log10(mean_conc_imps_x95),'o','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7]);
%     % plot(log10(mean_x_imps_conc95),log10(mean_conc_imps_conc95),'o','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0]);
%     % plot(log10(mean_x_imps_xconc95),log10(mean_conc_imps_xconc95),'o','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor',[0.8 0.6 0.7]);
%     % plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
%     % xlim([log10(min_imp) log10(max_imp)])
%     % ylim([log10(min_imp) log10(max_imp)])
%     % plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
%     % plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
%     % xlabel('Importance x (log10)')
%     % ylabel('Importance conc (log10)')
%     % title(['Prediction importance for all files'], 'Interpreter', 'none')
% 
%     mean_x_imps_x95=all_mean_x_imps(logical(above_thr_x_imps));
%     mean_x_imps_y95=all_mean_x_imps(logical(above_thr_y_imps));
%     mean_x_imps_conc95=all_mean_x_imps(logical(above_thr_conc_imps));
% 
%     mean_y_imps_x95=all_mean_y_imps(logical(above_thr_x_imps));
%     mean_y_imps_y95=all_mean_y_imps(logical(above_thr_y_imps));
%     mean_y_imps_conc95=all_mean_y_imps(logical(above_thr_conc_imps));
% 
%     mean_conc_imps_x95=all_mean_conc_imps(logical(above_thr_x_imps));
%     mean_conc_imps_y95=all_mean_conc_imps(logical(above_thr_y_imps));
%     mean_conc_imps_conc95=all_mean_conc_imps(logical(above_thr_conc_imps));
% 
%     mean_conc_imps_xconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));
%     mean_x_imps_xconc95=all_mean_x_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));
% 
%     mean_conc_imps_yconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));
%     mean_y_imps_yconc95=all_mean_y_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));
% 
%     mean_y_imps_xy95=all_mean_y_imps(logical(above_thr_x_imps)&logical(above_thr_y_imps));
%     mean_x_imps_xy95=all_mean_x_imps(logical(above_thr_x_imps)&logical(above_thr_y_imps));
% 
%     % Create a sample data matrix for PCA and t-SNE
%     % Keep track of which are >95% for the different dimensions
%     mean_x_imps=[];
%     mean_y_imps=[];
%     mean_conc_imps=[];
%     ii_95_class=[]; %1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s
%     ii_included=0;
% 
%     for ii=1:length(all_mean_x_imps)
%         include_ii=0;
% 
%         if above_thr_x_imps(ii)==1
%             this_ii_95_class=1;
%             include_ii=1;
%         end
% 
%         if above_thr_y_imps(ii)==1
%             this_ii_95_class=2;
%             include_ii=1;
%         end
% 
%         if above_thr_conc_imps(ii)==1
%             this_ii_95_class=3;
%             include_ii=1;
%         end
% 
%         if (above_thr_x_imps(ii)==1)&(above_thr_conc_imps(ii)==1)
%             this_ii_95_class=4;
%             include_ii=1;
%         end
% 
%         if (above_thr_y_imps(ii)==1)&(above_thr_conc_imps(ii)==1)
%             this_ii_95_class=5;
%             include_ii=1;
%         end
% 
%         if (above_thr_x_imps(ii)==1)&(above_thr_y_imps(ii)==1)
%             this_ii_95_class=6;
%             include_ii=1;
%         end
% 
%         if (above_thr_x_imps(ii)==1)&(above_thr_conc_imps(ii)==1)&(above_thr_y_imps(ii)==1)
%             this_ii_95_class=7;
%             include_ii=1;
%         end
% 
%         if include_ii==1
%             ii_included=ii_included+1;
%             ii_95_class(ii_included)=this_ii_95_class;
%             mean_x_imps(ii_included)=log10(all_mean_x_imps(ii));
%             mean_y_imps(ii_included)=log10(all_mean_y_imps(ii));
%             mean_conc_imps(ii_included)=log10(all_mean_conc_imps(ii));
%         end
% 
%     end
% 
%     data = [mean_x_imps; mean_y_imps; mean_conc_imps];
% 
%     % %This was the old way to do this
%     % mean_x_imps=log10([mean_x_imps_x95' mean_x_imps_y95' mean_x_imps_conc95']);
%     % mean_y_imps=log10([mean_y_imps_x95' mean_y_imps_y95' mean_y_imps_conc95']);
%     % mean_conc_imps=log10([mean_conc_imps_x95' mean_conc_imps_y95' mean_conc_imps_conc95']);
%     %
%     % ii_x_imps=[1:length(mean_x_imps_x95)];
%     % ii_y_imps=[length(mean_x_imps_x95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)];
%     % ii_conc_imps=[length(mean_x_imps_x95)+length(mean_y_imps_y95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)+length(mean_conc_imps_conc95)];
%     %
%     % data = [mean_x_imps; mean_y_imps; mean_conc_imps]; % Replace with your actual data
% 
%     %Now plot x vs conc for all above thr
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
%     plot(log10(mean_x_imps_x95),log10(mean_conc_imps_x95),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0]);
%     plot(log10(mean_x_imps_conc95),log10(mean_conc_imps_conc95),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5]);
%     plot(log10(mean_x_imps_xconc95),log10(mean_conc_imps_xconc95),'o','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor',[0.95 0.9 0.25]);
%     plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
%     xlim([log10(min_imp) log10(max_imp)])
%     ylim([log10(min_imp) log10(max_imp)])
%     plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
%     plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
%     xlabel('Importance x (log10)')
%     ylabel('Importance conc (log10)')
%     title(['Prediction importance above 95 percentile ' group_label{this_group}], 'Interpreter', 'none')
% 
%     %Now plot y vs conc for all above thr
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
%     plot(log10(mean_y_imps_y95),log10(mean_conc_imps_y95),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9]);
%     plot(log10(mean_y_imps_conc95),log10(mean_conc_imps_conc95),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5]);
%     plot(log10(mean_y_imps_yconc95),log10(mean_conc_imps_yconc95),'o','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7]);
%     plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
%     xlim([log10(min_imp) log10(max_imp)])
%     ylim([log10(min_imp) log10(max_imp)])
%     plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
%     plot([log10(min_imp) log10(max_imp)],[log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],'-k')
%     xlabel('Importance y (log10)')
%     ylabel('Importance conc (log10)')
%     title(['Prediction importance above 95 percentile ' group_label{this_group}], 'Interpreter', 'none')
% 
% 
% 
%     %Now plot x vs y for all above thr
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
% 
% 
%     plot(log10(mean_x_imps_x95),log10(mean_y_imps_x95),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0]);
%     plot(log10(mean_x_imps_y95),log10(mean_y_imps_y95),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9]);
%     plot(log10(mean_x_imps_xy95),log10(mean_y_imps_xy95),'o','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0]);
%     plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
%     xlim([log10(min_imp) log10(max_imp)])
%     ylim([log10(min_imp) log10(max_imp)])
%     plot([log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
%     plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
%     xlabel('Importance x (log10)')
%     ylabel('Importance y (log10)')
%     title(['Prediction importance above 95 percentile ' group_label{this_group}], 'Interpreter', 'none')
% 
%     % Step 1: Create a sample data matrix
%     % Assume you have 3 variables and 100 samples
%     %
%     % mean_x_imps=log10([mean_x_imps_x95' mean_x_imps_y95' mean_x_imps_conc95']);
%     % mean_y_imps=log10([mean_y_imps_x95' mean_y_imps_y95' mean_y_imps_conc95']);
%     % mean_conc_imps=log10([mean_conc_imps_x95' mean_conc_imps_y95' mean_conc_imps_conc95']);
%     %
%     % ii_x_imps=[1:length(mean_x_imps_x95)];
%     % ii_y_imps=[length(mean_x_imps_x95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)];
%     % ii_conc_imps=[length(mean_x_imps_x95)+length(mean_y_imps_y95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)+length(mean_conc_imps_conc95)];
%     %
%     % data = [mean_x_imps; mean_y_imps; mean_conc_imps]; % Replace with your actual data
%     % 
%     % % Step 2: Perform PCA
%     % [coeff, score, latent, tsquared, explained] = pca(data');
%     % 
%     % % Step 3: Plot the first two principal components
%     % figureNo=figureNo+1;
%     % try
%     %     close(figureNo)
%     % catch
%     % end
%     % hFig=figure(figureNo);
%     % hold on
%     % 
%     % ax=gca;ax.LineWidth=3;
%     % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
%     % scatter(score(:,1), score(:,2), 'filled');
%     % xlabel('1st Principal Component');
%     % ylabel('2nd Principal Component');
%     % title(['PCA of Data' group_label{this_group}]');
%     % grid on;
%     % 
%     % % Optionally display explained variance
%     % disp('Explained Variance by Each Principal Component:');
%     % disp(explained);
% 
% 
%     % Step 1 Transpose the data if necessary (N samples x 3 variables)
%     data = data';
% 
%     rng('default') % for reproducibility
% 
%     % Step 2: Run t-SNE
%     % The default parameters are usually sufficient, but you can adjust them
%     % For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
%     % mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
%     % mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
%     mappedX = tsne(data,'Algorithm','exact','Distance','cosine'); %works better
% 
%     % mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
%     % mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK
% 
%     % Step 3: Plot the results
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
%     for ii=1:size(data,1)
%         %1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s
%         switch ii_95_class(ii)
%             case 1
%                 %1 is >=x95
%                 plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0]);
%             case 2
%                 %2 is >=y95
%                 plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9]);
%             case 3
%                 %3 is >=conc95
%                 plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5]);
%             case 4
%                 %4 is >=x95 and conc95
%                 plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor',[0.95 0.9 0.25]);
%             case 5
%                 %5 is >=y95 and conc95
%                 plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7]);
%             case 6
%                 %6 is >=x95 and y95
%                 plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0]);
%             case 7
%                 %7 is >= all 95s
%                 plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor',[0.8 0.6 0.7]);
%         end
%     end
% 
%     these_ylim=ylim;
%     these_xlim=xlim;
%     text(these_xlim(1)+0.85*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.95*(these_ylim(2)-these_ylim(1)),'x','Color',[0.9 0.6 0],'FontWeight','bold','FontSize',16)
%     text(these_xlim(1)+0.20*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.92*(these_ylim(2)-these_ylim(1)),'y','Color',[0.35 0.7 0.9],'FontWeight','bold','FontSize',16)
%     text(these_xlim(1)+0.45*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.2*(these_ylim(2)-these_ylim(1)),'odor','Color',[0 0.6 0.5],'FontWeight','bold','FontSize',16)
%     xlabel('t-SNE Component 1');
%     ylabel('t-SNE Component 2');
%     title(['t-SNE Prediction Importance ' group_label{this_group}]);
%     pffft=1;
%     % grid on;
% end

fclose(fileID);

pffft=1;