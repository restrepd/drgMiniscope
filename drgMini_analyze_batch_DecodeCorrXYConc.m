%drgMini_analyze_batch_DecodeCorrXYConc
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
        % 
        % %Trained with all trials
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01062025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01062025.m';

        %Trained with hits only
         save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
        choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

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

%Exclude files with p value for R1_all_trials > 0.05
thr_rho=0.05;

%Now calculate p values for conc R1 all files
ii_run=1;
ii_for_corr=0;
R1_all_trials_pre=[];
P_rho_all_trials_pre=[];
these_fileNos_pre=[];


for fileNo=1:length(handles_conc.arena_file)
    
        op_all_trials=[];
        op_decod_all_trials=[];
        

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
        

        for trNo=1:trials.odor_trNo

            op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
            op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
            if op_predictedend>length(odor_plume_template)
                op_predictedend=length(odor_plume_template);
            end

            op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];


        end

        if ~isempty(op_all_trials)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(op_all_trials,op_decod_all_trials);
                R1_all_trials_pre(fileNo)=R1(1,2);
                P_rho_all_trials_pre(fileNo)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(op_all_trials,op_decod_all_trials,'Type','Spearman');
                R1_all_trials_pre(fileNo)=R1;
                P_rho_all_trials_pre(fileNo)=this_P_rho_all_trials;
            end
        else
            R1_all_trials_pre(fileNo)=NaN;
        end
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

try
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
                if (handles_conc.group(fileNo)==grNo)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
    fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1 cm from floor and bins before\n'])
    fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1 cm from floor and bins before\n']);
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
            if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
            if (sum(handles_XY.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
catch
end

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
        if (sum(handles_conc.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
        if (sum(handles_XY.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
            if (handles_conc.group(fileNo)==this_group)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
            if (sum(handles_conc.group(fileNo)==these_new_groups.gr(grNo).groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' odor in one or both spouts\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' odor in one or both spouts\n']);
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
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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


title(['R1 for prediction of odor hit, miss, shuffled'])
ylabel('R1')
ylim([-0.2 0.65])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for R1 vs trial type\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for R1 vs trial type\n']);
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
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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
R2_all_trials=[];
R2_hits=[];
R2_miss=[];
R2_all_trials_sh=[];
R2sh_all_trials=[];
R2sh_hits=[];
R2sh_miss=[];
R2sh_all_trials_sh=[];
Var_per_point_all=[];
Var_per_point_hits=[];
Var_per_point_miss=[];
Var_per_point_all_sh=[];
P_rho_all_trials=[];
these_fileNos=[];


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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

        if ~isempty(op_all_hits90)
            if is_pearson==1
                R1=corrcoef(op_all_hits90,op_decod_all_hits90);
                R1_hits90(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_hits90,op_decod_all_hits90,'Type','Spearman');
                R1_hits90(ii_for_corr)=R1;
            end
        else
            R1_hits90(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hitso)
            if is_pearson==1
                R1=corrcoef(op_all_hitso,op_decod_all_hitso);
                R1_hitso(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_hitso,op_decod_all_hitso,'Type','Spearman');
                R1_hitso(ii_for_corr)=R1;
            end
        else
            R1_hitso(ii_for_corr)=NaN;
        end

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

        if ~isempty(op_between_trials)
            if is_pearson==1
                R1=corrcoef(op_between_trials,op_decod_between_trials);
                R1_between(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_between_trials,op_decod_between_trials,'Type','Spearman');
                R1_between(ii_for_corr)=R1;
            end
        else
            R1_between(ii_for_corr)=NaN;
        end

        %Calculate R2, the fraction of variance explained
        if ~isempty(op_all_trials)
            this_R2=1-sum((op_decod_all_trials-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
            R2_all_trials(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_all(ii_for_corr)=sum((op_all_trials-mean(op_all_trials)).^2)/length(op_all_trials);
        else
            R2_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits)
            this_R2=1-sum((op_decod_all_hits-op_all_hits).^2)/(sum((op_all_hits-mean(op_all_hits)).^2));
            R2_hits(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_hits(ii_for_corr)=(sum((op_all_hits-mean(op_all_hits)).^2))/length(op_all_hits);
        else
            R2_hits(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_miss)
            this_R2=1-sum((op_decod_all_miss-op_all_miss).^2)/(sum((op_all_miss-mean(op_all_miss)).^2));
            R2_miss(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_miss(ii_for_corr)=(sum((op_all_miss-mean(op_all_miss)).^2))/length(op_all_miss);
        else
            R2_miss(ii_for_corr)=NaN;
        end

        %Shuffled
        if ~isempty(op_predicted_sh_conv)
            these_R2s=[];
            these_Vars=[];
            for ii_sh=1:size(op_predicted_sh_conv,2)
                this_op_all_trials_sh=op_decod_trials_sh.ii_sh(ii_sh).oppsh;
                this_R2=1-sum((this_op_all_trials_sh-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
                these_R2s=[these_R2s drgMini_rectifyR2(this_R2)];
                these_Vars=[these_Vars sum((op_all_trials-mean(op_all_trials)).^2)/length(op_all_trials)];
            end
            R2_all_trials_sh(ii_for_corr)=mean(these_R2s);
            Var_per_point_all_sh(ii_for_corr)=mean(these_Vars);
        else
            R2_all_trials_sh(ii_for_corr)=NaN;
        end

        %Calculate R2sh, the fraction of variance shared
        if ~isempty(op_all_trials)
            if sum((op_all_trials-mean(op_all_trials)).^2)>=sum((op_decod_all_trials-mean(op_decod_all_trials)).^2)
                this_R2sh=1-sum((op_decod_all_trials-op_all_trials).^2)/sum((op_all_trials-mean(op_all_trials)).^2);
            else
                this_R2sh=1-sum((op_decod_all_trials-op_all_trials).^2)/sum((op_decod_all_trials-mean(op_decod_all_trials)).^2);
            end
            R2sh_all_trials(ii_for_corr)=this_R2sh;
        else
            R2sh_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits)
            if (sum((op_all_hits-mean(op_all_hits)).^2))>=(sum((op_decod_all_hits-mean(op_decod_all_hits)).^2))
                this_R2sh=1-sum((op_decod_all_hits-op_all_hits).^2)/(sum((op_all_hits-mean(op_all_hits)).^2));
            else
                this_R2sh=1-sum((op_decod_all_hits-op_all_hits).^2)/(sum((op_decod_all_hits-mean(op_decod_all_hits)).^2));
            end
            R2sh_hits(ii_for_corr)=this_R2sh;
        else
            R2_hits(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_miss)
            if (sum((op_all_miss-mean(op_all_miss)).^2))>=(sum((op_decod_all_miss-mean(op_decod_all_miss)).^2))
                this_R2sh=1-sum((op_decod_all_miss-op_all_miss).^2)/(sum((op_all_miss-mean(op_all_miss)).^2));
            else
                this_R2sh=1-sum((op_decod_all_miss-op_all_miss).^2)/(sum((op_decod_all_miss-mean(op_decod_all_miss)).^2));
            end
            R2sh_miss(ii_for_corr)=this_R2sh;
        else
            R2sh_miss(ii_for_corr)=NaN;
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
            R2sh_all_trials_sh(ii_for_corr)=mean(these_R2shs);
        else
            R2sh_all_trials_sh(ii_for_corr)=NaN;
        end

    end
end

%Make R2sh zero for shuffled
R2sh_all_trials=(R2sh_all_trials-mean(R2sh_all_trials_sh))/(1-mean(R2sh_all_trials_sh));
R2sh_hits=(R2sh_hits-mean(R2sh_all_trials_sh))/(1-mean(R2sh_all_trials_sh));
R2sh_miss=(R2sh_miss-mean(R2sh_all_trials_sh))/(1-mean(R2sh_all_trials_sh));
R2sh_all_trials_sh=(R2sh_all_trials_sh-mean(R2sh_all_trials_sh))/(1-mean(R2sh_all_trials_sh));

%Let's plot an R2 bar graph for all, hits, miss and shuffled 
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

R2type=[];
R2type{1}='R2_all_trials';
R2type{2}='R2_hits';
R2type{3}='R2_miss';
% R2type{2}='R2mod_hits';
% R2type{3}='R2mod_miss';
R2type{4}='R2_all_trials_sh';


R2type_label=[];
R2type_label{1}='all';
R2type_label{2}='hits';
R2type_label{3}='miss';
R2type_label{4}='shuffled';

allR1s=[];
for ii_R2type=1:length(R2type)

    these_R2s=[];
    eval(['these_R2s=' R2type{ii_R2type} ';'])
   
    for ii=1:length(these_R2s)
        all_R1s.R2type(ii_R2type).R1(ii)=these_R2s(ii);
    end

    these_R2s=these_R2s(~isnan(these_R2s));

    %plot bar
    switch ii_R2type
        case 1
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, violin_x]=drgViolinPoint(these_R2s...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;
 
    ii_next=0;
    for ii=1:length(all_R1s.R2type(ii_R2type).R1)
        if ~isnan(all_R1s.R2type(ii_R2type).R1(ii))
            ii_next=ii_next+1;
            all_R1s.R2type(ii_R2type).violin_x(ii)=violin_x(ii_next);
        end
    end

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=these_R2s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=ii_R2type*ones(1,length(these_R2s));
    glm_r1_ii=glm_r1_ii+length(these_R2s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2s;
    input_r1_data(id_r1_ii).description=[R2type_label{ii_R2type}];

end

%Plot lines between points
for ii=1:length(all_R1s.R2type(ii_R2type).R1)
    these_R2s=[];
    these_violin_x=[];
    no_nans=1;
    for ii_R2type=1:length(R2type)
        if isnan(all_R1s.R2type(ii_R2type).R1(ii))
            no_nans=0;
        end
    end

    if no_nans==1
        for ii_R2type=1:length(R2type)
            these_R2s=[these_R2s all_R1s.R2type(ii_R2type).R1(ii)];
            these_violin_x=[these_violin_x all_R1s.R2type(ii_R2type).violin_x(ii)];
        end

        plot(these_violin_x,these_R2s,'-','Color',[0.7 0.7 0.7])
    end
end

% x_pos=3;
% text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 2 3])
xticklabels({'all','hits','miss','shuffled'})


title(['R2 for predcition of odor concentration'])
ylabel('R2')
ylim([-1 0.65])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R2 vs trial type\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R2 vs trial type\n']);
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


%Let's plot an R2sh bar graph for all, hits, miss and shuffled 
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

R2shtype=[];
R2shtype{1}='R2sh_all_trials';
R2shtype{2}='R2sh_hits';
R2shtype{3}='R2sh_miss';
% R2shtype{2}='R2shmod_hits';
% R2shtype{3}='R2shmod_miss';
R2shtype{4}='R2sh_all_trials_sh';


R2shtype_label=[];
R2shtype_label{1}='all';
R2shtype_label{2}='hits';
R2shtype_label{3}='miss';
R2shtype_label{4}='shuffled';

allR1s=[];
for ii_R2shtype=1:length(R2shtype)

    these_R2shs=[];
    eval(['these_R2shs=' R2shtype{ii_R2shtype} ';'])
   
    for ii=1:length(these_R2shs)
        all_R1s.R2shtype(ii_R2shtype).R1(ii)=these_R2shs(ii);
    end

    these_R2shs=these_R2shs(~isnan(these_R2shs));

    %plot bar
    switch ii_R2shtype
        case 1
            bar(bar_offset,mean(these_R2shs),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_R2shs),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_R2shs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_R2shs),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_R2shs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, violin_x]=drgViolinPoint(these_R2shs...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;
 
    ii_next=0;
    for ii=1:length(all_R1s.R2shtype(ii_R2shtype).R1)
        if ~isnan(all_R1s.R2shtype(ii_R2shtype).R1(ii))
            ii_next=ii_next+1;
            all_R1s.R2shtype(ii_R2shtype).violin_x(ii)=violin_x(ii_next);
        end
    end

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R2shs))=these_R2shs;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R2shs))=ii_R2shtype*ones(1,length(these_R2shs));
    glm_r1_ii=glm_r1_ii+length(these_R2shs);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2shs;
    input_r1_data(id_r1_ii).description=[R2shtype_label{ii_R2shtype}];

end

%Plot lines between points
for ii=1:length(all_R1s.R2shtype(ii_R2shtype).R1)
    these_R2shs=[];
    these_violin_x=[];
    no_nans=1;
    for ii_R2shtype=1:length(R2shtype)
        if isnan(all_R1s.R2shtype(ii_R2shtype).R1(ii))
            no_nans=0;
        end
    end

    if no_nans==1
        for ii_R2shtype=1:length(R2shtype)
            these_R2shs=[these_R2shs all_R1s.R2shtype(ii_R2shtype).R1(ii)];
            these_violin_x=[these_violin_x all_R1s.R2shtype(ii_R2shtype).violin_x(ii)];
        end

        plot(these_violin_x,these_R2shs,'-','Color',[0.7 0.7 0.7])
    end
end

% x_pos=3;
% text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 2 3])
xticklabels({'all','hits','miss','shuffled'})


title(['R2sh for predcition of odor concentration'])
ylabel('R2sh')
ylim([-0.5 0.7])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R2sh vs trial type\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R2sh vs trial type\n']);
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

%Let's plot an variance bar graph for all, hits, miss and shuffled 
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

Vartype=[];
Vartype{1}='Var_per_point_all';
Vartype{2}='Var_per_point_hits';
Vartype{3}='Var_per_point_miss';
% Vartype{2}='R2mod_hits';
% Vartype{3}='R2mod_miss';
Vartype{4}='Var_per_point_all_sh';

Vartype_label{1}='all';
Vartype_label{2}='hits';
Vartype_label{3}='miss';
% Vartype{2}='R2mod_hits';
% Vartype{3}='R2mod_miss';
Vartype_label{4}='shuffled';

allR1s=[];
for ii_Vartype=1:length(Vartype)

    these_Vars=[];
    eval(['these_Vars=' Vartype{ii_Vartype} ';'])
   
    for ii=1:length(these_Vars)
        all_R1s.Vartype(ii_Vartype).R1(ii)=these_Vars(ii);
    end

    these_Vars=these_Vars(~isnan(these_Vars));

    %plot bar
    switch ii_Vartype
        case 1
            bar(bar_offset,mean(these_Vars),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_Vars),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_Vars),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_Vars),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_Vars),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, violin_x]=drgViolinPoint(these_Vars...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;
 
    ii_next=0;
    for ii=1:length(all_R1s.Vartype(ii_Vartype).R1)
        if ~isnan(all_R1s.Vartype(ii_Vartype).R1(ii))
            ii_next=ii_next+1;
            all_R1s.Vartype(ii_Vartype).violin_x(ii)=violin_x(ii_next);
        end
    end

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=these_Vars;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=ii_Vartype*ones(1,length(these_Vars));
    glm_r1_ii=glm_r1_ii+length(these_Vars);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_Vars;
    input_r1_data(id_r1_ii).description=[Vartype_label{ii_Vartype}];

end

%Plot lines between points
for ii=1:length(all_R1s.Vartype(ii_Vartype).R1)
    these_Vars=[];
    these_violin_x=[];
    no_nans=1;
    for ii_Vartype=1:length(Vartype)
        if isnan(all_R1s.Vartype(ii_Vartype).R1(ii))
            no_nans=0;
        end
    end

    if no_nans==1
        for ii_Vartype=1:length(Vartype)
            these_Vars=[these_Vars all_R1s.Vartype(ii_Vartype).R1(ii)];
            these_violin_x=[these_violin_x all_R1s.Vartype(ii_Vartype).violin_x(ii)];
        end

        plot(these_violin_x,these_Vars,'-','Color',[0.7 0.7 0.7])
    end
end

% x_pos=3;
% text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 2 3])
xticklabels({'all','hits','miss','shuffled'})


title(['Variance per point for predcition of odor concentration'])
ylabel('Variance')
ylim([-1 20])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' variance vs trial type\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' variance vs trial type\n']);
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




%Let's plot R1 vs variance bar graph for hits and miss
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

Vartype=[];
Vartype{1}='Var_per_point_hits';
Vartype{2}='Var_per_point_miss';

type_label{1}='hits';
type_label{2}='miss';

R1type{1}='R1_hits';
R1type{2}='R1_miss';





for ii_Vartype=1:length(Vartype)

    these_Vars=[];
    eval(['these_Vars=' Vartype{ii_Vartype} ';'])

    these_R1s=[];
    eval(['these_R1s=' R1type{ii_Vartype} ';'])


    eval(['plot(these_Vars, these_R1s,''' marker_per_group{ii_Vartype} ''',''MarkerSize'',8,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_Vartype} ')'])


    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=these_R1s;
    glm_r1.variance(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=these_Vars;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=ii_Vartype*ones(1,length(these_Vars));
    glm_r1_ii=glm_r1_ii+length(these_Vars);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R1s;
    input_r1_data(id_r1_ii).description=[Vartype_label{ii_Vartype}];

end


title(['Odor R1 vs. variance per point for predcition of odor concentration'])
xlabel('Variance')
ylabel('R1')

ylim([-0.2 0.8])
xlim([0 13])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' odor R1 vs variance and hit vs miss (trial type)\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' odor R1 vs variance and hit vs miss (trial type)\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.variance',glm_r1.trial_type',...
    'VariableNames',{'R1','variance','trial_type'});
mdl = fitglm(tbl,'R1~variance+trial_type'...
    ,'CategoricalVars',[3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R1 vs variance and hit vs miss (trial type)\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R1 vs variance and hit vs miss (trial type)e\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);



%Let's plot R2 vs variance bar graph for hits and miss
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

Vartype=[];
Vartype{1}='Var_per_point_hits';
Vartype{2}='Var_per_point_miss';

type_label{1}='hits';
type_label{2}='miss';

R2type{1}='R2_hits';
R2type{2}='R2_miss';





for ii_Vartype=1:length(Vartype)

    these_Vars=[];
    eval(['these_Vars=' Vartype{ii_Vartype} ';'])

    these_R2s=[];
    eval(['these_R2s=' R2type{ii_Vartype} ';'])


    eval(['plot(these_Vars, these_R2s,''' marker_per_group{ii_Vartype} ''',''MarkerSize'',8,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_Vartype} ')'])


    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=these_R2s;
    glm_r1.variance(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=these_Vars;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=ii_Vartype*ones(1,length(these_Vars));
    glm_r1_ii=glm_r1_ii+length(these_Vars);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2s;
    input_r1_data(id_r1_ii).description=[Vartype_label{ii_Vartype}];

end


title(['Odor R2 vs. variance per point for predcition of odor concentration'])
xlabel('Variance')
ylabel('R2')

ylim([-1 0.5])
xlim([0 13])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' odor R2 vs variance and hit vs miss (trial type)\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' odor R2 vs variance and hit vs miss (trial type)\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.variance',glm_r1.trial_type',...
    'VariableNames',{'R2','variance','trial_type'});
mdl = fitglm(tbl,'R2~variance+trial_type'...
    ,'CategoricalVars',[3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2 vs variance and hit vs miss (trial type)\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2 vs variance and hit vs miss (trial type)e\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);



%Let's plot R2sh vs variance bar graph for hits and miss
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

Vartype=[];
Vartype{1}='Var_per_point_hits';
Vartype{2}='Var_per_point_miss';

type_label{1}='hits';
type_label{2}='miss';

R2shtype{1}='R2sh_hits';
R2shtype{2}='R2sh_miss';





for ii_Vartype=1:length(Vartype)

    these_Vars=[];
    eval(['these_Vars=' Vartype{ii_Vartype} ';'])

    these_R2shs=[];
    eval(['these_R2shs=' R2shtype{ii_Vartype} ';'])


    eval(['plot(these_Vars, these_R2shs,''' marker_per_group{ii_Vartype} ''',''MarkerSize'',8,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_Vartype} ')'])


    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=these_R2shs;
    glm_r1.variance(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=these_Vars;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_Vars))=ii_Vartype*ones(1,length(these_Vars));
    glm_r1_ii=glm_r1_ii+length(these_Vars);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2shs;
    input_r1_data(id_r1_ii).description=[Vartype_label{ii_Vartype}];

end


title(['Odor R2sh vs. variance per point for predcition of odor concentration'])
xlabel('Variance')
ylabel('R2sh')

ylim([-1 0.5])
xlim([0 13])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' odor R2sh vs variance and hit vs miss (trial type)\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' odor R2sh vs variance and hit vs miss (trial type)\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.variance',glm_r1.trial_type',...
    'VariableNames',{'R2sh','variance','trial_type'});
mdl = fitglm(tbl,'R2sh~variance+trial_type'...
    ,'CategoricalVars',[3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2sh vs variance and hit vs miss (trial type)\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2sh vs variance and hit vs miss (trial type)e\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Let's plot R2 vs R1 for hits and miss
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
R1type{1}='R1_hits';
R1type{2}='R1_miss';

type_label{1}='hits';
type_label{2}='miss';

R2type{1}='R2_hits';
R2type{2}='R2_miss';





for ii_R1type=1:length(R1type)

    these_R1s=[];
    eval(['these_R1s=' R1type{ii_R1type} ';'])

    these_R2s=[];
    eval(['these_R2s=' R2type{ii_R1type} ';'])


    eval(['plot(these_R1s, these_R2s,''' marker_per_group{ii_R1type} ''',''MarkerSize'',8,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_R1type} ')'])


    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R2s;
    glm_r1.R1(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_R1type*ones(1,length(these_R1s));
    glm_r1_ii=glm_r1_ii+length(these_R1s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2s;
    input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

end


title(['Odor R2 vs. R1 per point for predcition of odor concentration'])
xlabel('R1')
ylabel('R2')

ylim([-1 0.5])
xlim([-0.1 0.6])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' odor R2 vs R1 and hit vs miss (trial type)\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' odor R2 vs R1 and hit vs miss (trial type)\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.R1',glm_r1.trial_type',...
    'VariableNames',{'R2','R1','trial_type'});
mdl = fitglm(tbl,'R2~R1+trial_type'...
    ,'CategoricalVars',[3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2 vs R1 and hit vs miss (trial type)\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2 vs R1 and hit vs miss (trial type)e\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Let's plot R2sh vs R1 for hits and miss
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
R1type{1}='R1_hits';
R1type{2}='R1_miss';

type_label{1}='hits';
type_label{2}='miss';

R2shtype{1}='R2sh_hits';
R2shtype{2}='R2sh_miss';





for ii_R1type=1:length(R1type)

    these_R1s=[];
    eval(['these_R1s=' R1type{ii_R1type} ';'])

    these_R2shs=[];
    eval(['these_R2shs=' R2shtype{ii_R1type} ';'])


    eval(['plot(these_R1s, these_R2shs,''' marker_per_group{ii_R1type} ''',''MarkerSize'',8,''MarkerEdgeColor'',''none'',''MarkerFaceColor'', ' color_okabe_ito{ii_R1type} ')'])


    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R2shs;
    glm_r1.R1(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_R1type*ones(1,length(these_R1s));
    glm_r1_ii=glm_r1_ii+length(these_R1s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2shs;
    input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

end


title(['Odor R2sh vs. R1 per point for predcition of odor concentration'])
xlabel('R1')
ylabel('R2sh')

ylim([-0.6 0.8])
xlim([-0.1 0.6])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' odor R2sh vs R1 and hit vs miss (trial type)\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' odor R2sh vs R1 and hit vs miss (trial type)\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.R1',glm_r1.trial_type',...
    'VariableNames',{'R2sh','R1','trial_type'});
mdl = fitglm(tbl,'R2sh~R1+trial_type'...
    ,'CategoricalVars',[3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2sh vs R1 and hit vs miss (trial type)\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' odor R2sh vs R1 and hit vs miss (trial type)e\n']);

 
[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

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
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' trial type\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' trial type\n']);
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
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
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

            
        end

        ii_for_corr=ii_for_corr+1;
        if ~isempty(op_all_trials_bt)
            if is_pearson==1
                R1=corrcoef(op_all_trials_bt,op_decod_all_trials_bt);
                R1_all_trials_bt(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_trials_bt,op_decod_all_trials_bt,'Type','Spearman');
                R1_all_trials_bt(ii_for_corr)=R1;
            end
        else
            R1_all_trials_bt(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_trials_at)
            if is_pearson==1
                R1=corrcoef(op_all_trials_at,op_decod_all_trials_at);
                R1_all_trials_at(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_trials_at,op_decod_all_trials_at, 'Type','Spearman');
                R1_all_trials_at(ii_for_corr)=R1;
            end
        else
            R1_all_trials_at(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits_bt)
            if is_pearson==1
                R1=corrcoef(op_all_hits_bt,op_decod_all_hits_bt);
                R1_hits_bt(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_hits_bt,op_decod_all_hits_bt,'Type','Spearman');
                R1_hits_bt(ii_for_corr)=R1;
            end
        else
            R1_hits_bt(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_hits_at)
            if is_pearson==1
                R1=corrcoef(op_all_hits_at,op_decod_all_hits_at);
                R1_hits_at(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_hits_at,op_decod_all_hits_at,'Type','Spearman');
                R1_hits_at(ii_for_corr)=R1;
            end
        else
            R1_hits_at(ii_for_corr)=NaN;
        end

         if ~isempty(op_all_miss_bt)
             if is_pearson==1
                 R1=corrcoef(op_all_miss_bt,op_decod_all_miss_bt);
                 R1_miss_bt(ii_for_corr)=R1(1,2);
             else
                 R1=corr(op_all_miss_bt,op_decod_all_miss_bt,'Type','Spearman');
                 R1_miss_bt(ii_for_corr)=R1;
             end
        else
            R1_miss_bt(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_miss_at)
            if is_pearson==1
                R1=corrcoef(op_all_miss_at,op_decod_all_miss_at);
                R1_miss_at(ii_for_corr)=R1(1,2);
            else
                R1=corr(op_all_miss_at,op_decod_all_miss_at,'Type','Spearman');
                R1_miss_at(ii_for_corr)=R1;
            end
        else
            R1_miss_at(ii_for_corr)=NaN;
        end

        these_R1s=[];
        for ii_sh=1:handles_conc.n_shuffle
            if ~isempty(op_all_trials_bt_sh.ii_sh(ii_sh).allt)
                if is_pearson==1
                    R1=corrcoef(op_all_trials_bt_sh.ii_sh(ii_sh).allt,op_decod_all_trials_bt_sh.ii_sh(ii_sh).allt);
                    these_R1s=[these_R1s R1(1,2)];
                else
                    R1=corr(op_all_trials_bt_sh.ii_sh(ii_sh).allt,op_decod_all_trials_bt_sh.ii_sh(ii_sh).allt,'Type','Spearman');
                    these_R1s=[these_R1s R1];
                end

            else
                these_R1s=[these_R1s NaN];
            end
        end
        R1_all_trials_bt_sh(ii_for_corr)=mean(these_R1s);

         these_R1s=[];
         for ii_sh=1:handles_conc.n_shuffle
             if ~isempty(op_all_trials_at_sh.ii_sh(ii_sh).allt)
                 if is_pearson==1
                     R1=corrcoef(op_all_trials_at_sh.ii_sh(ii_sh).allt,op_decod_all_trials_at_sh.ii_sh(ii_sh).allt);
                     these_R1s=[these_R1s R1(1,2)];
                 else
                     R1=corr(op_all_trials_at_sh.ii_sh(ii_sh).allt,op_decod_all_trials_at_sh.ii_sh(ii_sh).allt,'Type','Spearman');
                     these_R1s=[these_R1s R1];
                end
            else
                these_R1s=[these_R1s NaN];
            end
        end
        R1_all_trials_at_sh(ii_for_corr)=mean(these_R1s);
 
        %Now calculate R2
        if ~isempty(op_all_trials)
            this_R2=1-sum((op_decod_all_trials_bt-op_all_trials_bt).^2)/sum((op_all_trials_bt-mean(op_all_trials_bt)).^2);
            R2_all_trials_bt(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_all_bt(ii_for_corr)=sum((op_all_trials-mean(op_all_trials)).^2)/length(op_all_trials);
        else
            R2_all_trials_bt(ii_for_corr)=NaN;
        end

        % if ~isempty(op_all_hits)
        %     this_R2=1-sum((op_decod_all_hits-op_all_hits).^2)/( ((length(op_all_hits)/length(op_all_trials))...
        %         *sum((op_all_trials-mean(op_all_trials)).^2)));
        %     R2mod_hits(ii_for_corr)=drgMini_rectifyR2(this_R2);
        % else
        %     R2mod_hits(ii_for_corr)=NaN;
        % end
        % 
        % if ~isempty(op_all_miss)
        %     this_R2=1-sum((op_decod_all_miss-op_all_miss).^2)/( ((length(op_all_miss)/length(op_all_trials))...
        %         *sum((op_all_trials-mean(op_all_trials)).^2)));
        %     R2mod_miss(ii_for_corr)=drgMini_rectifyR2(this_R2);
        % else
        %     R2mod_miss(ii_for_corr)=NaN;
        % end

        if ~isempty(op_all_hits)
            this_R2=1-sum((op_decod_all_hits_bt-op_all_hits_bt).^2)/(sum((op_all_hits_bt-mean(op_all_hits_bt)).^2));
            R2_hits_bt(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_hits_bt(ii_for_corr)=(sum((op_all_hits_bt-mean(op_all_hits_bt)).^2))/length(op_all_hits_bt);
        else
            R1_hits(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_miss)
            this_R2=1-sum((op_decod_all_miss_bt-op_all_miss_bt).^2)/(sum((op_all_miss_bt-mean(op_all_miss_bt)).^2));
            R2_miss_bt(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_miss_bt(ii_for_corr)=(sum((op_all_miss_bt-mean(op_all_miss_bt)).^2))/length(op_all_miss_bt);
        else
            R2_miss(ii_for_corr)=NaN;
        end

        %Shuffled
        if ~isempty(op_decod_all_trials_bt_sh)
            these_R2s=[];
            these_Vars=[];
            for ii_sh=1:length(op_decod_all_trials_bt_sh.ii_sh)
                this_op_all_trials_sh_bt=op_decod_all_trials_bt_sh.ii_sh(ii_sh).allt;
                this_R2=1-sum((this_op_all_trials_sh_bt-op_all_trials_bt).^2)/sum((op_all_trials_bt-mean(op_all_trials_bt)).^2);
                these_R2s=[these_R2s drgMini_rectifyR2(this_R2)];
                these_Vars=[these_Vars sum((op_all_trials_bt-mean(op_all_trials_bt)).^2)/length(op_all_trials_bt)];
            end
            R2_all_trials_sh_bt(ii_for_corr)=mean(these_R2s);
            Var_per_point_all_sh_bt(ii_for_corr)=mean(these_Vars);
        else
            R2_all_trials_sh_bt(ii_for_corr)=NaN;
        end

         if ~isempty(op_all_trials)
            this_R2=1-sum((op_decod_all_trials_at-op_all_trials_at).^2)/sum((op_all_trials_at-mean(op_all_trials_at)).^2);
            R2_all_trials_at(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_all_at(ii_for_corr)=sum((op_all_trials_at-mean(op_all_trials_at)).^2)/length(op_all_trials_at);
        else
            R2_all_trials_at(ii_for_corr)=NaN;
        end

        % if ~isempty(op_all_hits)
        %     this_R2=1-sum((op_decod_all_hits-op_all_hits).^2)/( ((length(op_all_hits)/length(op_all_trials))...
        %         *sum((op_all_trials-mean(op_all_trials)).^2)));
        %     R2mod_hits(ii_for_corr)=drgMini_rectifyR2(this_R2);
        % else
        %     R2mod_hits(ii_for_corr)=NaN;
        % end
        % 
        % if ~isempty(op_all_miss)
        %     this_R2=1-sum((op_decod_all_miss-op_all_miss).^2)/( ((length(op_all_miss)/length(op_all_trials))...
        %         *sum((op_all_trials-mean(op_all_trials)).^2)));
        %     R2mod_miss(ii_for_corr)=drgMini_rectifyR2(this_R2);
        % else
        %     R2mod_miss(ii_for_corr)=NaN;
        % end

        if ~isempty(op_all_hits)
            this_R2=1-sum((op_decod_all_hits_at-op_all_hits_at).^2)/(sum((op_all_hits_at-mean(op_all_hits_at)).^2));
            R2_hits_at(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_hits_at(ii_for_corr)=(sum((op_all_hits_at-mean(op_all_hits_at)).^2))/length(op_all_hits_at);
        else
            R1_hits(ii_for_corr)=NaN;
        end

        if ~isempty(op_all_miss)
            this_R2=1-sum((op_decod_all_miss_at-op_all_miss_at).^2)/(sum((op_all_miss_at-mean(op_all_miss_at)).^2));
            R2_miss_at(ii_for_corr)=drgMini_rectifyR2(this_R2);
            Var_per_point_miss_at(ii_for_corr)=(sum((op_all_miss_at-mean(op_all_miss_at)).^2))/length(op_all_miss_at);
        else
            R2_miss(ii_for_corr)=NaN;
        end

        %Shuffled
      if ~isempty(op_decod_all_trials_at_sh)
            these_R2s=[];
            these_Vars=[];
            for ii_sh=1:length(op_decod_all_trials_at_sh.ii_sh)
                this_op_all_trials_sh_at=op_decod_all_trials_at_sh.ii_sh(ii_sh).allt;
                this_R2=1-sum((this_op_all_trials_sh_at-op_all_trials_at).^2)/sum((op_all_trials_at-mean(op_all_trials_at)).^2);
                these_R2s=[these_R2s drgMini_rectifyR2(this_R2)];
                these_Vars=[these_Vars sum((op_all_trials_at-mean(op_all_trials_at)).^2)/length(op_all_trials_at)];
            end
            R2_all_trials_sh_at(ii_for_corr)=mean(these_R2s);
            Var_per_point_all_sh_at(ii_for_corr)=mean(these_Vars);
        else
            R2_all_trials_sh_at(ii_for_corr)=NaN;
        end

    end
end


%Let's plot R1 within vs, within before and after last turn
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
xlim([-1 3])
 
%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1 before after last turn\n'])
fprintf(fileID, ['\nglm for Fig 1 ' num2str(figureNo) ' before after last turn\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 before after last turn\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 before after last turn\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

 

%Let's plot R2 within vs, within before and after last turn
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

R2type=[];
R2type{1}='R2_all_trials';
R2type{2}='R2_all_trials_bt';
R2type{3}='R2_all_trials_at';


R2type_label=[];
R2type_label{1}='within';
R2type_label{2}='before';
R2type_label{3}='after';

allR1s=[];
for ii_R2type=1:length(R2type)

    these_R2s=[];
    eval(['these_R2s=' R2type{ii_R2type} ';'])
   
    for ii=1:length(these_R2s)
        all_R1s.R2type(ii_R2type).R1(ii)=these_R2s(ii);
    end

    these_R2s=these_R2s(~isnan(these_R2s));

    %plot bar
    switch ii_R2type
        case 1
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, violin_x]=drgViolinPoint(these_R2s...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;
 
    ii_next=0;
    for ii=1:length(all_R1s.R2type(ii_R2type).R1)
        if ~isnan(all_R1s.R2type(ii_R2type).R1(ii))
            ii_next=ii_next+1;
            all_R1s.R2type(ii_R2type).violin_x(ii)=violin_x(ii_next);
        end
    end

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=these_R2s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=ii_R2type*ones(1,length(these_R2s));
    glm_r1_ii=glm_r1_ii+length(these_R2s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2s;
    input_r1_data(id_r1_ii).description=[R2type_label{ii_R2type}];

end

%Plot lines between points
for ii=1:length(all_R1s.R2type(ii_R2type).R1)
    these_R2s=[];
    these_violin_x=[];
    no_nans=1;
    for ii_R2type=1:length(R2type)
        if isnan(all_R1s.R2type(ii_R2type).R1(ii))
            no_nans=0;
        end
    end
 
    if no_nans==1
        for ii_R2type=1:length(R2type)
            these_R2s=[these_R2s all_R1s.R2type(ii_R2type).R1(ii)];
            these_violin_x=[these_violin_x all_R1s.R2type(ii_R2type).violin_x(ii)];
        end

        plot(these_violin_x,these_R2s,'-','Color',[0.7 0.7 0.7])
    end
end

% x_pos=3;
% text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 2 3])
xticklabels({'within','before','after'})


title(['R2 for prediction of odor concentration before and after last turn'])
ylabel('R2')
ylim([-0.5 0.5])
xlim([-1 3])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R2 before after last turn\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R2 before after last turn\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 before after last turn\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 before after last turn\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

 
%Do R1 bar graph for 1 and 2 cm all trials vs hit, miss, after last turn
ii_run=1;
R1type=[];
R1type{1}='R1_all_trials_at';
R1type{2}='R1_hits_at';
R1type{3}='R1_miss_at';
R1type{4}='R1_all_trials_at_sh';

R1type_label=[];
R1type_label{1}='within';
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


title(['R1 for prediction of odor hit, miss, shuffled, after last turn'])
ylabel('R1')
ylim([-0.2 0.8])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' after turn R1\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' after turn R1\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' after turn R1\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' after turn R1\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

 
%Do R2 bar graph for 1 and 2 cm all trials vs hit, miss, after last turn
ii_run=1;
R2type=[];
R2type{1}='R2_all_trials_at';
R2type{2}='R2_hits_at';
R2type{3}='R2_miss_at';
R2type{4}='R2_all_trials_sh_at';

R2type_label=[];
R2type_label{1}='within';
R2type_label{2}='hits';
R2type_label{3}='miss';
R2type_label{4}='shuffled';

iiR2type_reorder=[2 1 3 4];

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

for ii_R2type=1:length(R2type)

    these_R1s=[];

    eval(['these_R2s=' R2type{ii_R2type} ';'])

    all_R1s.R2type(ii_R2type).R1=these_R2s;

    %plot bar
    switch ii_R2type
        case 1
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, all_R1s.R2type(ii_R2type).violin_x]=drgViolinPoint(these_R2s...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=these_R2s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=iiR2type_reorder(ii_R2type)*ones(1,length(these_R2s));
    glm_r1_ii=glm_r1_ii+length(these_R2s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2s;
    input_r1_data(id_r1_ii).description=[R2type_label{ii_R2type}];

end
bar_offset=bar_offset+1;

%Plot lines between points
for ii=1:length(these_R2s)
    these_R2s=[];
    these_x=[];
    for ii_R2type=1:length(R2type)
        these_R2s=[these_R2s all_R1s.R2type(ii_R2type).R1(ii)];
        these_x=[these_x all_R1s.R2type(ii_R2type).violin_x(ii)];
    end
    plot(these_x,these_R2s,'-','Color',[0.7 0.7 0.7])
end


xticks([0 1 2 3])
xticklabels({'within','hit','miss','shuffled'})


title(['R2 for prediction of odor hit, miss, shuffled, after last turn'])
ylabel('R2')
ylim([-1 0.5])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' after turn R2\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' after turn R2\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' after turn R1\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' after turn R1\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Do variance bar graph for 1 and 2 cm all trials vs hit, miss, after last turn
ii_run=1;
R2type=[];
R2type{1}='Var_per_point_all_at';
R2type{2}='Var_per_point_hits_at';
R2type{3}='Var_per_point_miss_at';
R2type{4}='Var_per_point_all_sh_at';

R2type_label=[];
R2type_label{1}='within';
R2type_label{2}='hits';
R2type_label{3}='miss';
R2type_label{4}='shuffled';

iiR2type_reorder=[2 1 3 4];

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

for ii_R2type=1:length(R2type)

    these_R1s=[];

    eval(['these_R2s=' R2type{ii_R2type} ';'])

    all_R1s.R2type(ii_R2type).R1=these_R2s;

    %plot bar
    switch ii_R2type
        case 1
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
        case 2
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
        case 3
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
        case 4
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
        case 5
            bar(bar_offset,mean(these_R2s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
    end

    %Violin plot
    [mean_out, CIout, all_R1s.R2type(ii_R2type).violin_x]=drgViolinPoint(these_R2s...
        ,edges,bar_offset,rand_offset,'k','k',4);
    bar_offset=bar_offset+1;

    glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=these_R2s;
    glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R2s))=iiR2type_reorder(ii_R2type)*ones(1,length(these_R2s));
    glm_r1_ii=glm_r1_ii+length(these_R2s);

    id_r1_ii=id_r1_ii+1;
    input_r1_data(id_r1_ii).data=these_R2s;
    input_r1_data(id_r1_ii).description=[R2type_label{ii_R2type}];

end
bar_offset=bar_offset+1;

%Plot lines between points
for ii=1:length(these_R2s)
    these_R2s=[];
    these_x=[];
    for ii_R2type=1:length(R2type)
        these_R2s=[these_R2s all_R1s.R2type(ii_R2type).R1(ii)];
        these_x=[these_x all_R1s.R2type(ii_R2type).violin_x(ii)];
    end
    plot(these_x,these_R2s,'-','Color',[0.7 0.7 0.7])
end


xticks([0 1 2 3])
xticklabels({'within','hit','miss','shuffled'})


title(['Variance for prediction of odor hit, miss, shuffled, after last turn'])
ylabel('Variance')
ylim([0 12])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' after turn R1\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' after turn R1\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' after turn R1\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' after turn R1\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

%Now calculate the hits et al for XY
these_groups=[1 5];
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
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)

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
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
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

        x_predicted_sh_conv=zeros(size(x_predicted_sh,1),size(x_predicted_sh,2))
        y_predicted_sh_conv=zeros(size(y_predicted_sh,1),size(y_predicted_sh,2))

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

        ii_for_corr=ii_for_corr+1;
        these_XYfileNos(ii_for_corr)=fileNo;

        if ~isempty(x_all_trials)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(x_all_trials,x_decod_all_trials);
                R1_x_all_trials(ii_for_corr)=R1(1,2);
                P_rho_x_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(x_all_trials,x_decod_all_trials,'Type','Spearman');
                R1_x_all_trials(ii_for_corr)=R1;
                P_rho_x_all_trials(ii_for_corr)=this_P_rho_all_trials;
            end
        else
            R1_x_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_trials)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(y_all_trials,y_decod_all_trials);
                R1_y_all_trials(ii_for_corr)=R1(1,2);
                P_rho_y_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(y_all_trials,y_decod_all_trials,'Type','Spearman');
                R1_y_all_trials(ii_for_corr)=R1;
                P_rho_y_all_trials(ii_for_corr)=this_P_rho_all_trials;
            end
        else
            R1_y_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials)
              
            R1XY=corr2([x_all_trials y_all_trials],[x_decod_all_trials y_decod_all_trials]);
            R1_XY_all_trials(ii_for_corr)=R1XY;
       
        else
            R1_XY_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_hits)
            if is_pearson==1
                R1=corrcoef(x_all_hits,x_decod_all_hits);
                R1_x_hits(ii_for_corr)=R1(1,2);
            else
                R1=corr(x_all_hits,x_decod_all_hits,'Type','Spearman');
                R1_x_hits(ii_for_corr)=R1;
            end
        else
            R1_x_hits(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_hits)
            if is_pearson==1
                R1=corrcoef(y_all_hits,y_decod_all_hits);
                R1_y_hits(ii_for_corr)=R1(1,2);
            else
                R1=corr(y_all_hits,y_decod_all_hits,'Type','Spearman');
                R1_y_hits(ii_for_corr)=R1;
            end
        else
            R1_y_hits(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials)
            R1XY=corr2([x_all_hits y_all_hits],[x_decod_all_hits y_decod_all_hits]);
            R1_XY_hits(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_hits(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_miss)
            if is_pearson==1
                R1=corrcoef(x_all_miss,x_decod_all_miss);
                R1_x_miss(ii_for_corr)=R1(1,2);
            else
                R1=corr(x_all_miss,x_decod_all_miss,'Type','Spearman');
                R1_x_miss(ii_for_corr)=R1;
            end
        else
            R1_x_miss(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_miss)
            if is_pearson==1
                R1=corrcoef(y_all_miss,y_decod_all_miss);
                R1_y_miss(ii_for_corr)=R1(1,2);
            else
                R1=corr(y_all_miss,y_decod_all_miss,'Type','Spearman');
                R1_y_miss(ii_for_corr)=R1;
            end
        else
            R1_y_miss(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials)
            R1XY=corr2([x_all_miss y_all_miss],[x_decod_all_miss y_decod_all_miss]);
            R1_XY_miss(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_miss(ii_for_corr)=NaN;
        end


        if ~isempty(x_all_trials)
            these_R1s=[];
            for ii_sh=1:size(x_predicted_sh,2)
                if is_pearson==1
                    [R1,this_P_rho_all_trials]=corrcoef(x_all_trials,x_decod_all_trials_sh.ii_sh(ii_sh).xdec);
                    these_R1s=[these_R1s R1(1,2)];
                else
                    [R1,this_P_rho_all_trials]=corr(x_all_trials,x_decod_all_trials_sh.ii_sh(ii_sh).xdec,'Type','Spearman');
                    these_R1s=[these_R1s R1];
                end
            end
            R1_x_all_trials_sh(ii_for_corr)=mean(these_R1s);
        else
            R1_x_all_trials_sh(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_trials)
            these_R1s=[];
            for ii_sh=1:size(y_predicted_sh,2)
                if is_pearson==1
                    [R1,this_P_rho_all_trials]=corrcoef(y_all_trials,y_decod_all_trials_sh.ii_sh(ii_sh).ydec);
                    these_R1s=[these_R1s R1(1,2)];
                else
                    [R1,this_P_rho_all_trials]=corr(y_all_trials,y_decod_all_trials_sh.ii_sh(ii_sh).ydec,'Type','Spearman');
                    these_R1s=[these_R1s R1];
                end
            end
            R1_y_all_trials_sh(ii_for_corr)=mean(these_R1s);
        else
            R1_y_all_trials_sh(ii_for_corr)=NaN;
        end



        if ~isempty(x_all_trials)
            these_R1XYs=[];
            for ii_sh=1:size(y_predicted_sh,2)
                these_R1XYs=[these_R1XYs corr2([x_all_trials y_all_trials],[x_decod_all_trials_sh.ii_sh(ii_sh).xdec y_decod_all_trials_sh.ii_sh(ii_sh).ydec])];
            end
            R1_XY_all_trials_sh(ii_for_corr)=mean(these_R1XYs);
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_all_trials_sh(ii_for_corr)=NaN;
        end

    end
end


%Do the R1_XY bar graph for 1 and 2 cm all trials vs hit, miss
ii_run=1;
R1type=[];
R1type{1}='R1_XY_all_trials';
R1type{2}='R1_XY_hits';
R1type{3}='R1_XY_miss';
R1type{4}='R1_XY_all_trials_sh';

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


title(['R1 XY for prediction of odor hit, miss, shuffled'])
ylabel('R1')
ylim([-0.2 1])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1_XY hit, miss, within, shuffled\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1_XY hit, miss, within, shuffled\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1_XY hit, miss, within, shuffled\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1_XY hit, miss, within, shuffled\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Do the R1_x and R1_y bar graph for 1 and 2 cm all trials vs hit, miss
ii_run=1;
R1t=[];

R1t.xy(1).R1type{1}='R1_x_all_trials';
R1t.xy(1).R1type{2}='R1_x_hits';
R1t.xy(1).R1type{3}='R1_x_miss';
R1t.xy(1).R1type{4}='R1_x_all_trials_sh';

R1t.xy(2).R1type{1}='R1_y_all_trials';
R1t.xy(2).R1type{2}='R1_y_hits';
R1t.xy(2).R1type{3}='R1_y_miss';
R1t.xy(2).R1type{4}='R1_y_all_trials_sh';

R1type_label=[];
R1type_label{1}='within';
R1type_label{2}='hits';
R1type_label{3}='miss';
R1type_label{4}='shuffled';

xy_label{1}='x';
xy_label{2}='y';

iiR1type_reorder=[2 3 1 4];

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

for ii_xy=1:2
    all_R1s=[];
    for ii_R1type=1:length(R1type)

        these_R1s=[];

        eval(['these_R1s=' R1t.xy(ii_xy).R1type{ii_R1type} ';'])

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
        glm_r1.xy(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_xy*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type} ' ' xy_label{ii_xy}];

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

end
xticks([0 1 2 3 5 6 7 8])
xticklabels({'x within','x hit','x miss','x shuffled','y within','y hit','y miss','y shuffled'})


title(['R1 x and y for prediction of odor hit, miss, shuffled'])
ylabel('R1')
ylim([-0.2 1])
xlim([-1 9])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1 x and y, hit, miss, within, shuffled\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1 x an dy, hit, miss, within, shuffled\n']);
fprintf(fileID, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.trial_type',glm_r1.xy',...
    'VariableNames',{'R1','trial_type','xy'});
mdl = fitglm(tbl,'R1~trial_type+xy+trial_type*xy'...
    ,'CategoricalVars',[2])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 x and y, hit, miss, within, shuffled\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 x and y, hit, miss, within, shuffled\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

%Now re-calculate R1 within before and after last turn
ii_run=1;
these_groups=[1 5];
ii_for_corr=0;


R1_XY_all_trials_bt=[];
R1_XY_hits_bt=[];
R1_XY_miss_bt=[];

R1_XY_all_trials_at=[];
R1_XY_hits_at=[];
R1_XY_miss_at=[];

R1_x_all_trials_bt=[];
R1_x_hits_bt=[];
R1_x_miss_bt=[];

R1_x_all_trials_at=[];
R1_x_hits_at=[];
R1_x_miss_at=[];

R1_y_all_trials_bt=[];
R1_y_hits_bt=[];
R1_y_miss_bt=[];

R1_y_all_trials_at=[];
R1_y_hits_at=[];
R1_y_miss_at=[];

R1_y_all_trials_at_sh=[];
R1_x_all_trials_at_sh=[];
R1_XY_all_trials_at_sh=[];



these_XYfileNos=[];


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
        
        x_all_trials_bt=[];
        x_decod_all_trials_bt=[];
        x_decod_all_trials_sh_bt=[];
        x_all_hits_bt=[];
        x_decod_all_hits_bt=[];
        x_all_miss_bt=[];
        x_decod_all_miss_bt=[];
        P_rho_x_all_trials_bt=[];

        y_all_trials_bt=[];
        y_decod_all_trials_bt=[];
        y_decod_all_trials_sh_bt=[];
        y_all_hits_bt=[];
        y_decod_all_hits_bt=[];
        y_all_miss_bt=[];
        y_decod_all_miss_bt=[];
        P_rho_x_all_trials_bt=[];

          x_all_trials_at=[];
        x_decod_all_trials_at=[];
        x_decod_all_trials_sh_at=[];
        x_all_hits_at=[];
        x_decod_all_hits_at=[];
        x_all_miss_at=[];
        x_decod_all_miss_at=[];
        P_rho_x_all_trials_at=[];

        y_all_trials_at=[];
        y_decod_all_trials_at=[];
        y_decod_all_trials_sh_at=[];
        y_all_hits_at=[];
        y_decod_all_hits_at=[];
        y_all_miss_at=[];
        y_decod_all_miss_at=[];
        P_rho_x_all_trials_at=[];

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
        arena_file=handles_XY.arena_file{fileNo};
        %load the ouptut file
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        x_predicted=handles_out.x_predicted;
        y_predicted=handles_out.y_predicted;
        x_predicted_sh=handles_out.x_predicted_sh(:,1);
        y_predicted_sh=handles_out.y_predicted_sh(:,1); %Note: all sh are identical
        XYtest=handles_out.XYtest;

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

        x_predicted_sh_conv=conv(x_predicted_sh,conv_win_gauss,'same');
        y_predicted_sh_conv=conv(y_predicted_sh,conv_win_gauss,'same');

        minY1=min(XYtest(:,1));
        x_predicted_sh_conv(x_predicted_sh_conv<minY1)=minY1;
        maxY1=max(XYtest(:,1));
        x_predicted_sh_conv(x_predicted_sh_conv>maxY1)=maxY1;

        minY2=min(XYtest(:,2));
        y_predicted_sh_conv(y_predicted_sh_conv<minY2)=minY2;
        maxY2=max(XYtest(:,2));
        y_predicted_sh_conv(y_predicted_sh_conv>maxY2)=maxY2;

        for trNo=1:trials.odor_trNo



            %Find the last turn
            ii_turn=find(handles_out_angle.angles.trial(trNo).delta_x>100,1,'last');
            if ~isempty(ii_turn)
                ii_predictedstart_bt=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                ii_predictedend_bt=ii_predictedstart_bt+handles_out_angle.angles.trial(trNo).ii_turns(ii_turn)-1;
                if ii_predictedend_bt>size(XYtest,1)
                    ii_predictedend_bt=size(XYtest,1);
                end

 
                ii_predictedstart_at=ii_predictedstart_bt+handles_out_angle.angles.trial(trNo).ii_turns(ii_turn);
                ii_predictedend_at=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if ii_predictedend_at>size(XYtest,1)
                    ii_predictedend_at=size(XYtest,1);
                end

                x_all_trials_bt=[x_all_trials_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,1)];
                x_decod_all_trials_bt=[x_decod_all_trials_bt; x_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                y_all_trials_bt=[y_all_trials_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,2)];
                y_decod_all_trials_bt=[y_decod_all_trials_bt; y_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                x_decod_all_trials_sh_bt=[x_decod_all_trials_sh_bt; x_predicted_sh_conv(ii_predictedstart_bt:ii_predictedend_bt)];
                y_decod_all_trials_sh_bt=[y_decod_all_trials_sh_bt; y_predicted_sh_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                x_all_trials_at=[x_all_trials_at; XYtest(ii_predictedstart_at:ii_predictedend_at,1)];
                x_decod_all_trials_at=[x_decod_all_trials_at; x_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                y_all_trials_at=[y_all_trials_at; XYtest(ii_predictedstart_at:ii_predictedend_at,2)];
                y_decod_all_trials_at=[y_decod_all_trials_at; y_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                x_decod_all_trials_sh_at=[x_decod_all_trials_sh_at; x_predicted_sh_conv(ii_predictedstart_at:ii_predictedend_at)];
                y_decod_all_trials_sh_at=[y_decod_all_trials_sh_at; y_predicted_sh_conv(ii_predictedstart_at:ii_predictedend_at)];


                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
                        x_all_hits_bt=[x_all_hits_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,1)];
                        x_decod_all_hits_bt=[x_decod_all_hits_bt; x_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        y_all_hits_bt=[y_all_hits_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,2)];
                        y_decod_all_hits_bt=[y_decod_all_hits_bt; y_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        x_all_hits_at=[x_all_hits_at; XYtest(ii_predictedstart_at:ii_predictedend_at,1)];
                        x_decod_all_hits_at=[x_decod_all_hits_at; x_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                        y_all_hits_at=[y_all_hits_at; XYtest(ii_predictedstart_at:ii_predictedend_at,2)];
                        y_decod_all_hits_at=[y_decod_all_hits_at; y_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                    case 2
                        %Lane 1 miss orange
                        x_all_miss_bt=[x_all_miss_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,1)];
                        x_decod_all_miss_bt=[x_decod_all_miss_bt; x_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        y_all_miss_bt=[y_all_miss_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,2)];
                        y_decod_all_miss_bt=[y_decod_all_miss_bt; y_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        x_all_miss_at=[x_all_miss_at; XYtest(ii_predictedstart_at:ii_predictedend_at,1)];
                        x_decod_all_miss_at=[x_decod_all_miss_at; x_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                        y_all_miss_at=[y_all_miss_at; XYtest(ii_predictedstart_at:ii_predictedend_at,2)];
                        y_decod_all_miss_at=[y_decod_all_miss_at; y_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                    case 3
                        %Lane 4 hit blue
                        x_all_hits_bt=[x_all_hits_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,1)];
                        x_decod_all_hits_bt=[x_decod_all_hits_bt; x_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        y_all_hits_bt=[y_all_hits_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,2)];
                        y_decod_all_hits_bt=[y_decod_all_hits_bt; y_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        x_all_hits_at=[x_all_hits_at; XYtest(ii_predictedstart_at:ii_predictedend_at,1)];
                        x_decod_all_hits_at=[x_decod_all_hits_at; x_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                        y_all_hits_at=[y_all_hits_at; XYtest(ii_predictedstart_at:ii_predictedend_at,2)];
                        y_decod_all_hits_at=[y_decod_all_hits_at; y_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];
                    case 4
                        %Lane 4 miss sky blue
                        x_all_miss_bt=[x_all_miss_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,1)];
                        x_decod_all_miss_bt=[x_decod_all_miss_bt; x_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        y_all_miss_bt=[y_all_miss_bt; XYtest(ii_predictedstart_bt:ii_predictedend_bt,2)];
                        y_decod_all_miss_bt=[y_decod_all_miss_bt; y_predicted_conv(ii_predictedstart_bt:ii_predictedend_bt)];

                        x_all_miss_at=[x_all_miss_at; XYtest(ii_predictedstart_at:ii_predictedend_at,1)];
                        x_decod_all_miss_at=[x_decod_all_miss_at; x_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];

                        y_all_miss_at=[y_all_miss_at; XYtest(ii_predictedstart_at:ii_predictedend_at,2)];
                        y_decod_all_miss_at=[y_decod_all_miss_at; y_predicted_conv(ii_predictedstart_at:ii_predictedend_at)];
                end
            end

        end

        ii_for_corr=ii_for_corr+1;

        if ~isempty(x_all_trials_bt)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(x_all_trials_bt,x_decod_all_trials_bt);
                R1_x_all_trials_bt(ii_for_corr)=R1(1,2);
                P_rho_x_all_trials_bt(ii_for_corr)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(x_all_trials_bt,x_decod_all_trials_bt,'Type','Spearman');
                R1_x_all_trials_bt(ii_for_corr)=R1;
                P_rho_x_all_trials_bt(ii_for_corr)=this_P_rho_all_trials;
            end
        else
            R1_x_all_trials_bt(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_trials_bt)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(y_all_trials_bt,y_decod_all_trials_bt);
                R1_y_all_trials_bt(ii_for_corr)=R1(1,2);
                P_rho_y_all_trials_bt(ii_for_corr)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(y_all_trials_bt,y_decod_all_trials_bt,'Type','Spearman');
                R1_y_all_trials_bt(ii_for_corr)=R1;
                P_rho_y_all_trials_bt(ii_for_corr)=this_P_rho_all_trials;
            end
        else
            R1_y_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_bt)

            R1XY=corr2([x_all_trials_bt y_all_trials_bt],[x_decod_all_trials_bt y_decod_all_trials_bt]);
            R1_XY_all_trials_bt(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_all_trials_bt(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_hits_bt)
            if is_pearson==1
                R1=corrcoef(x_all_hits_bt,x_decod_all_hits_bt);
                R1_x_hits_bt(ii_for_corr)=R1(1,2);
            else
                R1=corr(x_all_hits_bt,x_decod_all_hits_bt,'Type','Spearman');
                R1_x_hits_bt(ii_for_corr)=R1;
            end
        else
            R1_x_hits_bt(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_hits_bt)
            if is_pearson==1
                R1=corrcoef(y_all_hits_bt,y_decod_all_hits_bt);
                R1_y_hits_bt(ii_for_corr)=R1(1,2);
            else
                R1=corr(y_all_hits_bt,y_decod_all_hits_bt,'Type','Spearman');
                R1_y_hits_bt(ii_for_corr)=R1;
            end
        else
            R1_y_hits_bt(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_bt)
            R1XY=corr2([x_all_hits_bt y_all_hits_bt],[x_decod_all_hits_bt y_decod_all_hits_bt]);
            R1_XY_hits_bt(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_hits_bt(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_miss_bt)
            if is_pearson==1
                R1=corrcoef(x_all_miss_bt,x_decod_all_miss_bt);
                R1_x_miss_bt(ii_for_corr)=R1(1,2);
            else
                R1=corr(x_all_miss_bt,x_decod_all_miss_bt,'Type','Spearman');
                R1_x_miss_bt(ii_for_corr)=R1;
            end
        else
            R1_x_miss_bt(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_miss_bt)
            if is_pearson==1
                R1=corrcoef(y_all_miss_bt,y_decod_all_miss_bt);
                R1_y_miss_bt(ii_for_corr)=R1(1,2);
            else
                R1=corr(y_all_miss_bt,y_decod_all_miss_bt,'Type','Spearman');
                R1_y_miss_bt(ii_for_corr)=R1;
            end
        else
            R1_y_miss_bt(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_bt)
            R1XY=corr2([x_all_miss_bt y_all_miss_bt],[x_decod_all_miss_bt y_decod_all_miss_bt]);
            R1_XY_miss_bt(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_miss_bt(ii_for_corr)=NaN;
        end



        if ~isempty(x_all_trials_at)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(x_all_trials_at,x_decod_all_trials_at);
                R1_x_all_trials_at(ii_for_corr)=R1(1,2);
                P_rho_x_all_trials_at(ii_for_corr)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(x_all_trials_at,x_decod_all_trials_at,'Type','Spearman');
                R1_x_all_trials_at(ii_for_corr)=R1;
                P_rho_x_all_trials_at(ii_for_corr)=this_P_rho_all_trials;
            end
        else
            R1_x_all_trials_at(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_trials_at)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(y_all_trials_at,y_decod_all_trials_at);
                R1_y_all_trials_at(ii_for_corr)=R1(1,2);
                P_rho_y_all_trials_at(ii_for_corr)=this_P_rho_all_trials(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(y_all_trials_at,y_decod_all_trials_at,'Type','Spearman');
                R1_y_all_trials_at(ii_for_corr)=R1;
                P_rho_y_all_trials_at(ii_for_corr)=this_P_rho_all_trials;
            end
        else
            R1_y_all_trials(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_at)
            R1XY=corr2([x_all_trials_at y_all_trials_at],[x_decod_all_trials_at y_decod_all_trials_at]);
            R1_XY_all_trials_at(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_all_trials_at(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_hits_at)
            if is_pearson==1
                R1=corrcoef(x_all_hits_at,x_decod_all_hits_at);
                R1_x_hits_at(ii_for_corr)=R1(1,2);
            else
                R1=corr(x_all_hits_at,x_decod_all_hits_at,'Type','Spearman');
                R1_x_hits_at(ii_for_corr)=R1;
            end
        else
            R1_x_hits_at(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_hits_at)
            if is_pearson==1
                R1=corrcoef(y_all_hits_at,y_decod_all_hits_at);
                R1_y_hits_at(ii_for_corr)=R1(1,2);
            else
                R1=corr(y_all_hits_at,y_decod_all_hits_at,'Type','Spearman');
                R1_y_hits_at(ii_for_corr)=R1;
            end
        else
            R1_y_hits_at(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_at)
            R1XY=corr2([x_all_hits_at y_all_hits_at],[x_decod_all_hits_at y_decod_all_hits_at]);
            R1_XY_hits_at(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_hits_at(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_miss_at)
            if is_pearson==1
                R1=corrcoef(x_all_miss_at,x_decod_all_miss_at);
                R1_x_miss_at(ii_for_corr)=R1(1,2);
            else
                R1=corr(x_all_miss_at,x_decod_all_miss_at,'Type','Spearman');
                R1_x_miss_at(ii_for_corr)=R1;
            end
        else
            R1_x_miss_at(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_miss_at)
            if is_pearson==1
                R1=corrcoef(y_all_miss_at,y_decod_all_miss_at);
                R1_y_miss_at(ii_for_corr)=R1(1,2);
            else
                R1=corr(y_all_miss_at,y_decod_all_miss_at,'Type','Spearman');
                R1_y_miss_at(ii_for_corr)=R1;
            end
        else
            R1_y_miss_at(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_at)
            R1XY=corr2([x_all_miss_at y_all_miss_at],[x_decod_all_miss_at y_decod_all_miss_at]);
            R1_XY_miss_at(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_miss_at(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_at)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(x_all_trials_at,x_decod_all_trials_sh_at);
                R1_x_all_trials_at_sh(ii_for_corr)=R1(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(x_all_trials_at,x_decod_all_trials_sh_at,'Type','Spearman');
                R1_x_all_trials_at_sh(ii_for_corr)=R1;
            end
            % P_rho_x_all_trials_bt(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_x_all_trials_at_sh(ii_for_corr)=NaN;
        end

        if ~isempty(y_all_trials_at)
            if is_pearson==1
                [R1,this_P_rho_all_trials]=corrcoef(y_all_trials_at,y_decod_all_trials_sh_at);
                R1_y_all_trials_at_sh(ii_for_corr)=R1(1,2);
            else
                [R1,this_P_rho_all_trials]=corr(y_all_trials_at,y_decod_all_trials_sh_at,'Type','Spearman');
                R1_y_all_trials_at_sh(ii_for_corr)=R1;
            end
            % P_rho_y_all_trials_bt(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_y_all_trials_at_sh(ii_for_corr)=NaN;
        end

        if ~isempty(x_all_trials_at)
            R1XY=corr2([x_all_trials_at y_all_trials_at],[x_decod_all_trials_sh_at y_decod_all_trials_sh_at]);
            R1_XY_all_trials_at_sh(ii_for_corr)=R1XY;
            % P_rho_XY_all_trials(ii_for_corr)=this_P_rho_all_trials(1,2);
        else
            R1_XY_all_trials_at_sh(ii_for_corr)=NaN;
        end



    end
end



%Let's plot R1_XY within vs, within before and after last turn
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
R1type{1}='R1_XY_all_trials';
R1type{2}='R1_XY_all_trials_bt';
R1type{3}='R1_XY_all_trials_at';


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


title(['R1 XY for prediction of odor concentration before and after last turn'])
ylabel('R1')
ylim([-0.2 1])
xlim([-1 3])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1_XY before after last turn\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1_XY before after last turn\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1_XY before after last turn\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1_XY before after last turn\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

 
%Do R1 bar graph for 1 and 2 cm all trials vs hit, miss, after last turn
ii_run=1;
R1type=[];
R1type{1}='R1_XY_all_trials_at';
R1type{2}='R1_XY_hits_at';
R1type{3}='R1_XY_miss_at';
R1type{4}='R1_XY_all_trials_at_sh';

R1type_label=[];
R1type_label{1}='within';
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


title(['R1 XY for prediction of odor hit, miss, shuffled, after last turn'])
ylabel('R1')
ylim([-0.2 1])
xlim([-1 4])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1 XY after turn\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1 XY after turn\n']);
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
fprintf(1, ['\n\nRanksum or t-test p values for R1 XY Fig 19 after turn\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for R1 XY Fig 19 after turn\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


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

R1XYtype{1}='R1_XY_hits_at';
R1XYtype{2}='R1_XY_miss_at';

R1type{1}='R1_hits_at';
R1type{2}='R1_miss_at';

for ii_R1type=1:2
    %R1 for conc all trials

    R1_conc=[];
    eval(['R1_conc=' R1type{ii_R1type} ';'])

    R1_XY=[];
    eval(['R1_XY=' R1XYtype{ii_R1type} ';'])

       switch ii_R1type
        case 1
        
            plot(R1_XY,R1_conc,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',8)
        case 2
            plot(R1_XY,R1_conc,'o','MarkerEdgeColor','none','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',8) 
    end


end

plot([-0.1 1],[-0.1 1],'-k')
text(0.5,0.5,'Hit','Color','k')
text(0.5,0.4,'Miss','Color',[0.7 0.7 0.7])
title('Within trial R1 after turn')
ylim([-0.1 1])
xlim([-0.1 1])
xlabel('R1 XY')
ylabel('R1 odor')



%Let's plot R1 x and y within vs, within before and after last turn
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

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];


%Plot the different R1s
these_groups=[1 5]; %1,2 cm all tirals

R1t=[];
R1t.xy(1).R1type{1}='R1_x_all_trials';
R1t.xy(1).R1type{2}='R1_x_all_trials_bt';
R1t.xy(1).R1type{3}='R1_x_all_trials_at';

R1t.xy(2).R1type{1}='R1_y_all_trials';
R1t.xy(2).R1type{2}='R1_y_all_trials_bt';
R1t.xy(2).R1type{3}='R1_y_all_trials_at';

R1type_label=[];
R1type_label{1}='within';
R1type_label{2}='before';
R1type_label{3}='after';

for ii_xy=1:2
    allR1s=[];
    for ii_R1type=1:length(R1t.xy(ii_xy).R1type)

        these_R1s=[];
        eval(['these_R1s=' R1t.xy(ii_xy).R1type{ii_R1type} ';'])

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
        glm_r1.xy(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_xy*ones(1,length(these_R1s));
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
        for ii_R1type=1:length(R1t.xy(ii_xy).R1type)
            if isnan(all_R1s.R1type(ii_R1type).R1(ii))
                no_nans=0;
            end
        end

        if no_nans==1
            for ii_R1type=1:length(R1t.xy(ii_xy).R1type)
                these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)];
                these_violin_x=[these_violin_x all_R1s.R1type(ii_R1type).violin_x(ii)];
            end

            plot(these_violin_x,these_R1s,'-','Color',[0.7 0.7 0.7])
        end
    end
    bar_offset=bar_offset+1;
end

% x_pos=3;
% text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 2  4 5 6])
xticklabels({'within x','before x','after x','within y','before y','after y'})


title(['R1 x,y for prediction of odor concentration before and after last turn'])
ylabel('R1')
ylim([-0.2 1])
xlim([-1 7])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1 x,y before after last turn\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1 x,y before after last turn\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.trial_type',glm_r1.xy',...
    'VariableNames',{'R1','turn','xy'});
mdl = fitglm(tbl,'R1~turn+xy+xy*turn'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 x,y before after last turn\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for Fig ' num2str(figureNo) ' R1 x,y before after last turn\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

  
%Do R1 x,y bar graph for 1 and 2 cm all trials vs hit, miss, after last turn
ii_run=1;
R1t=[];
R1t.xy(1).R1type{1}='R1_x_all_trials_at';
R1t.xy(1).R1type{2}='R1_x_hits_at';
R1t.xy(1).R1type{3}='R1_x_miss_at';
R1t.xy(1).R1type{4}='R1_x_all_trials_at_sh';

R1t.xy(2).R1type{1}='R1_y_all_trials_at';
R1t.xy(2).R1type{2}='R1_y_hits_at';
R1t.xy(2).R1type{3}='R1_y_miss_at';
R1t.xy(2).R1type{4}='R1_y_all_trials_at_sh';

R1type_label=[];
R1type_label{1}='within';
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
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

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
for ii_xy=1:2
    for ii_R1type=1:length(R1t.xy(ii_xy).R1type)

        these_R1s=[];

        eval(['these_R1s=' R1t.xy(ii_xy).R1type{ii_R1type} ';'])

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
        glm_r1.xy(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_xy*ones(1,length(these_R1s));
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

end
xticks([0 1 2 3 5 6 7 8])
xticklabels({'within x','hit x','miss x','shuffled x','within y','hit y','miss y','shuffled y'})


title(['R1 x,y for prediction of odor hit, miss, shuffled, after last turn'])
ylabel('R1')
ylim([-0.2 1])
xlim([-1 9])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' R1 x,y after turn\n'])
fprintf(1, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' R1 x,y after turn\n']);
fprintf(fileID, ['\n1 Hit, 2 Miss, 3 All trials, 4 Shuffled\n'])
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.trial_type',glm_r1.xy',...
    'VariableNames',{'R1','trial_type','xy'});
mdl = fitglm(tbl,'R1~trial_type+xy+trial_type*xy'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for R1 x,y Fig 22 after turn\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for R1 x,y Fig 22 after turn\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

pffft=1;


fclose(fileID);

pffft=1;