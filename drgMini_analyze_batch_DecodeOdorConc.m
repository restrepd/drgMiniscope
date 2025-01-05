%drgMini_analyze_batch_DecodeOdorConc
close all
clear all

is_sphgpu=1;

if is_sphgpu==1
    addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
    addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
    addpath('/home/restrepd/Documents/MATLAB/drgMaster')
    addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end

fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');

% handles.algo=[1 2 3 4 5 3 3 3 3 3 3];
%1=fitrnet
%2=fitrgp
%3=fitrtree
%4=fitglm
%5=fitrsvm

% handles.weber_fechner=[1 1 1 1 1 0 0 0 1 1 1];
%0 is Stevens Law, R proportional to C^alpha
%1 is Weber-Flechner law R proportional to log(C)
%See Copelli et al DOI: 10.1103/PhysRevE.65.060901

% handles.alpha=[1 1 1 1 1 1 0.5 2 1 1 1];

%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4, 2 cm from floor
%Group 2 is rewarded, with odor lane 4, no odor in lane 1
%Group 3 is rewarded, with odor lane 1, no odor in lane 4
%Group 4 is rewarded, with no odor in lane 1 and lane 4
%Group 5 is rewarded, with ISO1 in both lane 1 and lane 4, 1 cm from floor

% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

group_order=[1 4 5 2 3];



figureNo=0;

display_each_figure=1;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgOdorConcChoices_*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgMini_batch_DecodeOdorConc run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
handles.no_files=length(handles.dFF_file);


first_file=handles.first_file;
first_run=handles.first_run;

%Check whether all output files exist
all_files_present=1;


%Check each output file
for fileNo=1:handles.no_files
    for ii_run=1:length(handles.weber_fechner)
        %does the file exist?
        arena_file=handles.arena_file{fileNo};
        if ~(exist([handles.save_path arena_file(1:end-4) handles.save_tag{ii_run} '.mat'],'file')==2)
            all_files_present=0;
            fprintf(1, ['File number %d run %d does not exist\n'],fileNo, ii_run);
        end
    end
end

%Analyze batch output
if all_files_present==1

    %Compare groups with fittree, weber frechner and only current data
    %point
    ii_run=3;

    groups=unique(handles.group);

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
            for fileNo=1:handles.no_files
                arena_file=handles.arena_file{fileNo};
                if handles.group(fileNo)==this_group
                    %load the ouptut file
                    load([handles.save_path arena_file(1:end-4) handles.save_tag{ii_run} '.mat'])
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


    title(['R1, correlation between decoded odor concentraiton and odor plume mean'])
    ylabel('R1')
    ylim([-0.2 0.65])
    xlim([-1 17])

    %Do the R1 bar graph for 2 cm all trials vs hit, miss
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
    this_group=1; %2 cm all tirals
    for ii_R1type=1:length(R1type)

        these_R1s=[];
        for fileNo=1:handles.no_files
            arena_file=handles.arena_file{fileNo};
            if handles.group(fileNo)==this_group
                %load the ouptut file
                load([handles.save_path arena_file(1:end-4) handles.save_tag{ii_run} '.mat'])
                eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
                % these_R1s=[these_R1s handles_out.R1.all_trials];
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
        [mean_out, CIout]=drgViolinPoint(these_R1s...
            ,edges,bar_offset,rand_offset,'k','k',4);
        bar_offset=bar_offset+1;

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
        glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_R1type*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

    end
    bar_offset=bar_offset+1;

    %Plot lines between points
    for ii=1:length(these_R1s)
        these_R1s=[];
        for ii_R1type=1:length(R1type)
            these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)]
        end
        plot([0:length(R1type)-1],these_R1s,'-','Color',[0.7 0.7 0.7])
    end

    % x_pos=3;
    % text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
    % text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
    % text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
    % text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

    xticks([0 1 2 3])
    xticklabels({'within','hit','miss','shuffled'})


    title(['R1 for predcition of odor concentration'])
    ylabel('R1')
    ylim([-0.2 0.65])
    xlim([-1 4])

    %Perform the glm  for errors
    fprintf(1, ['\nglm for trial type\n'])
    fprintf(fileID, ['\nglm for trial type\n']);
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


    %Do the R1 bar graph for 2 cm Weber-Fechner vs. Stevens
    R1type=[];
    R1type{1}='all_trials';
    R1type{2}='all_hits';
    R1type{3}='all_miss';
    R1type{4}='all_trials_sh';

    R1type_label=[];
    R1type_label{1}='wf';
    R1type_label{2}='s0p5';
    R1type_label{3}='s1';
    R1type_label{4}='s2';

    ii_run_order=[3 7 6 8]; %All 2 cm

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

    ii_R1type=1;
    this_group=1;
    all_R1s=[];
    for ii_run=1:length(ii_run_order)

        this_ii_run=ii_run_order(ii_run);
        these_R1s=[];
        for fileNo=1:handles.no_files
            arena_file=handles.arena_file{fileNo};
            if handles.group(fileNo)==this_group
                %load the ouptut file
                load([handles.save_path arena_file(1:end-4) handles.save_tag{this_ii_run} '.mat'])
                eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
                % these_R1s=[these_R1s handles_out.R1.all_trials];
            end

        end
        all_R1s.run(ii_run).R1=these_R1s;
        %plot bar
        switch ii_run
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

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
        glm_r1.group(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_run*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

    end
    bar_offset=bar_offset+1;

    %Plot lines between points
    for ii=1:length(these_R1s)
        these_R1s=[];
        for ii_run=1:length(ii_run_order)
            these_R1s=[these_R1s all_R1s.run(ii_run).R1(ii)]
        end
        plot([0:length(ii_run_order)-1],these_R1s,'-','Color',[0.7 0.7 0.7])
    end

    % x_pos=3;
    % text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
    % text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
    % text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
    % text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

    xticks([0 1 2 3])
    xticklabels({'wf','s0p5','s1','s2'})


    title(['R1 for predcition of odor concentration (1 cm)'])
    ylabel('R1')
    ylim([-0.2 0.65])
    xlim([-1 4])

    %Perform the glm  for errors
    fprintf(1, ['\nglm for wf vs s\n'])
    fprintf(fileID, ['\nglm for wf vs s\n']);
    %
    % fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
    % fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

    tbl = table(glm_r1.data',glm_r1.group',...
        'VariableNames',{'R1','wfs'});
    mdl = fitglm(tbl,'R1~wfs'...
        ,'CategoricalVars',[2])

    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);


    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for R1 vs weber frechner or stevens\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for R1 vs. frechner or stevens\n']);


    [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

    %Do the R1 bar graph for 2 cm vs. number of past points
    R1type=[];
    R1type{1}='all_trials';
    R1type{2}='all_hits';
    R1type{3}='all_miss';
    R1type{4}='all_trials_sh';

    R1type_label=[];
    R1type_label{1}='0';
    R1type_label{2}='1';
    R1type_label{3}='2';
    R1type_label{4}='4';

    ii_run_order=[3 9 10 11]; %All 2 cm

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

    ii_R1type=1;
    this_group=1;

    all_R1s=[];

    for ii_run=1:length(ii_run_order)

        this_ii_run=ii_run_order(ii_run);
        these_R1s=[];
        for fileNo=1:handles.no_files
            arena_file=handles.arena_file{fileNo};
            if handles.group(fileNo)==this_group
                %load the ouptut file
                load([handles.save_path arena_file(1:end-4) handles.save_tag{this_ii_run} '.mat'])
                eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
                % these_R1s=[these_R1s handles_out.R1.all_trials];
            end

        end

        all_R1s.run(ii_run).R1=these_R1s;
        %plot bar
        switch ii_run
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

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
        glm_r1.group(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_run*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

    end
    bar_offset=bar_offset+1;

    %Plot lines between points
    for ii=1:length(these_R1s)
        these_R1s=[];
        for ii_run=1:length(ii_run_order)
            these_R1s=[these_R1s all_R1s.run(ii_run).R1(ii)]
        end
        plot([0:length(ii_run_order)-1],these_R1s,'-','Color',[0.7 0.7 0.7])
    end

    % x_pos=3;
    % text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
    % text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
    % text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
    % text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

    xticks([0 1 2 3])
    xticklabels({'0','1','2','4'})
    xlabel('Number of previous points included')


    title(['R1 for predcition of odor concentration (1 cm)'])
    ylabel('R1')
    ylim([-0.2 0.65])
    xlim([-1 4])

    %Perform the glm  for errors
    fprintf(1, ['\nglm for wf vs s\n'])
    fprintf(fileID, ['\nglm for wf vs s\n']);
    %
    % fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
    % fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

    tbl = table(glm_r1.data',glm_r1.group',...
        'VariableNames',{'R1','wfs'});
    mdl = fitglm(tbl,'R1~wfs'...
        ,'CategoricalVars',[2])

    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);


    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for R1 vs weber frechner or stevens\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for R1 vs. frechner or stevens\n']);


    [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);

    %Do the R1 bar graph for 2 cm different algorithms
    R1type=[];
    R1type{1}='all_trials';
    R1type{2}='all_hits';
    R1type{3}='all_miss';
    R1type{4}='all_trials_sh';

    R1type_label=[];
    R1type_label{1}='0';
    R1type_label{2}='1';
    R1type_label{3}='2';
    R1type_label{4}='4';

    ii_run_order=[1 4 3 2 5]; %All 2 cm

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

    ii_R1type=1;
    this_group=1;
    all_R1s=[];
    for ii_run=1:length(ii_run_order)

        this_ii_run=ii_run_order(ii_run);
        these_R1s=[];
        for fileNo=1:handles.no_files
            arena_file=handles.arena_file{fileNo};
            if handles.group(fileNo)==this_group
                %load the ouptut file
                load([handles.save_path arena_file(1:end-4) handles.save_tag{this_ii_run} '.mat'])
                eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
                % these_R1s=[these_R1s handles_out.R1.all_trials];
            end

        end
        all_R1s.run(ii_run).R1=these_R1s;
        %plot bar
        switch ii_run
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

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
        glm_r1.group(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=ii_run*ones(1,length(these_R1s));
        glm_r1_ii=glm_r1_ii+length(these_R1s);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_R1s;
        input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];

    end
    bar_offset=bar_offset+1;

    %Plot lines between points
    for ii=1:length(these_R1s)
        these_R1s=[];
        for ii_run=1:length(ii_run_order)
            these_R1s=[these_R1s all_R1s.run(ii_run).R1(ii)]
        end
        plot([0:length(ii_run_order)-1],these_R1s,'-','Color',[0.7 0.7 0.7])
    end

    % x_pos=3;
    % text(x_pos,0.53,'ann','Color',[230/255 159/255 0/255])
    % text(x_pos,0.49,'gp','Color',[86/255 180/255 233/255])
    % text(x_pos,0.45,'tree','Color',[0/255 158/255 115/255])
    % text(x_pos,0.41,'glm','Color',[240/255 228/255 66/255])
    % text(x_pos,0.37,'svz','Color',[0/255 114/255 178/255])

    xticks([0 1 2 3 4])
    xticklabels({'ann','glm','tree','gp','svm'})

    title(['R1 for predcition of odor concentration (1 cm)'])
    ylabel('R1')
    ylim([-0.2 0.65])
    xlim([-1 5])

    %Perform the glm  for errors
    fprintf(1, ['\nglm for wf vs s\n'])
    fprintf(fileID, ['\nglm for wf vs s\n']);
    %
    % fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
    % fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

    tbl = table(glm_r1.data',glm_r1.group',...
        'VariableNames',{'R1','wfs'});
    mdl = fitglm(tbl,'R1~wfs'...
        ,'CategoricalVars',[2])

    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);


    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for R1 vs weber frechner or stevens\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for R1 vs. frechner or stevens\n']);


    [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);



end


fclose(fileID);

pffft=1;