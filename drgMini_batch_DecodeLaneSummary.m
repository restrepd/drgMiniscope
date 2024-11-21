%drgMini_batch_DecodeLaneSummary
close all
clear all
[choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniLanePredChoices_*.m'},'Select the .m file with all the choices for analysis');

fileNo=1; %This is the file that will be processed

fprintf(1, ['\ndrgMini_batch_dFFPrediction run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
handles.no_files=length(handles.dFF_file);


first_file=handles.first_file;

%Parallel batch processing for each file
all_files_present=1;

for filNum=first_file:handles.no_files


    %Make sure that all the files exist
    dFF_file=handles.dFF_file{filNum};

    this_path=handles.this_path{filNum};
    % else
    %     this_path=handles.this_path;
    % end

    if exist([this_path dFF_file])==0
        fprintf(1, ['Program will be terminated because file No %d, ' dFF_file ' does not exist\n'],filNum);
        all_files_present=0;
    end

    arena_file=handles.arena_file{filNum};

    if exist([this_path arena_file])==0
        fprintf(1, ['Program will be terminated because file No %d, ' arena_file ' does not exist\n'],filNum);
        all_files_present=0;
    end
end

%Now run the decoder for all files and all choices of algorithms
tic
if all_files_present==1

    %Calculate around start, around end and in the middle
    dt_window=2; %seconds
    dt_per_window=1;
    start_end_dt_span=3;

    all_sta_accuracy=[];
    all_sta_sh_accuracy=[];
    all_end_accuracy=[];
    all_end_sh_accuracy=[];

    ii_ROIst=0;
    ii_ROIend=0;

    fileNo_per_ROI=[];
    iiROI_per_ROI=[];



        %Load lane decode output file
        dFF_file=handles.dFF_file{fileNo};
        this_path=handles.this_path{fileNo};

        load([this_path dFF_file(1:end-4) handles.suffix])

        no_neurons=length(all_handles.run);
        trial_dt=all_handles.run(1).handles_out.trial_dt;

        %Calculate accuracy aligned to start
        for this_ROI=1:no_neurons

            %Now calculate the sliding window averages
            ii_ROIst=ii_ROIst+1;
            fileNo_per_ROI=[fileNo_per_ROI fileNo];
            iiROI_per_ROI=[iiROI_per_ROI this_ROI];
            time_span=all_handles.run(this_ROI).handles_out.start_alligned.time_span;
            accuracy_all_trials=mean(all_handles.run(this_ROI).handles_out.start_alligned.accuracy_all_trials,2);
            accuracy_sh_all_trials=mean(all_handles.run(this_ROI).handles_out.start_alligned.accuracy_sh_all_trials,2);
            ii_dt=0;

            %Aligned to start
            for ii=1:(2*start_end_dt_span)+1
                ii_dt=ii_dt+1;
                t=-start_end_dt_span+(ii-1)*dt_per_window;
                these_times=(time_span>t-dt_window/2)&(time_span<=t+dt_window/2);
                all_sta_accuracy(ii_ROIst,ii_dt)=mean(accuracy_all_trials(these_times));
                all_sta_sh_accuracy(ii_ROIst,ii_dt)=mean(accuracy_sh_all_trials(these_times));
            end

            %Middle
            ii_dt=ii_dt+1;
            t_middle=mean(trial_dt/2);
            these_times=(time_span>t_middle-dt_window/2)&(time_span<=t_middle+dt_window/2);
            all_sta_accuracy(ii_ROIst,ii_dt)=mean(accuracy_all_trials(these_times));
            all_sta_sh_accuracy(ii_ROIst,ii_dt)=mean(accuracy_sh_all_trials(these_times));

            %Aligned to end
            for ii=1:(2*start_end_dt_span)+1
                ii_dt=ii_dt+1;
                t=-start_end_dt_span+(ii-1)*dt_per_window+mean(trial_dt);
                these_times=(time_span>t-dt_window/2)&(time_span<=t+dt_window/2);
                all_sta_accuracy(ii_ROIst,ii_dt)=mean(accuracy_all_trials(these_times));
                all_sta_sh_accuracy(ii_ROIst,ii_dt)=mean(accuracy_sh_all_trials(these_times));
            end

        end


        %Calculate accuracy aligned to end
        for this_ROI=1:no_neurons

            %Now calculate the sliding window averages
            ii_ROIend=ii_ROIend+1;
            time_span=all_handles.run(this_ROI).handles_out.end_alligned.time_span;
            accuracy_all_trials=mean(all_handles.run(this_ROI).handles_out.end_alligned.accuracy_all_trials,2);
            accuracy_sh_all_trials=mean(all_handles.run(this_ROI).handles_out.end_alligned.accuracy_sh_all_trials,2);
            ii_dt=0;

            %Aligned to start
            for ii=1:(2*start_end_dt_span)+1
                ii_dt=ii_dt+1;
                t=-start_end_dt_span+(ii-1)*dt_per_window;
                these_times=(time_span>t-dt_window/2)&(time_span<=t+dt_window/2);
                all_end_accuracy(ii_ROIend,ii_dt)=mean(accuracy_all_trials(these_times));
                all_end_sh_accuracy(ii_ROIend,ii_dt)=mean(accuracy_sh_all_trials(these_times));
            end

            %Middle
            ii_dt=ii_dt+1;
            t_middle=mean(trial_dt/2);
            these_times=(time_span>t_middle-dt_window/2)&(time_span<=t_middle+dt_window/2);
            all_end_accuracy(ii_ROIend,ii_dt)=mean(accuracy_all_trials(these_times));
            all_end_sh_accuracy(ii_ROIend,ii_dt)=mean(accuracy_sh_all_trials(these_times));

            %Aligned to end
            for ii=1:(2*start_end_dt_span)+1
                ii_dt=ii_dt+1;
                t=-start_end_dt_span+(ii-1)*dt_per_window+mean(trial_dt);
                these_times=(time_span>t-dt_window/2)&(time_span<=t+dt_window/2);
                all_end_accuracy(ii_ROIend,ii_dt)=mean(accuracy_all_trials(these_times));
                all_end_sh_accuracy(ii_ROIend,ii_dt)=mean(accuracy_sh_all_trials(these_times));
            end

        end
    

    %Plot bar graph aligned to start
    figNo=0;
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig2 = figure(figNo);
    set(hFig2, 'units','normalized','position',[.1 .1 .6 .3])
    hold on

    bar_offset=0;
    edges=[0:0.05:1];
    rand_offset=0.8;
    for ii=1:ii_dt

        %Plot accuracy
        bar_offset=bar_offset+1;
        these_accuracies=zeros(1,size(all_sta_accuracy,1));
        these_accuracies(1,:)=all_sta_accuracy(:,ii);
        bar(bar_offset,mean(these_accuracies),'LineWidth', 3,'EdgeColor','none','FaceColor',[150/255 150/255 150/255])
        [mean_out, CIout]=drgViolinPoint(these_accuracies,edges,bar_offset,rand_offset,'k','k',3);

        %Plot shifted accuracy
        bar_offset=bar_offset+1;
        these_accuracies=zeros(1,size(all_sta_sh_accuracy,1));
        these_accuracies(1,:)=all_sta_sh_accuracy(:,ii);
        bar(bar_offset,mean(these_accuracies),'LineWidth', 3,'EdgeColor','none','FaceColor',[255/255 150/255 150/255])
        [mean_out, CIout]=drgViolinPoint(these_accuracies,edges,bar_offset,rand_offset,'k','k',3);

        bar_offset=bar_offset+1;

    end

    xticks([1.5 4.5 7.5 10.5 13.5 16.5 19.5 22.5 25.5 28.5 31.5 34.5 37.5 40.5 43.5])
    xticklabels({'st-3','st-2','st-1','st','st+1','st+2','st+3','mid','end-3','end-2','end-1','end','end+1','end+2','end+3'})
    ylabel('Accuracy')
    ylim([0 1])
    if handles.align_training_start==1
        title('Accuracy, trained at start, aligned to start')
    else
        title('Accuracy, trained at end, aligned to start')
    end

    %Plot bar graph aligned to end
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig2 = figure(figNo);
    set(hFig2, 'units','normalized','position',[.1 .1 .6 .3])
    hold on

    bar_offset=0;
    edges=[0:0.05:1];
    rand_offset=0.8;
    for ii=1:ii_dt

        %Plot accuracy
        bar_offset=bar_offset+1;
        these_accuracies=zeros(1,size(all_end_accuracy,1));
        these_accuracies(1,:)=all_end_accuracy(:,ii);
        bar(bar_offset,mean(these_accuracies),'LineWidth', 3,'EdgeColor','none','FaceColor',[150/255 150/255 150/255])
        [mean_out, CIout]=drgViolinPoint(these_accuracies,edges,bar_offset,rand_offset,'k','k',3);

        %Plot shifted accuracy
        bar_offset=bar_offset+1;
        these_accuracies=zeros(1,size(all_end_sh_accuracy,1));
        these_accuracies(1,:)=all_end_sh_accuracy(:,ii);
        bar(bar_offset,mean(these_accuracies),'LineWidth', 3,'EdgeColor','none','FaceColor',[255/255 150/255 150/255])
        [mean_out, CIout]=drgViolinPoint(these_accuracies,edges,bar_offset,rand_offset,'k','k',3);

        bar_offset=bar_offset+1;

    end

    xticks([1.5 4.5 7.5 10.5 13.5 16.5 19.5 22.5 25.5 28.5 31.5 34.5 37.5 40.5 43.5])
    xticklabels({'st-3','st-2','st-1','st','st+1','st+2','st+3','mid','end-3','end-2','end-1','end','end+1','end+2','end+3'})
    ylabel('Accuracy')
    ylim([0 1])
    if handles.align_training_start==1
        title('Accuracy, trained at start, aligned to end')
    else
        title('Accuracy, trained at end, aligned to end')
    end


    %Now let's take only the accuracies above 95 percentile
    ii_sta_99p_accuracy=zeros(1,size(all_end_accuracy,2));
    all_sta_99p_accuracy=[];
    ii_ROIs_sta_99p=[];
    ii_end_99p_accuracy=zeros(1,size(all_end_accuracy,2));
    all_end_99p_accuracy=[];
    ii_ROIs_end_99p=[];
    start_ROIs=[];
    end_ROIs=[];
    all_ROIs=[1:size(all_sta_accuracy,1)];
    all_max_end_99p_accuracy=-1*ones(1,size(all_end_accuracy,1));
    all_max_sta_99p_accuracy=-1*ones(1,size(all_end_accuracy,1));

    for ii=1:size(all_end_accuracy,2)

        these_accuracies=zeros(1,size(all_sta_accuracy,1));
        these_accuracies(1,:)=all_end_accuracy(:,ii);
        these_sh_accuracies=zeros(1,size(all_end_sh_accuracy,1));
        these_sh_accuracies(1,:)=all_end_sh_accuracy(:,ii);
        these_99p_accuracies=[];
        these_99p_accuracies=these_accuracies(these_accuracies>prctile(these_sh_accuracies,99));
        ii_end_99p_accuracy(ii)=length(these_99p_accuracies);
        all_end_99p_accuracy(ii,1:ii_end_99p_accuracy(ii))=these_99p_accuracies;
        these_ROIs=all_ROIs(these_accuracies>prctile(these_sh_accuracies,99));
        for jj_ROI=1:length(these_ROIs)
            if these_99p_accuracies(jj_ROI)>all_max_end_99p_accuracy(these_ROIs(jj_ROI))
                all_max_end_99p_accuracy(these_ROIs(jj_ROI))=these_99p_accuracies(jj_ROI);
            end
        end
        end_ROIs=[end_ROIs these_ROIs];

        these_accuracies=zeros(1,size(all_end_accuracy,1));
        these_accuracies(1,:)=all_sta_accuracy(:,ii);
        these_sh_accuracies=zeros(1,size(all_sta_sh_accuracy,1));
        these_sh_accuracies(1,:)=all_sta_sh_accuracy(:,ii);
        these_99p_accuracies=[];
        these_99p_accuracies=these_accuracies(these_accuracies>prctile(these_sh_accuracies,99));
        ii_sta_99p_accuracy(ii)=length(these_99p_accuracies);
        all_sta_99p_accuracy(ii,1:ii_sta_99p_accuracy(ii))=these_99p_accuracies;
        these_ROIs=all_ROIs(these_accuracies>prctile(these_sh_accuracies,99));
        for jj_ROI=1:length(these_ROIs)
            if these_99p_accuracies(jj_ROI)>all_max_sta_99p_accuracy(these_ROIs(jj_ROI))
                all_max_sta_99p_accuracy(these_ROIs(jj_ROI))=these_99p_accuracies(jj_ROI);
            end
        end
        start_ROIs=[start_ROIs these_ROIs];
    end

    end_ROIs=unique(end_ROIs);
    start_ROIs=unique(start_ROIs);

    all_max_end_99p_accuracy=all_max_end_99p_accuracy(all_max_end_99p_accuracy>-1);
    all_max_sta_99p_accuracy=all_max_sta_99p_accuracy(all_max_sta_99p_accuracy>-1);


    %Plot bar graph aligned to start

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig2 = figure(figNo);
    set(hFig2, 'units','normalized','position',[.1 .1 .6 .3])
    hold on

    bar_offset=0;
    edges=[0:0.05:1];
    rand_offset=0.8;
    for ii=1:ii_dt

        %Plot accuracy
        bar_offset=bar_offset+1;
        these_accuracies=zeros(1,ii_sta_99p_accuracy(ii));
        these_accuracies(1,:)=all_sta_99p_accuracy(ii,1:ii_sta_99p_accuracy(ii));
        bar(bar_offset,mean(these_accuracies),'LineWidth', 3,'EdgeColor','none','FaceColor',[150/255 150/255 150/255])
        [mean_out, CIout]=drgViolinPoint(these_accuracies,edges,bar_offset,rand_offset,'k','k',3);

        bar_offset=bar_offset+1;

    end

    xticks([1 3 5 7 9 11 13 15 17 19 21 23 25 27 29])
    xticklabels({'st-3','st-2','st-1','st','st+1','st+2','st+3','mid','end-3','end-2','end-1','end','end+1','end+2','end+3'})
    ylabel('Accuracy')
    ylim([0.5 0.8])
    if handles.align_training_start==1
        title('Significant ccuracy, trained at start, aligned to start')
    else
        title('Significant accuracy, trained at end, aligned to start')
    end

    %Plot bar graph aligned to end
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig2 = figure(figNo);
    set(hFig2, 'units','normalized','position',[.1 .1 .6 .3])
    hold on

    bar_offset=0;
    edges=[0:0.05:1];
    rand_offset=0.8;
    for ii=1:ii_dt

        %Plot accuracy
        bar_offset=bar_offset+1;
        these_accuracies=zeros(1,ii_end_99p_accuracy(ii));
        these_accuracies(1,:)=all_end_99p_accuracy(ii,1:ii_end_99p_accuracy(ii));
        bar(bar_offset,mean(these_accuracies),'LineWidth', 3,'EdgeColor','none','FaceColor',[150/255 150/255 150/255])
        [mean_out, CIout]=drgViolinPoint(these_accuracies,edges,bar_offset,rand_offset,'k','k',3);

    
        bar_offset=bar_offset+1;

    end

    xticks([1 3 5 7 9 11 13 15 17 19 21 23 25 27 29])
    xticklabels({'st-3','st-2','st-1','st','st+1','st+2','st+3','mid','end-3','end-2','end-1','end','end+1','end+2','end+3'})
    ylabel('Accuracy')
    ylim([0.5 0.8])
    if handles.align_training_start==1
        title('Significant accuracy, trained at start, aligned to end')
    else
        title('Significant accuracy, trained at end, aligned to end')
    end
 

    %Sort all the trials based on their time course
    %Normalize time to start to end
     
    norm_time=[-0.2:0.05:1.3];

    %Aligned to end
    all_norm_accuracies=zeros(length(norm_time),length(end_ROIs));
    iiROIs_included=0;
    for iiROI=end_ROIs

        fNo=fileNo_per_ROI(iiROI);
        iROI=iiROI_per_ROI(iiROI);
        time_span=all_handles.run(iROI).handles_out.end_alligned.time_span;
        accuracy_all_trials=mean(all_handles.run(iROI).handles_out.end_alligned.accuracy_all_trials,2);
        iiROIs_included=iiROIs_included+1;
        trial_dt=all_handles.run(iROI).handles_out.trial_dt;
        these_norm_times=norm_time*mean(trial_dt);
        time_span=all_handles.run(iROI).handles_out.end_alligned.time_span;
        for ii_t=1:length(these_norm_times)
            if sum(these_norm_times(ii_t)==time_span)==1
                all_norm_accuracies(ii_t,iiROIs_included)=accuracy_all_trials(these_norm_times(ii_t)==time_span);
            else
                ii_after=find(these_norm_times(ii_t)<time_span,1,'first');
                ii_before=find(these_norm_times(ii_t)>time_span,1,'last');
                all_norm_accuracies(ii_t,iiROIs_included)=accuracy_all_trials(ii_before)+...
                    (these_norm_times(ii_t)-time_span(ii_before))*...
                    (accuracy_all_trials(ii_after)-accuracy_all_trials(ii_before))/...
                    (time_span(ii_after)-time_span(ii_before));
            end
        end
    end

    %Do clustering

    croscorr_traces=corrcoef(all_norm_accuracies);

    %Set autocorrelations to zero
    for ii=1:size(croscorr_traces,1)
        croscorr_traces(ii,ii)=0;
    end
    Z = linkage(croscorr_traces,'complete','correlation');
    no_clusters=2;
    clusters = cluster(Z,'Maxclust',no_clusters);
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    %Do cutoff for 4 clusters
    cutoff = median([Z(end-(no_clusters-1),3) Z(end-(no_clusters-2),3)]);
    [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
    set(H,'LineWidth',2)
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

    %re-sort the matrix
    for ii=1:size(croscorr_traces,1)
        for jj=1:size(croscorr_traces,1)
            perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
        end
    end



    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
    hold on
    pcolor(perm_croscorr_traces)
    colormap fire
    shading flat

    caxis([-1  1])
    title(['Cross correlations for all ROIs'])
    % xlim([1 handles_out.all_div_ii_dFF])
    % ylim([1 handles_out.all_div_ii_dFF])

    %Plot rainbow
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


    prain=[0:0.6/99:0.6];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    %Plot timecourses for all ROIs


    %S+
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
    hold on



    sorted_handles_out.all_norm_accuracies=[];
    for ii=1:size(all_norm_accuracies,2)
        sorted_handles_out.all_norm_accuracies(:,ii)=all_norm_accuracies(:,outperm(ii));
    end

    time_span_mat=repmat(norm_time,ii,1)';
    ROI_mat=repmat(1:ii,length(norm_time),1);

    pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_norm_accuracies)
    colormap fire
    shading flat

    caxis([prctile(all_norm_accuracies(:),1) prctile(all_norm_accuracies(:),99)])


    plot([0 0],[0 size(all_norm_accuracies,2)],'-r')
    plot([1 1],[0 size(all_norm_accuracies,2)],'-r')


    xlim([-0.2 1.3])
    ylim([1 ii])
    title(['Accuracy aligned to end'])
    xlabel('Time (sec)')
    ylabel('ROI number')

    %Aligned to start
    all_norm_accuracies=zeros(length(norm_time),length(end_ROIs));
    iiROIs_included=0;
    for iiROI=end_ROIs

        fNo=fileNo_per_ROI(iiROI);
        iROI=iiROI_per_ROI(iiROI);
        time_span=all_handles.run(iROI).handles_out.start_alligned.time_span;
        accuracy_all_trials=mean(all_handles.run(iROI).handles_out.start_alligned.accuracy_all_trials,2);
        iiROIs_included=iiROIs_included+1;
        trial_dt=all_handles.run(iROI).handles_out.trial_dt;
        these_norm_times=norm_time*mean(trial_dt);
        time_span=all_handles.run(iROI).handles_out.start_alligned.time_span;
        for ii_t=1:length(these_norm_times)
            if sum(these_norm_times(ii_t)==time_span)==1
                all_norm_accuracies(ii_t,iiROIs_included)=accuracy_all_trials(these_norm_times(ii_t)==time_span);
            else
                ii_after=find(these_norm_times(ii_t)<time_span,1,'first');
                ii_before=find(these_norm_times(ii_t)>time_span,1,'last');
                all_norm_accuracies(ii_t,iiROIs_included)=accuracy_all_trials(ii_before)+...
                    (these_norm_times(ii_t)-time_span(ii_before))*...
                    (accuracy_all_trials(ii_after)-accuracy_all_trials(ii_before))/...
                    (time_span(ii_after)-time_span(ii_before));
            end
        end
    end

    %Do clustering

    croscorr_traces=corrcoef(all_norm_accuracies);

    %Set autocorrelations to zero
    for ii=1:size(croscorr_traces,1)
        croscorr_traces(ii,ii)=0;
    end
    Z = linkage(croscorr_traces,'complete','correlation');
    no_clusters=2;
    clusters = cluster(Z,'Maxclust',no_clusters);
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    %Do cutoff for 4 clusters
    cutoff = median([Z(end-(no_clusters-1),3) Z(end-(no_clusters-2),3)]);
    [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
    set(H,'LineWidth',2)
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

    %re-sort the matrix
    for ii=1:size(croscorr_traces,1)
        for jj=1:size(croscorr_traces,1)
            perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
        end
    end



    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
    hold on
    pcolor(perm_croscorr_traces)
    colormap fire
    shading flat

    caxis([-1  1])
    title(['Cross correlations for all ROIs'])
    % xlim([1 handles_out.all_div_ii_dFF])
    % ylim([1 handles_out.all_div_ii_dFF])

    %Plot rainbow
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


    prain=[0:0.6/99:0.6];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    %Plot timecourses for all ROIs


    %S+
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
    hold on



    sorted_handles_out.all_norm_accuracies=[];
    for ii=1:size(all_norm_accuracies,2)
        sorted_handles_out.all_norm_accuracies(:,ii)=all_norm_accuracies(:,outperm(ii));
    end

    time_span_mat=repmat(norm_time,ii,1)';
    ROI_mat=repmat(1:ii,length(norm_time),1);

    pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_norm_accuracies)
    colormap fire
    shading flat

    caxis([prctile(all_norm_accuracies(:),1) prctile(all_norm_accuracies(:),99)])


    plot([0 0],[0 size(all_norm_accuracies,2)],'-r')
    plot([1 1],[0 size(all_norm_accuracies,2)],'-r')


    xlim([-0.2 1.3])
    ylim([1 ii])
    title(['Accuracy aligned to start'])
    xlabel('Time (sec)')
    ylabel('ROI number')

    %Now compare to Moser analysis
    pred_file=handles.pred_file{fileNo};
    load([this_path pred_file])
    pffft=1;

    %Aligned to end
    %Plot max accuracy vs spatial correlation
    spatial_rhol1l4=handles_out.spatial_rhol1l4;
    max_all_end_accuracy=zeros(1,no_neurons);
    for ii_neuron=1:no_neurons
        these_accuracies=zeros(1,size(all_end_accuracy,2));
        these_accuracies(1,:)=all_end_accuracy(ii_neuron,:);
        max_all_end_accuracy(ii_neuron)=max(these_accuracies);
    end

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(spatial_rhol1l4,max_all_end_accuracy,'ob')
    xlabel('Spatial rho')
    ylabel('Accuracy')
    title('Rho vs. accuracy aligned to end')

    %Plot max accuracy vs center of mass
    delta_center_of_mass=handles_out.delta_center_of_mass;

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(delta_center_of_mass,max_all_end_accuracy,'ob')
    xlabel('Delta center of mass')
    ylabel('Accuracy')
    title('dCOM vs. accuracy aligned to end')

    %Now plot only those above 99 percentile
    
    %Plot max accuracy vs spatial correlation


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(spatial_rhol1l4(end_ROIs),all_max_end_99p_accuracy,'ob')
    xlabel('Spatial rho')
    ylabel('Accuracy')
    title('Rho vs. significant accuracy aligned to end')

    %Plot max accuracy vs center of mass
  

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(delta_center_of_mass(end_ROIs),all_max_end_99p_accuracy,'ob')
    xlabel('Delta center of mass')
    ylabel('Accuracy')
    title('dCOM vs. accuracy aligned to end')

    %Aligned to start

    %Plot max accuracy vs spatial correlation
  
    max_all_sta_accuracy=zeros(1,no_neurons);
    for ii_neuron=1:no_neurons
        these_accuracies=zeros(1,size(all_sta_accuracy,2));
        these_accuracies(1,:)=all_sta_accuracy(ii_neuron,:);
        max_all_sta_accuracy(ii_neuron)=max(these_accuracies);
    end

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(spatial_rhol1l4,max_all_sta_accuracy,'ob')
    xlabel('Spatial rho')
    ylabel('Accuracy')
    title('Rho vs. accuracy aligned to start')

    %Plot max accuracy vs center of mass
    

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(delta_center_of_mass,max_all_sta_accuracy,'ob')
    xlabel('Delta center of mass')
    ylabel('Accuracy')
    title('dCOM vs. accuracy aligned to start')

    %Now plot only those above 99 percentile
    
    %Plot max accuracy vs spatial correlation


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(spatial_rhol1l4(start_ROIs),all_max_sta_99p_accuracy,'ob')
    xlabel('Spatial rho')
    ylabel('Accuracy')
    title('Rho vs. significant accuracy aligned to start')

    %Plot max accuracy vs center of mass
  

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
    hold on

    plot(delta_center_of_mass(start_ROIs),all_max_sta_99p_accuracy,'ob')
    xlabel('Delta center of mass')
    ylabel('Accuracy')
    title('dCOM vs. accuracy aligned to start')

    pffft=1;

end
