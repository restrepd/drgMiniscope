
% drgMini_display_differential_response_per_ROI_hit_vs_miss_v2

%Now show the divergent responses for after
input_file='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/glm_div_end.mat';

load(input_file)

handles_out_div=[];
figureNo=0;

try
    colormap fire
catch
    loadedData=load('/data2/DRMatlab/drgMiniscope/fireColormap2023b.mat');
    fire_map=loadedData.fireMap;
    colormap(fire_map)
end
this_cmap=colormap;
this_cmap(1,:)=[0.3 0.3 0.3];

these_adiv_names{1}='all_div_hit';
these_adiv_names{2}='all_div_miss';

these_adiv_names_le1{1}='all_div_hit_le1';
these_adiv_names_le1{2}='all_div_miss_le1';

these_adiv_names_le0{1}='all_div_hit_le0';
these_adiv_names_le0{2}='all_div_miss_le0';


trial_type_labels{1}='Hit';
trial_type_labels{2}='Miss';

%Reject small cahnges
SD_mult=4;

no_clusters=6;
no_clusters_le1=3;
no_clusters_le0=3;

%Get all divergent ROIs
all_div_hit=[];
all_div_miss=[];
all_div_hit_end=[];

%And all divergent with large_enough==1
all_div_hit_le1=[];
all_div_miss_le1=[];
all_div_hit_end_le1=[];

%And all divergent with large_enough==0
all_div_hit_le0=[];
all_div_miss_le0=[];
all_div_hit_end_le0=[];

all_div_to_sort=[];
all_div_to_sort_le0=[];
all_div_to_sort_le1=[];
 
for fileNo=1:length(handles_out2.file)
    if ~isempty(handles_out2.file(fileNo).ROI)
        odor_aligned_time_bins=handles_out2.file(fileNo).odor_aligned_time_bins;
        these_significant_aft=handles_out2.file(fileNo).significant_aft;
        these_significant_bef=handles_out2.file(fileNo).significant_bef;
        no_neurons=length(these_significant_bef);
 
        for ii_ROI=1:no_neurons
            these_bef=[handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF(odor_aligned_time_bins<0)...
                handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF(odor_aligned_time_bins<0)];
            sd_these_bef=std(these_bef);
            abs_delta_aft_hit=abs(handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF(odor_aligned_time_bins>=0)...
                -mean(handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF(odor_aligned_time_bins<0)));
            abs_delta_aft_miss=abs(handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF(odor_aligned_time_bins>=0)...
                -mean(handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF(odor_aligned_time_bins<0)));
            if (sum(abs_delta_aft_hit>=SD_mult*sd_these_bef)>0)||(sum(abs_delta_aft_miss>=SD_mult*sd_these_bef)>0)
                large_enough=1;
                handles_out_div.file(fileNo).large_enough(ii_ROI)=1;
            else
                large_enough=0;
                handles_out_div.file(fileNo).large_enough(ii_ROI)=0;
            end
            if (these_significant_aft(ii_ROI)==1)&(these_significant_bef(ii_ROI)==0)
                handles_out_div.file(fileNo).sig_div(ii_ROI)=1;
                 %aligned to odor encounter
                all_div_hit=[all_div_hit; handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF...
                    /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];
                all_div_miss=[all_div_miss; handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF...
                    /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];

                %aligned to reward
                all_div_hit_end=[all_div_hit_end; handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF_end...
                    /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF_end],95)];
                if (large_enough==1)
                    
                    %aligned to odor encounter
                    all_div_hit_le1=[all_div_hit_le1; handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF...
                        /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];
                    all_div_miss_le1=[all_div_miss_le1; handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF...
                        /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];

                    %aligned to reward
                    all_div_hit_end_le1=[all_div_hit_end_le1; handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF_end...
                        /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF_end],95)];
                else
                    
                    
                     %aligned to odor encounter
                    all_div_hit_le0=[all_div_hit_le0; handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF...
                        /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];
                    all_div_miss_le0=[all_div_miss_le0; handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF...
                        /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];

                    %aligned to reward
                    all_div_hit_end_le0=[all_div_hit_end_le0; handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF_end...
                        /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF_end],95)];
                end
            else
                handles_out_div.file(fileNo).sig_div(ii_ROI)=0;
            end
        end

        %Sort using the two ii_types that were used in the comparison
        all_div_to_sort=[all_div_hit all_div_miss];
        all_div_to_sort_le1=[all_div_hit_le1 all_div_miss_le1];
        all_div_to_sort_le0=[all_div_hit_le0 all_div_miss_le0];
      
    end
end

save([input_file(1:end-4) '_sig.mat'],'handles_out_div')

all_div=[all_div_hit;all_div_miss];
all_div_le1=[all_div_hit_le1;all_div_miss_le1];
all_div_le0=[all_div_hit_le0;all_div_miss_le0];

%Do the crosscorrelation/linkage analysis for all the divergent ROIs
croscorr_traces=corrcoef(all_div_to_sort');

Z = linkage(croscorr_traces,'complete','correlation');


handles_out2.clusters = cluster(Z,'Maxclust',no_clusters);
figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);


n = size(Z,1);              % number of merges = N-1 if you had N points

if no_clusters <= 1 || no_clusters > n+1
    error('no_clusters must be between 2 and %d', n+1);
end

% Take the last (no_clusters-1) linkage heights in column 3
% last_heights = Z(n-(no_clusters-2):n, 3);
% cutoff = median(last_heights);
% 
% [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);

[H,T,outperm]=dendrogram(Z,0,'Orientation','left','ClusterIndices',handles_out2.clusters);



set(H,'LineWidth',2)
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

clusters=cluster(Z,'Maxclust',no_clusters);

%I am forcing the clusters
% clusters(36)=1;
% clusters(8)=1;
% clusters(7)=1;
% clusters(5)=1;
% clusters(4)=1;
% clusters(9)=1;
% clusters(29)=1;
% clusters(10)=1;
% clusters(2)=1;

%re-sort the matrix
for ii=1:size(croscorr_traces,1)
    for jj=1:size(croscorr_traces,1)
        perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
    end
end


%Print the crosscorrelation matrix
figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
hold on
pcolor(perm_croscorr_traces)
colormap(this_cmap)
shading flat

% caxis([-1  1])
clim([-1 1])
% eval(['title([''Cross correlations for all ROIs ' trial_type_comp_labels{ii_comp}  '''])'])
title('Cross correlations for all ROIs for hit vs miss')

ROIs_included=size(all_div_to_sort,1);
xlim([1 ROIs_included])
ylim([1 ROIs_included])

%Print a pseudocolor of the dFF time courses aligned to odorant encounter
for ii_type=1:2
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
    hold on

    eval(['sorted_handles_out.' these_adiv_names{ii_type} '=[];'])
    eval(['sorted_handles_out.' these_adiv_names{ii_type} '= ' these_adiv_names{ii_type} '(outperm,:);'])

    eval(['ii_included=size(sorted_handles_out.' these_adiv_names{ii_type} ',1);'])
    
    time_span_mat=repmat(odor_aligned_time_bins,ii_included,1);
    ROI_mat=repmat(1:ii_included,length(odor_aligned_time_bins),1)';

    eval(['pcolor(time_span_mat,ROI_mat,sorted_handles_out.' these_adiv_names{ii_type} ')'])
    
    colormap(this_cmap)
    shading flat

    caxis([prctile(all_div(:),1) prctile(all_div(:),99.9)])

    plot([0 0],[0 ii_included],'-k','LineWidth',2)


    xlim([-3 3])
    ylim([1 ii_included])
    eval(['title([''' trial_type_labels{ii_type} '''])'])
    xlabel('Time (sec)')
    ylabel('ROI number')
end


%Print a pseudocolor of the dFF time courses for
% hits aligned to reward

figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
hold on

sorted_handles_out.all_div_hit_end=[];
sorted_handles_out.all_div_hit_end=all_div_hit_end(outperm,:);

pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_hit_end);

colormap(this_cmap)
shading flat

caxis([prctile(all_div_hit_end(:),1) prctile(all_div_hit_end(:),99.9)])

plot([0 0],[0 ii_included],'-k','LineWidth',2)

xlim([-3 3])
ylim([1 ii_included])
title('dFF timecourse for hits aligned to reward')
xlabel('Time (sec)')
ylabel('ROI number')

%Report the number of ROIs with high divergence
fprintf(1, ['\nFor hit vs miss there were ' num2str(length(clusters)) ' divergent ROIs\n'])

%Plot the average timecourses per cluster
do_std=0;
these_ylim=[];
for clus=1:no_clusters
    %Report the number of ROIs per cluster
    fprintf(1, ['\nFor hit vs miss cluster No '...
        num2str(clus) ' has ' num2str(sum(clusters==clus)) ' ROIs\n'])

    %Plot the average of this cluster for the type of trials compared
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .2 .18 .18])


    hold on

    %get the dF/F


    %Miss
    this_cluster_miss=zeros(sum(clusters==clus),size(all_div_miss,2));
    ii_included=0;
    no_ROIs_this_clus= size(all_div_miss,1);
    for ii=1:no_ROIs_this_clus
        if clusters(ii)==clus
            ii_included=ii_included+1;
            this_cluster_miss(ii_included,:) = all_div_miss(ii,:);
        end
    end


    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_miss);
        meanpv=mean(this_cluster_miss,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_miss);
        meanpv=mean(this_cluster_miss,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end


    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(this_cluster_miss), CIpv','cmap',[86/255 180/255 233/255]);

    %Hits
    this_cluster_hit=zeros(sum(clusters==clus),size(all_div_miss,2));

    ii_included=0;
    no_ROIs_this_clus= size(all_div_hit,1);
    for ii=1:no_ROIs_this_clus
        if clusters(ii)==clus
            ii_included=ii_included+1;
            this_cluster_hit(ii_included,:) = all_div_hit(ii,:);
        end
    end

    % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_hit);
        meanpv=mean(this_cluster_hit,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_hit);
        meanpv=mean(this_cluster_hit,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end


    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins, mean(this_cluster_hit), CIpv','cmap',[230/255 159/255 0/255]);


    plot(odor_aligned_time_bins',mean(this_cluster_miss)','Color',[86/255 180/255 233/255],'LineWidth',1.5);
    plot(odor_aligned_time_bins',mean(this_cluster_hit)','Color',[230/255 159/255 0/255],'LineWidth',1.5);



    xlim([-3 3])
    these_ylim=[these_ylim; ylim];
    xlabel('Time(sec)')
    ylabel('dFF')
    title(['dFF mean for hit vs miss cluster No ' num2str(clus)])
    
end


fig_minus=0;
for clus=1:no_clusters

    figNo=figureNo+fig_minus;
    figure(figNo)
    ylim([min(these_ylim(:)) max(these_ylim(:))])
    this_ylim=ylim;

    %Odor on markers
    plot([0 0],this_ylim,'-k')

    text(-2.5,this_ylim(1)+0.9*(this_ylim(2)-this_ylim(1)),'Miss','Color',[86/255 180/255 233/255])
    text(-2.5,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),'Hit','Color',[230/255 159/255 0/255])


    %Odor on markers
    plot([0 0],this_ylim,'-k')

    fig_minus=fig_minus-1;

end

%Do the crosscorrelation/linkage analysis for large divergence ROIs 
% (large_enough == 1)
croscorr_traces_le1=corrcoef(all_div_to_sort_le1');

Z_le1 = linkage(croscorr_traces_le1,'complete','correlation');


handles_out2.clusters_le1 = cluster(Z_le1,'Maxclust',no_clusters_le1);

figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);


n_le1 = size(Z_le1,1);              % number of merges = N-1 if you had N points

if no_clusters_le1 <= 1 || no_clusters_le1 > n_le1+1
    error('no_clusters must be between 2 and %d', n_le1+1);
end

% Take the last (no_clusters-1) linkage heights in column 3
% last_heights_le1 = Z_le1(n_le1-(no_clusters_le1-2):n_le1, 3);
% cutoff_le1 = median(last_heights_le1);
% 
% [H_le1,T_le1,outperm_le1]=dendrogram(Z_le1,0,'Orientation','left','ColorThreshold',cutoff_le1);

[H_le1,T_le1,outperm_le1]=dendrogram(Z_le1,0,'Orientation','left','ClusterIndices',handles_out2.clusters_le1);


set(H_le1,'LineWidth',2)
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

clusters_le1=cluster(Z_le1,'Maxclust',no_clusters_le1);

%I am forcing the clusters
% clusters_le1(36)=1;
% clusters_le1(8)=1;
% clusters_le1(7)=1;
% clusters_le1(5)=1;
% clusters_le1(4)=1;
% clusters_le1(9)=1;
% clusters_le1(29)=1;
% clusters_le1(10)=1;
% clusters_le1(2)=1;

%re-sort the matrix
for ii=1:size(croscorr_traces_le1,1)
    for jj=1:size(croscorr_traces_le1,1)
        perm_croscorr_traces_le1(ii,jj)=croscorr_traces_le1(outperm_le1(ii),outperm_le1(jj));
    end
end


%Print the crosscorrelation matrix
figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
hold on
pcolor(perm_croscorr_traces_le1)
colormap(this_cmap)
shading flat

clim([-1 1])

title('Cross correlations for all ROIs for hit vs miss for large divergence')

ROIs_included_le1=size(all_div_to_sort_le1,1);
xlim([1 ROIs_included_le1])
ylim([1 ROIs_included_le1])

%Plot rainbow
figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


prain=[-1:0.025:1];
pcolor(repmat([1:10],length(prain),1)',repmat(prain,10,1),repmat(prain,10,1))
%             colormap jet
colormap(this_cmap)
shading interp
ax=gca;
set(ax,'XTickLabel','')

%Print a pseudocolor of the dFF time courses aligned to odorant encounter
for ii_type=1:2
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
    hold on

    eval(['sorted_handles_out.' these_adiv_names_le1{ii_type} '=[];'])
    eval(['sorted_handles_out.' these_adiv_names_le1{ii_type} '= ' these_adiv_names_le1{ii_type} '(outperm_le1,:);'])

    eval(['ii_included_le1=size(sorted_handles_out.' these_adiv_names_le1{ii_type} ',1);'])
    
    time_span_mat=repmat(odor_aligned_time_bins,ii_included_le1,1);
    ROI_mat_le1=repmat(1:ii_included_le1,length(odor_aligned_time_bins),1)';

    eval(['pcolor(time_span_mat,ROI_mat_le1,sorted_handles_out.' these_adiv_names_le1{ii_type} ')'])
    
    colormap(this_cmap)
    shading flat

    caxis([prctile(all_div_le1(:),1) prctile(all_div_le1(:),99.9)])

    plot([0 0],[0 ii_included_le1],'-k','LineWidth',2)


    xlim([-3 3])
    ylim([1 ii_included_le1])
    eval(['title([''' trial_type_labels{ii_type} ' high divergence' ''' ])' ])
    xlabel('Time (sec)') 
    ylabel('ROI number')
end


%Print a pseudocolor of the dFF time courses for
% hits aligned to reward

figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
hold on

sorted_handles_out.all_div_hit_end_le1=[];
sorted_handles_out.all_div_hit_end_le1=all_div_hit_end_le1(outperm_le1,:);

pcolor(time_span_mat,ROI_mat_le1,sorted_handles_out.all_div_hit_end_le1);

colormap(this_cmap)
shading flat

caxis([prctile(all_div_hit_end_le1(:),1) prctile(all_div_hit_end_le1(:),99.9)])

plot([0 0],[0 ii_included_le1],'-k','LineWidth',2)

xlim([-3 3])
ylim([1 ii_included_le1])
title('Hits aligned to reward, high divergence')
xlabel('Time (sec)')
ylabel('ROI number')

%Report the number of ROIs with high divergence
fprintf(1, ['\nFor hit vs miss high divergence there were ' num2str(length(clusters_le1)) ' divergent ROIs\n'])

%Plot the average timecourses per cluster
do_std=0;
these_ylim=[];

for clus=1:no_clusters_le1
    %Report the number of ROIs per cluster
    fprintf(1, ['\nFor hit vs miss high divergence cluster No '...
        num2str(clus) ' has ' num2str(sum(clusters_le1==clus)) ' ROIs\n'])

    %Plot the average of this cluster for the type of trials compared
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .2 .18 .18])


    hold on

    %get the dF/F


    %Miss
    this_cluster_miss_le1=zeros(sum(clusters_le1==clus),size(all_div_miss_le1,2));
    ii_included_le1=0;
    no_ROIs_this_clus_le1= size(all_div_miss_le1,1);
    for ii=1:no_ROIs_this_clus_le1
        if clusters_le1(ii)==clus
            ii_included_le1=ii_included_le1+1;
            this_cluster_miss_le1(ii_included_le1,:) = all_div_miss_le1(ii,:);
        end
    end


    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_miss_le1);
        meanpv=mean(this_cluster_miss_le1,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_miss_le1);
        meanpv=mean(this_cluster_miss_le1,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end


    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(this_cluster_miss_le1), CIpv','cmap',[86/255 180/255 233/255]);

    %Hits
    this_cluster_hit_le1=zeros(sum(clusters_le1==clus),size(all_div_miss_le1,2));

    ii_included_le1=0;
    no_ROIs_this_clus_le1= size(all_div_hit_le1,1);
    for ii=1:no_ROIs_this_clus_le1
        if clusters_le1(ii)==clus
            ii_included_le1=ii_included_le1+1;
            this_cluster_hit_le1(ii_included_le1,:) = all_div_hit_le1(ii,:);
        end
    end

    % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_hit_le1);
        meanpv=mean(this_cluster_hit_le1,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_hit_le1);
        meanpv=mean(this_cluster_hit_le1,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end



    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins, mean(this_cluster_hit_le1), CIpv','cmap',[230/255 159/255 0/255]);


    plot(odor_aligned_time_bins',mean(this_cluster_miss_le1)','Color',[86/255 180/255 233/255],'LineWidth',1.5);
    plot(odor_aligned_time_bins',mean(this_cluster_hit_le1)','Color',[230/255 159/255 0/255],'LineWidth',1.5);


    % ylim([-0.5 1.5])
    xlim([-3 3])
    these_ylim=[these_ylim; ylim];




    xlabel('Time(sec)')
    ylabel('dFF')
    title(['dFF high divergence clNo ' num2str(clus)])
    
end


fig_minus=0;
for clus=1:no_clusters_le1

    figNo=figureNo+fig_minus;
    figure(figNo)
    ylim([min(these_ylim(:)) max(these_ylim(:))])
    this_ylim=ylim;

    %Odor on markers
    plot([0 0],this_ylim,'-k')

    text(-2.5,this_ylim(1)+0.9*(this_ylim(2)-this_ylim(1)),'Miss','Color',[86/255 180/255 233/255])
    text(-2.5,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),'Hit','Color',[230/255 159/255 0/255])


    %Odor on markers
    plot([0 0],this_ylim,'-k')

    fig_minus=fig_minus-1;

end


%Do the crosscorrelation/linkage analysis for low divergence ROIs 
% (large_enough == 0)
croscorr_traces_le0=corrcoef(all_div_to_sort_le0');

Z_le0 = linkage(croscorr_traces_le0,'complete','correlation');

% flat clusters
handles_out2.clusters_le0 = cluster(Z_le0,'Maxclust',no_clusters_le0);

figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);


n_le0 = size(Z_le0,1);              % number of merges = N-1 if you had N points

if no_clusters_le0 <= 1 || no_clusters_le0 > n_le0+1
    error('no_clusters must be between 2 and %d', n_le0+1);
end

% Take the last (no_clusters-1) linkage heights in column 3
% last_heights_le0 = Z_le0(n_le0-(no_clusters_le0-2):n_le0, 3);
% cutoff_le0 = median(last_heights_le0);
% 
% [H_le0,T_le0,outperm_le0]=dendrogram(Z_le0,0,'Orientation','left','ColorThreshold',cutoff_le0);


[H_le0,T_le0,outperm_le0]=dendrogram(Z_le0,0,'Orientation','left','ClusterIndices',handles_out2.clusters_le0);


set(H_le0,'LineWidth',2)
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

clusters_le0=cluster(Z_le0,'Maxclust',no_clusters_le0);

%I am forcing the clusters
% clusters_le0(36)=1;
% clusters_le0(8)=1;
% clusters_le0(7)=1;
% clusters_le0(5)=1;
% clusters_le0(4)=1;
% clusters_le0(9)=1;
% clusters_le0(29)=1;
% clusters_le0(10)=1;
% clusters_le0(2)=1;

%re-sort the matrix
for ii=1:size(croscorr_traces_le0,1)
    for jj=1:size(croscorr_traces_le0,1)
        perm_croscorr_traces_le0(ii,jj)=croscorr_traces_le0(outperm_le0(ii),outperm_le0(jj));
    end
end


%Print the crosscorrelation matrix
figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
hold on
pcolor(perm_croscorr_traces_le0)
colormap(this_cmap)
shading flat

clim([-1 1])

title('Cross correlations for all ROIs for hit vs miss for large divergence')

ROIs_included_le0=size(all_div_to_sort_le0,1);
xlim([1 ROIs_included_le0])
ylim([1 ROIs_included_le0])

%Print a pseudocolor of the dFF time courses aligned to odorant encounter
for ii_type=1:2
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
    hold on

    eval(['sorted_handles_out.' these_adiv_names_le0{ii_type} '=[];'])
    eval(['sorted_handles_out.' these_adiv_names_le0{ii_type} '= ' these_adiv_names_le0{ii_type} '(outperm_le0,:);'])

    eval(['ii_included_le0=size(sorted_handles_out.' these_adiv_names_le0{ii_type} ',1);'])
    
    time_span_mat=repmat(odor_aligned_time_bins,ii_included_le0,1);
    ROI_mat_le0=repmat(1:ii_included_le0,length(odor_aligned_time_bins),1)';

    eval(['pcolor(time_span_mat,ROI_mat_le0,sorted_handles_out.' these_adiv_names_le0{ii_type} ')'])
    
    colormap(this_cmap)
    shading flat

    caxis([prctile(all_div_le0(:),1) prctile(all_div_le0(:),99.9)])

    plot([0 0],[0 ii_included_le0],'-k','LineWidth',2)


    xlim([-3 3])
    ylim([1 ii_included_le0])
    eval(['title([''' trial_type_labels{ii_type} ' low divergence' ''' ])' ])
    xlabel('Time (sec)') 
    ylabel('ROI number')
end


%Print a pseudocolor of the dFF time courses for
% hits aligned to reward

figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
hold on

sorted_handles_out.all_div_hit_end_le0=[];
sorted_handles_out.all_div_hit_end_le0=all_div_hit_end_le0(outperm_le0,:);

pcolor(time_span_mat,ROI_mat_le0,sorted_handles_out.all_div_hit_end_le0);

colormap(this_cmap)
shading flat

caxis([prctile(all_div_hit_end_le0(:),1) prctile(all_div_hit_end_le0(:),99.9)])

plot([0 0],[0 ii_included_le0],'-k','LineWidth',2)

xlim([-3 3])
ylim([1 ii_included_le0])
title('Hits aligned to reward, low divergence')
xlabel('Time (sec)')
ylabel('ROI number')

%Report the number of ROIs with high divergence
fprintf(1, ['\nFor hit vs miss low divergence there were ' num2str(length(clusters_le0)) ' divergent ROIs\n'])

%Plot the average timecourses per cluster
do_std=0;
these_ylim=[];
for clus=1:no_clusters_le0
    %Report the number of ROIs per cluster
    fprintf(1, ['\nFor hit vs miss low divergence cluster No '...
        num2str(clus) ' has ' num2str(sum(clusters_le0==clus)) ' ROIs\n'])

    %Plot the average of this cluster for the type of trials compared
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .2 .18 .18])


    hold on

    %get the dF/F


    %Miss
    this_cluster_miss_le0=zeros(sum(clusters_le0==clus),size(all_div_miss_le0,2));
    ii_included_le0=0;
    no_ROIs_this_clus_le0= size(all_div_miss_le0,1);
    for ii=1:no_ROIs_this_clus_le0
        if clusters_le0(ii)==clus
            ii_included_le0=ii_included_le0+1;
            this_cluster_miss_le0(ii_included_le0,:) = all_div_miss_le0(ii,:);
        end
    end


    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_miss_le0);
        meanpv=mean(this_cluster_miss_le0,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_miss_le0);
        meanpv=mean(this_cluster_miss_le0,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end


    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(this_cluster_miss_le0), CIpv','cmap',[86/255 180/255 233/255]);

    %Hits
    this_cluster_hit_le0=zeros(sum(clusters_le0==clus),size(all_div_miss_le0,2));

    ii_included_le0=0;
    no_ROIs_this_clus_le0= size(all_div_hit_le0,1);
    for ii=1:no_ROIs_this_clus_le0
        if clusters_le0(ii)==clus
            ii_included_le0=ii_included_le0+1;
            this_cluster_hit_le0(ii_included_le0,:) = all_div_hit_le0(ii,:);
        end
    end

    % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_hit_le0);
        meanpv=mean(this_cluster_hit_le0,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_hit_le0);
        meanpv=mean(this_cluster_hit_le0,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end



    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins, mean(this_cluster_hit_le0), CIpv','cmap',[230/255 159/255 0/255]);


    plot(odor_aligned_time_bins',mean(this_cluster_miss_le0)','Color',[86/255 180/255 233/255],'LineWidth',1.5);
    plot(odor_aligned_time_bins',mean(this_cluster_hit_le0)','Color',[230/255 159/255 0/255],'LineWidth',1.5);


    % ylim([-0.5 1.5])
    xlim([-3 3])
    these_ylim=[these_ylim; ylim];




    xlabel('Time(sec)')
    ylabel('dFF')
    title(['dFF low divergence clNo ' num2str(clus)])
    
end


fig_minus=0;
for clus=1:no_clusters_le0

    figNo=figureNo+fig_minus;
    figure(figNo)
    ylim([min(these_ylim(:)) max(these_ylim(:))])
    this_ylim=ylim;

    %Odor on markers
    plot([0 0],this_ylim,'-k','LineWidth', 2)

    text(-2.5,this_ylim(1)+0.9*(this_ylim(2)-this_ylim(1)),'Miss','Color',[86/255 180/255 233/255])
    text(-2.5,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),'Hit','Color',[230/255 159/255 0/255])


    % %Odor on markers
    % plot([0 0],this_ylim,'-k')

    fig_minus=fig_minus-1;

end

pffft=1;

