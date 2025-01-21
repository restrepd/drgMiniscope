%drgMini_analyze_batch_ImpMoserXYConc
close all
clear all

is_sphgpu=0;

switch is_sphgpu
    case 0

        % %There was a bug and the shuffled runs were the same for all shuffled runs
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';
        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        %For the files below the shuffled runs should be different
        %Trained with all trials
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01062025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01062025.m';

         %Trained with hits only
         save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
         choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');



    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end


percentile_stringency=95; %mean_imps at this percentile stringency are included as cells of high prediction importance

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


% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])
eval(['handles_Moser=' choiceMoserFileName(1:end-2) ';'])

figureNo=0;

colormap fire
this_cmap=colormap;
this_cmap(1,:)=[0.3 0.3 0.3];

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);


%Let's look at predictive importance
these_groups=[1 5];
ii_run=1;
all_mean_conc_imps=[];
above_thr_conc_imps=[];
all_mean_x_imps=[];
above_thr_x_imps=[];
all_mean_y_imps=[];
above_thr_y_imps=[];
thr_x_imps=zeros(1,length(handles_conc.arena_file));
thr_y_imps=zeros(1,length(handles_conc.arena_file));
thr_conc_imps=zeros(1,length(handles_conc.arena_file));
all_spatial_rhol1l4=[];
all_delta_center_of_mass=[];
all_information_content=[];
all_sparsity=[];
conc_spatial_rhol1l4=[];
conc_delta_center_of_mass=[];
conc_information_content=[];
conc_sparsity=[];
x_spatial_rhol1l4=[];
x_delta_center_of_mass=[];
x_information_content=[];
x_sparsity=[];
y_spatial_rhol1l4=[];
y_delta_center_of_mass=[];
y_information_content=[];
y_sparsity=[];

file_numbers=[];


figureNo=figureNo+1;
these_R1s=[];
imps=[];

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
        arena_file=handles_conc.arena_file{fileNo};

        %Get predictor importance for conc decoding
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        all_conc_imps=[];
        for trNo=1:length(handles_out.imp.trial)
            these_imps=handles_out.imp.trial(trNo).imp;
            all_conc_imps=[all_conc_imps these_imps'];
        end
        sum_all_conc_imps=sum(all_conc_imps);
        mat_sum_all_conc_imps=repmat(sum_all_conc_imps,size(all_conc_imps,1),1);
        all_conc_imps=all_conc_imps./mat_sum_all_conc_imps; %Normalize to all added imps=1
        mean_conc_imps=mean(all_conc_imps,2);
        thr_conc_imps(fileNo)=prctile(mean_conc_imps,percentile_stringency);
        % imps.file(fileNo).conc.mean_conc_imps=mean_conc_imps;
        % imps.file(fileNo).conc.thr_conc_imps=prctile(mean_conc_imps,percentile_stringency);
        all_mean_conc_imps=[all_mean_conc_imps; mean_conc_imps];
        above_thr_conc_imps=[above_thr_conc_imps; mean_conc_imps>=thr_conc_imps(fileNo)];
        file_numbers=[file_numbers; fileNo*ones(length(mean_conc_imps),1)];
        these_ROIs=[1:length(mean_conc_imps)];
        imps.file(fileNo).ROIs_conc=these_ROIs(mean_conc_imps>=thr_conc_imps(fileNo));

        %Get predictor importance for x and y decoding
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        all_x_imps=[];
        all_y_imps=[];
        for trNo=1:length(handles_out.imp.trial)
            these_x_imps=handles_out.imp.trial(trNo).impx;
            all_x_imps=[all_x_imps these_x_imps'];
            these_y_imps=handles_out.imp.trial(trNo).impy;
            all_y_imps=[all_y_imps these_y_imps'];
        end
        sum_all_x_imps=sum(all_x_imps);
        mat_sum_all_x_imps=repmat(sum_all_x_imps,size(all_x_imps,1),1);
        all_x_imps=all_x_imps./mat_sum_all_x_imps; %Normalize to all added imps=1
        mean_x_imps=mean(all_x_imps,2);
        thr_x_imps(fileNo)=prctile(mean_x_imps,percentile_stringency);
        all_mean_x_imps=[all_mean_x_imps; mean_x_imps];
        above_thr_x_imps=[above_thr_x_imps; mean_x_imps>=thr_x_imps(fileNo)];
        imps.file(fileNo).ROIs_x=these_ROIs(mean_x_imps>=thr_x_imps(fileNo));


        sum_all_y_imps=sum(all_y_imps);
        mat_sum_all_y_imps=repmat(sum_all_y_imps,size(all_y_imps,1),1);
        all_y_imps=all_y_imps./mat_sum_all_y_imps; %Normalize to all added imps=1
        mean_y_imps=mean(all_y_imps,2);
        thr_y_imps(fileNo)=prctile(mean_y_imps,percentile_stringency);
        all_mean_y_imps=[all_mean_y_imps; mean_y_imps];
        above_thr_y_imps=[above_thr_y_imps; mean_y_imps>=thr_y_imps(fileNo)];
        imps.file(fileNo).ROIs_y=these_ROIs(mean_y_imps>=thr_y_imps(fileNo));

        %Get Moser analysis results
        load([save_PathMoser arena_file(1:end-4) handles_Moser.save_tag '.mat'])
        % imps.file(fileNo).spatial_rhol1l4=handles_out.spatial_rhol1l4;
        % imps.file(fileNo).delta_center_of_mass=handles_out.delta_center_of_mass;
        % imps.file(fileNo).information_content=handles_out.information_content;
        % imps.file(fileNo).sparsity=handles_out.sparsity;

        all_spatial_rhol1l4=[all_spatial_rhol1l4; handles_out.spatial_rhol1l4];
        all_delta_center_of_mass=[all_delta_center_of_mass; handles_out.delta_center_of_mass];
        all_information_content=[all_information_content; handles_out.information_content];
        all_sparsity=[all_sparsity; handles_out.sparsity];

        imps.file(fileNo).spatial_rhol1l4=handles_out.spatial_rhol1l4;
        imps.file(fileNo).delta_center_of_mass=handles_out.delta_center_of_mass;
        imps.file(fileNo).information_content=handles_out.information_content;
        imps.file(fileNo).information_contentl1=handles_out.information_contentl1;
        imps.file(fileNo).information_contentl4=handles_out.information_contentl4;
        imps.file(fileNo).sparsity=handles_out.sparsity;

        conc_spatial_rhol1l4=[conc_spatial_rhol1l4; handles_out.spatial_rhol1l4((mean_conc_imps>=thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];
        conc_delta_center_of_mass=[conc_delta_center_of_mass; handles_out.delta_center_of_mass((mean_conc_imps>=thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];
        conc_information_content=[conc_information_content; handles_out.information_content((mean_conc_imps>=thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];
        conc_sparsity=[conc_sparsity; handles_out.sparsity((mean_conc_imps>=thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];

        x_spatial_rhol1l4=[x_spatial_rhol1l4; handles_out.spatial_rhol1l4((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps>=thr_x_imps(fileNo)))];
        x_delta_center_of_mass=[x_delta_center_of_mass; handles_out.delta_center_of_mass((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps>=thr_x_imps(fileNo)))];
        x_information_content=[x_information_content; handles_out.information_content((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps>=thr_x_imps(fileNo)))];
        x_sparsity=[x_sparsity; handles_out.sparsity((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps<thr_y_imps(fileNo))&(mean_x_imps>=thr_x_imps(fileNo)))];

        y_spatial_rhol1l4=[y_spatial_rhol1l4; handles_out.spatial_rhol1l4((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps>=thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];
        y_delta_center_of_mass=[y_delta_center_of_mass; handles_out.delta_center_of_mass((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps>=thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];
        y_information_content=[y_information_content; handles_out.information_content((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps>=thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];
        y_sparsity=[y_sparsity; handles_out.sparsity((mean_conc_imps<thr_conc_imps(fileNo))&(mean_y_imps>=thr_y_imps(fileNo))&(mean_x_imps<thr_x_imps(fileNo)))];

        pfft=1;

        %Let's look at importance relationships between the different decoders

        %x vs conc
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

        max_imp=max([max(mean_x_imps) max(mean_conc_imps)]);
        min_imp_non_zero=min([min(mean_x_imps(mean_x_imps~=0)) min(mean_conc_imps(mean_conc_imps~=0))]);
        min_imp=0.9*min_imp_non_zero;
        mean_x_imps(mean_x_imps==0)=min_imp;
        mean_conc_imps(mean_conc_imps==0)=min_imp;

        plot(log10(mean_x_imps),log10(mean_conc_imps),'ob','MarkerFaceColor','b');
        plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
        xlim([log10(min_imp) log10(max_imp)])
        ylim([log10(min_imp) log10(max_imp)])
        plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
        plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
        xlabel('Importance x (log10)')
        ylabel('Importance conc (log10)')
        title(['Prediction importance for ' arena_file(1:14)], 'Interpreter', 'none')

        %y vs conc
        try
            close(figureNo+1)
        catch
        end
        hFig=figure(figureNo+1);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

        max_imp=max([max(mean_y_imps) max(mean_conc_imps)]);
        min_imp_non_zero=min([min(mean_y_imps(mean_y_imps~=0)) min(mean_conc_imps(mean_conc_imps~=0))]);
        min_imp=0.9*min_imp_non_zero;
        mean_y_imps(mean_y_imps==0)=min_imp;
        mean_conc_imps(mean_conc_imps==0)=min_imp;

        plot(log10(mean_y_imps),log10(mean_conc_imps),'ob','MarkerFaceColor','b')
        plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
        xlim([log10(min_imp) log10(max_imp)])
        ylim([log10(min_imp) log10(max_imp)])
        plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
        plot([log10(min_imp) log10(max_imp)],[log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],'-k')

        xlabel('Importance y (log10)')
        ylabel('Importance conc (log10)')
        title(['Prediction importance for ' arena_file(1:14)], 'Interpreter', 'none')

        %x vs y
        try
            close(figureNo+2)
        catch
        end
        hFig=figure(figureNo+2);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

        max_imp=max([max(mean_y_imps) max(mean_x_imps)]);
        min_imp_non_zero=min([min(mean_y_imps(mean_y_imps~=0)) min(mean_x_imps(mean_x_imps~=0))]);
        min_imp=0.9*min_imp_non_zero;
        mean_y_imps(mean_y_imps==0)=min_imp;
        mean_x_imps(mean_x_imps==0)=min_imp;

        plot(log10(mean_x_imps),log10(mean_y_imps),'ob','MarkerFaceColor','b')
        plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
        xlim([log10(min_imp) log10(max_imp)])
        ylim([log10(min_imp) log10(max_imp)])
        plot([log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
        plot([log10(min_imp) log10(max_imp)],[log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],'-k')
        xlabel('Importance x (log10)')
        ylabel('Importance y (log10)')
        title(['Prediction importance for ' arena_file(1:14)], 'Interpreter', 'none')

        pffft=1;
    end
end
figureNo=figureNo+2;

mean_x_imps_x95=all_mean_x_imps(logical(above_thr_x_imps));
mean_x_imps_y95=all_mean_x_imps(logical(above_thr_y_imps));
mean_x_imps_conc95=all_mean_x_imps(logical(above_thr_conc_imps));

mean_y_imps_x95=all_mean_y_imps(logical(above_thr_x_imps));
mean_y_imps_y95=all_mean_y_imps(logical(above_thr_y_imps));
mean_y_imps_conc95=all_mean_y_imps(logical(above_thr_conc_imps));

mean_conc_imps_x95=all_mean_conc_imps(logical(above_thr_x_imps));
mean_conc_imps_y95=all_mean_conc_imps(logical(above_thr_y_imps));
mean_conc_imps_conc95=all_mean_conc_imps(logical(above_thr_conc_imps));

mean_conc_imps_xconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));
mean_x_imps_xconc95=all_mean_x_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));

mean_conc_imps_yconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));
mean_y_imps_yconc95=all_mean_y_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));

mean_y_imps_xy95=all_mean_y_imps(logical(above_thr_x_imps)&logical(above_thr_y_imps));
mean_x_imps_xy95=all_mean_x_imps(logical(above_thr_x_imps)&logical(above_thr_y_imps));

% Create a sample data matrix for PCA and t-SNE
% Keep track of which are >95% for the different dimensions
mean_x_imps=[];
mean_y_imps=[];
mean_conc_imps=[];
ii_95_class=[]; %1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s
ii_included=0;

for ii=1:length(all_mean_x_imps)
    include_ii=0;

    if above_thr_x_imps(ii)==1
        this_ii_95_class=1;
        include_ii=1;
    end

    if above_thr_y_imps(ii)==1
        this_ii_95_class=2;
        include_ii=1;
    end

    if above_thr_conc_imps(ii)==1
        this_ii_95_class=3;
        include_ii=1;
    end

    if (above_thr_x_imps(ii)==1)&(above_thr_conc_imps(ii)==1)
        this_ii_95_class=4;
        include_ii=1;
    end

    if (above_thr_y_imps(ii)==1)&(above_thr_conc_imps(ii)==1)
        this_ii_95_class=5;
        include_ii=1;
    end

    if (above_thr_x_imps(ii)==1)&(above_thr_y_imps(ii)==1)
        this_ii_95_class=6;
        include_ii=1;
    end

    if (above_thr_x_imps(ii)==1)&(above_thr_conc_imps(ii)==1)&(above_thr_y_imps(ii)==1)
        this_ii_95_class=7;
        include_ii=1;
    end

    if isinf(log10(all_mean_x_imps(ii)))
        include_ii=0;
    end

    if isinf(log10(all_mean_y_imps(ii)))
        include_ii=0;
    end

    if isinf(log10(all_mean_conc_imps(ii)))
        include_ii=0;
    end

    if include_ii==1
        ii_included=ii_included+1;
        ii_95_class(ii_included)=this_ii_95_class;
        mean_x_imps(ii_included)=log10(all_mean_x_imps(ii));
        mean_y_imps(ii_included)=log10(all_mean_y_imps(ii));
        mean_conc_imps(ii_included)=log10(all_mean_conc_imps(ii));
    end

end

data = [mean_x_imps; mean_y_imps; mean_conc_imps];

%Now plot x vs conc for all above thr
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

plot(log10(mean_x_imps_x95),log10(mean_conc_imps_x95),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0]);
plot(log10(mean_x_imps_conc95),log10(mean_conc_imps_conc95),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5]);
plot(log10(mean_x_imps_xconc95),log10(mean_conc_imps_xconc95),'o','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor',[0.95 0.9 0.25]);
plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
xlim([log10(min_imp) log10(max_imp)])
ylim([log10(min_imp) log10(max_imp)])
plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
xlabel('Importance x (log10)')
ylabel('Importance conc (log10)')
title(['Prediction importance above 95 percentile '], 'Interpreter', 'none')

%Now plot y vs conc for all above thr
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

plot(log10(mean_y_imps_y95),log10(mean_conc_imps_y95),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9]);
plot(log10(mean_y_imps_conc95),log10(mean_conc_imps_conc95),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5]);
plot(log10(mean_y_imps_yconc95),log10(mean_conc_imps_yconc95),'o','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7]);
plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
xlim([log10(min_imp) log10(max_imp)])
ylim([log10(min_imp) log10(max_imp)])
plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
plot([log10(min_imp) log10(max_imp)],[log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],'-k')
xlabel('Importance y (log10)')
ylabel('Importance conc (log10)')
title(['Prediction importance above 95 percentile ' ], 'Interpreter', 'none')



%Now plot x vs y for all above thr
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])



plot(log10(mean_x_imps_x95),log10(mean_y_imps_x95),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0]);
plot(log10(mean_x_imps_y95),log10(mean_y_imps_y95),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9]);
plot(log10(mean_x_imps_xy95),log10(mean_y_imps_xy95),'o','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0]);
plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
xlim([log10(min_imp) log10(max_imp)])
ylim([log10(min_imp) log10(max_imp)])
plot([log10(thr_y_imps(fileNo)) log10(thr_y_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
xlabel('Importance x (log10)')
ylabel('Importance y (log10)')
title(['Prediction importance above 95 percentile ' ], 'Interpreter', 'none')




% Step 1 Transpose the data if necessary (N samples x 3 variables)
data = data';

rng('default') % for reproducibility

% Step 2: Run t-SNE
% The default parameters are usually sufficient, but you can adjust them
% For example, set 'NumPCAComponents' to reduce dimensionality before t-SNE
% mappedX = tsne(data, 'NumPCAComponents', 2, 'Perplexity', 30);
% mappedX = tsne(data,'Algorithm','exact','Distance','mahalanobis'); %Works well
mappedX = tsne(data,'Algorithm','exact','Distance','cosine'); %works better

% mappedX = tsne(data,'Algorithm','exact','Distance','chebychev'); %Ok
% mappedX = tsne(data,'Algorithm','exact','Distance','euclidean'); %OK

% Step 3: Plot the results
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
for ii=1:size(data,1)
    %1 is >=x95, 2 is >=y95, 3 is >=conc95, 4 is >=x95 and conc95, 5 is >=y95 and conc95, 6 is >=x95 and y95, 7 is >= all 95s
    switch ii_95_class(ii)
        case 1
            %1 is >=x95
            plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.9 0.6 0],'MarkerEdgeColor',[0.9 0.6 0]);
        case 2
            %2 is >=y95
            plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.35 0.7 0.9],'MarkerEdgeColor',[0.35 0.7 0.9]);
        case 3
            %3 is >=conc95
            plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0 0.6 0.5],'MarkerEdgeColor',[0 0.6 0.5]);
        case 4
            %4 is >=x95 and conc95
            plot(mappedX(ii,1),mappedX(ii,2),'s','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor','k');
        case 5
            %5 is >=y95 and conc95
            plot(mappedX(ii,1),mappedX(ii,2),'s','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor','k');
        case 6
            %6 is >=x95 and y95
            plot(mappedX(ii,1),mappedX(ii,2),'s','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor','k');
        case 7
            %7 is >= all 95s
            plot(mappedX(ii,1),mappedX(ii,2),'s','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor','k');
    end
end

these_ylim=ylim;
these_xlim=xlim;
text(these_xlim(1)+0.10*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.95*(these_ylim(2)-these_ylim(1)),'x','Color',[0.9 0.6 0],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.80*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.92*(these_ylim(2)-these_ylim(1)),'y','Color',[0.35 0.7 0.9],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.25*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.1*(these_ylim(2)-these_ylim(1)),'odor','Color',[0 0.6 0.5],'FontWeight','bold','FontSize',16)
xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
title(['t-SNE Prediction Importance ' ]);
pffft=1;
% grid on;

%Let's look at the Moser/spatial information parameters
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

[f_aic,x_aic] = drg_ecdf(conc_spatial_rhol1l4);
plot(x_aic,f_aic,'Color',[0 0.6 0.5],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(x_spatial_rhol1l4);
plot(x_aic,f_aic,'Color',[0.9 0.6 0],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(y_spatial_rhol1l4);
plot(x_aic,f_aic,'Color',[0.35 0.7 0.9],'LineWidth',3)

title('Spatial correlation')
xlabel('Correlation')

these_ylim=ylim;
these_xlim=xlim;
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.3*(these_ylim(2)-these_ylim(1)),'odor','Color',[0 0.6 0.5],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.2*(these_ylim(2)-these_ylim(1)),'x','Color',[0.9 0.6 0],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.1*(these_ylim(2)-these_ylim(1)),'y','Color',[0.35 0.7 0.9],'FontWeight','bold','FontSize',16)


figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

[f_aic,x_aic] = drg_ecdf(conc_delta_center_of_mass);
plot(x_aic,f_aic,'Color',[0 0.6 0.5],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(x_delta_center_of_mass);
plot(x_aic,f_aic,'Color',[0.9 0.6 0],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(y_delta_center_of_mass);
plot(x_aic,f_aic,'Color',[0.35 0.7 0.9],'LineWidth',3)

title('Delta center of mass')
xlabel('Delta center of mass')

these_ylim=ylim;
these_xlim=xlim;
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.3*(these_ylim(2)-these_ylim(1)),'odor','Color',[0 0.6 0.5],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.2*(these_ylim(2)-these_ylim(1)),'x','Color',[0.9 0.6 0],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.1*(these_ylim(2)-these_ylim(1)),'y','Color',[0.35 0.7 0.9],'FontWeight','bold','FontSize',16)


figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

[f_aic,x_aic] = drg_ecdf(conc_information_content);
plot(x_aic,f_aic,'Color',[0 0.6 0.5],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(x_information_content);
plot(x_aic,f_aic,'Color',[0.9 0.6 0],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(y_information_content);
plot(x_aic,f_aic,'Color',[0.35 0.7 0.9],'LineWidth',3)

title('SSI')
xlabel('SSI')

these_ylim=ylim;
these_xlim=xlim;
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.3*(these_ylim(2)-these_ylim(1)),'odor','Color',[0 0.6 0.5],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.2*(these_ylim(2)-these_ylim(1)),'x','Color',[0.9 0.6 0],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.1*(these_ylim(2)-these_ylim(1)),'y','Color',[0.35 0.7 0.9],'FontWeight','bold','FontSize',16)



figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

[f_aic,x_aic] = drg_ecdf(conc_sparsity);
plot(x_aic,f_aic,'Color',[0 0.6 0.5],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(x_sparsity);
plot(x_aic,f_aic,'Color',[0.9 0.6 0],'LineWidth',3)

[f_aic,x_aic] = drg_ecdf(y_sparsity);
plot(x_aic,f_aic,'Color',[0.35 0.7 0.9],'LineWidth',3)

title('Sparsity y vs all (blue)')
xlabel('Sparsity')

these_ylim=ylim;
these_xlim=xlim;
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.3*(these_ylim(2)-these_ylim(1)),'odor','Color',[0 0.6 0.5],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.2*(these_ylim(2)-these_ylim(1)),'x','Color',[0.9 0.6 0],'FontWeight','bold','FontSize',16)
text(these_xlim(1)+0.50*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.1*(these_ylim(2)-these_ylim(1)),'y','Color',[0.35 0.7 0.9],'FontWeight','bold','FontSize',16)

%Now plot the space activity maps for each of the high prediction
%importance ROIs
figNo=figureNo+1;
x=25:50:475;
y=24:48:456;
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Get XY and dFF per trial
        arena_file=handles_XY.arena_file{fileNo};
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        no_neurons=handles_out.no_neurons;


        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        angles=handles_out.angles;

        these_important_ROIs=unique([imps.file(fileNo).ROIs_conc imps.file(fileNo).ROIs_x imps.file(fileNo).ROIs_y]);

        %Now show the spatial activity maps
        no_place_cells=0;
        place_cells=[];
        no_lane_trial_cells=0;
        lane_trial_cells=[];
        for ii_ROI=1:length(these_important_ROIs)

            this_ROI=these_important_ROIs(ii_ROI);
            spatial_rhol1l4=imps.file(fileNo).spatial_rhol1l4;
            delta_center_of_mass=imps.file(fileNo).delta_center_of_mass;
            SSI=imps.file(fileNo).information_content;
            SSIl1=imps.file(fileNo).information_contentl1;
            SSIl4=imps.file(fileNo).information_contentl4;
            sparsity=imps.file(fileNo).sparsity;

            %Initialize variables
            this_dFF_activity=zeros(10,10);
            this_dFF_activity_n=zeros(10,10);
            sum_dFF_activity=0;

            this_dFFl1_activity=zeros(10,10);
            this_dFFl1_activity_n=zeros(10,10);
            sum_dFFl1_activity=0;

            this_dFFl4_activity=zeros(10,10);
            this_dFFl4_activity_n=zeros(10,10);
            sum_dFFl4_activity=0;


            for trNo=1:trials.odor_trNo
                these_x=trials.trial(trNo).XYtest(:,1);
                these_y=trials.trial(trNo).XYtest(:,2);
                these_dFF=trials.trial(trNo).XdFFtest(:,these_important_ROIs(ii_ROI));

                for ii_t=1:length(these_x)
                    this_x_ii=ceil(these_x(ii_t)/50);
                    if this_x_ii==11
                        this_x_ii=10;
                    end

                    this_y_ii=ceil(these_y(ii_t)/48);
                    if this_y_ii==11
                        this_y_ii=10;
                    end

                    this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                    this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
                    sum_dFF_activity=sum_dFF_activity+these_dFF(ii_t);

                    if trials.lane_per_trial(trNo)==1
                        this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
                        sum_dFFl1_activity=sum_dFFl1_activity+these_dFF(ii_t);
                    else
                        this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_dFF(ii_t);
                        this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
                        sum_dFFl4_activity=sum_dFFl4_activity+these_dFF(ii_t);
                    end
                end
            end

            for ii_x=1:10
                for ii_y=1:10
                    if this_dFF_activity_n(ii_x,ii_y)~=0
                        this_dFF_activity(ii_x,ii_y)=this_dFF_activity(ii_x,ii_y)/this_dFF_activity_n(ii_x,ii_y);
                    end
                end
            end

            for ii_x=1:10
                for ii_y=1:10
                    if this_dFFl1_activity_n(ii_x,ii_y)~=0
                        this_dFFl1_activity(ii_x,ii_y)=this_dFFl1_activity(ii_x,ii_y)/this_dFFl1_activity_n(ii_x,ii_y);
                    end
                end
            end

            for ii_x=1:10
                for ii_y=1:10
                    if this_dFFl4_activity_n(ii_x,ii_y)~=0
                        this_dFFl4_activity(ii_x,ii_y)=this_dFFl4_activity(ii_x,ii_y)/this_dFFl4_activity_n(ii_x,ii_y);
                    end
                end
            end

            %Plot the space activity maps
            try
                close(figNo)
            catch
            end


            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .75 .75])



            % %Plot dFF timecourses
            % subplot(2,3,1:3)
            % hold on
            % ii_plot=0;

            %Find turn points and time vectors
            all_ii_turns=zeros(1,trials.odor_trNo);
            all_ii_ends=zeros(1,trials.odor_trNo);
            include_trial=zeros(1,trials.odor_trNo);
            for trNo=1:trials.odor_trNo
                this_ii_last_turn=find(angles.trial(trNo).delta_x>=100,1,'last');
                if ~isempty(this_ii_last_turn)
                    all_ii_turns(trNo)=angles.trial(trNo).ii_turns(this_ii_last_turn);
                    all_ii_ends(trNo)= size(trials.trial(trNo).XYtest,1);
                    include_trial(trNo)=1;
                else
                    pffft=1;
                end
            end

            delta_below_zero_ii=max(all_ii_turns(include_trial==1));
            delta_above_zero_ii=max(all_ii_ends(include_trial==1)-all_ii_turns(include_trial==1));
            time_bins=handles_XY.dt*([1:delta_below_zero_ii+delta_above_zero_ii]-delta_below_zero_ii);

            hit1_dFF=[];
            ii_hit1=0;
            miss1_dFF=[];
            ii_miss1=0;
            hit4_dFF=[];
            ii_hit4=0;
            miss4_dFF=[];
            ii_miss4=0;

            for trNo=1:trials.odor_trNo
                if include_trial(trNo)==1
                    these_dFF=trials.trial(trNo).XdFFtest(:,these_important_ROIs(ii_ROI));
                    %Okabe_Ito colors
                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            ii_hit1=ii_hit1+1;
                            hit1_dFF(ii_hit1,1:length(time_bins))=min(these_dFF);
                            hit1_dFF(ii_hit1,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                        case 2
                            %Lane 1 miss orange
                            ii_miss1=ii_miss1+1;
                            miss1_dFF(ii_miss1,1:length(time_bins))=min(these_dFF);
                            miss1_dFF(ii_miss1,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                        case 3
                            %Lane 4 hit blue
                            ii_hit4=ii_hit4+1;
                            hit4_dFF(ii_hit4,1:length(time_bins))=min(these_dFF);
                            hit4_dFF(ii_hit4,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                        case 4
                            %Lane 4 miss sky blue
                            ii_miss4=ii_miss4+1;
                            miss4_dFF(ii_miss4,1:length(time_bins))=min(these_dFF);
                            miss4_dFF(ii_miss4,delta_below_zero_ii-all_ii_turns(trNo)+1:delta_below_zero_ii+(all_ii_ends(trNo)-all_ii_turns(trNo)))=these_dFF;
                    end
                end
            end

            %Plot dFFs for lane 1 and 4
            y_gap=2;
            hit_miss_1=zeros(size(hit1_dFF,1)+size(miss1_dFF,1)+y_gap,delta_below_zero_ii+delta_above_zero_ii);
            hit_miss_1(size(miss1_dFF,1)+y_gap+1:size(hit1_dFF,1)+size(miss1_dFF,1)+y_gap,:)=hit1_dFF;
            hit_miss_1(1:size(miss1_dFF,1),:)=miss1_dFF;

            y_trials1=[1:size(hit_miss_1,1)];

            hit_miss_4=zeros(size(hit4_dFF,1)+size(miss4_dFF,1)+y_gap,delta_below_zero_ii+delta_above_zero_ii);
            hit_miss_4(size(miss4_dFF,1)+y_gap+1:size(hit4_dFF,1)+size(miss4_dFF,1)+y_gap,:)=hit4_dFF;
            hit_miss_4(1:size(miss4_dFF,1),:)=miss4_dFF;

            y_trials4=[1:size(hit_miss_4,1)];

            y_trials_end=max([y_trials1(end) y_trials4(end)])

            %Plot dFF for lane 1
            subplot(2, 6, [1 2 3]);
            hold on
            
            drg_pcolor(repmat(time_bins,length(y_trials1),1),repmat(y_trials1,length(time_bins),1)',hit_miss_1)
            colormap(this_cmap)
            clim([0 max(hit_miss_1(:))])
            shading flat
            plot([0 0],[y_trials1(1) y_trials1(end)+1],'-w','LineWidth',3)

            %y_trials1=1 to 2 is the first miss
            %size(miss1_dFF,1) to size(miss1_dFF,1) + 1 is the last miss,
            %size(miss1_dFF,1)+1 to size(miss1_dFF,1)+1+y_gap is the in between miss
            %and hit
            %size(miss1_dFF,1)+1+y_gap to size(miss1_dFF,1)+1+y_gap+size(hit1_dFF,1) is
            %hit

            rectangle('Position', [time_bins(1), size(miss1_dFF,1)+1, 1.03*(time_bins(end)-time_bins(1)), y_gap], ... % [x, y, width, height]
                'FaceColor', 'white', ...    % Fill color (white)
                'EdgeColor', 'none');        % No border for the rectangle

            ylim([0 y_trials_end+2])
            xlim([time_bins(1)-0.05*(time_bins(end)-time_bins(1)) time_bins(end)+0.05*(time_bins(end)-time_bins(1)) ])
            xlabel('Time (sec)')
            yticks([1+(size(miss1_dFF,1)/2) size(miss1_dFF,1)+1+y_gap+(size(hit1_dFF,1)/2)])
            yticklabels({'Misses','Hits'})
            title(['Lane 1 dFF'])

            %Plot dFF for lane 4
            subplot(2, 6, [4 5 6]);
            hold on
            
            drg_pcolor(repmat(time_bins,length(y_trials4),1),repmat(y_trials4,length(time_bins),1)',hit_miss_4)
            colormap(this_cmap)
            clim([0 max(hit_miss_4(:))])
            shading flat
            plot([0 0],[y_trials4(1) y_trials4(end)+1],'-w','LineWidth',3)

            rectangle('Position', [time_bins(1), size(miss4_dFF,1)+1, 1.03*(time_bins(end)-time_bins(1)), y_gap], ... % [x, y, width, height]
                'FaceColor', 'white', ...    % Fill color (white)
                'EdgeColor', 'none');        % No border for the rectangle
            ylim([0 y_trials_end+2])
            xlim([time_bins(1)-0.05*(time_bins(end)-time_bins(1)) time_bins(end)+0.05*(time_bins(end)-time_bins(1)) ])
            xlabel('Time (sec)')
            yticks([1+(size(miss4_dFF,1)/2) size(miss4_dFF,1)+1+y_gap+(size(hit4_dFF,1)/2)])
            yticklabels({'Misses','Hits'})
            title(['Lane 4 dFF'])

            %Plot space activity maps as in Moser https://doi.org/10.1126/science.1114037
            max_this_dFF_activity=max(this_dFF_activity(:));
            max_this_dFFl4_activity=max(this_dFFl4_activity(:));
            max_this_dFFl1_activity=max(this_dFFl1_activity(:));

            colormap fire
            this_cmap=colormap;
            this_cmap(1,:)=[0.3 0.3 0.3];

            if max_this_dFFl1_activity>max_this_dFFl4_activity
                %Lane 1
                subplot(2, 6, [7 8]);
                max_activity=max_this_dFFl1_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_dFFl1_activity=this_dFFl1_activity;
                this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks(0:48:480)
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 1, SSI= ' num2str(SSIl1(this_ROI))])

                %Lane 4 normalized to lane 1
                subplot(2, 6, [9 10]);
                max_activity=max_this_dFFl1_activity;
                this_masked_dFFl4_activity=this_dFFl4_activity;
                this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks(0:48:480)
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 4'])

                %Lane 4 normalized to lane 4
                subplot(2, 6, [11 12]);
                max_activity=max_this_dFFl4_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_dFFl4_activity=this_dFFl4_activity;
                this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks(0:48:480)
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 4, SSI= ' num2str(SSIl4(this_ROI))])

            else


                %Lane 4 normalized to lane 4
                subplot(2, 6, [7 8]);
                max_activity=max_this_dFFl4_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.000001;
                end
                this_masked_dFFl4_activity=this_dFFl4_activity;
                this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks(0:48:480)
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 4, SSI= ' num2str(SSIl4(this_ROI))])

                %Lane 1 normalized to lane 4
                subplot(2, 6, [9 10]);
                max_activity=max_this_dFFl4_activity;
                this_masked_dFFl1_activity=this_dFFl1_activity;
                this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks(0:48:480)
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 1'])

                %Lane 1 normalized to lane 1
                subplot(2, 6, [11 12]);
                max_activity=max_this_dFFl1_activity;
                delta_ac=max_activity/255;
                if delta_ac==0
                    delta_ac=0.00001;
                end
                this_masked_dFFl1_activity=this_dFFl1_activity;
                this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
                drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
                colormap(this_cmap)
                clim([-1.5*delta_ac max_activity])          
                shading interp
                set(gca, 'YDir', 'reverse');

                yticks(0:48:480)
                xticks(0:50:500)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['Lane 1, SSI= ' num2str(SSIl1(this_ROI))])
            end


            sgt_legend=['dFF map file  No ' num2str(fileNo) ' ROI No ' num2str(this_ROI) ', rho ' ...
                num2str(spatial_rhol1l4(this_ROI)) ', dcm ' num2str(delta_center_of_mass(this_ROI)),...
                ', SSI ' num2str(SSI(this_ROI)),', PI = '];

            if sum(this_ROI==imps.file(fileNo).ROIs_conc)>0
                sgt_legend=[sgt_legend ' odor '];
            end

            if sum(this_ROI==imps.file(fileNo).ROIs_x)>0
                sgt_legend=[sgt_legend ' x '];
            end

            if sum(this_ROI==imps.file(fileNo).ROIs_y)>0
                sgt_legend=[sgt_legend ' y '];
            end

            sgtitle(sgt_legend)
   
            pffft=1;

        end
    end
end
fclose(fileID);

pffft=1;

