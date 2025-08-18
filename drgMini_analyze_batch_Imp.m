%drgMini_analyze_batch_Imp
%
%Gathers all the importance information into outputPredictionImportance.mat
%

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

         % %Trained with hits only
         % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
         % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m';

          %Trained with hits only taking on account when mouse detects the odor
        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc04192024/';
        choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_04192024.m'

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';
        
        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        save_PathPredImp='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FilePredImp='outputPredictionImportancev2.mat';
        % 
        % %The output of drgMini_information_contentv2 is saved here
        % save_PathIC='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        % save_FileIC='outputPerROIInformationContent.mat';

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
all_imps_ROI=[];
above_thr_x_imps=[];
all_mean_y_imps=[];
above_thr_y_imps=[];
thr_x_imps=zeros(1,length(handles_conc.arena_file));
thr_y_imps=zeros(1,length(handles_conc.arena_file));
thr_conc_imps=zeros(1,length(handles_conc.arena_file));
all_spatial_rhol1l4=[];
all_delta_center_of_mass=[];
all_information_content=[];
all_information_contentl1=[];
all_information_contentl4=[];
all_information_content_sh=[];
all_information_contentl1_sh=[];
all_information_contentl4_sh=[];
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

conc_ssi_all_info_lane1=[];
conc_ssi_all_info_lane4=[];
conc_ssi_all_info_both_lanes=[];
conc_ssi_all_info_op_bin=[];
conc_idxssi=[];

x_ssi_all_info_lane1=[];
x_ssi_all_info_lane4=[];
x_ssi_all_info_both_lanes=[];
x_ssi_all_info_op_bin=[];
x_idxssi=[];

y_ssi_all_info_lane1=[];
y_ssi_all_info_lane4=[];
y_ssi_all_info_both_lanes=[];
y_ssi_all_info_op_bin=[];
y_idxssi=[];


 
file_numbers=[];


figureNo=figureNo+1;
these_R1s=[];
imps=[];


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
        fprintf(1,['File No ' num2str(fileNo) '\n'])
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
        all_imps_ROI=[all_imps_ROI these_ROIs];

        %Save the data in imps
        imps.file(fileNo).ROIs_conc=these_ROIs(mean_conc_imps>=thr_conc_imps(fileNo));
        imps.thr_conc_imps=thr_conc_imps;
        imps.all_mean_conc_imps=all_mean_conc_imps;
        imps.all_imps_ROI=all_imps_ROI;
        imps.file_numbers=file_numbers;
 
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

        imps.thr_x_imps=thr_x_imps;
        imps.all_mean_x_imps=all_mean_x_imps;

        sum_all_y_imps=sum(all_y_imps);
        mat_sum_all_y_imps=repmat(sum_all_y_imps,size(all_y_imps,1),1);
        all_y_imps=all_y_imps./mat_sum_all_y_imps; %Normalize to all added imps=1
        mean_y_imps=mean(all_y_imps,2);
        thr_y_imps(fileNo)=prctile(mean_y_imps,percentile_stringency);
        all_mean_y_imps=[all_mean_y_imps; mean_y_imps];
        above_thr_y_imps=[above_thr_y_imps; mean_y_imps>=thr_y_imps(fileNo)];
        imps.file(fileNo).ROIs_y=these_ROIs(mean_y_imps>=thr_y_imps(fileNo));

        imps.thr_y_imps=thr_y_imps;
        imps.all_mean_y_imps=all_mean_y_imps;

        %Get Moser analysis results
        load([save_PathMoser arena_file(1:end-4) handles_Moser.save_tag '.mat'])
        % imps.file(fileNo).spatial_rhol1l4=handles_out.spatial_rhol1l4;
        % imps.file(fileNo).delta_center_of_mass=handles_out.delta_center_of_mass;
        % imps.file(fileNo).information_content=handles_out.information_content;
        % imps.file(fileNo).sparsity=handles_out.sparsity;

        all_spatial_rhol1l4=[all_spatial_rhol1l4; handles_out.spatial_rhol1l4];
        all_delta_center_of_mass=[all_delta_center_of_mass; handles_out.delta_center_of_mass];
        all_information_content=[all_information_content; handles_out.information_content];
        all_information_contentl1=[all_information_contentl1; handles_out.information_contentl1];
        all_information_contentl4=[all_information_contentl4; handles_out.information_contentl4];
        all_information_content_sh=[all_information_content_sh; handles_out.sh_information_content(:)];
        all_information_contentl1_sh=[all_information_contentl1_sh; handles_out.sh_information_contentl1(:)];
        all_information_contentl4_sh=[all_information_contentl4_sh; handles_out.sh_information_contentl4(:)];
        all_sparsity=[all_sparsity; handles_out.sparsity];

        imps.file(fileNo).spatial_rhol1l4=handles_out.spatial_rhol1l4;
        imps.file(fileNo).delta_center_of_mass=handles_out.delta_center_of_mass;
        imps.file(fileNo).information_content=handles_out.information_content;
        imps.file(fileNo).information_contentl1=handles_out.information_contentl1;
        imps.file(fileNo).information_contentl4=handles_out.information_contentl4;
        imps.file(fileNo).sh_information_content=handles_out.sh_information_content;
        imps.file(fileNo).sh_information_contentl1=handles_out.sh_information_contentl1;
        imps.file(fileNo).sh_information_contentl4=handles_out.sh_information_contentl4;
        imps.file(fileNo).sparsity=handles_out.sparsity;
        imps.all_spatial_rhol1l4=all_spatial_rhol1l4;
        imps.all_delta_center_of_mass=all_delta_center_of_mass;
        imps.all_Fusi_SSI=all_information_content;
        imps.all_Fusi_SSIl4=all_information_contentl4;
        imps.all_Fusi_SSIl1=all_information_contentl1;
        imps.all_Fusi_SSI_sh=all_information_content_sh;
        imps.all_Fusi_SSIl4_sh=all_information_contentl4_sh;
        imps.all_Fusi_SSIl1_sh=all_information_contentl1_sh;


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
mean_imp_fileNo=[];
mean_imp_ROI=[];
sig_pred_imp_x=[];
sig_pred_imp_y=[];
sig_pred_imp_conc=[];


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
        mean_imp_fileNo(ii_included)=file_numbers(ii);
        mean_imp_ROI(ii_included)=all_imps_ROI(ii);
        if (this_ii_95_class==1)||(this_ii_95_class==4)||(this_ii_95_class==6)||(this_ii_95_class==7)
            sig_pred_imp_x(ii_included)=1;
        else
            sig_pred_imp_x(ii_included)=0;
        end
        if (this_ii_95_class==2)||(this_ii_95_class==5)||(this_ii_95_class==6)||(this_ii_95_class==7)
            sig_pred_imp_y(ii_included)=1;
        else
            sig_pred_imp_y(ii_included)=0;
        end
        if (this_ii_95_class==3)||(this_ii_95_class==3)||(this_ii_95_class==5)||(this_ii_95_class==7)
            sig_pred_imp_conc(ii_included)=1;
        else
            sig_pred_imp_conc(ii_included)=0;
        end
    end

end
 
imps.ii_95_class=ii_95_class;
imps.mean_x_imps=mean_x_imps;
imps.mean_y_imps=mean_y_imps;
imps.mean_conc_imps=mean_conc_imps;
imps.mean_imp_fileNo=mean_imp_fileNo;
imps.mean_imp_ROI=mean_imp_ROI;
imps.sig_pred_imp_x=sig_pred_imp_x;
imps.sig_pred_imp_y=sig_pred_imp_y;
imps.sig_pred_imp_conc=sig_pred_imp_conc;



save([save_PathPredImp save_FilePredImp],'imps')


pffft=1;

