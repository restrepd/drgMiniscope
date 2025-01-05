%drgMini_analyze_batch_ImpMoserXYConc
close all
clear all

is_sphgpu=0;

switch is_sphgpu
    case 0

        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';

        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        choiceAngleFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end




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

figureNo=0;
% 
% %Let's look at dependence of correlation on bins before for all 2 cm files
% figureNo = figureNo + 1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% this_group=1;
% edges=[0:0.05:1];
% rand_offset=1;
% 
% %R1 for conc all trials
% meanR1_conc=[];
% bins_before=handles_conc.bins_before;
% 
% for ii_run=1:length(handles_conc.bins_before)
%     these_R1s=[];
%     for fileNo=1:length(handles_conc.arena_file)
%         if handles_conc.group(fileNo)==this_group
%             arena_file=handles_conc.arena_file{fileNo};
%             %load the ouptut file
%             load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
%             these_R1s=[these_R1s handles_out.R1.all_trials];
%         end
%     end
% 
%     [mean_out, CIout]=drgViolinPoint(these_R1s...
%         ,edges,bins_before(ii_run),rand_offset,[0.9 0.6 0],[0.9 0.6 0],4);
%     meanR1_conc=[meanR1_conc mean(these_R1s)];
% end
% 
% plot(bins_before,meanR1_conc,'o-','MarkerSize',10,'MarkerEdgeColor',[0.9 0.6 0],'MarkerFaceColor',[0.9 0.6 0],'Color',[0.9 0.6 0])
% 
% %R1 for XY all trials
% meanR1_XY=[];
% bins_before=handles_XY.bins_before;
% 
% for ii_run=1:length(handles_XY.bins_before)
%     these_R1s=[];
%     for fileNo=1:length(handles_XY.arena_file)
%         if handles_XY.group(fileNo)==this_group
%             arena_file=handles_conc.arena_file{fileNo};
%             %load the ouptut file
%             try
%                 load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
%                 these_R1s=[these_R1s handles_out.R1.within_XY];
%             catch
%             end
%         end
%     end
% 
%     [mean_out, CIout]=drgViolinPoint(these_R1s...
%         ,edges,bins_before(ii_run),rand_offset,[0.35 0.7 0.9],[0.35 0.7 0.9],4);
%     meanR1_XY=[meanR1_XY mean(these_R1s)];
% end
% 
% plot(bins_before,meanR1_XY,'o-','MarkerSize',10,'MarkerEdgeColor',[0.35 0.7 0.9],'MarkerFaceColor',[0.35 0.7 0.9],'Color',[0.35 0.7 0.9])
% title('Within trial R1')
% xlabel('Bins before')
% ylabel('R1')
% 
% %For 8 bins before plot R1 for conc vs XY
% ii_run=4;
% 
% figureNo = figureNo + 1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% this_group=1;
% edges=[0:0.05:1];
% rand_offset=1;
% 
% %R1 for conc all trials
% R1_conc=[];
% 
% for fileNo=1:length(handles_conc.arena_file)
%     if handles_conc.group(fileNo)==this_group
%         arena_file=handles_conc.arena_file{fileNo};
%         %load the ouptut file
%         load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
%         R1_conc=[R1_conc handles_out.R1.all_trials];
%     end
% end
% 
% 
% %R1 for XY all trials
% R1_XY=[];
% 
% 
% 
% for fileNo=1:length(handles_XY.arena_file)
%     if handles_XY.group(fileNo)==this_group
%         arena_file=handles_conc.arena_file{fileNo};
%         %load the ouptut file
%         load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
%         R1_XY=[R1_XY handles_out.R1.within_XY];
%     end
% end
% 
% plot(R1_conc,R1_XY,'ok','MarkerFaceColor','k')
% title('Within trial R1')
% xlabel('R1 conc')
% ylabel('R1 XY')
% 
% 
% %Compare groups with fittree, weber frechner and only current data
% %point
% groups=unique(handles_conc.group);
% 
% group_order=[1 4 5 2 3];
% 
% R1type=[];
% R1type{1}='all_trials';
% R1type{2}='between_trials';
% R1type{3}='all_hits';
% R1type{4}='all_miss';
% R1type{5}='all_trials_sh';
% 
% R1type_label=[];
% R1type_label{1}='within';
% R1type_label{2}='between';
% R1type_label{3}='hits';
% R1type_label{4}='miss';
% R1type_label{5}='shuffled';
% 
% %Do the R1 bar graph
% figureNo = figureNo + 1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .6 .3])
% 
% bar_offset=0;
% 
% edges=[0:0.05:1];
% rand_offset=0.5;
% 
% %Plot the different R1s
% for ii_group=1:3
%     for ii_R1type=1:length(R1type)
% 
%         this_group=group_order(ii_group);
% 
%         these_R1s=[];
%         for fileNo=1:length(handles_conc.arena_file)
%             arena_file=handles_conc.arena_file{fileNo};
%             if handles_conc.group(fileNo)==this_group
%                 %load the ouptut file
%                 load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
%                 eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
%                 % these_R1s=[these_R1s handles_out.R1.all_trials];
%             end
% 
%         end
% 
%         %plot bar
%         switch ii_R1type
%             case 1
%                 bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%             case 2
%                 bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%             case 3
%                 bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%             case 4
%                 bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
%             case 5
%                 bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
%         end
% 
%         %Violin plot
%         [mean_out, CIout]=drgViolinPoint(these_R1s...
%             ,edges,bar_offset,rand_offset,'k','k',4);
%         bar_offset=bar_offset+1;
%     end
%     bar_offset=bar_offset+1;
% end
% 
% x_pos=15.3;
% text(x_pos,0.57,'within','Color',[230/255 159/255 0/255])
% text(x_pos,0.53,'between','Color',[86/255 180/255 233/255])
% text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])
% 
% xticks([2 8 14])
% xticklabels({'2 cm','No odor','1 cm','Pre','Odor','Reinf'})
% 
% 
% title(['R1, correlation between decoded odor concentraiton and odor plume mean'])
% ylabel('R1')
% ylim([-0.2 0.65])
% xlim([-1 17])
% 
% %Do the R1 bar graph for 2 cm all trials vs hit, miss
% R1type=[];
% R1type{1}='all_trials';
% R1type{2}='all_hits';
% R1type{3}='all_miss';
% R1type{4}='all_trials_sh';
% 
% R1type_label=[];
% R1type_label{1}='within';
% R1type_label{2}='hits';
% R1type_label{3}='miss';
% R1type_label{4}='shuffled';
% 
% iiR1type_reorder=[2 3 1 4];
% 
% figureNo = figureNo + 1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% bar_offset=0;
% 
% edges=[0:0.05:1];
% rand_offset=0.5;
% 
% glm_r1=[];
% glm_r1_ii=0;
% 
% id_r1_ii=0;
% input_r1_data=[];
% 
% all_R1s=[];
% %Plot the different R1s
% these_groups=[1 5]; % 1 and 2 cm all tirals
% for ii_R1type=1:length(R1type)
% 
%     these_R1s=[];
%     for fileNo=1:length(handles_conc.arena_file)
%         arena_file=handles_conc.arena_file{fileNo};
%         if sum(handles_conc.group(fileNo)==these_groups)>0
%             %load the ouptut file
%             load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
%             eval(['these_R1s=[these_R1s handles_out.R1.' R1type{ii_R1type} '];'])
%             % these_R1s=[these_R1s handles_out.R1.all_trials];
%         end
% 
%     end
%     all_R1s.R1type(ii_R1type).R1=these_R1s;
%     %plot bar
%     switch ii_R1type
%         case 1
%             bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
%         case 2
%             bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
%         case 3
%             bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
%         case 4
%             bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])
%         case 5
%             bar(bar_offset,mean(these_R1s),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255])
%     end
% 
%     %Violin plot
%     [mean_out, CIout]=drgViolinPoint(these_R1s...
%         ,edges,bar_offset,rand_offset,'k','k',4);
%     bar_offset=bar_offset+1;
% 
%     glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=these_R1s;
%     glm_r1.trial_type(glm_r1_ii+1:glm_r1_ii+length(these_R1s))=iiR1type_reorder(ii_R1type)*ones(1,length(these_R1s));
%     glm_r1_ii=glm_r1_ii+length(these_R1s);
% 
%     id_r1_ii=id_r1_ii+1;
%     input_r1_data(id_r1_ii).data=these_R1s;
%     input_r1_data(id_r1_ii).description=[R1type_label{ii_R1type}];
% 
% end
% bar_offset=bar_offset+1;
% 
% %Plot lines between points
% for ii=1:length(these_R1s)
%     these_R1s=[];
%     for ii_R1type=1:length(R1type)
%         these_R1s=[these_R1s all_R1s.R1type(ii_R1type).R1(ii)]
%     end
%     plot([0:length(R1type)-1],these_R1s,'-','Color',[0.7 0.7 0.7])
% end
% 
% % x_pos=3;
% % text(x_pos,0.53,'within','Color',[230/255 159/255 0/255])
% % text(x_pos,0.49,'hit','Color',[0/255 158/255 115/255])
% % text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% % text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])
% 
% xticks([0 1 2 3])
% xticklabels({'within','hit','miss','shuffled'})
% 
% 
% title(['R1 for predcition of odor concentration'])
% ylabel('R1')
% ylim([-0.2 0.65])
% xlim([-1 4])
% 
% %Perform the glm  for errors
% fprintf(1, ['\nglm for trial type\n'])
% fprintf(fileID, ['\nglm for trial type\n']);
% %
% % fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% % fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);
% 
% tbl = table(glm_r1.data',glm_r1.trial_type',...
%     'VariableNames',{'R1','trial_type'});
% mdl = fitglm(tbl,'R1~trial_type'...
%     ,'CategoricalVars',[2])
% 
% txt = evalc('mdl');
% txt=regexp(txt,'<strong>','split');
% txt=cell2mat(txt);
% txt=regexp(txt,'</strong>','split');
% txt=cell2mat(txt);
% 
% fprintf(fileID,'%s\n', txt);
% 
% 
% %Do the ranksum/t-test
% fprintf(1, ['\n\nRanksum or t-test p values for R1 vs trial typ\n'])
% fprintf(fileID, ['\n\nRanksum or t-test p values for R1 vs. trial type\n']);
% 
% 
% [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);
% 


%Let's look at imp
for this_group=1:5
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
    file_numbers=[];

    
    figureNo=figureNo+1;
    these_R1s=[];
    for fileNo=1:length(handles_conc.arena_file)
        if handles_conc.group(fileNo)==this_group
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
            thr_conc_imps(fileNo)=prctile(mean_conc_imps,95);
            all_mean_conc_imps=[all_mean_conc_imps; mean_conc_imps];
            above_thr_conc_imps=[above_thr_conc_imps; mean_conc_imps>=thr_conc_imps(fileNo)];
            file_numbers=[file_numbers; fileNo*ones(length(mean_conc_imps),1)];

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
            thr_x_imps(fileNo)=prctile(mean_x_imps,95);
            all_mean_x_imps=[all_mean_x_imps; mean_x_imps];
            above_thr_x_imps=[above_thr_x_imps; mean_x_imps>=thr_x_imps(fileNo)];


            sum_all_y_imps=sum(all_y_imps);
            mat_sum_all_y_imps=repmat(sum_all_y_imps,size(all_y_imps,1),1);
            all_y_imps=all_y_imps./mat_sum_all_y_imps; %Normalize to all added imps=1
            mean_y_imps=mean(all_y_imps,2);
            thr_y_imps(fileNo)=prctile(mean_y_imps,95);
            all_mean_y_imps=[all_mean_y_imps; mean_y_imps];
            above_thr_y_imps=[above_thr_y_imps; mean_y_imps>=thr_y_imps(fileNo)];

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
    %
    % %Now plot x vs conc for all above thr
    % try
    %     close(figureNo)
    % catch
    % end
    % hFig=figure(figureNo);
    % hold on
    %
    % ax=gca;ax.LineWidth=3;
    % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    %
    % mean_x_imps_x95=all_mean_x_imps(logical(above_thr_x_imps));
    % mean_conc_imps_x95=all_mean_conc_imps(logical(above_thr_x_imps));
    %
    % mean_conc_imps_conc95=all_mean_conc_imps(logical(above_thr_conc_imps));
    % mean_x_imps_conc95=all_mean_x_imps(logical(above_thr_conc_imps));
    %
    % mean_conc_imps_xconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));
    % mean_x_imps_xconc95=all_mean_x_imps(logical(above_thr_conc_imps)&logical(above_thr_x_imps));
    %
    % mean_conc_imps_yconc95=all_mean_conc_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));
    % mean_y_imps_yconc95=all_mean_y_imps(logical(above_thr_conc_imps)&logical(above_thr_y_imps));
    %
    % plot(log10(mean_x_imps_x95),log10(mean_conc_imps_x95),'o','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7]);
    % plot(log10(mean_x_imps_conc95),log10(mean_conc_imps_conc95),'o','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0]);
    % plot(log10(mean_x_imps_xconc95),log10(mean_conc_imps_xconc95),'o','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor',[0.8 0.6 0.7]);
    % plot([log10(min_imp) log10(max_imp)],[log10(min_imp) log10(max_imp)],'-k','LineWidth',2)
    % xlim([log10(min_imp) log10(max_imp)])
    % ylim([log10(min_imp) log10(max_imp)])
    % plot([log10(thr_conc_imps(fileNo)) log10(thr_conc_imps(fileNo))],[log10(min_imp) log10(max_imp)],'-k')
    % plot([log10(min_imp) log10(max_imp)],[log10(thr_x_imps(fileNo)) log10(thr_x_imps(fileNo))],'-k')
    % xlabel('Importance x (log10)')
    % ylabel('Importance conc (log10)')
    % title(['Prediction importance for all files'], 'Interpreter', 'none')

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

        if include_ii==1
            ii_included=ii_included+1;
            ii_95_class(ii_included)=this_ii_95_class;
            mean_x_imps(ii_included)=log10(all_mean_x_imps(ii));
            mean_y_imps(ii_included)=log10(all_mean_y_imps(ii));
            mean_conc_imps(ii_included)=log10(all_mean_conc_imps(ii));
        end

    end

    data = [mean_x_imps; mean_y_imps; mean_conc_imps];

    % %This was the old way to do this
    % mean_x_imps=log10([mean_x_imps_x95' mean_x_imps_y95' mean_x_imps_conc95']);
    % mean_y_imps=log10([mean_y_imps_x95' mean_y_imps_y95' mean_y_imps_conc95']);
    % mean_conc_imps=log10([mean_conc_imps_x95' mean_conc_imps_y95' mean_conc_imps_conc95']);
    %
    % ii_x_imps=[1:length(mean_x_imps_x95)];
    % ii_y_imps=[length(mean_x_imps_x95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)];
    % ii_conc_imps=[length(mean_x_imps_x95)+length(mean_y_imps_y95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)+length(mean_conc_imps_conc95)];
    %
    % data = [mean_x_imps; mean_y_imps; mean_conc_imps]; % Replace with your actual data

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
    title(['Prediction importance above 95 percentile ' group_label{this_group}], 'Interpreter', 'none')

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
    title(['Prediction importance above 95 percentile ' group_label{this_group}], 'Interpreter', 'none')



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
    title(['Prediction importance above 95 percentile ' group_label{this_group}], 'Interpreter', 'none')

    % Step 1: Create a sample data matrix
    % Assume you have 3 variables and 100 samples
    %
    % mean_x_imps=log10([mean_x_imps_x95' mean_x_imps_y95' mean_x_imps_conc95']);
    % mean_y_imps=log10([mean_y_imps_x95' mean_y_imps_y95' mean_y_imps_conc95']);
    % mean_conc_imps=log10([mean_conc_imps_x95' mean_conc_imps_y95' mean_conc_imps_conc95']);
    %
    % ii_x_imps=[1:length(mean_x_imps_x95)];
    % ii_y_imps=[length(mean_x_imps_x95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)];
    % ii_conc_imps=[length(mean_x_imps_x95)+length(mean_y_imps_y95)+1:length(mean_x_imps_x95)+length(mean_y_imps_y95)+length(mean_conc_imps_conc95)];
    %
    % data = [mean_x_imps; mean_y_imps; mean_conc_imps]; % Replace with your actual data
    % 
    % % Step 2: Perform PCA
    % [coeff, score, latent, tsquared, explained] = pca(data');
    % 
    % % Step 3: Plot the first two principal components
    % figureNo=figureNo+1;
    % try
    %     close(figureNo)
    % catch
    % end
    % hFig=figure(figureNo);
    % hold on
    % 
    % ax=gca;ax.LineWidth=3;
    % set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    % scatter(score(:,1), score(:,2), 'filled');
    % xlabel('1st Principal Component');
    % ylabel('2nd Principal Component');
    % title(['PCA of Data' group_label{this_group}]');
    % grid on;
    % 
    % % Optionally display explained variance
    % disp('Explained Variance by Each Principal Component:');
    % disp(explained);


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
                plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.95 0.9 0.25],'MarkerEdgeColor',[0.95 0.9 0.25]);
            case 5
                %5 is >=y95 and conc95
                plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0 0.45 0.7],'MarkerEdgeColor',[0 0.45 0.7]);
            case 6
                %6 is >=x95 and y95
                plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.8 0.4 0],'MarkerEdgeColor',[0.8 0.4 0]);
            case 7
                %7 is >= all 95s
                plot(mappedX(ii,1),mappedX(ii,2),'o','MarkerFaceColor',[0.8 0.6 0.7],'MarkerEdgeColor',[0.8 0.6 0.7]);
        end
    end

    these_ylim=ylim;
    these_xlim=xlim;
    text(these_xlim(1)+0.85*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.95*(these_ylim(2)-these_ylim(1)),'x','Color',[0.9 0.6 0],'FontWeight','bold','FontSize',16)
    text(these_xlim(1)+0.20*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.92*(these_ylim(2)-these_ylim(1)),'y','Color',[0.35 0.7 0.9],'FontWeight','bold','FontSize',16)
    text(these_xlim(1)+0.45*(these_xlim(2)-these_xlim(1)),these_ylim(1)+0.2*(these_ylim(2)-these_ylim(1)),'odor','Color',[0 0.6 0.5],'FontWeight','bold','FontSize',16)
    xlabel('t-SNE Component 1');
    ylabel('t-SNE Component 2');
    title(['t-SNE Prediction Importance ' group_label{this_group}]);
    pffft=1;
    % grid on;
end

fclose(fileID);

pffft=1;