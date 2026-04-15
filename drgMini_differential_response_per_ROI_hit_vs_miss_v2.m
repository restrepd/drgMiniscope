%drgMini_differential_response_per_ROI_hit_vs_miss_v2
close all
clear all

is_sphgpu=1;

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
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m';

        %Trained with all trials and mult=1 taking on account when mouse detects the odor
        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc_1_11222025/';
        choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_1_11222025.m';

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        %This is not used here, all the Moser infor is input through PredImp
        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';


        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        save_PathPredImp='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FilePredImp='outputPredictionImportance.mat';


    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/data2/DRMatlab/drgMiniscope')
        addpath('/data2/DRMatlab/m new/Chi Squared')
        addpath('/data2/DRMatlab/drgMaster')
        addpath(genpath('/data2/DRMatlab/m new/kakearney-boundedline-pkg-32f2a1f'))

        %Trained with all trials and mult=1 taking on account when mouse detects the odor
        save_PathConc='/data2/SFTP/DecodeDynOdorConc_1_11222025/';
        choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_1_11222025.m';

        save_PathAngle='/data2/SFTP/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        %This one has the dFF per trial
        save_PathXY='/data2/SFTP/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        choiceBatchPathName='/data2/SFTP/PreProcessed/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

end

save_path_file='/data2/SFTP/glm_div.mat';

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
% eval(['handles_Moser=' choiceMoserFileName(1:end-2) ';'])

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

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

these_groups=[1 5];
ii_run=1;

n_shuffle_SI=100;


x=25:50:475;
y=24:48:456;

trial_type_labels{1}='Hit';
trial_type_labels{2}='Miss';
% trial_type_labels{3}='Hit4';
% trial_type_labels{4}='Miss4';

these_adiv_names{1}='all_div_hit';
these_adiv_names{2}='all_div_miss';
% these_adiv_names{3}='all_div_hit4';
% these_adiv_names{4}='all_div_miss4';

%Labels for comparisons
ii_comp=0;
ii_type1=1;
ii_type2=2;
ii_comp=ii_comp+1;
trial_type_comp_labels{ii_comp}=[trial_type_labels{ii_type1} ' vs ' trial_type_labels{ii_type2}];
trial_type_comp_ii_type1(ii_comp)=ii_type1;
trial_type_comp_ii_type2(ii_comp)=ii_type2;
% these_types=[1:2];
% these_types=these_types((these_types~=ii_type1)&(these_types~=ii_type2))
% trial_type_comp_ii_type3(ii_comp)=these_types(1);
% trial_type_comp_ii_type4(ii_comp)=these_types(2);
% trial_type_no_comp_labels{ii_comp}=[trial_type_labels{these_types(1)} ' vs ' trial_type_labels{these_types(2)}];


%Save data for before divergence
p_values_bef=[];
ii_all_pvalues_bef=0;
all_ROI_files_bef=[];
all_iiROIs_bef=[];
pFDRs_per_file_bef=[];

%Save data for after divergence
p_values_aft=[];
ii_all_pvalues_aft=0;
all_ROI_files_aft=[];
all_iiROIs_aft=[];
pFDRs_per_file_aft=[];

handles_out2=[];

%Show a subset of the spatial maps
figureNo=figureNo+1;


handles_out2.no_neurons=0;
handles_out2.no_significant_bef=0;
handles_out2.no_significant_aft=0;

handles_out2.no_files=0;
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Get trials
        arena_file=handles_conc.arena_file{fileNo};
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        % no_neurons=handles_out.no_neurons

        %Get dFF
        this_path=handles_conc.this_path;
        dFF_file=handles_conc.dFF_file;
        arena_file=handles_conc.arena_file;

        %Get dFF
        dFF=[];
        this_dFF_file=dFF_file{fileNo};
        this_this_path=this_path{fileNo};
        if strcmp(this_dFF_file(end-3:end),'.mat')
            %This reads the extract file
            load([this_this_path this_dFF_file])
            dFF=zeros(size(output.temporal_weights,1),size(output.temporal_weights,2));
            for traceNo=1:size(output.temporal_weights,2)
                dFF(:,traceNo)=output.temporal_weights(:,traceNo);
            end
        else
            if strcmp(dFF_file(end-4:end),'.hdf5')
                %This is an hdf5 generated by CaImAn
                dFF=h5read([this_this_path this_dFF_file],'/estimates/F_dff' );
            else
                %This is a csv file created from ImageJ
                dFF=readmatrix([this_this_path this_dFF_file]);
            end
        end

        % load([this_path{fileNo} arena_file{fileNo}])
        % no_time_points=length(arena.xsync);

        no_time_points=size(dFF,1);
        dt_miniscope=handles_conc.dt_miniscope;
        dt=handles_conc.dt;

        dFF_times=[1:no_time_points]*dt_miniscope;

        no_neurons=size(dFF,2);
        no_time_bins=round(dFF_times(end)/dt);
        time_binned=[1:no_time_bins]*dt-dt/2;
        neural_data=zeros(no_time_bins,no_neurons);
        % pos_binned=zeros(no_time_bins,2);

        for ii_time_bin=1:no_time_bins
            time_from=time_binned(ii_time_bin)-dt/2;
            time_to=time_binned(ii_time_bin)+dt/2;
            % pos_binned(ii_time_bin,:)=mean(pos((dFF_times>=time_from)&(dFF_times<time_to),:),1);
            for ii_neuron=1:no_neurons
                neural_data(ii_time_bin,ii_neuron)=mean(dFF((dFF_times>=time_from)&(dFF_times<time_to),ii_neuron),1);
            end
        end

        these_p_values_per_ROI_bef=zeros(no_neurons,1);
        these_p_values_per_ROI_aft=zeros(no_neurons,1);

        for ii_ROI=1:no_neurons


            %These encompass times from the earliest start time to the latest end times in all trials
            hit_dFF=[];
            ii_hit=0;
            miss_dFF=[];
            ii_miss=0;


            %Let's do the glm and ZETA stats from -3 to 3 sec in 1 sec bins
            time_range=[-3 3];
            ii_for_time_bins=(time_range(2)-time_range(1))/handles_XY.dt;
            odor_aligned_time_bins=[0:ii_for_time_bins]*handles_XY.dt+time_range(1);


            %Let's find ROIs whose dFF does not differ before odor and
            %diverges after odor application
            glm_div_bef_ii=0;
            glm_div_bef=[];

            glm_div_aft_ii=0;
            glm_div_aft=[];

            %GLM will be calculated on a 0.5 sec bins
            trimmed_dt=0.5;
            trimmed_time_bins=[time_range(1)+(trimmed_dt/2):trimmed_dt:time_range(2)-(trimmed_dt/2)];

            for trNo=1:trials.odor_trNo

                this_ii_odor_encounter=trials.odor_encounter_ii(trNo);

                if ((this_ii_odor_encounter-ii_for_time_bins/2) > 0) && (((this_ii_odor_encounter+ii_for_time_bins/2) < no_time_bins))

                    these_dFF=zeros(1,length(odor_aligned_time_bins));
                    these_dFF(1,:)=neural_data(this_ii_odor_encounter-ii_for_time_bins/2:this_ii_odor_encounter+ii_for_time_bins/2,ii_ROI);

                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            ii_hit=ii_hit+1;
                            hit_dFF(ii_hit,1:length(odor_aligned_time_bins))=these_dFF;

                            for ii_tr=1:length(trimmed_time_bins)
                                these_dFF_trimmed=these_dFF( (odor_aligned_time_bins>=(trimmed_time_bins(ii_tr)-(trimmed_dt/2)))...
                                    &(odor_aligned_time_bins<(trimmed_time_bins(ii_tr)+(trimmed_dt/2))) );
                                if trimmed_time_bins(ii_tr)>0
                                    glm_div_aft.data(glm_div_aft_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_aft.trial_type(glm_div_aft_ii+1)=1;
                                    glm_div_aft.time(glm_div_aft_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_aft_ii=glm_div_aft_ii+1;
                                end
                                if trimmed_time_bins(ii_tr)>0
                                    glm_div_bef.data(glm_div_bef_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_bef.trial_type(glm_div_bef_ii+1)=1;
                                    glm_div_bef.time(glm_div_bef_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_bef_ii=glm_div_bef_ii+1;
                                end
                            end
                        case 2
                            %Lane 1 miss orange
                            ii_miss=ii_miss+1;
                            miss_dFF(ii_miss,1:length(odor_aligned_time_bins))=these_dFF;

                            for ii_tr=1:length(trimmed_time_bins)
                                these_dFF_trimmed=these_dFF( (odor_aligned_time_bins>=(trimmed_time_bins(ii_tr)-(trimmed_dt/2)))...
                                    &(odor_aligned_time_bins<(trimmed_time_bins(ii_tr)+(trimmed_dt/2))) );
                                if trimmed_time_bins(ii_tr)>0
                                    glm_div_aft.data(glm_div_aft_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_aft.trial_type(glm_div_aft_ii+1)=2;
                                    glm_div_aft.time(glm_div_aft_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_aft_ii=glm_div_aft_ii+1;
                                end
                                if trimmed_time_bins(ii_tr)>0
                                    glm_div_bef.data(glm_div_bef_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_bef.trial_type(glm_div_bef_ii+1)=2;
                                    glm_div_bef.time(glm_div_bef_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_bef_ii=glm_div_bef_ii+1;
                                end
                            end
                        case 3
                            %Lane 4 hit blue
                            ii_hit=ii_hit+1;
                            hit_dFF(ii_hit,1:length(odor_aligned_time_bins))=these_dFF;

                            for ii_tr=1:length(trimmed_time_bins)
                                these_dFF_trimmed=these_dFF( (odor_aligned_time_bins>=(trimmed_time_bins(ii_tr)-(trimmed_dt/2)))...
                                    &(odor_aligned_time_bins<(trimmed_time_bins(ii_tr)+(trimmed_dt/2))) );
                                if trimmed_time_bins(ii_tr)>0
                                    glm_div_aft.data(glm_div_aft_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_aft.trial_type(glm_div_aft_ii+1)=1;
                                    glm_div_aft.time(glm_div_aft_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_aft_ii=glm_div_aft_ii+1;
                                end
                                if trimmed_time_bins(ii_tr)>0
                                    glm_div_bef.data(glm_div_bef_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_bef.trial_type(glm_div_bef_ii+1)=1;
                                    glm_div_bef.time(glm_div_bef_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_bef_ii=glm_div_bef_ii+1;
                                end
                            end
                        case 4
                            %Lane 4 miss sky blue
                            ii_miss=ii_miss+1;
                            miss_dFF(ii_miss,1:length(odor_aligned_time_bins))=these_dFF;

                            for ii_tr=1:length(trimmed_time_bins)
                                these_dFF_trimmed=these_dFF( (odor_aligned_time_bins>=(trimmed_time_bins(ii_tr)-(trimmed_dt/2)))...
                                    &(odor_aligned_time_bins<(trimmed_time_bins(ii_tr)+(trimmed_dt/2))) );
                                if trimmed_time_bins(ii_tr)>0
                                    glm_div_aft.data(glm_div_aft_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_aft.trial_type(glm_div_aft_ii+1)=2;
                                    glm_div_aft.time(glm_div_aft_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_aft_ii=glm_div_aft_ii+1;
                                end
                                if trimmed_time_bins(ii_tr)<=0
                                    glm_div_bef.data(glm_div_bef_ii+1)=mean(these_dFF_trimmed);
                                    glm_div_bef.trial_type(glm_div_bef_ii+1)=2;
                                    glm_div_bef.time(glm_div_bef_ii+1)=trimmed_time_bins(ii_tr);
                                    glm_div_bef_ii=glm_div_bef_ii+1;
                                end
                            end
                    end
                    pffft=1;
                end
            end

            handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF=mean(hit_dFF,1);
            handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF=mean(miss_dFF,1);
            handles_out2.file(fileNo).odor_aligned_time_bins=odor_aligned_time_bins;


            %Plot the dFF timecourses
            try
                close(figureNo)
            catch
            end


            hFig = figure(figureNo);
            set(hFig, 'units','normalized','position',[.1 .1 .75 .75])

            %Plot dFFs for lane 1 and 4
            y_gap=2;
            hit_miss=zeros(size(hit_dFF,1)+size(miss_dFF,1)+y_gap,length(odor_aligned_time_bins));
            hit_miss(size(miss_dFF,1)+y_gap+1:size(hit_dFF,1)+size(miss_dFF,1)+y_gap,:)=hit_dFF;
            hit_miss(1:size(miss_dFF,1),:)=miss_dFF;

            y_trials=[1:size(hit_miss,1)];


            y_trials_end=y_trials(end);
            all_hit_miss=hit_miss(:);
            c_percentile=prctile(all_hit_miss(all_hit_miss>0),99);

            %Plot dFF for lane 1
            % subplot(2, 6, [1 2 3]);
            hold on

            drg_pcolor(repmat(odor_aligned_time_bins,length(y_trials),1),repmat(y_trials,length(odor_aligned_time_bins),1)',hit_miss)
            colormap(this_cmap)
            if max(max(max(hit_miss(:))))>0
                clim([0 c_percentile]);
            end

            shading flat
            plot([0 0],[y_trials(1) y_trials(end)+1],'-w','LineWidth',3)



            rectangle('Position', [odor_aligned_time_bins(1), size(miss_dFF,1)+1, 1.03*(odor_aligned_time_bins(end)-odor_aligned_time_bins(1)), y_gap], ... % [x, y, width, height]
                'FaceColor', 'white', ...    % Fill color (white)
                'EdgeColor', 'none');        % No border for the rectangle

            ylim([0 y_trials_end+2])
            xlim([odor_aligned_time_bins(1)-0.05*(odor_aligned_time_bins(end)-odor_aligned_time_bins(1)) odor_aligned_time_bins(end)+0.05*(odor_aligned_time_bins(end)-odor_aligned_time_bins(1)) ])
            xlabel('Time (sec)')
            yticks([1+(size(miss_dFF,1)/2) size(miss_dFF,1)+1+y_gap+(size(hit_dFF,1)/2)])
            yticklabels({'Misses','Hits'})
            title(['dFF'])

            title('dFF for hit vs miss')

            %Now do glm before
            tbl = table(glm_div_bef.data',glm_div_bef.trial_type',glm_div_bef.time',...
                'VariableNames',{'dFF','trial_type','time'});
            mdl = fitglm(tbl,'dFF~trial_type+time'...
                ,'CategoricalVars',[2])

            %Save the p values
            ii_all_pvalues_bef=ii_all_pvalues_bef+1;
            all_ROI_files_bef(ii_all_pvalues_bef)=fileNo;
            all_iiROIs_bef(ii_all_pvalues_bef)=ii_ROI;


            if isnan(coefTest(mdl,[0 1 0]))
                p_values_bef(ii_all_pvalues_bef,1)=1;
            else
                p_values_bef(ii_all_pvalues_bef,1)=coefTest(mdl,[0 1 0]);
            end


            these_p_values_per_ROI_bef(ii_ROI,1)=coefTest(mdl,[0 1 0]);

            fprintf(1, '\n')
            fprintf(1, ['\nFor hit vs. miss p value before: '...
                num2str(p_values_bef(ii_all_pvalues_bef,1)) '\n'])

            %Now do glm after
            tbl = table(glm_div_aft.data',glm_div_aft.trial_type',glm_div_aft.time',...
                'VariableNames',{'dFF','trial_type','time'});
            mdl = fitglm(tbl,'dFF~trial_type+time'...
                ,'CategoricalVars',[2])

            %Save the p values
            ii_all_pvalues_aft=ii_all_pvalues_aft+1;
            all_ROI_files_aft(ii_all_pvalues_aft)=fileNo;
            all_iiROIs_aft(ii_all_pvalues_aft)=ii_ROI;

            if isnan(coefTest(mdl,[0 1 0]))
                p_values_aft(ii_all_pvalues_aft,1)=1;
            else
                p_values_aft(ii_all_pvalues_aft,1)=coefTest(mdl,[0 1 0]);
            end


            these_p_values_per_ROI_aft(ii_ROI,1)=coefTest(mdl,[0 1 0]);
            fprintf(1, '\n')
            fprintf(1, ['\nFor hit vs. miss p value after: '...
                num2str(p_values_aft(ii_all_pvalues_aft,1)) '\n'])




            pfft=1;

            if (fileNo==7)&(ii_ROI==154)
                pffft=1;
            end

            if (fileNo==13)&(ii_ROI==124)
                pffft=1;
            end

            if (fileNo==14)&(ii_ROI==8)
                pffft=1;
            end

        end

        %Get FDR p values
        fprintf(1, ['\n\n\nFor file No ' num2str(fileNo) ' the total number of neurons is ' num2str(no_neurons) '\n\n'])

        pFDRs_per_file_bef(fileNo)=drsFDRpval(these_p_values_per_ROI_bef(:));
        handles_out2.file(fileNo).pFDR_bef=pFDRs_per_file_bef(fileNo);
        fprintf(1, ['\n p FDR before: ' num2str(pFDRs_per_file_bef(fileNo)) '\n'])

        pFDRs_per_file_aft(fileNo)=drsFDRpval(these_p_values_per_ROI_aft(:));
        handles_out2.file(fileNo).pFDR_aft=pFDRs_per_file_aft(fileNo);
        fprintf(1, ['\n p FDR after: ' num2str(pFDRs_per_file_aft(fileNo)) '\n'])

        handles_out2.no_files=handles_out2.no_files+1;

        these_p_values_bef=zeros(no_neurons,1);
        these_p_values_bef(:,1)=these_p_values_per_ROI_bef(:,1);
        handles_out2.file(fileNo).these_p_values_bef=these_p_values_bef;
        handles_out2.file(fileNo).significant_bef=these_p_values_bef<=pFDRs_per_file_bef(fileNo);
        fprintf(1, ['\nFor hit vs miss before number significant ' num2str(sum(these_p_values_bef<=pFDRs_per_file_bef(fileNo))) '\n'])
        handles_out2.no_significant_bef=handles_out2.no_significant_bef+sum(handles_out2.file(fileNo).significant_bef);


        these_p_values_aft=zeros(no_neurons,1);
        these_p_values_aft(:,1)=these_p_values_per_ROI_aft(:,1);
        handles_out2.file(fileNo).these_p_values_aft=these_p_values_aft;
        handles_out2.file(fileNo).significant_aft=these_p_values_aft<=pFDRs_per_file_aft(fileNo);
        fprintf(1, ['\nFor hit vs miss before number significant ' num2str(sum(these_p_values_aft<=pFDRs_per_file_aft(fileNo))) '\n'])
        handles_out2.no_significant_aft=handles_out2.no_significant_aft+sum(handles_out2.file(fileNo).significant_aft);

        handles_out2.no_neurons=handles_out2.no_neurons+no_neurons;
    end
end

%Report the overall fraction of significant divergent responses
fprintf(1, ['\n\nPercent of cells showing divergent responses between hit and miss before odor in ' num2str(handles_out2.no_files)  ' sessions ' ...
    num2str(100*handles_out2.no_significant_bef/handles_out2.no_neurons) '\n'])
fprintf(1, ['\n\nPercent of cells showing divergent responses between hit and miss after odor in ' num2str(handles_out2.no_files)  ' sessions ' ...
    num2str(100*handles_out2.no_significant_aft/handles_out2.no_neurons) '\n'])
fprintf(1, ['\nTotal number of cells ' num2str(handles_out2.no_neurons)  '\n\n'])



%Now show the divergent responses for after


%Get all divergent ROIs
all_div_hit=[];
all_div_miss=[];
all_div_to_sort=[];


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
        these_significant_aft=handles_out2.file(fileNo).significant_aft;
        these_significant_bef=handles_out2.file(fileNo).significant_bef;
        no_neurons=length(these_significant_bef);
        for ii_ROI=1:no_neurons
            if (these_significant_aft(ii_ROI)==1)&(these_significant_bef(ii_ROI)==0)
                all_div_hit=[all_div_hit; handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF...
                    /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];
                all_div_miss=[all_div_miss; handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF...
                    /prctile([handles_out2.file(fileNo).ROI(ii_ROI).hit_dFF handles_out2.file(fileNo).ROI(ii_ROI).miss_dFF],95)];
            end
        end
        %Sort using the two ii_types that were used in the comparison
        all_div_to_sort=[all_div_hit all_div_miss];
        % eval(['all_div_to_sort=[' these_adiv_names{trial_type_comp_ii_type1(ii_comp)} ' '...
        %     these_adiv_names{trial_type_comp_ii_type2(ii_comp)} '];'])
    end
end

all_div=[all_div_hit;all_div_miss];


croscorr_traces=corrcoef(all_div_to_sort');

Z = linkage(croscorr_traces,'complete','correlation');

no_clusters=3;
handles_out2.clusters = cluster(Z,'Maxclust',no_clusters);
figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

%Do cutoff for no_clusters
cutoff = median([Z(end-(no_clusters-1),no_clusters) Z(end-(no_clusters-2),no_clusters)]);
[H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
set(H,'LineWidth',2)
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

%re-sort the matrix
for ii=1:size(croscorr_traces,1)
    for jj=1:size(croscorr_traces,1)
        perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
    end
end



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

    % sorted_handles_out.all_div_hit1=[];
    % sorted_handles_out.all_div_hit1=all_div_hit1(outperm,:);

    eval(['ii_included=size(sorted_handles_out.' these_adiv_names{ii_type} ',1);'])

    % ii_included=size(sorted_handles_out.all_div_hit1,1);

    time_span_mat=repmat(odor_aligned_time_bins,ii_included,1);
    ROI_mat=repmat(1:ii_included,length(odor_aligned_time_bins),1)';

    eval(['pcolor(time_span_mat,ROI_mat,sorted_handles_out.' these_adiv_names{ii_type} ')'])
    % pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_hit1)
    colormap(this_cmap)
    shading flat

    % caxis([prctile(sorted_handles_out.all_div_hit1(:),1) prctile(sorted_handles_out.all_div_hit1(:),99.9)])

    caxis([prctile(all_div(:),1) prctile(all_div(:),99.9)])

    plot([0 0],[0 ii_included],'-r')


    xlim([-3 3])
    ylim([1 ii_included])
    eval(['title([''' trial_type_comp_labels{ii_comp} ': ' trial_type_labels{ii_type} '''])'])
    xlabel('Time (sec)')
    ylabel('ROI number')
end




%Plot the average timecourses per cluster
do_std=0;
these_ylim=[];
for clus=1:no_clusters
    %Report the number of ROIs per cluster
    fprintf(1, ['\nFor '  trial_type_labels{trial_type_comp_ii_type1(ii_comp)}...
        ' vs ' trial_type_labels{trial_type_comp_ii_type2(ii_comp)} ' cluster No '...
        num2str(clus) ' has ' num2str(sum(handles_out2.clusters==clus)) ' ROIs\n'])

    clusters=handles_out2.clusters;

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

    % these_adiv_names{1}='all_div_hit1';
    % these_adiv_names{2}='all_div_miss1';
    % these_adiv_names{3}='all_div_hit4';
    % these_adiv_names{4}='all_div_miss1';

    % this_adiv_name1=these_adiv_names{trial_type_comp_ii_type1(ii_comp)};
    % this_adiv_name2=these_adiv_names{trial_type_comp_ii_type2(ii_comp)};


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

    % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

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


    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(this_cluster_miss), CIpv','cmap',[158/255 31/255 99/255]);

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



    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins, mean(this_cluster_hit), CIpv','cmap',[0 114/255 178/255]);


    plot(odor_aligned_time_bins',mean(this_cluster_miss)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
    plot(odor_aligned_time_bins',mean(this_cluster_hit)','Color',[0 114/255 178/255],'LineWidth',1.5);


    % ylim([-0.5 1.5])
    xlim([-3 3])
    these_ylim=[these_ylim; ylim];




    xlabel('Time(sec)')
    ylabel('dFF')
    title(['dFF mean for hit vs miss cluster No ' num2str(clus)])
    % eval(['title([''' trial_type_comp_labels{ii_comp} ''' '' significantly diff, cluster no ''  num2str(clus) ])'])

end


fig_minus=0;
for clus=1:no_clusters

    figNo=figureNo+fig_minus;
    figure(figNo)
    ylim([min(these_ylim(:)) max(these_ylim(:))])
    this_ylim=ylim;

    %Odor on markers
    plot([0 0],this_ylim,'-k')

    text(-2.5,this_ylim(1)+0.9*(this_ylim(2)-this_ylim(1)),'Hit','Color',[158/255 31/255 99/255])
    text(-2.5,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),'Miss','Color',[0 114/255 178/255])


    %Odor on markers
    plot([0 0],this_ylim,'-k')

    fig_minus=fig_minus-1;

end

pffft=1;


save(save_path_file,'handles_out2')


fclose(fileID);

pffft=1;

