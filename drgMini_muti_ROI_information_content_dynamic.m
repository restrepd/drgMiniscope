%drgMini_multi_ROI_information_content_dynamic
close all
clear all

is_sphgpu=0; %0=Mac, 1=sphgpu, 2=Alpine
handles.is_sphgpu=is_sphgpu;

overwrite=0; %If 0 do not overwrite
handles.overwrite=overwrite;
display_figures=1;
handles.display_figures=display_figures;
is_gpu=0; %1=GPU, 0=CPU
handles.is_gpu=is_gpu;

handles_choices.trial_start_offset=-15; %This was -10
handles_choices.trial_end_offset=15;

%Time bins for decoding in seconds
dt=0.1;
handles_choices.dt=dt;

%Speed for background air flow in mm/sec
air_flow_speed=50; %5 cm/sec = 50 mm/sec
handles.air_flow_speed=air_flow_speed;

switch is_sphgpu
    case 0
        %Mac

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
         % 
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

        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        %This is not used here, all the Moser infor is input through PredImp
        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';


        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        %This file is saved by drgMini_analyze_batch_ImpMoserXYConc
        save_PathPredImp='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FilePredImp='outputPredictionImportance.mat';

    case 1
        %sphgpu

        % fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        % addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        % addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))


         %Trained with hits only
         save_PathConc='/data2/SFTP//DecodeOdorConc01122025/';
         choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m';

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/data2/SFTP//OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        save_PathAngle='/data2/SFTP//Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        %This is not used here, all the Moser infor is input through PredImp
        save_PathMoser='/data2/SFTP//Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';

% 
%         choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
%         fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        %This file is saved by drgMini_analyze_batch_ImpMoserXYConc
        save_PathPredImp='/data2/SFTP/PreProcessed/';
        save_FilePredImp='outputPredictionImportance.mat';

        choiceBatchPathName='/data2/SFTP/PreProcessed/';
        
    case 2
        %Alpine
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/drgMiniscope')
%         addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/projects/drestrepo@xsede.org/software/DR_matlab/drgMaster')
%         addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))


         %Trained with hits only
         save_PathConc='/pl/active/restrepo-lab/DecodeOdorConc01122025/';
         choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m';

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/pl/active/restrepo-lab/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        save_PathAngle='/pl/active/restrepo-lab/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        %This is not used here, all the Moser infor is input through PredImp
        save_PathMoser='/pl/active/restrepo-lab/Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';

% 
%         choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
%         fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        %This file is saved by drgMini_analyze_batch_ImpMoserXYConc
        save_PathPredImp='/pl/active/restrepo-lab/PreProcessed/';
        save_FilePredImp='outputPredictionImportance.mat';

        choiceBatchPathName='/pl/active/restrepo-lab/PreProcessed/';
        
      
end


outputFile='multi_ROI_info_content_dynamic.mat';


tic

if overwrite==0
    load([choiceBatchPathName outputFile])
else
    handles_out2=[];
end

all_n_bits=[1 3 6 12 24 2000]; %Number of bits (ROIs) for most and least important prediction ROIs
                            %And for best

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
 
these_groups=[1 5];
ii_run=1;

n_shuffle_SI=100;

%Load the odor plumes
switch is_sphgpu
    case 0
        %Mac
        load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odor_plume_patternsDR.mat')
    case 1
        %sphgpu
        load('/data2/SFTP/Odor Arena Plumes/odor_plume_patternsDR.mat')
    case 2
        %Alpine
        load('/pl/active/restrepo-lab/PreProcessed/Odor Arena Plumes/odor_plume_patternsDR.mat')     
end
op_threshold=-8;
for cm_from_floor=1:2
    mean_plume_l1=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1';
    mean_plume_l4=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4';
    x_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume;
    y_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume;

    %Generate binary plumes
    binary_plumel1=zeros(10,10);
    binary_plumel4=zeros(10,10);

    for ii_x=1:10
        for ii_y=1:10
            these_opl1=[];
            these_opl1=mean_plume_l1( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
            this_mean_opl1=mean(these_opl1(:));
            if this_mean_opl1<op_threshold
                binary_plumel1(ii_x,ii_y)=0;
            else
                binary_plumel1(ii_x,ii_y)=1;
            end

             these_opl4=[];
            these_opl4=mean_plume_l4( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
            this_mean_opl4=mean(these_opl4(:));
            if this_mean_opl4<op_threshold
                binary_plumel4(ii_x,ii_y)=0;
            else
                binary_plumel4(ii_x,ii_y)=1;
            end
        end
    end

    odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1=binary_plumel1;
    odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4=binary_plumel4;

end
 
%Select ROIs for best MI between multiple ROI binary dFF and binary op
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Get XY and dFF per trial
        arena_file=handles_XY.arena_file{fileNo};
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        no_neurons=handles_out.no_neurons;
        
        handles_out2.file(fileNo).no_neurons=no_neurons;

        neural_data=[];
        binary_odor_plume=[];

       
        

        %These are used for shuffled calculations
        rng(fileNo)
        these_rnd=rand(1,2*n_shuffle_SI*no_neurons);
        ii_rnd=0;
        shuffled_neural_data=[];

        for ii_ROI=1:no_neurons
            all_bin_dFF=[];
            all_bin_op=[];
            hits=[];
            
            for trNo=1:trials.odor_trNo
                these_x=trials.trial(trNo).XYtest(:,1);
                these_y=trials.trial(trNo).XYtest(:,2);

                these_dFF=trials.trial(trNo).XdFFtest(:,ii_ROI);

                ii_t_odor_on=1-handles_choices.trial_start_offset;

                these_bin_dFF=(these_dFF>0);
                all_bin_dFF=[all_bin_dFF;these_bin_dFF];

                these_bin_op=zeros(length(these_x),1);
                these_hits=zeros(length(these_x),1);
                for ii_t=1:length(these_x)
                    this_x_ii=ceil(these_x(ii_t)/50);
                    if this_x_ii==11
                        this_x_ii=10;
                    end

                    this_y_ii=ceil(these_y(ii_t)/48);
                    if this_y_ii==11
                        this_y_ii=10;
                    end

                  
                    %Keep track of binary op
                    if ii_t<=ii_t_odor_on
                        %Before odor is turned on this is zero
                        these_bin_op(ii_t)=0;
                    else
                        x_on=(ii_t-ii_t_odor_on)*dt*air_flow_speed;
                        if these_x(ii_t)>x_on
                            these_bin_op(ii_t)=0;
                        else
                            if trials.lane_per_trial(trNo)==1

                                switch handles_conc.group(fileNo)
                                    case 1
                                        %2 cm from the floor
                                        cm_from_floor=2;
                                        these_bin_op(ii_t)=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                                    case 5
                                        %1 cm from the floor
                                        cm_from_floor=1;
                                        these_bin_op(ii_t)=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1(this_x_ii,this_y_ii);
                                end
                            else
                                switch handles_conc.group(fileNo)
                                    case 1
                                        %2 cm from the floor
                                        cm_from_floor=2;
                                        these_bin_op(ii_t)=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);
                                    case 5
                                        %1 cm from the floor
                                        cm_from_floor=1;
                                        these_bin_op(ii_t)=odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4(this_x_ii,this_y_ii);
                                end
                            end
                        end
                    end
                    %Also keep track of hits and misses
                    %trials.odor_trial_type(trNo)
                    %1 Hit Lane 1
                    %2 Miss Lane 1
                    %3 Hit Lane 4
                    %4 Miss Lane 1
                    if (trials.odor_trial_type(trNo)==1)||(trials.odor_trial_type(trNo)==3)
                        these_hits(ii_t)=1;
                    else
                        these_hits(ii_t)=0;
                    end

                end
                all_bin_op=[all_bin_op;these_bin_op];
                hits=[hits;these_hits];
            end
            neural_data(:,ii_ROI)=all_bin_dFF;
            binary_odor_plume=all_bin_op;

            for ii_sh=1:n_shuffle_SI
                %With this one you flip and roll dFF
                ii_rnd=ii_rnd+1;
                offset_ii=floor(these_rnd(ii_rnd)*length(all_bin_dFF));
                all_bin_dFF_reversed=zeros(1,length(all_bin_dFF));

                for ii_trl=1:length(all_bin_dFF)
                    this_ii_trl=ii_trl+offset_ii;
                    if this_ii_trl>length(all_bin_dFF)
                        offset_ii=-ii_trl+1;
                        this_ii_trl=ii_trl+offset_ii;
                    end
                    all_bin_dFF_reversed(1,length(all_bin_dFF)-ii_trl+1)=all_bin_dFF(this_ii_trl);
                end
                shuffled_neural_data.sh(ii_sh).ROI(ii_ROI).neural_data=all_bin_dFF_reversed;
            end
        end

        ii=0;
        for ii_neurons=all_n_bits
            ii=ii+1;
            if ii_neurons>size(neural_data,2)
                no_neurons=size(neural_data,2);
            else
                no_neurons=ii_neurons;
            end
            if overwrite==1
                proceed=1;
            else
                proceed=0;
                if fileNo==length(handles_out2.best.file)
                    if ii_neurons>length(handles_out2.best.file(fileNo).ii_neuron)
                        proceed=1;
                    end
                end
                if fileNo>length(handles_out2.best.file)
                    proceed=1;
                end
            end
            if proceed==1
                start_toc=toc;
                fprintf(1,['Starting best file No ' num2str(fileNo) ' n bits= ' num2str(ii_neurons) '\n\n'] )
                if is_gpu==1
                    [handles_out2.best.file(fileNo).ii_neuron(ii).selected_best_neurons,...
                        handles_out2.best.file(fileNo).ii_neuron(ii).best_mi,...
                        handles_out2.best.file(fileNo).ii_neuron(ii).selected_worst_neurons,...
                        handles_out2.best.file(fileNo).ii_neuron(ii).worst_mi]...
                        = best_mi_selection_gpu(neural_data, binary_odor_plume, ii_neurons);
                else
                    [handles_out2.best.file(fileNo).ii_neuron(ii).selected_neurons,...
                        handles_out2.best.file(fileNo).ii_neuron(ii).best_mi,...
                        handles_out2.best.file(fileNo).ii_neuron(ii).selected_worst_neurons,...
                        handles_out2.best.file(fileNo).ii_neuron(ii).worst_mi]...
                        = best_mi_selection_cpu(neural_data, binary_odor_plume, ii_neurons);
                end
                fprintf(1,['Duration (mins) ' num2str((toc-start_toc)/60) '\n\n'] )
                fprintf(1,['Done with file No ' num2str(fileNo) ' n bits= ' num2str(no_neurons) ' best mi= '...
                    num2str(handles_out2.best.file(fileNo).ii_neuron(ii).best_mi) '\n\n'] )
                fprintf(1,['Done with file No ' num2str(fileNo) ' n bits= ' num2str(no_neurons) ' worst mi= '...
                    num2str(handles_out2.best.file(fileNo).ii_neuron(ii).worst_mi) '\n\n'] )
                handles_out2.all_n_bits=all_n_bits;
                save([choiceBatchPathName outputFile],'handles','handles_out2')
            end
        end


        %Now calculate the predictive importance multi ROI mutual information
        %Load information on prediction importance
        load([save_PathPredImp save_FilePredImp])

        imp_ROIs=imps.all_imps_ROI(imps.file_numbers==fileNo);
        all_imps=imps.all_mean_conc_imps(imps.file_numbers==fileNo);

        ii=0;
        for ii_neurons=all_n_bits
            ii=ii+1;
            if ii_neurons>size(neural_data,2)
                no_neurons=size(neural_data,2);
            else
                no_neurons=ii_neurons;
            end
            if overwrite==1
                proceed=1;
            else
                if ~isfield(handles_out2,'imps')
                    proceed=1;
                else
                    proceed=0;
                    if fileNo==length(handles_out2.best.file)
                        if fileNo>length(handles_out2.imps.file)
                            proceed=1;
                        else
                            if ii_neurons>length(handles_out2.imps.file(fileNo).ii_neuron)
                                proceed=1;
                            end
                        end
                    end
                    if fileNo>length(handles_out2.best.file)
                        proceed=1;
                    end
                end
            end
            if proceed==1
                start_toc=toc;
                fprintf(1,['Starting predictive importance file No ' num2str(fileNo) ' n bits= ' num2str(ii_neurons) '\n\n'] )

                [handles_out2.imps.file(fileNo).ii_neuron(ii).selected_most_imp_neurons,...
                    handles_out2.imps.file(fileNo).ii_neuron(ii).most_imp_mi,...
                    handles_out2.imps.file(fileNo).ii_neuron(ii).selected_least_imp_neurons,...
                    handles_out2.imps.file(fileNo).ii_neuron(ii).least_imp_mi]...
                    = imps_mi_selection_cpu(neural_data, binary_odor_plume, ii_neurons,imp_ROIs,all_imps);

                fprintf(1,['Duration (mins) ' num2str((toc-start_toc)/60) '\n\n'] )
                fprintf(1,['Done with file No ' num2str(fileNo) ' n bits= ' num2str(no_neurons) ' most important mi= '...
                    num2str(handles_out2.imps.file(fileNo).ii_neuron(ii).most_imp_mi) '\n\n'] )
                fprintf(1,['Done with file No ' num2str(fileNo) ' n bits= ' num2str(no_neurons) ' least important mi= '...
                    num2str(handles_out2.imps.file(fileNo).ii_neuron(ii).least_imp_mi) '\n\n'] )
                handles_out2.all_n_bits=all_n_bits;
                save([choiceBatchPathName outputFile],'handles','handles_out2')
            end
        end

        %Now calculate MI for hits and misses for all ROIs
        pfft=1;
        no_neurons=size(neural_data,2);

        %Hits
        hit_neural_data=zeros(sum(hits),no_neurons);
        hit_neural_data=neural_data(logical(hits),:);
        hit_binary_odor_plume=zeros(sum(hits),1);
        hit_binary_odor_plume(:,1)=binary_odor_plume(logical(hits),1);
        handles_out2.hit_miss.file(fileNo).hit_mi=drg_compute_mutual_information_cpu(hit_neural_data,hit_binary_odor_plume);
        fprintf(1,['For file No ' num2str(fileNo) ' n bits= ' num2str(no_neurons) ' mi for hits= '...
                    num2str(handles_out2.hit_miss.file(fileNo).hit_mi) '\n\n'] )

        %Misses
        miss_neural_data=zeros(sum(~hits),no_neurons);
        miss_neural_data=neural_data(~logical(hits),:);
        miss_binary_odor_plume=zeros(sum(~hits),1);
        miss_binary_odor_plume(:,1)=binary_odor_plume(~logical(hits),1);
        handles_out2.hit_miss.file(fileNo).miss_mi=drg_compute_mutual_information_cpu(miss_neural_data,miss_binary_odor_plume);
        fprintf(1,['For file No ' num2str(fileNo) ' n bits= ' num2str(no_neurons) ' mi for misses= '...
                    num2str(handles_out2.hit_miss.file(fileNo).miss_mi) '\n\n'] )

        save([choiceBatchPathName outputFile],'handles','handles_out2')

        %Now calculate shuffled MIs
        fprintf(1,['Starting shuffled mi calculation for file No ' num2str(fileNo)  '\n\n'] )
        start_toc=toc;
        for ii=1:length(all_n_bits)

            these_best_ROIs=handles_out2.best.file(fileNo).ii_neuron(ii).selected_neurons;
            these_worst_ROIs=handles_out2.best.file(fileNo).ii_neuron(ii).selected_worst_neurons;
            these_most_imp_ROIs=handles_out2.imps.file(fileNo).ii_neuron(ii).selected_most_imp_neurons;
            these_least_imp_ROIs=handles_out2.imps.file(fileNo).ii_neuron(ii).selected_least_imp_neurons;
            for ii_sh=1:n_shuffle_SI
                these_best_neural_data=zeros(size(neural_data,1),length(these_best_ROIs));
                these_worst_neural_data=zeros(size(neural_data,1),length(these_best_ROIs));
                these_most_imp_neural_data=zeros(size(neural_data,1),length(these_best_ROIs));
                these_least_imp_neural_data=zeros(size(neural_data,1),length(these_best_ROIs));
                for jj_ROI=1:length(these_best_ROIs)
                    ii_ROI=these_best_ROIs(jj_ROI);
                    these_best_neural_data(:,jj_ROI)=shuffled_neural_data.sh(ii_sh).ROI(ii_ROI).neural_data;
                    ii_ROI=these_worst_ROIs(jj_ROI);
                    these_worst_neural_data(:,jj_ROI)=shuffled_neural_data.sh(ii_sh).ROI(ii_ROI).neural_data;
                    ii_ROI=these_most_imp_ROIs(jj_ROI);
                    these_most_imp_neural_data(:,jj_ROI)=shuffled_neural_data.sh(ii_sh).ROI(ii_ROI).neural_data;
                    ii_ROI=these_least_imp_ROIs(jj_ROI);
                    these_least_imp_neural_data(:,jj_ROI)=shuffled_neural_data.sh(ii_sh).ROI(ii_ROI).neural_data;
                end

                this_mi = drg_compute_mutual_information_cpu(these_best_neural_data, binary_odor_plume);
                handles_out2.best.file(fileNo).ii_neuron(ii).sh_best_mi(ii_sh)=this_mi;

                this_mi = drg_compute_mutual_information_cpu(these_worst_neural_data, binary_odor_plume);
                handles_out2.best.file(fileNo).ii_neuron(ii).sh_worst_mi(ii_sh)=this_mi;

                this_mi = drg_compute_mutual_information_cpu(these_most_imp_neural_data, binary_odor_plume);
                handles_out2.imps.file(fileNo).ii_neuron(ii).sh_most_imp_mi(ii_sh)=this_mi;

                this_mi = drg_compute_mutual_information_cpu(these_least_imp_neural_data, binary_odor_plume);
                handles_out2.imps.file(fileNo).ii_neuron(ii).sh_least_imp_mi(ii_sh)=this_mi;
            end

        end
        fprintf(1,['Duration (mins) ' num2str((toc-start_toc)/60) '\n\n'] )
        fprintf(1,['Done with shuffled mi calculation file No ' num2str(fileNo) '\n\n'] )

        save([choiceBatchPathName outputFile],'handles','handles_out2')





    end
end

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

pffft=1;


function [selected_best_neurons,best_mi,selected_worst_neurons,worst_mi] = best_mi_selection_cpu(neural_data, odor_plume, n_subset)
    % Inputs:
    % neural_data: matrix of size [num_samples x num_neurons], where each column represents a neuron's activity.
    % odor_plume: vector of size [num_samples x 1], representing the odor plume signal.
    % max_neurons: maximum number of neurons to select.
    %
    % Output:
    % selected_neurons: indices of the selected neurons.

    num_neurons = size(neural_data, 2);
    num_combinations=500000; %Number of random combinations if nchooseek crashes

    if n_subset>num_neurons
        n_subset=num_neurons;
    end
    try
        all_combinations = nchoosek(1:num_neurons, n_subset);
    catch
        % Initialize the output matrix
        all_combinations = [];

        % Generate random combinations
        for i = 1:num_combinations
            all_combinations(i, :) = sort(randperm(num_neurons, n_subset));
        end
        % Remove any duplicate combinations
        all_combinations = unique(all_combinations, 'rows');

        pffft=1;
    end

    

    % best_mi = -inf;

    mi_combined=zeros(1,size(all_combinations,1));
    parfor ii_combs = 1:size(all_combinations,1)
        current_subset = all_combinations(ii_combs,:); % Add current neuron to the subset
        combined_activity = neural_data(:, current_subset); % Combine selected neuron activities

        % Compute mutual information between combined activity and odor plume
        mi_combined(ii_combs) = drg_compute_mutual_information_cpu(combined_activity, odor_plume);

        % if mi_combined > best_mi
        %     best_mi = mi_combined;
        %     selected_neurons=current_subset;
        % end
    end
    [best_mi,ii_best]=max(mi_combined);
    selected_best_neurons=zeros(1,size(all_combinations,2));
    selected_best_neurons(1,:)=all_combinations(ii_best,:);
    [worst_mi,ii_worst]=min(mi_combined);
    selected_worst_neurons=zeros(1,size(all_combinations,2));
    selected_worst_neurons(1,:)=all_combinations(ii_best,:);

end

function [selected_neurons_best,best_mi,selected_neurons_worst,worst_mi] = best_mi_selection_gpu(neural_data, odor_plume, n_subset)
    % Inputs:
    % neural_data: matrix of size [num_samples x num_neurons], where each column represents a neuron's activity.
    % odor_plume: vector of size [num_samples x 1], representing the odor plume signal.
    % max_neurons: maximum number of neurons to select.
    %
    % Output:
    % selected_neurons: indices of the selected neurons.

    num_neurons = size(neural_data, 2);
    num_combinations=500000; %Number of random combinations if nchooseek crashes

    if n_subset>num_neurons
        n_subset=num_neurons;
    end
    try
        all_combinations = nchoosek(1:num_neurons, n_subset);
    catch
        % Initialize the output matrix
        all_combinations = [];

        % Generate random combinations
        for i = 1:num_combinations
            all_combinations(i, :) = sort(randperm(num_neurons, n_subset));
        end
        % Remove any duplicate combinations
        all_combinations = unique(all_combinations, 'rows');

        pffft=1;
    end

    selected_neurons_best = []; % Initialize empty set of selected neurons
    best_mi = -inf;
    selected_neurons_worst = []; % Initialize empty set of selected neurons
    worst_mi = inf;

    for ii_combs = 1:size(all_combinations,1)
        current_subset = all_combinations(ii_combs,:); % Add current neuron to the subset
        combined_activity = neural_data(:, current_subset); % Combine selected neuron activities

        % Compute mutual information between combined activity and odor plume
        mi_combined = drg_compute_mutual_information_gpu(combined_activity, odor_plume);

        if mi_combined > best_mi
            best_mi = mi_combined;
            selected_neurons_best=current_subset;
        end
        if mi_combined < best_mi
            worst_mi = mi_combined;
            selected_neurons_worst=current_subset;
        end
    end

end

function [most_subset,mi_most_imp,least_subset,mi_least_imp]...
    = imps_mi_selection_cpu(neural_data, odor_plume, n_subset,imp_ROIs,imps)
    % Inputs:
    % neural_data: matrix of size [num_samples x num_neurons], where each column represents a neuron's activity.
    % odor_plume: vector of size [num_samples x 1], representing the odor plume signal.
    % max_neurons: maximum number of neurons to select.
    %
    % Output:
    % selected_neurons: indices of the selected neurons.

    num_neurons = size(neural_data, 2);
    num_combinations=500000; %Number of random combinations if nchooseek crashes

    if n_subset>num_neurons
        n_subset=num_neurons;
    end

    %Sort ROIs according to predictive importance
    to_sort=[imps';imp_ROIs]';
    sorted_rows=sortrows(to_sort);

    %most
    most_subset=sorted_rows(end-n_subset+1:end,2);
    combined_activity = neural_data(:, most_subset); % Combine selected neuron activities

    % Compute mutual information between combined activity and odor plume
    mi_most_imp = drg_compute_mutual_information_cpu(combined_activity, odor_plume);

    %least
    least_subset=sorted_rows(1:n_subset,2);
    combined_activity = neural_data(:, least_subset); % Combine selected neuron activities

    % Compute mutual information between combined activity and odor plume
    mi_least_imp = drg_compute_mutual_information_cpu(combined_activity, odor_plume);


end


function mi = drg_compute_mutual_information_gpu(data, target)
    % Move data to GPU
    data = gpuArray(data);
    target = gpuArray(target);
    
    [n_timepoints, n_bits] = size(data);
    
    % Precompute indices
    target_indices = target + 1;
    data_indices = data + 1;
    
    % Compute cumulative counts using GPU operations
    cum_bindFF = zeros(n_bits, 2, 'gpuArray');
    cum_op_bin = zeros(1, 2, 'gpuArray');
    cum_op_bindFF_bin = zeros(n_bits, 2, 2, 'gpuArray');
    
    % Use accumarray for faster counting
    cum_op_bin = accumarray(target_indices, ones(n_timepoints, 1, 'gpuArray'), [2, 1]);
    
    for jj_ROI = 1:n_bits
        cum_bindFF(jj_ROI, :) = accumarray(data_indices(:, jj_ROI), ones(n_timepoints, 1, 'gpuArray'), [2, 1]);
        cum_op_bindFF_bin(jj_ROI, :, :) = accumarray([data_indices(:, jj_ROI), target_indices], ones(n_timepoints, 1, 'gpuArray'), [2, 2]);
    end
    
    % Calculate probabilities
    p_bindFF = cum_bindFF / n_timepoints;
    p_op_bin = cum_op_bin / n_timepoints;
    p_op_bindFF_bin = cum_op_bindFF_bin / n_timepoints;
    

    % Calculate mutual information
    mi=0;

    %Calculate info for op_bin
    for ii_op=1:2
        for ii_bits=1:n_bits
            for ii_bin_dFF=1:2
                if p_op_bindFF_bin(ii_bits,ii_bin_dFF,ii_op)~=0
                    mi=mi+p_op_bindFF_bin(ii_bits,ii_bin_dFF,ii_op)*...
                        log2(p_op_bindFF_bin(ii_bits,ii_bin_dFF,ii_op)/(p_op_bin(ii_op)*p_bindFF(ii_bits,ii_bin_dFF)));
                end
            end
        end
    end

    % Move result back to CPU
    mi = gather(mi);
    gpuDevice(1).reset();
end


function mi = drg_compute_mutual_information_cpu(data, target)
    %Computations with CPU
    % data = gpuArray(data);
    % target = gpuArray(target);
    
    [n_timepoints, n_bits] = size(data);
    
    % Precompute indices
    target_indices = target + 1;
    data_indices = data + 1;
    
    % Compute cumulative counts using CPU operations
    cum_bindFF = zeros(n_bits, 2);
    cum_op_bin = zeros(1, 2);
    cum_op_bindFF_bin = zeros(n_bits, 2, 2);
    
    % Use accumarray for faster counting
    cum_op_bin = accumarray(target_indices, ones(n_timepoints, 1), [2, 1]);
    
    for jj_ROI = 1:n_bits
        cum_bindFF(jj_ROI, :) = accumarray(data_indices(:, jj_ROI), ones(n_timepoints, 1), [2, 1]);
        cum_op_bindFF_bin(jj_ROI, :, :) = accumarray([data_indices(:, jj_ROI), target_indices], ones(n_timepoints, 1), [2, 2]);
    end
    
    % Calculate probabilities
    p_bindFF = cum_bindFF / n_timepoints;
    p_op_bin = cum_op_bin / n_timepoints;
    p_op_bindFF_bin = cum_op_bindFF_bin / n_timepoints;
    

    % Calculate mutual information
    mi=0;

    %Calculate info for op_bin
    for ii_op=1:2
        for ii_bits=1:n_bits
            for ii_bin_dFF=1:2
                if p_op_bindFF_bin(ii_bits,ii_bin_dFF,ii_op)~=0
                    mi=mi+p_op_bindFF_bin(ii_bits,ii_bin_dFF,ii_op)*...
                        log2(p_op_bindFF_bin(ii_bits,ii_bin_dFF,ii_op)/(p_op_bin(ii_op)*p_bindFF(ii_bits,ii_bin_dFF)));
                end
            end
        end
    end

 
end




pffft=1;