%drgMini_batch_DecodeOdorPlume
close all
clear all

run_slurm=1;

algo_name{1}='SVZ';
algo_name{2}='NN';
algo_name{3}='Tree';
algo_name{4}='Bayes';
algo_name{5}='glm';
algo_name{6}='LD';
algo_name{7}='NNopt';
 
training_labels{1}='all';
training_labels{2}='hits';
training_labels{3}='misses';

if run_slurm==1
    
else
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniOPPredChoices_*.m'},'Select the .m file with all the choices for analysis');
end


fprintf(1, ['\ndrgMini_batch_dFFPrediction run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
handles.no_files=length(handles.dFF_file);

if handles.is_sphgpu==1
    addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
    addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
    addpath('/home/restrepd/Documents/MATLAB/drgMaster')
    addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end


first_file=1;

%Parallel batch processing for each file

%Please note that the file names have to be in the same sequence as handles_conc
addpath(choiceBatchPathName)
eval(['handles_conc=' handles.choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_Angle=' handles.choiceAngleFileName(1:end-2) ';'])
eval(['handles_XY=' handles.choiceXYFileName(1:end-2) ';'])
%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,handles.save_PathAngle, handles_conc, handles.save_PathConc);

all_files_present=1;
these_groups=[1 5];
for fileNo=1:length(handles.dFF_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)&(~isempty(handles.dFF_file{fileNo}))

        %Make sure that all the files exist
        dFF_file=handles.dFF_file{fileNo};

        this_path=handles.this_path{fileNo};
        % else
        %     this_path=handles.this_path;
        % end

        if exist([this_path dFF_file])==0
            fprintf(1, ['Program will be terminated because file No %d, ' dFF_file ' does not exist\n'],fileNo);
            all_files_present=0;
        end

        arena_file=handles.arena_file{fileNo};

        if exist([this_path arena_file])==0
            fprintf(1, ['Program will be terminated because file No %d, ' arena_file ' does not exist\n'],fileNo);
            all_files_present=0;
        end
    end
end

%Now run the decoder for all files and all choices of algorithms
tic
ii_run=1;
if all_files_present==1

    if ~exist(handles.save_path(1:end-1),'dir')
        mkdir(handles.save_path(1:end-1))
    end

    if handles.resume_processing==1
        load([handles.save_path handles.save_file])
        first_file=length(all_handles.file); 
        resume_processing=1;
    else
        all_handles=[];
        resume_processing=0;
    end

    %Process each file separately
    for fileNo=first_file:length(handles.dFF_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)&(~isempty(handles.dFF_file{fileNo}))
            %Load prediction file and define place and odor cells
            % load([handles.this_path{fileNo} handles.pred_file{fileNo}],'handles_out')

            proceed=1;
            first_ml_algo=1;
            first_training_range=1;
            no_ROIs_per_subset=handles.no_ROIs_per_subset;

            if (resume_processing==1)&(fileNo==first_file)
                %Determine where to start
                proceed=0;
                if ~isfield(all_handles,'file')
                    proceed=1;
                else
                    if size(all_handles.file,2)<first_file
                        proceed=1;
                    else
                        if ~isfield(all_handles.file(fileNo),'ml_algo')
                            proceed=1;
                        else
                            last_ml_algo=size(all_handles.file(fileNo).ml_algo,2);
                            if ~isfield(all_handles.file(fileNo).ml_algo(last_ml_algo),'training')
                                proceed=1;
                                first_ml_algo=last_ml_algo;
                            else
                                last_training=size(all_handles.file(fileNo).ml_algo(last_ml_algo).training,2);
                                if ~isfield(all_handles.file(fileNo).ml_algo(last_ml_algo).training,'subset')
                                    proceed=1;
                                    first_ml_algo=last_ml_algo;
                                    first_training_range=last_training;
                                else
                                    last_subset=size(all_handles.file(fileNo).ml_algo(last_ml_algo).training(last_training).subset,2);
                                    if last_subset==length(handles.no_ROIs_per_subset)
                                        %All subsets were processed
                                        proceed=1;
                                        first_ml_algo=last_ml_algo;
                                        first_training_range=last_training+1;
                                        no_ROIs_per_subset=handles.no_ROIs_per_subset;
                                    else
                                        %Some subsets are missing in this
                                        %training set
                                        proceed=1;
                                        first_ml_algo=last_ml_algo;
                                        first_training_range=last_training;
                                        no_ROIs_per_subset=handles.no_ROIs_per_subset(last_subset+1:end);
                                    end
                                end

                            end
                        end
                    end
                end
            end

            dFF_file=handles.dFF_file{fileNo};

            this_path=handles.this_path{fileNo};

            %Get XY and dFF per trial
            arena_file=handles_XY.arena_file{fileNo};
            load([handles.save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
            trials=handles_out.trials;
            no_neurons=handles_out.no_neurons;

            all_handles.file(fileNo).no_neurons=no_neurons;



            %These are used for nchoosek calculations
            rng(fileNo)

            figNo=0;
            for which_ml_algo=first_ml_algo:3
                %1 SVZ
                %2 nn
                %3 tree
                %4 bayesian
                %5 glm
                %6 linear
                %7 NN parameter optimized
                for which_training_range=first_training_range:3
                    %1 all trials
                    %hits
                    %misses
                    
                    for ii_neurons=no_ROIs_per_subset
                        ii_subset=find(ii_neurons==handles.no_ROIs_per_subset);
                        if ii_neurons>no_neurons
                            n_subset=no_neurons;
                        else
                            n_subset=ii_neurons;
                        end

                        try
                            all_combinations = nchoosek(1:no_neurons, n_subset);
                            if size(all_combinations,1)>handles.max_nchoosek
                                % Initialize the output matrix
                                all_combinations = [];

                                % Generate random combinations
                                for i = 1:handles.max_nchoosek
                                    all_combinations(i, :) = sort(randperm(no_neurons, n_subset));
                                end
                                % Remove any duplicate combinations
                                all_combinations = unique(all_combinations, 'rows');
                            end
                        catch
                            % Initialize the output matrix
                            all_combinations = [];

                            % Generate random combinations
                            for i = 1:handles.max_nchoosek
                                all_combinations(i, :) = sort(randperm(no_neurons, n_subset));
                            end
                            % Remove any duplicate combinations
                            all_combinations = unique(all_combinations, 'rows');
                        end
                        fprintf(1, ['Started processing file number %d, algo %d,training range %d, subset %d\n']...
                            ,fileNo, which_ml_algo, which_training_range, ii_subset);
                        first_toc=toc;
                        all_mean_hit_after=[];
                        all_mean_miss_after=[];
                        for ii_comb=1:size(all_combinations,1)

                            handles_choices.this_path=handles.this_path{fileNo};
                            handles_choices.dFF_file=handles.dFF_file{fileNo};
                            handles_choices.arena_file=handles.arena_file{fileNo};
                            arena_file=handles_Angle.arena_file{fileNo};
                            handles_choices.angle_file=[arena_file(1:end-4) handles_Angle.save_tag '.mat'];
                            handles_choices.save_PathAngle=handles.save_PathAngle;
                            handles_choices.xy_arena_file=[arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'];
                            handles_choices.save_PathXY=handles.save_PathXY;
                            handles_choices.dt=handles.dt;
                            handles_choices.n_shuffle=handles.n_shuffle;
                            handles_choices.which_training_range=which_training_range;
                            handles_choices.which_ml_algo=which_ml_algo;
                            handles_choices.ii_cost=handles.ii_cost;
                            handles_choices.display_figures=0;
                            this_combination=zeros(1,size(all_combinations,2));
                            this_combination(1,:)=all_combinations(ii_comb,:);
                            handles_choices.neurons_included=this_combination;
                            handles_choices.is_gpu=handles.is_gpu;
                            handles_choices.handles_conc=handles_conc;

                            handles_out2=drgMini_DecodeOdorPlume(handles_choices);

                            all_handles.file(fileNo).ml_algo(which_ml_algo).training(which_training_range).subset(ii_subset).ROI_combination(ii_comb).handles_out2=handles_out2;
                            all_mean_hit_after=[all_mean_hit_after handles_out2.mean_accuracy_hit_after];
                            all_mean_miss_after=[all_mean_miss_after handles_out2.mean_accuracy_miss_after];

                            save([handles.save_path handles.save_file],'all_handles')
                        end
                        [max_hit, max_ii]=max(all_mean_hit_after);

                        fprintf(1, ['Data processed fileNo %d, ' algo_name{which_ml_algo} ' ' training_labels{which_training_range}...
                            ' ROI subset %d\n'],fileNo,ii_subset);
                        fprintf(1, ['Max hits accuracy %d, miss %d\n'],...
                            max_hit,all_mean_miss_after(max_ii));

                        fprintf(1,['Processing time %d mins\n\n'],(toc-first_toc)/(60));

                        figNo=figNo+1;
                        try
                            close(figNo)
                        catch
                        end

                        hFig = figure(figNo);

                        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


                        hold on
                        plot(all_mean_miss_after,all_mean_hit_after,'ob')
                        plot([0 1],[0.5 0.5],'-k')
                        plot([0.5 0.5],[0 1],'-k')
                        title(['FileNo ' num2str(fileNo) ', ' algo_name{which_ml_algo} ' ' training_labels{which_training_range}...
                            ' subset ' num2str(ii_subset)])
                        xlabel('accuracy (miss)')
                        ylabel('accuracy (hits)')
                        xlim([0 1])
                        ylim([0 1])

                        pffft=1;
                    end
                end
            end
        end
    end

end



pffft=1;