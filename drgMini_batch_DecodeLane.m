%drgMini_batch_DecodeLane
close all
clear all

[choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniLanePredChoices_*.m'},'Select the .m file with all the choices for analysis');


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



    %Process each file separately
    for fileNo=first_file:handles.no_files
        %Load prediction file and define place and odor cells
        % load([handles.this_path{fileNo} handles.pred_file{fileNo}],'handles_out')

        all_handles=[];

        dFF_file=handles.dFF_file{fileNo};

        this_path=handles.this_path{fileNo};

        if handles.resume_processing==1
            load([this_path dFF_file(1:end-4) handles.suffix])
            first_file=length(all_handles.file);
        end

        
        dFF=[];
        if strcmp(dFF_file(end-3:end),'.mat')
            %This reads the extract file
            load([this_path dFF_file])
            dFF=zeros(size(output.temporal_weights,1),size(output.temporal_weights,2));
            for traceNo=1:size(output.temporal_weights,2)
                dFF(:,traceNo)=output.temporal_weights(:,traceNo);
            end
        else
            if strcmp(dFF_file(end-4:end),'.hdf5')
                %This is an hdf5 generated by CaImAn
                dFF=h5read([this_path dFF_file],'/estimates/F_dff' );
            else
                %This is a csv file created from ImageJ
                dFF=readmatrix([this_path dFF_file]);
            end
        end


        no_neurons=size(dFF,2);
        dFF=[];


        %
        % spatial_rhol1l4=handles_out.spatial_rhol1l4;
        % no_neurons=length(spatial_rhol1l4);

        if handles.resume_processing==1
            first_ROI=length(all_handles.file(fileNo).run);
            handles.resume_processing=0;
        else
            first_ROI=1;
        end

        for this_ROI=first_ROI:no_neurons

            first_toc=toc;
            handles_choices.this_path=handles.this_path{fileNo};
            handles_choices.dFF_file=handles.dFF_file{fileNo};
            handles_choices.arena_file=handles.arena_file{fileNo};
            try
                handles_choices.pred_file=handles.pred_file{fileNo}; %This is only nescessary for which_ROIs=2 or 3
            catch
            end
            handles_choices.process_these_ROIs=[this_ROI];
            handles_choices.which_ROIs=handles.which_ROIs; %when this function is called the user has to specify handles_choices2.process_these_ROIs
            handles_choices.bins_before=handles.bins_before;
            handles_choices.bins_current=handles.bins_current;
            handles_choices.bins_after=handles.bins_after;
            handles_choices.no_repeats=handles.no_repeats;
            handles_choices.dt=handles.dt;
            handles_choices.dt_miniscope=handles.dt_miniscope;
            handles_choices.n_shuffle=handles.n_shuffle;
            handles_choices.which_ml_algo=handles.which_ml_algo;
            handles_choices.align_training_start=handles.align_training_start;
            handles_choices.dt_training_start=handles.dt_training_start; %seconds from start alignment -5
            handles_choices.dt_training_end=handles.dt_training_end; %seconds from end alignment 5
            handles_choices.which_shuffle=handles.which_shuffle;

            handles_choices.is_sphgpu=handles.is_sphgpu;

            fprintf(1, ['Started processing file number %d, ROI No %d\n'],fileNo, this_ROI);

            all_handles.run(this_ROI).handles_out=drgMini_DecodeLaneOdorArenav6(handles_choices);

            save([this_path dFF_file(1:end-4) handles.suffix],'all_handles')

            fprintf(1, ['Data processed for file number %d, ROI No %d\n'],fileNo,this_ROI);


            fprintf(1,'\n\nProcessing time for file No %d is %d hours\n\n',fileNo,(toc-first_toc)/(60*60));

        end
    end
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end



pffft=1;