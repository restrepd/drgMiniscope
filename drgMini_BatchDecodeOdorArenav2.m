%drgBatchDecodeOdorArenav2
[choiceFileName,choiceBatchPathName] = uigetfile({'drgOdorArenaChoices_*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgMini_BatchDecodeOdorArenav2 run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];
 
addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
handles.no_files=length(handles.dFF_file);


first_file=handles.first_file;
first_run=handles.first_run;
handles_choices.save_path=handles.save_path;

%Parallel batch processing for each file
all_files_present=1;
for filNum=first_file:handles.no_files
       
    
    %Make sure that all the files exist
    dFF_file=handles.dFF_file{filNum};
    if iscell(handles.this_path)
        this_path=handles.this_path{filNum};
    else
        this_path=handles.this_path;
    end
     
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


%Process each file
tic
if all_files_present==1


    %Process each file separately
    for fileNo=first_file:handles.no_files
        if sum(handles.process_these_groups==handles.group(fileNo))>0

            first_toc=toc;
            handles_choices.this_path=handles.this_path{fileNo};
            handles_choices.dFF_file=handles.dFF_file{fileNo};
            handles_choices.arena_file=handles.arena_file{fileNo};
%             handles_choices.training_fraction=handles.training_fraction;
            handles_choices.bins_after=handles.bins_after;
            handles_choices.bins_current=handles.bins_current;
%             handles_choices.resume_processing=handles.resume_processing;
            handles_choices.dt=handles.dt;
            handles_choices.dt_miniscope=handles.dt_miniscope;
            handles_choices.n_shuffle=handles.n_shuffle;
            handles_choices.is_sphgpu=handles.is_sphgpu;
            handles_choices.z=handles.z;
            handles_choices.trial_start_offset=handles.trial_start_offset; %This was -10
            handles_choices.trial_end_offset=handles.trial_end_offset;
            handles_choices.group=handles.group(fileNo);
            handles_choices.save_results=handles.save_results;
            handles_choices.algo=handles.algo;

%             handles_choices.which_training_algorithm=algo;

            for ii_run=first_run:length(handles.bins_before)
                handles_choices.bins_before=handles.bins_before(ii_run);
                handles_choices.save_tag=handles.save_tag{ii_run};
                
                  %Run if overwrite is on or overwrite off and the file does not
                %exist
                if handles.overwrite_output==1
                    process_data=1;
                else
                    %Overwrite off, does the file exist?
                    if exist([handles_choices.save_path handles_choices.arena_file(1:end-4) handles_choices.save_tag '.mat'],'file')==2
                        process_data=0;
                    else
                        process_data=1;
                    end
                end
 
                if sum((logical(fileNo==handles.skip_these_files))&(logical(ii_run==handles.skip_these_runs)))>0
                    process_data=0;
                end

                if process_data==1
                    fprintf(1, ['Starting data fprocessing for file number %d, run %d\n'],fileNo,ii_run);

                    handles_out=drgMini_DecodeOdorArena_gpuv2(handles_choices);

                    fprintf(1, ['Data processed for file number %d, run %d\n'],fileNo,ii_run);

                else
                    fprintf(1, ['Data not processed for file number %d, run %d because output exists or was skipped\n'],fileNo, ii_run);
                end
                
            end
            first_run=1;
        end
    end
%     fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end



pffft=1;