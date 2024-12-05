%drgMini_batch_DecodeOdorConc
[choiceFileName,choiceBatchPathName] = uigetfile({'drgOdorConcChoices_*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgMini_batch_DecodeOdorConc run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];
 
addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
handles.no_files=length(handles.dFF_file);


first_file=handles.first_file;
first_run=handles.first_run;

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
handles_choices.save_path='/data2/SFTP/DecodeOdorConc/';
if all_files_present==1

    fprintf(1, ['Started processing odor concentration decoding\n']);
    %Process each file separately
    for fileNo=first_file:handles.no_files
        

            first_toc=toc;
            handles_choices.this_path=handles.this_path{fileNo};
            handles_choices.dFF_file=handles.dFF_file{fileNo};
            handles_choices.arena_file=handles.arena_file{fileNo};
            handles_choices.group=handles.group(fileNo);
    
            handles_choices.multiplier=handles.multiplier;
            handles_choices.lowest_conc=handles.lowest_conc;
            
            %Hill transform
            handles_choices.hill=handles.hill; %0=no Hill transform, 1=Hill transform
            handles_choices.k_half=handles.k_half; %Hill equation K1/2
            % handles_choices.actual_maxC=(handles_choices.k_half^handles_choices.n_hill); %Measured from the simulated data 0.0043

            handles_choices.maxC=handles.maxC;%maxC=0.0043 maximum of simulated odor plume
            handles_choices.n_hill=handles.n_hill; %Hill coefficient
            % handles_choices.max_overlap=handles.max_overlap; %Maximum overlap of the shuffled segments

            handles_choices.group=handles.group(fileNo);
            handles_choices.dt=handles.dt;
            handles_choices.dt_miniscope=handles.dt_miniscope;
            handles_choices.displayFigures=handles.displayFigures;
            % handles_choices.dt_decoding_op=handles.dt_decoding_op;
            % handles_choices.dt_decoding_xy=handles.dt_decoding_xy;
            
            handles_choices.cm_from_floor=handles.cm_from_floor(fileNo);
            handles_choices.resume_processing=handles.resume_processing;
            % handles_choices.n_shuffle_SI=handles.n_shuffle_SI;
            handles_choices.z=handles.z;
            handles_choices.n_shuffle=handles.n_shuffle;
            % handles_choices.max_overlap=handles.max_overlap;
            handles_choices.is_sphgpu=handles.is_sphgpu;
            % handles_choices.prctile_thr=handles.prctile_thr;
            % handles_choices.binary_dFF=handles.binary_dFF;
            handles_choices.trial_start_offset=handles.trial_start_offset; %This was -10
            handles_choices.trial_end_offset=handles.trial_end_offset;

            handles_choices.bins_after=handles.bins_after;
            handles_choices.bins_current=handles.bins_current;

            for ii_run=first_run:length(handles.weber_fechner)
                handles_choices.alpha=handles.alpha(ii_run);
                handles_choices.weber_fechner=handles.weber_fechner(ii_run);
                handles_choices.bins_before=handles.bins_before(ii_run);
                handles_choices.save_tag=handles.save_tag{ii_run};
                handles_choices.which_training_algorithm=handles.algo(ii_run); %1 ann, 2 glm

                fprintf(1, ['Started processing file number %d run %d\n'],fileNo, ii_run);

                %Added this try/catch here because of "Error using
                %parallel.Future/fetchOutputs"
                is_done=0;
                ii_failures=0;
                while is_done==0
                    try
                        handles_out=drgMini_DecodeOdorConc(handles_choices);
                        is_done=1;
                    catch
                        ii_failures=ii_failures+1;
                        fprintf(1, ['Failure number %d for file number %d run %d\n'],ii_failures,fileNo, ii_run);
                    end
                end


                fprintf(1, ['Data processed for file number %d run %d\n'],fileNo, ii_run);
            end

            first_run=1;

            fprintf(1,'\n\nProcessing time for file No %d is %d hours\n\n',fileNo,(toc-first_toc)/(60*60));

      
    end
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end



pffft=1;