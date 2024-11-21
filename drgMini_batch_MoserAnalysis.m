%drgMini_batch_MoserAnalysis
[choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniMoserChoices_*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgMini_batch_dFFPrediction run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];
 
addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
handles.no_files=length(handles.dFF_file);


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

    fprintf(1, ['Started processing Moser Analysis\n']);
    %Process each file separately
    for fileNo=first_file:handles.no_files
        

            first_toc=toc;
            handles_choices.this_path=handles.this_path{fileNo};
            handles_choices.dFF_file=handles.dFF_file{fileNo};
            handles_choices.arena_file=handles.arena_file{fileNo};
            handles_choices.group=handles.group(fileNo);
            handles_choices.alpha=handles.alpha;
            handles_choices.weber_fechner=handles.weber_fechner;
            handles_choices.multiplier=handles.multiplier;
            handles_choices.algo=handles.algo; %1 ann, 2 glm
            handles_choices.max_overlap=handles.max_overlap; %Maximum overlap of the shuffled segments
            handles_choices.group=handles.group(fileNo);
            handles_choices.dt=handles.dt;
            handles_choices.dt_miniscope=handles.dt_miniscope;
            handles_choices.displayFigures=handles.displayFigures;
            handles_choices.dt_decoding_op=handles.dt_decoding_op;
            handles_choices.dt_decoding_xy=handles.dt_decoding_xy;
            handles_choices.save_tag=handles.save_tag;
            handles_choices.cm_from_floor=handles.cm_from_floor(fileNo);
            handles_choices.resume_processing=handles.resume_processing;
            handles_choices.n_shuffle_SI=handles.n_shuffle_SI;
            handles_choices.z=handles.z;
            handles_choices.sh_repeats=handles.sh_repeats;
            handles_choices.max_overlap=handles.max_overlap;
            handles_choices.is_sphgpu=handles.is_sphgpu;
            handles_choices.prctile_thr=handles.prctile_thr;
            handles_choices.binary_dFF=handles.binary_dFF;


            fprintf(1, ['Started processing file number %d\n'],fileNo);
           
            handles_out=drgMini_MoserAnalysis(handles_choices);
              

            fprintf(1, ['Data processed for file number %d\n'],fileNo);


            fprintf(1,'\n\nProcessing time for file No %d is %d hours\n\n',fileNo,(toc-first_toc)/(60*60));

      
    end
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end



pffft=1;