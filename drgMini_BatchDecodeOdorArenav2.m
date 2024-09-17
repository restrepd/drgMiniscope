%drgBatchDecodeOdorArenav2
[choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniChoices_*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgBatchDecodeOdorArenav2 run for ' choiceFileName '\n\n']);

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

%Now run the decoder for all files and all choices of algorithms
pffft=1;

handles_out=[];

%Process each file
tic
if all_files_present==1


    %Process each file separately
    for fileNo=first_file:handles.no_files
        for algo=handles.which_training_algorithm

            first_toc=toc;
            handles_choices.this_path=handles.this_path;
            handles_choices.dFF_file=handles.dFF_file{fileNo};
            handles_choices.arena_file=handles.arena_file{fileNo};
            handles_choices.training_fraction=handles.training_fraction;
            handles_choices.bins_before=handles.bins_before;
            handles_choices.bins_current=handles.bins_current;
            handles_choices.bins_after=handles.bins_after;
            handles_choices.dt=handles.dt;
            handles_choices.dt_miniscope=handles.dt_miniscope;
            handles_choices.n_shuffle=handles.n_shuffle;
            handles_choices.which_training_algorithm=algo;

            fprintf(1, ['Started processing file number %d, algorithm number %d\n'],fileNo,algo);

            handles_out.file(fileNo).algo(algo).handles_out=drgDecodeOdorArenav2(handles_choices);

            fprintf(1, ['Data processed for file number %d, algorithm number %d\n'],fileNo,algo);


            fprintf(1,'\n\nProcessing time for file No %d is %d hours\n\n',fileNo,(toc-first_toc)/(60*60));

            %Save output file
            handles_out.last_file_processed=fileNo;
            handles_out.handles=handles;
            
            save([handles.this_path handles.outputFile],'handles_out','-v7.3')
        end
    end
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end



pffft=1;