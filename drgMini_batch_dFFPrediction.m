%drgMini_batch_dFFPrediction
[choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniPredChoices_*.m'},'Select the .m file with all the choices for analysis');


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

%Process each file


tic
if all_files_present==1


    %Process each file separately
    for fileNo=first_file:handles.no_files
        for ii_wf=1:length(handles.alpha)

            first_toc=toc;
            handles_choices.this_path=handles.this_path;
            handles_choices.dFF_file=handles.dFF_file{fileNo};
            handles_choices.arena_file=handles.arena_file{fileNo};
            handles_choices.training_fraction=handles.training_fraction;
            handles_choices.alpha=handles.alpha(ii_wf);
            handles_choices.weber_fechner=handles.weber_fechner(ii_wf);
            handles_choices.repeats=1;
            handles_choices.sh_repeats=5;
            handles_choices.algo=2; %1 ann, 2 glm
            handles_choices.max_overlap=3; %Maximum overlap of the shuffled segments
            handles_choices.group=handles.group(fileNo);
       

            %Note: The data brought into the Kording lab jupyter notebbok seems to be
            %binned in 200 msec bins
            dt=0.2;
            dt_miniscope=1/30;
            n_shuffle=5; %Note that n_shuffle is changed to a maximum of ii_n_training

            handles_choices.dt=dt;
            handles_choices.dt_miniscope=dt_miniscope;
            handles_choices.n_shuffle=n_shuffle;

            %The user can define what time period to use spikes from (with respect to the output).
            bins_before=10; %How many bins of neural data prior to the output are used for decoding, 4
            bins_current=1; %Whether to use concurrent time bin of neural data, 1
            bins_after=0; %How many bins of neural data after the output are used for decoding, 10
            handles_choices.bins_before=bins_before;
            handles_choices.bins_current=bins_current;
            handles_choices.bins_after=bins_after;

 
            handles_choices.op_path=handles.op_path;
            handles_choices.op_file=handles.op_file;

            handles_choices.displayFigures=1;

            fprintf(1, ['Started processing file number %d, weber_fechner group %d\n'],fileNo, ii_wf);

            handles_out=drgMini_dFFprediction(handles_choices)

            fprintf(1, ['Data processed for file number %d, weber_fechner group %d\n'],fileNo,ii_wf);


            fprintf(1,'\n\nProcessing time for file No %d is %d hours\n\n',fileNo,(toc-first_toc)/(60*60));

            % %Save output file
            % handles_out.last_file_processed=fileNo;
            % handles_out.handles=handles;
            % 
            % save([handles.this_path handles.outputFile],'handles_out','-v7.3')
        end
    end
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end



pffft=1;