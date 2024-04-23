%drgMini_batch_analysis_dFFprediction
close all
clear all

[choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniPredChoices_*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgBatchDecodeOdorArenav2 run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];
 
addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
handles.no_files=length(handles.dFF_file);

group_legend{1}='Odor plume';
group_legend{2}='Spatial';

figNo=0;

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

    %Setup histo figure
    edges=[0:0.05:1];
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .7 .7])
    groups=unique(handles.group);
    no_subplots=max(groups)*2*length(handles.alpha);
    ii_subplot=0;

    %Process each file separately

    for ii_wf=1:length(handles.alpha)



        R_per_group=[];
        for grNo=groups
            R_per_group(grNo).Rops=[];
            R_per_group(grNo).Rops_above_95=[];
            R_per_group(grNo).Rxys=[];
            R_per_group(grNo).Rxys_above_95=[];
            R_per_group(grNo).noROIs=0;
            R_per_group(grNo).noRops_above_95=0;
            R_per_group(grNo).noRxys_above_95=0;
        end

        for fileNo=first_file:handles.no_files

            this_path=handles.this_path;
            arena_file=handles.arena_file{fileNo};
            dFF_file=handles.dFF_file{fileNo};
            handles_choices.algo=2;
            handles_choices.alpha=handles.alpha(ii_wf);
            handles_choices.weber_fechner=handles.weber_fechner(ii_wf);
            handles_choices.group=handles.group(fileNo);

            load([this_path arena_file(1:end-4) 'decdFFa' num2str(handles_choices.algo) ...
                'wf' num2str(handles_choices.weber_fechner) 'g' num2str(handles_choices.group) ...
                'a' num2str(handles_choices.alpha) '.mat'])

            %Calculate 95th percentile
            noROIs=length(handles_out.per_ROI);
            if isfield(handles_out.per_ROI(1),'sh_repeats')
                noshRepeats=length(handles_out.per_ROI(1).sh_repeats);
            else
                noshRepeats=length(handles_out.per_ROI(1).repeats);
            end

            shRops=[];
            shRxys=[];
            for ii_ROI=1:noROIs
                for ii_sh_repeats=1:noshRepeats
                    shRops=[shRops handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rop];
                    shRxys=[shRxys handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rxy];
                end
            end

            % Bootstrap resampling op
            bootstrap_samples_op = bootstrp(1000, @prctile, shRops, 95);

            % Compute the 95th percentile of the bootstrap samples
            Ropsh_bootstrap_95th_percentile = prctile(bootstrap_samples_op, 95);

            % Bootstrap resampling xy
            bootstrap_samples_xy = bootstrp(1000, @prctile, shRxys, 95);

            % Compute the 95th percentile of the bootstrap samples
            Rxysh_bootstrap_95th_percentile = prctile(bootstrap_samples_xy, 95);

            %Find the Rops that are above 95th percentile
            noRepeats=length(handles_out.per_ROI(1).repeats);

            Rops=[];
            Rxys=[];
            Rops_above_95=[];
            Rxys_above_95=[];
            for ii_ROI=1:noROIs
                for ii_repeats=1:noRepeats
                    Rops=[Rops handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rop];
                    Rxys=[Rxys handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rxy];
                    if handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rop>Ropsh_bootstrap_95th_percentile
                        Rops_above_95=[Rops_above_95 handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rop];
                    end
                    if handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rxy>Rxysh_bootstrap_95th_percentile
                        Rxys_above_95=[Rxys_above_95 handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rxy];
                    end
                end
            end
 
            grNo= handles_choices.group;
            R_per_group(grNo).Rops_above_95=[R_per_group(grNo).Rops_above_95 Rops_above_95];
            R_per_group(grNo).Rxys_above_95=[ R_per_group(grNo).Rxys_above_95 Rxys_above_95];
            R_per_group(grNo).noROIs=R_per_group(grNo).noROIs+length(Rops);
            R_per_group(grNo).noRops_above_95=R_per_group(grNo).noRops_above_95+length(Rops_above_95);
            R_per_group(grNo).noRxys_above_95=R_per_group(grNo).noRxys_above_95+length(Rxys_above_95);

            pffft=1;
        end
        pffft=1;


        %Display histograms

        for grNo=groups
            %Display the histogram for xy

            subplot(length(handles.alpha), max(groups)*2, ii_subplot+2*(grNo-1)+1)

            histogram(R_per_group(grNo).Rxys_above_95,edges)
            xlabel('Rho xy')
            ylabel('No ROIs')
            xlim([0 0.5])
            title(['Rho xy ' group_legend{grNo} ' log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha)])

            %Display the histogram for op
            subplot(length(handles.alpha), max(groups)*2, ii_subplot+2*(grNo-1)+2)

            histogram(R_per_group(grNo).Rops_above_95,edges)
            xlabel('Rho xy')
            ylabel('No ROIs')
            xlim([0 0.5])
            title(['Rho op ' group_legend{grNo} ' log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha)])

        end
        ii_subplot=ii_subplot+max(groups)*2;
    end
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end



pffft=1;