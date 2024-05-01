%drgMini_drocnsampler
close all
clear all

[choiceFileName,choiceBatchPathName] = uigetfile({'drgMiniDropcnsampler*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgMini_drocnsamplerrun for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles_choice=' choiceFileName(1:end-2) ';'])
handles_choice.choiceFileName=choiceFileName;
handles_choice.choiceBatchPathName=choiceBatchPathName;

if ischar(handles_choice.PathName)
    cd(handles_choice.PathName)
end

%Parallel batch processing for each file
all_files_present=1;
for fileNo=handles_choice.first_file:handles_choice.no_files


    if exist([handles_choice.PathName{fileNo} handles_choice.spmFileName{fileNo}])==0
        fprintf(1, ['Program will be terminated because file No %d, ' handles_choice.spmFileName{fileNo} ' does not exist\n'],fileNo);
        all_files_present=0;
    end

    %    if exist([handles_choice.PathName{filNum} handles_choice.spmFileName{filNum}])==0
    %     fprintf(1, ['Program will be terminated because file No %d, ' dFF_file ' does not exist\n'],filNum);
    %     all_files_present=0;
    % end


end


if all_files_present==1

    for fileNo=handles_choice.first_file:handles_choice.no_files

        figNo=0;
        
        mini_rate=handles_choice.freq_miniscope(fileNo);

        % fileNo=1;
        % handles_choice.PathName{fileNo}='/Users/restrepd/Documents/Projects/Basal_Forebrain/290424/';
        % handles_choice.spmFileName{fileNo}='KO_urine_042924_120240429T140816.mat';
        % handles_choice.rhdFileName{fileNo}='ko_urine_2_240429_140804.rhd';

        %Read the dropc file
        handles=[];
        load([handles_choice.PathName{fileNo} handles_choice.spmFileName{fileNo}])

        %Read the rhd file
        adc_in=[];
        digital_in=[];
        acq_rate=[];
        [adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file([handles_choice.PathName{fileNo} handles_choice.rhdFileName{fileNo}],1);

        %Plot digital
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .7 .2])

        plot(digital_in)
        title('Digital output')

        %Plot camera TTL
         figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .35 .7 .2])
        plot(adc_in)
        title('Camera TTL')

        %Find first and last TTL
        first_ttl=find(adc_in>1,1,'first');
        last_ttl=find(adc_in>1,1,'last');

        delta_t_movie=(last_ttl-first_ttl)/acq_rate;
        no_frames=delta_t_movie*mini_rate;

        %Now find the trials
        %4=Final valve
        %16=Odorant
        %After 16 2^Odor No, 2, 4, 8, 16, etc

        at_end=0;
        ii=0;

        handles_out.which_odor=[];
        handles_out.trialNo=0;
        handles_out.FVon_ii=[];
        handles_out.OdorOn_ii=[];
        handles_out.OdorOff_ii=[];

        while at_end==0
            %Find final valve (4) and exit if at end
            if isempty(find(digital_in(ii+1:end)==4,1,'first'))
                at_end=1;
            else
                delta_ii_next_FV=find(digital_in(ii+1:end)==4,1,'first');


                delta_ii_next_odor=find(digital_in(ii+delta_ii_next_FV:end)>4,1,'first');
                delta_ii_end_odor=find(digital_in(ii+delta_ii_next_FV+delta_ii_next_odor:end)<16,1,'first');
                delta_ii_which_odor=find(digital_in(ii+delta_ii_next_FV+delta_ii_next_odor+delta_ii_end_odor:end)>0,1,'first');

                handles_out.trialNo=handles_out.trialNo+1;
                handles_out.which_odor(handles_out.trialNo)=digital_in(ii+delta_ii_next_FV+delta_ii_next_odor+delta_ii_end_odor+delta_ii_which_odor);
                handles_out.FVon_ii(handles_out.trialNo)=ii+delta_ii_next_FV;
                handles_out.OdorOn_ii(handles_out.trialNo)=ii+delta_ii_next_FV+delta_ii_next_odor;
                handles_out.OdorOff_ii(handles_out.trialNo)=ii+delta_ii_next_FV+delta_ii_next_odor+delta_ii_end_odor;

                delta_ii_which_odor_end=find(digital_in(ii+delta_ii_next_FV+delta_ii_next_odor+delta_ii_end_odor+delta_ii_which_odor:end)==0,1,'first');

                ii=ii+delta_ii_next_FV+delta_ii_next_odor+delta_ii_end_odor+delta_ii_which_odor+delta_ii_which_odor_end;
            end

        end


        %First read the traces from the EXTRACT output
        this_filename=handles_choice.extFileName{fileNo};

        if strcmp(this_filename(end-3:end),'.mat')
            %This reads the extract file
            load([handles_choice.PathName{fileNo} handles_choice.extFileName{fileNo}])
            traces=zeros(size(output.temporal_weights,1),size(output.temporal_weights,2));
            for traceNo=1:size(output.temporal_weights,2)
                traces(:,traceNo)=output.temporal_weights(:,traceNo);
            end
        else
            if strcmp(this_filename(end-4:end),'.hdf5')
                %This is an hdf5 generated by CaImAn
                traces=h5read([handles_choice.PathName{fileNo} handles_choice.extFileName{fileNo}],'/estimates/F_dff' );
            else
                %This is a csv file created from ImageJ
                traces=readmatrix([handles_choice.PathName{fileNo} handles_choice.extFileName{fileNo}]);
            end
        end

        traces=traces';

        fnameca=handles_choice.extFileName{fileNo};

        sz_traces=size(traces);
        no_traces=sz_traces(1);
        no_images=sz_traces(2);

        %Now plot aligned traces and trials
          figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .85 .8])


        hold on

        % Determine the y spacing of the traces
        y_shift=4*(prctile(traces(:),95)-prctile(traces(:),5));

           %Plot the traces
           dt=1/mini_rate;
        time=[1:no_images]*dt;
        for trNo=1:no_traces
            % for trNo=1:20
            plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
        end

        %Now add aligned trials
        first_ttl_ii=find(adc_in>1,1,'first');

            % Plot the traces
        these_lines{1}='-b';
        these_lines{2}='-r';
        these_lines{3}='-m';
        these_lines{8}='-g';
        these_lines{5}='-y';
        these_lines{6}='-k';
        these_lines{7}='-c';
        these_lines{4}='-k';

        odors=unique(handles_out.which_odor);
        for trialNo=1:handles_out.trialNo

            %first plot odor on
            delta_ii_ttl_start_to_odor_on=handles_out.OdorOn_ii(trialNo)-first_ttl_ii;
            delta_t_ttl_start_to_odor_on=delta_ii_ttl_start_to_odor_on/acq_rate;
            delta_movie_time=delta_t_ttl_start_to_odor_on;

            this_color=find(handles_out.which_odor(trialNo)==odors);

            plot([delta_movie_time delta_movie_time], [0 (no_traces+2)*y_shift],...
                these_lines{this_color},'LineWidth',1)

            %now plot odor off
            delta_ii_ttl_start_to_odor_off=handles_out.OdorOff_ii(trialNo)-first_ttl_ii;
            delta_t_ttl_start_to_odor_off=delta_ii_ttl_start_to_odor_off/acq_rate;
            delta_movie_time_off=delta_t_ttl_start_to_odor_off;

            plot([delta_movie_time_off delta_movie_time_off], [0 (no_traces+2)*y_shift],...
                these_lines{this_color},'LineWidth',1)

        end

        fprintf(1, ['\ndrgCaImAn_dropc run for ' handles_choice.extFileName{fileNo} '\n\n']);

        this_rhd_file_name=handles_choice.rhdFileName{fileNo};
        save([handles_choice.op_path this_rhd_file_name(1:end-4) '_proc.mat'],'handles','handles_out','handles_choices','-v7.3')

    end
end

