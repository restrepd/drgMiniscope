%drgMini_DisplayDecodeOdorConcv3
clear all
close all

%Trained with hits only
save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

%This one has the dFF per trial
save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

% addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])


%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

%Now plot pseudocolor locations for last turn angles, start and end points,
%etc

%First make point maps
ii_run=1;
these_groups=[1 5];
percent_correct=[];
file_dates=[];

%Initialize 
for fileNo=1:length(handles_conc.arena_file)
    percent_correct.file(fileNo).percent_correct=0;
    percent_correct.file(fileNo).percent_correct1=0;
    percent_correct.file(fileNo).percent_correct4=0;
    percent_correct.file(fileNo).date=datetime('04012024', 'InputFormat', 'MMddyyyy');
    percent_correct.file(fileNo).mouse=0;
end

all_percents=[];

figNo=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        %Load conc data
        arena_file=handles_conc.arena_file{fileNo};
        %load the ouptut file
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])

        trials=handles_out.trials;
        odor_plume_template=handles_out.odor_plume_template;
        op_predicted=handles_out.op_predicted;
        op_predicted_sh=handles_out.op_predicted_sh;
        n_shuffle=handles_choices.n_shuffle;
        dt=handles_choices.dt;

        %calculate percent correct
        percent_correct.file(fileNo).percent_correct=100*(sum(trials.hit4)+sum(trials.hit1))/trials.odor_trNo;
        percent_correct.file(fileNo).percent_correct1=100*sum(trials.hit1)/(sum(trials.hit1)+sum(trials.miss1));
        percent_correct.file(fileNo).percent_correct4=100*sum(trials.hit4)/(sum(trials.hit4)+sum(trials.miss4));
        eval(['this_date=handles_conc.date{' num2str(fileNo) '};'])
        percent_correct.file(fileNo).date=datetime(this_date, 'InputFormat', 'MMddyyyy');
        percent_correct.file(fileNo).mouse=handles_conc.mouse(fileNo);
        all_percents=[all_percents 100*(sum(trials.hit4)+sum(trials.hit1))/trials.odor_trNo];

        %Load XY data
        arena_file=handles_XY.arena_file{fileNo};
        %load the ouptut file
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        x_predicted=handles_out.x_predicted;
        y_predicted=handles_out.y_predicted;
        x_predicted_sh=handles_out.x_predicted_sh(:,1);
        y_predicted_sh=handles_out.y_predicted_sh(:,1); %Note: all sh are identical
        XYtest=handles_out.XYtest;
        no_neurons=handles_out.no_neurons;


        no_conv_points=11;
        % conv_win=ones(1,no_conv_points)/no_conv_points;
        conv_win_gauss = gausswin(no_conv_points);
        conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

        op_predicted_conv=conv(op_predicted,conv_win_gauss,'same');

        %Now limit the x and y to max and min
        minop=min(odor_plume_template);
        op_predicted_conv(op_predicted_conv<minop)=minop;
        maxop=max(odor_plume_template);
        op_predicted_conv(op_predicted_conv>maxop)=maxop;


        op_predicted_sh_conv=zeros(size(op_predicted_sh,1),size(op_predicted_sh,2));

        for ii_sh=1:n_shuffle
            this_op_predicted_sh=zeros(length(op_predicted),1);
            this_op_predicted_sh(:,1)=op_predicted_sh(:,ii_sh);

            this_op_predicted_sh_conv=[];

            this_op_predicted_sh_conv=conv(this_op_predicted_sh,conv_win_gauss,'same');


            %Now limit the x and y to max and min
            minop=min(odor_plume_template);
            this_op_predicted_sh_conv(this_op_predicted_sh_conv<minop)=minop;
            maxop=max(odor_plume_template);
            this_op_predicted_sh_conv(this_op_predicted_sh_conv>maxop)=maxop;

            op_predicted_sh_conv(:,ii_sh)=this_op_predicted_sh_conv;
        end

        spout4_xy=[0 70];
        spout1_xy=[0 450];


        %Keep track of the per trial decoding
        op_all_trials=[];
        op_decod_all_trials=[];
        op_all_trials_sh=[];
        op_all_hits=[];
        op_decod_all_hits=[];
        op_all_hits_sh=[];
        op_decod_all_hits_sh=[];
        op_all_miss=[];
        op_decod_all_miss=[];
        op_all_miss_sh=[];
        op_decod_all_miss_sh=[];
        op_decod_all_trials_sh=[];
        op_between_trials=[];
        op_decod_between_trials=[];
        op_between_trials_sh=[];
        op_decod_between_trials_sh=[];

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


        hold on
        ii_start=0;

        last_op_predictedend=1;
        wall_threshold=50;

        for trNo=1:trials.odor_trNo

            op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
            op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
            if op_predictedend>length(op_predicted)
                op_predictedend=length(op_predicted);
            end

            op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];

            op_between_trials=[op_between_trials; odor_plume_template(last_op_predictedend:op_predictedstart)'];
            op_decod_between_trials=[op_decod_between_trials; op_predicted_conv(last_op_predictedend:op_predictedstart)];
            last_op_predictedend=op_predictedend;

            ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;



            %Okabe_Ito colors
            switch trials.odor_trial_type(trNo)
                case 1
                    %Lane 1 hits orange
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[230/255 159/255 0/255],'LineWidth',3)
                    plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];
                case 2
                    %Lane 1 miss sky blue
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[86/255 180/255 233/255],'LineWidth',3)
                    plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];

                case 3
                    %Lane 4 hit vermillion
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[213/255 94/255 0/255],'LineWidth',3)
                    plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_hits=[op_decod_all_hits; op_predicted_conv(op_predictedstart:op_predictedend)];
                case 4
                    %Lane 4 miss bluish green
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[0/255 158/255 115/255],'LineWidth',3)
                    plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss=[op_decod_all_miss; op_predicted_conv(op_predictedstart:op_predictedend)];
            end
            ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
        end

        minC=min(op_all_trials)-0.1*(max(op_all_trials)-min(op_all_trials));
        maxC=max(op_all_trials)+0.1*(max(op_all_trials)-min(op_all_trials));
        ylim([minC maxC])

        if figNo==22
            pffft=1;
        end

        this_xlim=xlim;
        divisor=6;
        text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1))/divisor,minC+0.95*(maxC-minC),'Hit lane 1','Color',[230/255 159/255 0/255],'FontSize',16,'FontWeight','bold')
        text(this_xlim(1)+1.7*(this_xlim(2)-this_xlim(1))/divisor,minC+0.95*(maxC-minC),'Miss lane 1','Color',[86/255 180/255 233/255],'FontSize',16,'FontWeight','bold')
        text(this_xlim(1)+2.7*(this_xlim(2)-this_xlim(1))/divisor,minC+0.95*(maxC-minC),'Hit lane 4','Color',[213/255 94/255 0/255],'FontSize',16,'FontWeight','bold')
        text(this_xlim(1)+3.7*(this_xlim(2)-this_xlim(1))/divisor,minC+0.95*(maxC-minC),'Miss lane 4','Color',[0/255 158/255 115/255],'FontSize',16,'FontWeight','bold')
        text(this_xlim(1)+4.7*(this_xlim(2)-this_xlim(1))/divisor,minC+0.95*(maxC-minC),'Predicted','Color',[0/255 0/255 0/255],'FontSize',16,'FontWeight','bold')


        title(['log10(odor) per trial, file No ' num2str(fileNo) ],'FontSize',16,'FontWeight','bold')
        ylabel('log10(odor) A.U.')
        xlabel('Time (s)')


     
        %Plot the per trial results for x with nn trained with permuted input
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


        hold on
        ii_start=0;
        last_op_predictedend=1;
        for trNo=1:trials.odor_trNo

            op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
            op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
            if op_predictedend>length(op_predicted)
                op_predictedend=length(op_predicted);
            end

            op_all_trials_sh=[op_all_trials_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_trials_sh=[op_decod_all_trials_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];

            op_between_trials_sh=[op_between_trials_sh; odor_plume_template(last_op_predictedend:op_predictedstart)'];
            op_decod_between_trials_sh=[op_decod_between_trials_sh; op_predicted_sh_conv(last_op_predictedend:op_predictedstart,1)];
            last_op_predictedend=op_predictedend;

            ii_end=ii_start+length(op_predicted(op_predictedstart:op_predictedend))-1;

            %Plot accuracy per trial for permuted training control
            CIsp = bootci(1000, @mean, op_predicted_sh_conv(op_predictedstart:op_predictedend,:)');
            meansp=mean(op_predicted_sh_conv(op_predictedstart:op_predictedend,:)',1);
            CIsp(1,:)=meansp-CIsp(1,:);
            CIsp(2,:)=CIsp(2,:)-meansp;

            [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(op_predicted_sh_conv(op_predictedstart:op_predictedend,:)',1)', CIsp', '-k');


            switch trials.odor_trial_type(trNo)
                case 1
                    %Lane 1 hits
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[213/255 94/255 0],'LineWidth',3)
                    %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_hits_sh=[op_all_hits_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_hits_sh=[op_decod_all_hits_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
                case 2
                    %Lane 1 miss
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[230/255 159/255 0],'LineWidth',3)
                    %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_miss_sh=[op_all_miss_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss_sh=[op_decod_all_miss_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
                case 3
                    %Lane 4 hit
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[0 114/255 178/255],'LineWidth',3)
                    %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_hits_sh=[op_all_hits_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_hits_sh=[op_decod_all_hits_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
                case 4
                    %Lane 4 hit
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[86/255 180/255 233/255],'LineWidth',3)
                    %             plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend,1),'-k','LineWidth',1)
                    plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
                    plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
                    op_all_miss_sh=[op_all_miss_sh; odor_plume_template(op_predictedstart:op_predictedend)'];
                    op_decod_all_miss_sh=[op_decod_all_miss_sh; op_predicted_sh_conv(op_predictedstart:op_predictedend,1)];
            end
            ii_start=ii_start+length(op_predicted(op_predictedstart:op_predictedend))+20;
        end

        title(['log10(odor) permuted, orange:hit1, blue:mis1, green:hit4, yellow:miss4 k:predicted, file No ' num2str(fileNo) ])



        %Plot odor conc vs decoded
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


        hold on

        plot(op_all_trials,op_decod_all_trials,'.b')
        xlabel('Odor concentration')
        ylabel('Predicted')

        title('odor concentration per trial (trained per trial)')

        [R1,p_val]=corrcoef(op_all_hits,op_decod_all_hits);
        fprintf(1, ['file No ' num2str(fileNo) ', Correlation coefficient nn conv odor conc hits %d\n'],R1(1,2));
        fprintf(1, ['file No ' num2str(fileNo) ', P value nn conv odor conc hits %d\n'],p_val(1,2));
        fprintf(1, ['file No ' num2str(fileNo) ', Number of neurons  %d\n'],no_neurons);
        fprintf(1, ['file No ' num2str(fileNo) ', Number of trials  %d\n\n'],trials.odor_trNo);

        pfft=1;

    end
end

%Calculate percent correct behavior
date_zero = datetime('01012020', 'InputFormat', 'MMddyyyy');
these_mice=unique(handles_conc.mouse);
for ii_mouse=these_mice
    pcorr.mouse(ii_mouse).pcorr=[];
    pcorr.mouse(ii_mouse).pcorr1=[];
    pcorr.mouse(ii_mouse).pcorr4=[];
    pcorr.mouse(ii_mouse).delta_date=[];
end
for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
        
        this_date=percent_correct.file(fileNo).date;
        ii_mouse=percent_correct.file(fileNo).mouse;
        pcorr.mouse(ii_mouse).delta_date=[pcorr.mouse(ii_mouse).delta_date days(this_date - date_zero)];

        this_pcorr=percent_correct.file(fileNo).percent_correct;
        pcorr.mouse(ii_mouse).pcorr=[pcorr.mouse(ii_mouse).pcorr this_pcorr];

         this_pcorr=percent_correct.file(fileNo).percent_correct1;
        pcorr.mouse(ii_mouse).pcorr1=[pcorr.mouse(ii_mouse).pcorr1 this_pcorr];

         this_pcorr=percent_correct.file(fileNo).percent_correct4;
        pcorr.mouse(ii_mouse).pcorr4=[pcorr.mouse(ii_mouse).pcorr4 this_pcorr];

        pffft=1;
    end
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

for ii_mouse=these_mice
    this_pc_array=zeros(length(pcorr.mouse(ii_mouse).pcorr),2);
    for ii_pcorr=1:length(pcorr.mouse(ii_mouse).pcorr)
        this_pc_array(ii_pcorr,2)=pcorr.mouse(ii_mouse).pcorr(ii_pcorr);
        this_pc_array(ii_pcorr,1)=pcorr.mouse(ii_mouse).delta_date(ii_pcorr);
    end
    sorted_this_pc_array=sortrows(this_pc_array);
    these_pcs=zeros(1,length(pcorr.mouse(ii_mouse).pcorr));
    these_pcs=sorted_this_pc_array(:,2);
    switch ii_mouse
        case 1
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[213/255 94/255 0])
        case 2
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[230/255 159/255 0])
        case 3
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[0 114/255 178/255])
        case 4
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[86/255 180/255 233/255])
    end
end

title('Percent correct')
xlabel('Session No')
ylabel('Percent correct')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

for ii_mouse=these_mice
    this_pc_array=zeros(length(pcorr.mouse(ii_mouse).pcorr),2);
    for ii_pcorr=1:length(pcorr.mouse(ii_mouse).pcorr)
        this_pc_array(ii_pcorr,2)=pcorr.mouse(ii_mouse).pcorr1(ii_pcorr);
        this_pc_array(ii_pcorr,1)=pcorr.mouse(ii_mouse).delta_date(ii_pcorr);
    end
    sorted_this_pc_array=sortrows(this_pc_array);
    these_pcs=zeros(1,length(pcorr.mouse(ii_mouse).pcorr));
    these_pcs=sorted_this_pc_array(:,2);
    switch ii_mouse
        case 1
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[213/255 94/255 0])
        case 2
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[230/255 159/255 0])
        case 3
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[0 114/255 178/255])
        case 4
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[86/255 180/255 233/255])
    end
end

title('Percent correct lane 1')
xlabel('Session No')
ylabel('Percent correct')


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

for ii_mouse=these_mice
    this_pc_array=zeros(length(pcorr.mouse(ii_mouse).pcorr),2);
    for ii_pcorr=1:length(pcorr.mouse(ii_mouse).pcorr)
        this_pc_array(ii_pcorr,2)=pcorr.mouse(ii_mouse).pcorr4(ii_pcorr);
        this_pc_array(ii_pcorr,1)=pcorr.mouse(ii_mouse).delta_date(ii_pcorr);
    end
    sorted_this_pc_array=sortrows(this_pc_array);
    these_pcs=zeros(1,length(pcorr.mouse(ii_mouse).pcorr));
    these_pcs=sorted_this_pc_array(:,2);
    switch ii_mouse
        case 1
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[213/255 94/255 0])
        case 2
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[230/255 159/255 0])
        case 3
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[0 114/255 178/255])
        case 4
            plot([1:length(pcorr.mouse(ii_mouse).pcorr)],these_pcs,'-o','Color',[86/255 180/255 233/255])
    end
end

title('Percent correct lane 4')
xlabel('Session No')
ylabel('Percent correct')

pfft=1;