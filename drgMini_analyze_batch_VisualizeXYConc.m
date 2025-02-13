%drgMini_analyze_batch_VisualizeXYConc
close all
clear all

is_sphgpu=0;

switch is_sphgpu
    case 0

        % %There was a bug and the shuffled runs were the same for all shuffled runs
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';
        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        %For the files below the shuffled runs should be different
        %Trained with all trials
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01062025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01062025.m';

         %Trained with hits only
         save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01122025/';
         choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01122025.m'

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');



    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end



%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4, 2 cm from floor
%Group 2 is rewarded, with odor lane 4, no odor in lane 1
%Group 3 is rewarded, with odor lane 1, no odor in lane 4
%Group 4 is rewarded, with no odor in lane 1 and lane 4
%Group 5 is rewarded, with ISO1 in both lane 1 and lane 4, 1 cm from floor

group_label{1}='Odor both lanes 2 cm';
group_label{2}='Odor lane 4';
group_label{3}='Odor lane 1';
group_label{4}='No odor';
group_label{5}='Odor both lanes 1 cm';


% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])
eval(['handles_Moser=' choiceMoserFileName(1:end-2) ';'])

figNo=0;

min_ii=10; %Entries with ii less than this nmber are shown as grey
d_xy=40;

%Load odor plume and calculate the odor plume for lane 1 and lane 4

%Use Aaron's simulation
%Plot odor plume
if is_sphgpu==1
    load('/data/SFTP/PreProcessedDR/Odor Arena Plumes/odorArenaPlumesDR.mat')
else
    load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumesDR2.mat')
end

%Now make mean_plume matrices for lanes 1 and 4
%Determine which lanes have odorants
%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4
%Group 2 is rewarded, with odor lane 4, no odor in lane 1
%Group 3 is rewarded, with odor lane 1, no odor in lane 4
%Group 4 is rewarded, with no odor in lane 1 and lane 4

%Note: This assumes taht the gorups chosen arre 1 and 5



%Groups 1 and 5 are rewarded, odor ISO1 in both lane 1 and lane 4
%They differ in cm from floor
lane1_odor_on=1;
lane4_odor_on=1;
odor_plume_patterns=[];

colormap(fire)
this_cmap=colormap;
this_cmap(1,:)=[0.7 0.7 0.7];

maxC_odor=-10000;
minC_odor=100000;

for cm_from_floor=1:2

    %Calculate mean for lane 4
    this_lane=4; %Note that I had mistakenly made lane 1 close to y=0, this is =1, not 4 even though I am calcualting lane 4
    for ii_source=1:length(odor_plumes.source)
        if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==cm_from_floor)
            this_source=ii_source;
        end
    end
    x_for_plume=10*odor_plumes.source(this_source).x;
    y_for_plume=10*(odor_plumes.source(this_source).y-min(odor_plumes.source(this_source).y));
    mean_plume_l4=odor_plumes.source(this_source).mean_plume;

    %Calculate mean lane 1
    this_lane=1; %Note that I had mistakenly made lane 4 close to y=480, this is =4, not 1 even though I am calcualting lane 4
    for ii_source=1:length(odor_plumes.source)
        if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==cm_from_floor)
            this_source=ii_source;
        end
    end
    mean_plume_l1=odor_plumes.source(this_source).mean_plume;

    mean_plume_l4=mean_plume_l4-min(mean_plume_l4(:));
    mean_plume_l1=mean_plume_l1-min(mean_plume_l1(:));

    min_nonzero=min([min(mean_plume_l4(mean_plume_l4~=0)) min(mean_plume_l1(mean_plume_l1~=0))]);

    if lane4_odor_on==1
        mean_plume_l4(mean_plume_l4==0)=min_nonzero;
    else
        mean_plume_l4(:,:)=min_nonzero;
    end

    %Hill transform?
    if handles_conc.hill==1
        mean_plume_l4=handles_conc.maxC*(mean_plume_l4.^handles_conc.n_hill)./...
            ((mean_plume_l4.^handles_conc.n_hill)+(handles_conc.k_half.^handles_conc.n_hill));
        mean_plume_l1=handles_conc.maxC*(mean_plume_l1.^handles_conc.n_hill)./...
            ((mean_plume_l1.^handles_conc.n_hill)+(handles_conc.k_half.^handles_conc.n_hill));
    end

    %Shift the plume to 7 cm (70 mm)
    if handles_conc.weber_fechner==0
        mean_plume_l4=handles_conc.multiplier*mean_plume_l4.^handles_conc.alpha;
    else
        %Weber-Frechner
        mean_plume_l4=handles_conc.multiplier*log10(mean_plume_l4);
    end

    if lane1_odor_on==1
        mean_plume_l1(mean_plume_l1==0)=min_nonzero;
    else
        mean_plume_l1(:,:)=min_nonzero;
    end

    if handles_conc.weber_fechner==0
        mean_plume_l1=handles_conc.multiplier*mean_plume_l1.^handles_conc.alpha;
    else
        %Weber-Frechner
        mean_plume_l1=handles_conc.multiplier*log10(mean_plume_l1);
    end

    % % minC=min([min(mean_plume_l4(:)) min(mean_plume_l1(:))]);
    % minC=prctile([mean_plume_l4(:); mean_plume_l1(:)],0.5);
    % if minC<handles_conc.lowest_conc
    %     minC=handles_conc.lowest_conc;
    % end
    % maxC=max([max(mean_plume_l4(:)) max(mean_plume_l1(:))]);
    % if minC==maxC
    %     maxC=minC+0.1;
    % end

    %Now shift to the actual dimensions of the odorant arena

    %There are slight differences between lane 1 and 4 in the low odorant areas
    %Here I will merge the data for the simulated lane 1 and lane 4 into a
    %lane14 and then I will use those data for the shifted lanes

    %The dimensions of the chamber are 50 cm for x and 48 cm for y
    %Cut the simulated plumes to 480 mm in the y dimension
    new_mean_plume_l1=mean_plume_l1(y_for_plume>=20,:);
    mean_plume_l1=[];
    mean_plume_l1=new_mean_plume_l1;
    new_mean_plume_l4=mean_plume_l4(y_for_plume<=480,:);
    mean_plume_l4=new_mean_plume_l4;

    new_y_for_plume=y_for_plume(:,y_for_plume<=480);
    y_for_plume=[];
    y_for_plume=new_y_for_plume;

    %Now make a merged lane 1 and 4

    %Reflect plume l1
    y_reflected_mean_plume_l1=zeros(size(mean_plume_l1,1),size(mean_plume_l1,2));
    for ii_y=1:size(mean_plume_l1,1)
        y_reflected_mean_plume_l1(ii_y,:)=mean_plume_l1(size(mean_plume_l1,1)-ii_y+1,:);
    end

    %Mean14 plume l4
    mean14_plume_l4=(y_reflected_mean_plume_l1+mean_plume_l4)/2;

    %Reflect plume l4
    mean14_plume_l1=zeros(size(mean_plume_l1,1),size(mean_plume_l1,2));
    for ii_y=1:size(mean_plume_l1,1)
        mean14_plume_l1(ii_y,:)=mean14_plume_l4(size(mean_plume_l1,1)-ii_y+1,:);
    end

    %Now shift the center of lane 4 to 70 mm
    new_mean14_plume_l4(y_for_plume>=20,:)=mean14_plume_l4(y_for_plume<=460,:);
    from_y_ii=find(y_for_plume<70,1,'first');
    new_mean14_plume_l4(length(y_for_plume):-1:from_y_ii,:)=mean14_plume_l4((y_for_plume>50)&(y_for_plume<=120),:);
    mean14_plume_l4=new_mean14_plume_l4;

    %Replace with mean
    mean_plume_l4=mean14_plume_l4;
    mean_plume_l1=mean14_plume_l1;

    %Note that I am replacing minC with prctile because there are a small
    %number of points with very low min
    if cm_from_floor==1
        minC=prctile([mean14_plume_l4(:); mean14_plume_l1(:)],0.5);
        maxC=min([max(mean14_plume_l4(:)) max(mean14_plume_l1(:))]);
        if maxC==minC
            maxC=minC+0.1;
        end
    end

    mean_plume_l4(mean_plume_l4<minC)=minC;
    mean_plume_l1(mean_plume_l1<minC)=minC;


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l4')
    colormap(this_cmap)
    shading interp
    caxis([minC-0.1 maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')
    yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})

    title(['Mean odor plume lane 4, ' num2str(cm_from_floor) ' cm from floor'])

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l1')
    colormap(this_cmap)
    shading interp
    caxis([minC-0.1 maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')
    yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})

    title(['Mean odor plume lane 1, ' num2str(cm_from_floor) ' cm from floor'])
 
    odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume=x_for_plume;
    odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume=y_for_plume;
    odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1=mean_plume_l1;
    odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4=mean_plume_l4;

    maxC_odor=max([maxC_odor max(mean_plume_l1(:)) max(mean_plume_l4(:))]);
    minC_odor=min([minC_odor min(mean_plume_l1(:)) min(mean_plume_l4(:))]);
end
  
if is_sphgpu==1
    save('/data/SFTP/PreProcessedDR/Odor Arena Plumes/odor_plume_patternsDR.mat','odor_plume_patterns','x_for_plume','y_for_plume','-v7.3')
else
    save('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odor_plume_patternsDR.mat','odor_plume_patterns','x_for_plume','y_for_plume','-v7.3')
end
  
%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

%Now plot pseudocolor locations for last turn angles, start and end points,
%etc

%First make point maps
ii_run=1;
these_groups=[1 5];
ii_for_corr=0;

%These are the size of the slides
y_length=480;
x_length=500;

 

%Note: We make this array d_xy+1 because pcolor does not show the last
%row/column

%Do this for three time windows: before trial start, during odor on, after
%odor off
window_legend{1}='before odor on';
window_legend{2}='during odor on';
window_legend{3}='after odor on';

%Determine actual and predicted odorant in these areas
odor_decoder=[];
odor_decoder.lane(1).areas(1).x=[150 350]; %Lane 1 corridor
odor_decoder.lane(1).areas(1).y=[400 480];
odor_decoder.lane(1).areas(2).x=[400 500]; %Lane 1 far corner
odor_decoder.lane(1).areas(2).y=[400 480];
odor_decoder.lane(1).areas(3).x=[150 350]; %Lane 1 opposite corridor
odor_decoder.lane(1).areas(3).y=[0 100];
odor_decoder.lane(1).areas(4).x=[400 500]; %Lane 1 opposite far corner
odor_decoder.lane(1).areas(4).y=[0 100];

odor_decoder.lane(4).areas(1).x=[150 350]; %Lane 4 corridor
odor_decoder.lane(4).areas(1).y=[0 100];
odor_decoder.lane(4).areas(2).x=[400 500]; %Lane 4 far corner
odor_decoder.lane(4).areas(2).y=[0 100];
odor_decoder.lane(4).areas(3).x=[150 350]; %Lane 4 opposite corridor
odor_decoder.lane(4).areas(3).y=[400 480];
odor_decoder.lane(4).areas(4).x=[400 500]; %Lane 4 opposite far corner
odor_decoder.lane(4).areas(4).y=[400 480];

for ii_time_window=1:3
    for ii_areas=1:4
        for ii_lane=[1 4]
            for ii_hit_miss=1:2
                for fileNo=1:length(handles_conc.arena_file)
                    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                        odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted=[];
                        odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual=[];
                    end
                end
            end
        end
    end
end

for ii_time_window=1:3
    y_values=(y_length/(2*d_xy)):y_length/d_xy:y_length-(y_length/(2*d_xy));
    x_values=(x_length/(2*d_xy)):x_length/d_xy:x_length-(x_length/(2*d_xy));

    lane1_hit_odor_predicted=zeros(d_xy,d_xy);
    lane1_miss_odor_predicted=zeros(d_xy,d_xy);

    lane4_hit_odor_predicted=zeros(d_xy,d_xy);
    lane4_miss_odor_predicted=zeros(d_xy,d_xy);

    lane1_hit_odor_actual=zeros(d_xy,d_xy);
    lane1_miss_odor_actual=zeros(d_xy,d_xy);

    lane4_hit_odor_actual=zeros(d_xy,d_xy);
    lane4_miss_odor_actual=zeros(d_xy,d_xy);

    lane1_hit_x_predicted=zeros(d_xy,d_xy);
    lane1_miss_x_predicted=zeros(d_xy,d_xy);

    lane4_hit_x_predicted=zeros(d_xy,d_xy);
    lane4_miss_x_predicted=zeros(d_xy,d_xy);

    lane1_hit_x_actual=zeros(d_xy,d_xy);
    lane1_miss_x_actual=zeros(d_xy,d_xy);

    lane4_hit_x_actual=zeros(d_xy,d_xy);
    lane4_miss_x_actual=zeros(d_xy,d_xy);

    lane1_hit_y_predicted=zeros(d_xy,d_xy);
    lane1_miss_y_predicted=zeros(d_xy,d_xy);

    lane4_hit_y_predicted=zeros(d_xy,d_xy);
    lane4_miss_y_predicted=zeros(d_xy,d_xy);

    lane1_hit_y_actual=zeros(d_xy,d_xy);
    lane1_miss_y_actual=zeros(d_xy,d_xy);

    lane4_hit_y_actual=zeros(d_xy,d_xy);
    lane4_miss_y_actual=zeros(d_xy,d_xy);





    lane1_hit_odor_predicted_ii=zeros(d_xy,d_xy);
    lane1_miss_odor_predicted_ii=zeros(d_xy,d_xy);

    lane4_hit_odor_predicted_ii=zeros(d_xy,d_xy);
    lane4_miss_odor_predicted_ii=zeros(d_xy,d_xy);

    lane1_hit_odor_actual_ii=zeros(d_xy,d_xy);
    lane1_miss_odor_actual_ii=zeros(d_xy,d_xy);

    lane4_hit_odor_actual_ii=zeros(d_xy,d_xy);
    lane4_miss_odor_actual_ii=zeros(d_xy,d_xy);

    lane1_hit_x_predicted_ii=zeros(d_xy,d_xy);
    lane1_miss_x_predicted_ii=zeros(d_xy,d_xy);

    lane4_hit_x_predicted_ii=zeros(d_xy,d_xy);
    lane4_miss_x_predicted_ii=zeros(d_xy,d_xy);

    lane1_hit_x_actual_ii=zeros(d_xy,d_xy);
    lane1_miss_x_actual_ii=zeros(d_xy,d_xy);

    lane4_hit_x_actual_ii=zeros(d_xy,d_xy);
    lane4_miss_x_actual_ii=zeros(d_xy,d_xy);

    lane1_hit_y_predicted_ii=zeros(d_xy,d_xy);
    lane1_miss_y_predicted_ii=zeros(d_xy,d_xy);

    lane4_hit_y_predicted_ii=zeros(d_xy,d_xy);
    lane4_miss_y_predicted_ii=zeros(d_xy,d_xy);

    lane1_hit_y_actual_ii=zeros(d_xy,d_xy);
    lane1_miss_y_actual_ii=zeros(d_xy,d_xy);

    lane4_hit_y_actual_ii=zeros(d_xy,d_xy);
    lane4_miss_y_actual_ii=zeros(d_xy,d_xy);






    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

            %Load angle file
            angle_file=handles_Angle.arena_file{fileNo};
            %load the ouptut file
            load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
            meanAngles=[];
            for trNo=1:length(handles_out.angles.trial)
                if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
                    meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
                else
                    meanAngles(trNo)=0; %I need to fix this
                end
            end

            handles_out_angle=handles_out;

            %Load conc data
            arena_file=handles_conc.arena_file{fileNo};
            %load the ouptut file
            load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            trials=handles_out.trials;
            odor_plume_template=handles_out.odor_plume_template;
            op_predicted=handles_out.op_predicted;

            min_odor_plume_actual=min(odor_plume_template);
            max_odor_plume_actual=max(odor_plume_template);


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

            for trNo=1:trials.odor_trNo



                %Find the last turn
                this_ii_turn=find(handles_out_angle.angles.trial(trNo).delta_x>100,1,'last');
                if ~isempty(this_ii_turn)

                    ii_turns=handles_out_angle.angles.trial(trNo).ii_turns(this_ii_turn);

                    %Choose time relevant to this time window
                    switch ii_time_window
                        case 1
                            ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                            ii_predictedend=trials.odor_ii_start(trNo)-1;
                        case 2
                            ii_predictedstart=trials.odor_ii_start(trNo);
                            ii_predictedend=trials.odor_ii_end(trNo);

                        case 3
                            ii_predictedstart=trials.odor_ii_end(trNo)+1;
                            ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                    end

                    if ii_predictedend>size(XYtest,1)
                        ii_predictedend=size(XYtest,1);
                    end
                    ii_end=trials.odor_ii_end(trNo);
                    if ii_end>size(XYtest,1)
                        ii_end=size(XYtest,1);
                    end

                    these_odor_predicted=op_predicted(ii_predictedstart:ii_predictedend);
                    these_odor_actual=odor_plume_template(ii_predictedstart:ii_predictedend);

                    these_x_predicted=x_predicted(ii_predictedstart:ii_predictedend);
                    these_x_actual=XYtest(ii_predictedstart:ii_predictedend,1);

                    these_y_predicted=y_predicted(ii_predictedstart:ii_predictedend);
                    these_y_actual=these_y_predicted;
                    these_y_actual=XYtest(ii_predictedstart:ii_predictedend,2);

                    for ii_t=1:length(these_odor_predicted)

                        %Find values for actual points
                        this_actual_x=these_x_actual(ii_t);
                        this_actual_y=these_y_actual(ii_t);
                        this_actual_odor=these_odor_actual(ii_t);


                        %Find values for predicted points
                        this_predicted_x=these_x_predicted(ii_t);
                        this_predicted_y=these_y_predicted(ii_t);
                        this_predicted_odor=these_odor_predicted(ii_t);

                        ii_x=ceil(this_actual_x/(x_length/d_xy));
                        if ii_x==0
                            ii_x=1;
                        end
                        if ii_x>d_xy
                            ii_x=d_xy;
                        end

                        ii_y=ceil(this_actual_y/(y_length/d_xy));
                        if ii_y==0
                            ii_y=1;
                        end
                        if ii_y>d_xy
                            ii_y=d_xy;
                        end


                        %Okabe_Ito colors
                        switch trials.odor_trial_type(trNo)
                            case 1
                                %Lane 1 hits vermillion
                                lane1_hit_odor_predicted(ii_x,ii_y)=lane1_hit_odor_predicted(ii_x,ii_y)+this_predicted_odor;
                                lane1_hit_odor_predicted_ii(ii_x,ii_y)=lane1_hit_odor_predicted_ii(ii_x,ii_y)+1;

                                lane1_hit_odor_actual(ii_x,ii_y)=lane1_hit_odor_actual(ii_x,ii_y)+this_actual_odor;
                                lane1_hit_odor_actual_ii(ii_x,ii_y)=lane1_hit_odor_actual_ii(ii_x,ii_y)+1;

                                lane1_hit_x_predicted(ii_x,ii_y)=lane1_hit_x_predicted(ii_x,ii_y)+this_predicted_x;
                                lane1_hit_x_predicted_ii(ii_x,ii_y)=lane1_hit_x_predicted_ii(ii_x,ii_y)+1;

                                lane1_hit_x_actual(ii_x,ii_y)=lane1_hit_x_actual(ii_x,ii_y)+this_actual_x;
                                lane1_hit_x_actual_ii(ii_x,ii_y)=lane1_hit_x_actual_ii(ii_x,ii_y)+1;

                                lane1_hit_y_predicted(ii_x,ii_y)=lane1_hit_y_predicted(ii_x,ii_y)+this_predicted_y;
                                lane1_hit_y_predicted_ii(ii_x,ii_y)=lane1_hit_y_predicted_ii(ii_x,ii_y)+1;

                                lane1_hit_y_actual(ii_x,ii_y)=lane1_hit_y_actual(ii_x,ii_y)+this_actual_y;
                                lane1_hit_y_actual_ii(ii_x,ii_y)=lane1_hit_y_actual_ii(ii_x,ii_y)+1;

                                %Keep track of odor in the quantificaiton areas
                                ii_lane=1;
                                ii_hit_miss=1;
                                for ii_areas=1:4
                                    if (this_actual_x>=odor_decoder.lane(ii_lane).areas(ii_areas).x(1))&(this_actual_x<odor_decoder.lane(ii_lane).areas(ii_areas).x(2))
                                        if (this_actual_y>=odor_decoder.lane(ii_lane).areas(ii_areas).y(1))&(this_actual_y<odor_decoder.lane(ii_lane).areas(ii_areas).y(2))
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted this_predicted_odor];
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual this_actual_odor];
                                        end
                                    end
                                end

                            case 2
                                %Lane 1 miss orange
                                lane1_miss_odor_predicted(ii_x,ii_y)=lane1_miss_odor_predicted(ii_x,ii_y)+this_predicted_odor;
                                lane1_miss_odor_predicted_ii(ii_x,ii_y)=lane1_miss_odor_predicted_ii(ii_x,ii_y)+1;

                                lane1_miss_odor_actual(ii_x,ii_y)=lane1_miss_odor_actual(ii_x,ii_y)+this_actual_odor;
                                lane1_miss_odor_actual_ii(ii_x,ii_y)=lane1_miss_odor_actual_ii(ii_x,ii_y)+1;

                                lane1_miss_x_predicted(ii_x,ii_y)=lane1_miss_x_predicted(ii_x,ii_y)+this_predicted_x;
                                lane1_miss_x_predicted_ii(ii_x,ii_y)=lane1_miss_x_predicted_ii(ii_x,ii_y)+1;

                                lane1_miss_x_actual(ii_x,ii_y)=lane1_miss_x_actual(ii_x,ii_y)+this_actual_x;
                                lane1_miss_x_actual_ii(ii_x,ii_y)=lane1_miss_x_actual_ii(ii_x,ii_y)+1;

                                lane1_miss_y_predicted(ii_x,ii_y)=lane1_miss_y_predicted(ii_x,ii_y)+this_predicted_y;
                                lane1_miss_y_predicted_ii(ii_x,ii_y)=lane1_miss_y_predicted_ii(ii_x,ii_y)+1;

                                lane1_miss_y_actual(ii_x,ii_y)=lane1_miss_y_actual(ii_x,ii_y)+this_actual_y;
                                lane1_miss_y_actual_ii(ii_x,ii_y)=lane1_miss_y_actual_ii(ii_x,ii_y)+1;

                                %Keep track of odor in the quantification areas
                                ii_lane=1;
                                ii_hit_miss=2;
                                for ii_areas=1:4
                                    if (this_actual_x>=odor_decoder.lane(ii_lane).areas(ii_areas).x(1))&(this_actual_x<odor_decoder.lane(ii_lane).areas(ii_areas).x(2))
                                        if (this_actual_y>=odor_decoder.lane(ii_lane).areas(ii_areas).y(1))&(this_actual_y<odor_decoder.lane(ii_lane).areas(ii_areas).y(2))
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted this_predicted_odor];
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual this_actual_odor];
                                        end
                                    end
                                end

                            case 3
                                %Lane 4 hit blue
                                lane4_hit_odor_predicted(ii_x,ii_y)=lane4_hit_odor_predicted(ii_x,ii_y)+this_predicted_odor;
                                lane4_hit_odor_predicted_ii(ii_x,ii_y)=lane4_hit_odor_predicted_ii(ii_x,ii_y)+1;

                                lane4_hit_odor_actual(ii_x,ii_y)=lane4_hit_odor_actual(ii_x,ii_y)+this_actual_odor;
                                lane4_hit_odor_actual_ii(ii_x,ii_y)=lane4_hit_odor_actual_ii(ii_x,ii_y)+1;

                                lane4_hit_x_predicted(ii_x,ii_y)=lane4_hit_x_predicted(ii_x,ii_y)+this_predicted_x;
                                lane4_hit_x_predicted_ii(ii_x,ii_y)=lane4_hit_x_predicted_ii(ii_x,ii_y)+1;

                                lane4_hit_x_actual(ii_x,ii_y)=lane4_hit_x_actual(ii_x,ii_y)+this_actual_x;
                                lane4_hit_x_actual_ii(ii_x,ii_y)=lane4_hit_x_actual_ii(ii_x,ii_y)+1;

                                lane4_hit_y_predicted(ii_x,ii_y)=lane4_hit_y_predicted(ii_x,ii_y)+this_predicted_y;
                                lane4_hit_y_predicted_ii(ii_x,ii_y)=lane4_hit_y_predicted_ii(ii_x,ii_y)+1;

                                lane4_hit_y_actual(ii_x,ii_y)=lane4_hit_y_actual(ii_x,ii_y)+this_actual_y;
                                lane4_hit_y_actual_ii(ii_x,ii_y)=lane4_hit_y_actual_ii(ii_x,ii_y)+1;


                                %Keep track of odor in the quantificaiton areas
                                ii_lane=4;
                                ii_hit_miss=1;
                                for ii_areas=1:4
                                    if (this_actual_x>=odor_decoder.lane(ii_lane).areas(ii_areas).x(1))&(this_actual_x<odor_decoder.lane(ii_lane).areas(ii_areas).x(2))
                                        if (this_actual_y>=odor_decoder.lane(ii_lane).areas(ii_areas).y(1))&(this_actual_y<odor_decoder.lane(ii_lane).areas(ii_areas).y(2))
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted this_predicted_odor];
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual this_actual_odor];
                                        end
                                    end
                                end
                                pfft=1;
                            case 4
                                %Lane 4 miss sky blue
                                lane4_miss_odor_predicted(ii_x,ii_y)=lane4_miss_odor_predicted(ii_x,ii_y)+this_predicted_odor;
                                lane4_miss_odor_predicted_ii(ii_x,ii_y)=lane4_miss_odor_predicted_ii(ii_x,ii_y)+1;

                                lane4_miss_odor_actual(ii_x,ii_y)=lane4_miss_odor_actual(ii_x,ii_y)+this_actual_odor;
                                lane4_miss_odor_actual_ii(ii_x,ii_y)=lane4_miss_odor_actual_ii(ii_x,ii_y)+1;

                                lane4_miss_x_predicted(ii_x,ii_y)=lane4_miss_x_predicted(ii_x,ii_y)+this_predicted_x;
                                lane4_miss_x_predicted_ii(ii_x,ii_y)=lane4_miss_x_predicted_ii(ii_x,ii_y)+1;

                                lane4_miss_x_actual(ii_x,ii_y)=lane4_miss_x_actual(ii_x,ii_y)+this_actual_x;
                                lane4_miss_x_actual_ii(ii_x,ii_y)=lane4_miss_x_actual_ii(ii_x,ii_y)+1;

                                lane4_miss_y_predicted(ii_x,ii_y)=lane4_miss_y_predicted(ii_x,ii_y)+this_predicted_y;
                                lane4_miss_y_predicted_ii(ii_x,ii_y)=lane4_miss_y_predicted_ii(ii_x,ii_y)+1;

                                lane4_miss_y_actual(ii_x,ii_y)=lane4_miss_y_actual(ii_x,ii_y)+this_actual_y;
                                lane4_miss_y_actual_ii(ii_x,ii_y)=lane4_miss_y_actual_ii(ii_x,ii_y)+1;

                                %Keep track of odor in the quantificaiton areas
                                ii_lane=4;
                                ii_hit_miss=2;
                                for ii_areas=1:4
                                    if (this_actual_x>=odor_decoder.lane(ii_lane).areas(ii_areas).x(1))&(this_actual_x<odor_decoder.lane(ii_lane).areas(ii_areas).x(2))
                                        if (this_actual_y>=odor_decoder.lane(ii_lane).areas(ii_areas).y(1))&(this_actual_y<odor_decoder.lane(ii_lane).areas(ii_areas).y(2))
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted this_predicted_odor];
                                            odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual=...
                                                [odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual this_actual_odor];
                                        end
                                    end
                                end
                        end
                    end
                end

            end

            pffft=1;
        end
    end


    % max_odor_pred=-1000;
    lane1_hit_odor_predicted=lane1_hit_odor_predicted./lane1_hit_odor_predicted_ii;
    % max_odor_pred=max([max_odor_pred max(lane1_hit_odor_predicted(~isnan(lane1_hit_odor_predicted)))]);
    % min_odor=min(lane1_hit_odor_predicted(~isnan(lane1_hit_odor_predicted)));
    % lane1_hit_odor_predicted(isnan(lane1_hit_odor_predicted))=min_odor-0.001;
    % lane1_hit_odor_predicted(lane1_hit_odor_predicted_ii<min_ii)=min_odor-0.001;
    % minC_lane1_hit_pred=min_odor;

    lane1_miss_odor_predicted=lane1_miss_odor_predicted./lane1_miss_odor_predicted_ii;
    % max_odor_pred=max([max_odor_pred max(lane1_miss_odor_predicted(~isnan(lane1_miss_odor_predicted)))]);
    % min_odor=min(lane1_miss_odor_predicted(~isnan(lane1_miss_odor_predicted)));
    % lane1_miss_odor_predicted(isnan(lane1_miss_odor_predicted))=min_odor-0.001;
    % lane1_miss_odor_predicted(lane1_miss_odor_predicted_ii<min_ii)=min_odor-0.001;
    % minC_lane1_miss_pred=min_odor;

    lane4_hit_odor_predicted=lane4_hit_odor_predicted./lane4_hit_odor_predicted_ii;
    % min_odor=min(lane4_hit_odor_predicted(~isnan(lane4_hit_odor_predicted)));
    % max_odor_pred=max([max_odor_pred max(lane4_hit_odor_predicted(~isnan(lane4_hit_odor_predicted)))]);
    % lane4_hit_odor_predicted(isnan(lane4_hit_odor_predicted))=min_odor-0.001;
    % lane4_hit_odor_predicted(lane4_hit_odor_predicted_ii<min_ii)=min_odor-0.001;
    % minC_lane4_hit_pred=min_odor;

    lane4_miss_odor_predicted=lane4_miss_odor_predicted./lane4_miss_odor_predicted_ii;
    % max_odor_pred=max([max_odor_pred max(lane4_miss_odor_predicted(~isnan(lane4_miss_odor_predicted)))]);
    % min_odor=min(lane4_miss_odor_predicted(~isnan(lane4_miss_odor_predicted)));
    % lane4_miss_odor_predicted(isnan(lane4_miss_odor_predicted))=min_odor-0.001;
    % lane4_miss_odor_predicted(lane4_miss_odor_predicted_ii<min_ii)=min_odor-0.001;
    % minC_lane4_miss_pred=min_odor;

    max_odor_act=-10000;
    lane1_hit_odor_actual=lane1_hit_odor_actual./lane1_hit_odor_actual_ii;
    % max_odor_act=max([max_odor_act max(lane1_hit_odor_predicted(~isnan(lane1_hit_odor_actual)))]);
    % min_odor=min(lane1_hit_odor_predicted(~isnan(lane1_hit_odor_actual)));
    % lane1_hit_odor_actual(isnan(lane1_hit_odor_actual))=min_odor-0.001;
    % lane1_hit_odor_actual(lane1_hit_odor_actual_ii<min_ii)=min_odor-0.001;
    % minC_lane1_hit_act=min_odor;

    lane1_miss_odor_actual=lane1_miss_odor_actual./lane1_miss_odor_actual_ii;
    % max_odor_act=max([max_odor_act max(lane1_miss_odor_predicted(~isnan(lane1_miss_odor_actual)))]);
    % min_odor=min(lane1_miss_odor_predicted(~isnan(lane1_miss_odor_actual)));
    % lane1_miss_odor_actual(isnan(lane1_miss_odor_actual))=min_odor-0.001;
    % lane1_miss_odor_actual(lane1_miss_odor_actual_ii<min_ii)=min_odor-0.001;
    % minC_lane1_miss_act=min_odor;

    lane4_hit_odor_actual=lane4_hit_odor_actual./lane4_hit_odor_actual_ii;
    % max_odor_act=max([max_odor_act max(lane4_hit_odor_predicted(~isnan(lane4_hit_odor_actual)))]);
    % min_odor=min(lane4_hit_odor_predicted(~isnan(lane4_hit_odor_actual)));
    % lane4_hit_odor_actual(isnan(lane4_hit_odor_actual))=min_odor-0.001;
    % lane4_hit_odor_actual(lane4_hit_odor_actual_ii<min_ii)=min_odor-0.001;
    % minC_lane4_hit_act=min_odor;

    lane4_miss_odor_actual=lane4_miss_odor_actual./lane4_miss_odor_actual_ii;
    % max_odor_act=max([max_odor_act max(lane4_miss_odor_predicted(~isnan(lane4_miss_odor_actual)))]);
    % min_odor=min(lane4_miss_odor_predicted(~isnan(lane4_miss_odor_actual)));
    % lane4_miss_odor_actual(isnan(lane4_miss_odor_actual))=min_odor-0.001;
    % lane4_miss_odor_actual(lane4_miss_odor_actual_ii<min_ii)=min_odor-0.001;
    % minC_lane4_miss_act=min_odor;


    %Plot predicted odor
    %Lane 1 hits
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane1_hit_odor_predicted)

    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

    yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})

    title(['Lane 1 hit predicted odor ' window_legend{ii_time_window}])
    

    %Lane 1 miss
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane1_miss_odor_predicted)
    
    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

      yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})

    title(['Lane 1 miss predicted odor ' window_legend{ii_time_window}])

    %Lane 4 hits
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane4_hit_odor_predicted)
    
    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

      yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})

    title(['Lane 4 hit predicted odor ' window_legend{ii_time_window}])

    %Lane 4 miss
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane4_miss_odor_predicted)
    
    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

      yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})


    title(['Lane 4 miss predicted odor ' window_legend{ii_time_window}])



    %Plot actual odor
    %Lane 1 hits
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane1_hit_odor_actual)
    
    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

      yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})
 
    title(['Lane 1 hit actual odor '  window_legend{ii_time_window}])

    %Lane 1 miss
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane1_miss_odor_actual)
    
    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

      yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})

    title(['Lane 1 miss actual odor '  window_legend{ii_time_window}])

    %Lane 4 hits
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane4_hit_odor_actual)
    
    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

      yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})

    title(['Lane 4 hit actual odor '  window_legend{ii_time_window}])

    %Lane 4 miss
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    drg_pcolor(repmat(x_values,d_xy,1)',repmat(y_values,d_xy,1),lane4_miss_odor_actual)
    
    colormap fire
    shading flat
    caxis([minC_odor maxC_odor]);
    set(gca, 'Color', [0.7, 0.7, 0.7]); % Background for NaNs
    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')

      yticks([70 100 200 300 400 430])
    yticklabels({'lane 4','100','200','300','400','lane 1'})


    title(['Lane 4 miss actual odor '  window_legend{ii_time_window}])

end
 
%Plot the bar graph for corridor predicted odor
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

window_labels{1}='before odor';
window_labels{2}='during odor';

lane_labels{1}='lane 1';
lane_labels{4}='lane 4';

%Plot the different R1s
all_R1s=[];

ii_areas=1;
for ii_lane=[1 4]
    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_time_window=1:2
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted];
                    these_odor_per_file.ii_type(ii_these_odor_per_file).file(fileNo).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted;
                    these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)]
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    %Orange
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

             these_CIs(ii_these_odor_per_file).CIout=CIout;

            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.window(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_time_window*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.window(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_time_window*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
            mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file(fileNo).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
        end
    end

     %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;


end





% x_pos=3;
text(2,-1.5,'Lane 1')
text(9,-1.5,'Lane 4')
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4 7 8 10 11])
xticklabels({'Hit bef','Miss bef','Hit','Miss','Hit bef','Miss bef','Hit','Miss'})


title(['Predicted odor for corridor area'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 12])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridor\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridor\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.lane',glm_r1.window',glm_r1.hit_miss',...
    'VariableNames',{'odor','lane','window','hit_miss'});
mdl = fitglm(tbl,'odor~lane+window+hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the corridor\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the corridor\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridor per file\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridorper file\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.lane',glm_per_file.window',glm_per_file.hit_miss',...
    'VariableNames',{'odor','lane','window','hit_miss'});
mdl = fitglm(tbl,'odor~lane+window+hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the corridor per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the corridorper file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);



%Plot the bar graph for far corner predicted odor
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

corner_labels{2}='same corner';
corner_labels{4}='opposite corner';

lane_labels{1}='lane 1';
lane_labels{4}='lane 4';

%Plot the different R1s
all_R1s=[];

ii_time_window=2;

for ii_lane=[1 4]
    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_areas=[2 4] %These are the lane and opposite far corners
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted];
                    these_odor_per_file.ii_type(ii_these_odor_per_file).file(fileNo).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted;
                    these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)]
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

            these_CIs(ii_these_odor_per_file).CIout=CIout;
           
            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.corner(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_areas*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.corner(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_areas*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' corner_labels{ii_time_window}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
            mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file(fileNo).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
        end
    end

    %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;


end

% x_pos=3;
text(2,-1.5,'Lane 1')
text(9,-1.5,'Lane 4')
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4 7 8 10 11])
xticklabels({'Hit','Miss','Hit opp','Miss opp','Hit','Miss','Hit opp','Miss opp'})


title(['Predicted odor for far corner'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 12])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.lane',glm_r1.corner',glm_r1.hit_miss',...
    'VariableNames',{'odor','lane','corner','hit_miss'});
mdl = fitglm(tbl,'odor~lane+corner+hit_miss+lane*corner*hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the far corner\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the far corner\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner per file\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner per file\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.lane',glm_per_file.corner',glm_per_file.hit_miss',...
    'VariableNames',{'odor','lane','corner','hit_miss'});
mdl = fitglm(tbl,'odor~lane+corner+hit_miss+lane*corner*hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the corridor per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the corridorper file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);


%Plot the bar graph for corridor actual odor
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

window_labels{1}='before odor';
window_labels{2}='during odor';

lane_labels{1}='lane 1';
lane_labels{4}='lane 4';

%Plot the different R1s
all_R1s=[];

ii_areas=1;
for ii_lane=[1 4]
    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_time_window=1:2
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual];
                    these_odor_per_file.ii_type(ii_these_odor_per_file).file(fileNo).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual;
                    these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)]
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

             these_CIs(ii_these_odor_per_file).CIout=CIout;

            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.window(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_time_window*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.window(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_time_window*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
            mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file(fileNo).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
        end
    end

     %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;


end





% x_pos=3;
text(2,-1.5,'Lane 1')
text(9,-1.5,'Lane 4')
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4 7 8 10 11])
xticklabels({'Hit bef','Miss bef','Hit','Miss','Hit bef','Miss bef','Hit','Miss'})


title(['Actual odor for corridor area'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 12])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridor\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridor\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.lane',glm_r1.window',glm_r1.hit_miss',...
    'VariableNames',{'odor','lane','window','hit_miss'});
mdl = fitglm(tbl,'odor~lane+window+hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the corridor\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the corridor\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridor per file\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridorper file\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.lane',glm_per_file.window',glm_per_file.hit_miss',...
    'VariableNames',{'odor','lane','window','hit_miss'});
mdl = fitglm(tbl,'odor~lane+window+hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the corridor per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the corridorper file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);



%Plot the bar graph for far corner actual odor
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

corner_labels{2}='same corner';
corner_labels{4}='opposite corner';

lane_labels{1}='lane 1';
lane_labels{4}='lane 4';

%Plot the different R1s
all_R1s=[];

ii_time_window=2;

for ii_lane=[1 4]
    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_areas=[2 4] %These are the lane and opposite far corners
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual];
                    these_odor_per_file.ii_type(ii_these_odor_per_file).file(fileNo).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual;
                    these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)]
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

            these_CIs(ii_these_odor_per_file).CIout=CIout;
           
            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.corner(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_areas*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.corner(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_areas*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' corner_labels{ii_time_window}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
            mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file(fileNo).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
        end
    end

    %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;


end

% x_pos=3;
text(2,-1.5,'Lane 1')
text(9,-1.5,'Lane 4')
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4 7 8 10 11])
xticklabels({'Hit','Miss','Hit opp','Miss opp','Hit','Miss','Hit opp','Miss opp'})


title(['Actual odor for far corner'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 12])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.lane',glm_r1.corner',glm_r1.hit_miss',...
    'VariableNames',{'odor','lane','corner','hit_miss'});
mdl = fitglm(tbl,'odor~lane+corner+hit_miss+lane*corner*hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the far corner\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the far corner\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner per file\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner per file\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.lane',glm_per_file.corner',glm_per_file.hit_miss',...
    'VariableNames',{'odor','lane','corner','hit_miss'});
mdl = fitglm(tbl,'odor~lane+corner+hit_miss+lane*corner*hit_miss'...
    ,'CategoricalVars',[2,3,4])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the corridor per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the corridorper file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);


%Now merge lanes 

%Plot the bar graph for corridor predicted odor merged lanes
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

window_labels{1}='before odor';
window_labels{2}='during odor';



%Plot the different R1s
all_R1s=[];

ii_areas=1;

these_odor_per_file=[];
these_CIs=[];
ii_these_odor_per_file=0;
for ii=1:4
    these_odor_per_file.ii_type(ii).ii_file_lane=0
end

for ii_time_window=1:2
    for ii_hit_miss=1:2
        ii_these_odor_per_file=ii_these_odor_per_file+1;
        these_odors=[];
        these_mean_odors_per_file=[];
        ii_file_lane=0;
        for fileNo=1:length(handles_conc.arena_file)
            if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                %load the ouptut file
                for ii_lane=[1 4]
                    ii_file_lane=ii_file_lane+1;
                    these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted];
                    these_odor_per_file.ii_type(ii_these_odor_per_file).file_lane(ii_file_lane).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted;
                    these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)];
                end
            end
        end

        %plot bar
        switch ii_hit_miss+2*(ii_time_window-1)
            case 1
                bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
            case 2
                bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
            case 3
                bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
            case 4
                bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
            case 5
                bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
        end

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_odors...
            ,edges,bar_offset,rand_offset,'k','k',1);

        these_CIs(ii_these_odor_per_file).CIout=CIout;

        bar_offset=bar_offset+1;

        glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
        glm_r1.window(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_time_window*ones(1,length(these_odors));
        glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
        glm_r1_ii=glm_r1_ii+length(these_odors);

        id_r1_ii=id_r1_ii+1;
        input_r1_data(id_r1_ii).data=these_odors;
        input_r1_data(id_r1_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

        glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
        glm_per_file.window(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_time_window*ones(1,length(these_mean_odors_per_file));
        glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
        glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

        id_per_file_ii=id_per_file_ii+1;
        input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
        input_per_file_data(id_per_file_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

    end
    bar_offset=bar_offset+1;
end

%Plot lines between file means
bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
for file_lane=1:ii_file_lane
        mean_odors_per_type=[];
        for ii_type=1:4
            this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file_lane(file_lane).odors);
            mean_odors_per_type=[mean_odors_per_type this_mean_odor];
            plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',5)
        end
        plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
    
end

%Plot CIs again
for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
end


% x_pos=3;
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4])
xticklabels({'Hit bef','Miss bef','Hit','Miss'})


title(['Predicted odor for corridor area'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 5])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridor merged lanes\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridor merged lanes\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.window',glm_r1.hit_miss',...
    'VariableNames',{'odor','window','hit_miss'});
mdl = fitglm(tbl,'odor~window+hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the corridor merged lanes\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the corridor merged lanes\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridor per file merged lanes\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the corridorper file merged lanes\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.window',glm_per_file.hit_miss',...
    'VariableNames',{'odor','window','hit_miss'});
mdl = fitglm(tbl,'odor~window+hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the corridor per file merged lanes\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the corridorper file merged lanes\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);

% 
% 
% %Plot the bar graph for far corner predicted odor merging lanes 1 and 4
% figNo = figNo + 1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% bar_offset=0;
% 
% edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
% rand_offset=0.8;
% 
% glm_r1=[];
% glm_r1_ii=0;
% 
% id_r1_ii=0;
% input_r1_data=[];
% 
% glm_per_file=[];
% glm_per_file_ii=0;
% 
% id_per_file_ii=0;
% input_per_file_data=[];
% 
% 
% 
% hit_miss_labels{1}='Hit';
% hit_miss_labels{2}='Miss';
% 
% corner_labels{2}='same corner';
% corner_labels{4}='opposite corner';
% 
% lane_labels{1}='lane 1';
% lane_labels{4}='lane 4';
% 
% %Plot the different R1s
% all_R1s=[];
% 
% ii_time_window=2;
% 
% 
%     these_odor_per_file=[];
%     these_CIs=[];
%     ii_these_odor_per_file=0;
%     for ii_areas=[2 4] %These are the lane and opposite far corners
%         for ii_hit_miss=1:2
%             ii_these_odor_per_file=ii_these_odor_per_file+1;
%             these_odors=[];
%             these_mean_odors_per_file=[];
%             ii_file_lane=0;
%             for fileNo=1:length(handles_conc.arena_file)
%                 if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
%                     %load the ouptut file
%                     for ii_lane=[1 4]
%                         ii_file_lane=ii_file_lane+1;
%                     these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted];
%                     these_odor_per_file.ii_type(ii_these_odor_per_file).file_lane(ii_file_lane).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted;
%                     these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)];
%                     end
%                 end
%             end
% 
%             %plot bar
%             switch ii_hit_miss+2*(ii_time_window-1)
%                 case 1
%                     bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
%                 case 2
%                     bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
%                 case 3
%                     bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
%                 case 4
%                     bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
%                 case 5
%                     bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
%             end
% 
%             %Violin plot
%             [mean_out, CIout]=drgViolinPoint(these_odors...
%                 ,edges,bar_offset,rand_offset,'k','k',1);
% 
%             these_CIs(ii_these_odor_per_file).CIout=CIout;
% 
%             bar_offset=bar_offset+1;
% 
%             glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
%             % glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
%             glm_r1.corner(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_areas*ones(1,length(these_odors));
%             glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
%             glm_r1_ii=glm_r1_ii+length(these_odors);
% 
%             id_r1_ii=id_r1_ii+1;
%             input_r1_data(id_r1_ii).data=these_odors;
%             input_r1_data(id_r1_ii).description=[hit_miss_labels{ii_hit_miss} ' ' corner_labels{ii_areas}];
% 
%             glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
%             % glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
%             glm_per_file.corner(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_areas*ones(1,length(these_mean_odors_per_file));
%             glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
%             glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);
% 
%             id_per_file_ii=id_per_file_ii+1;
%             input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
%             input_per_file_data(id_per_file_ii).description=[ hit_miss_labels{ii_hit_miss} ' ' corner_labels{ii_areas}];
% 
%         end
%         bar_offset=bar_offset+1;
%     end
% 
%     %Plot lines between file means
%     bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
%     for file_lane=1:ii_file_lane
%            mean_odors_per_type=[];
%             for ii_type=1:4
%                 this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file_lane(file_lane).odors);
%                 mean_odors_per_type=[mean_odors_per_type this_mean_odor];
%                 plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
%             end
%             plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
%     end
% 
%     %Plot CIs again
%     for ii=1:4
%     plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
%     end
% 
% 
%     bar_offset=bar_offset+1;
% 
% 
% 
% 
% % x_pos=3;
% % text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% % text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])
% 
% xticks([0 1 3 4])
% xticklabels({'Hit','Miss','Hit opposite','Miss opposite'})
% xtickangle(45);
% 
% title(['Predicted odor for far corner merging the two lanes'])
% ylabel('log10(odor)')
% ylim([-10 -1])
% xlim([-1 5])
% 
% %Perform the glm per data point
% fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner lanes merged\n'])
% fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner lanes merged\n']);
% %
% % fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% % fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);
% 
% tbl = table(glm_r1.data',glm_r1.corner',glm_r1.hit_miss',...
%     'VariableNames',{'odor','corner','hit_miss'});
% mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
%     ,'CategoricalVars',[2,3])
% 
% txt = evalc('mdl');
% txt=regexp(txt,'<strong>','split');
% txt=cell2mat(txt);
% txt=regexp(txt,'</strong>','split');
% txt=cell2mat(txt);
% 
% fprintf(fileID,'%s\n', txt);
% 
% 
% %Do the ranksum/t-test per data point
% fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the far corner lanes merged\n'])
% fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the far corner lanes merged\n']);
% 
% 
% [output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);
% 
% 
% %Perform the glm per file
% fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner per file lanes merged\n'])
% fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner per file lanes merged\n']);
% %
% % fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% % fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);
% 
% tbl = table(glm_per_file.data',glm_per_file.corner',glm_per_file.hit_miss',...
%     'VariableNames',{'odor','corner','hit_miss'});
% mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
%     ,'CategoricalVars',[2,3])
% 
% txt = evalc('mdl');
% txt=regexp(txt,'<strong>','split');
% txt=cell2mat(txt);
% txt=regexp(txt,'</strong>','split');
% txt=cell2mat(txt);
% 
% fprintf(fileID,'%s\n', txt);
% 
% 
% %Do the ranksum/t-test per file
% fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the far corner lanes merged per file\n'])
% fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the far corner lanes merged per file\n']);
% 
% 
% [output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);



%Plot the bar graph for far corner vs opposite corner predicted odor merging lanes 1 and 4
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

corner_labels{2}='same corner';
corner_labels{4}='opposite corner';

lane_labels{1}='lane 1';
lane_labels{4}='lane 4';

%Plot the different R1s
all_R1s=[];

ii_time_window=2;


    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_areas=[2 4] %These are the lane and opposite far corners
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            ii_file_lane=0;
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    for ii_lane=[1 4]
                        ii_file_lane=ii_file_lane+1;
                    these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted];
                    these_odor_per_file.ii_type(ii_these_odor_per_file).file_lane(ii_file_lane).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted;
                    these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)];
                    end
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

            these_CIs(ii_these_odor_per_file).CIout=CIout;
           
            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            % glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.corner(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_areas*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[hit_miss_labels{ii_hit_miss} ' ' corner_labels{ii_areas}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            % glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.corner(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_areas*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[ hit_miss_labels{ii_hit_miss} ' ' corner_labels{ii_areas}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for file_lane=1:ii_file_lane
           mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file_lane(file_lane).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
    end

    %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;




% x_pos=3;
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4])
xticklabels({'Hit','Miss','Hit opposite','Miss opposite'})
xtickangle(45);


title(['Predicted odor for far corner merging the two lanes'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 5])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner lanes merged\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner lanes merged\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.corner',glm_r1.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the far corner lanes merged\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the far corner lanes merged\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner per file lanes merged\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor in the far corner per file lanes merged\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.corner',glm_per_file.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the far corner lanes merged per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the far corner lanes merged per file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);




%Plot the bar graph comparing the corridor to far corner predicted odor merging lanes 1 and 4
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

area_labels{1}='corridor';
area_labels{2}='far corner';


%Plot the different R1s
all_R1s=[];

ii_time_window=2;


    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_areas=[1 2] %These are the lane and far corners
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            ii_file_lane=0;
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    for ii_lane=[1 4]
                        ii_file_lane=ii_file_lane+1;
                        these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted];
                        these_odor_per_file.ii_type(ii_these_odor_per_file).file_lane(ii_file_lane).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_predicted;
                        these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)];
                    end
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

            these_CIs(ii_these_odor_per_file).CIout=CIout;
           
            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            % glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.corner(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_areas*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[hit_miss_labels{ii_hit_miss} ' ' area_labels{ii_areas}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            % glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.corner(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_areas*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[ hit_miss_labels{ii_hit_miss} ' ' area_labels{ii_areas}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for file_lane=1:ii_file_lane
           mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file_lane(file_lane).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
    end

    %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;




% x_pos=3;
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4])
xticklabels({'Hit corridor','Miss corridor','Hit corner','Miss corner'})
xtickangle(45);

title(['Predicted odor for corridor vs. far corner merging the two lanes'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 5])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor for corridor vs far corner lanes merged\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor for corridor vs far corner lanes merged\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.corner',glm_r1.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor for corridor vs far corner lanes merged\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor for corridor vs far corner merged\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' predicted odor for corridor vs far corner per file lanes merged\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' predicted odor for corridor vs far corner per file lanes merged\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.corner',glm_per_file.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for predicted odor in the far corner lanes merged per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values predicted odor in the far corner lanes merged per file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);


%Plot the bar graph for corridor actual odor merged lanes
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

window_labels{1}='before odor';
window_labels{2}='during odor';

lane_labels{1}='lane 1';
lane_labels{4}='lane 4';

%Plot the different R1s
all_R1s=[];

ii_areas=1;

    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_time_window=1:2
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            ii_file_lane=0;
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    for ii_lane=[1 4]
                        ii_file_lane=ii_file_lane+1;
                        these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual];
                        these_odor_per_file.ii_type(ii_these_odor_per_file).file_lane(ii_file_lane).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual;
                        these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)];
                    end
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

             these_CIs(ii_these_odor_per_file).CIout=CIout;

            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            glm_r1.window(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_time_window*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            glm_per_file.window(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_time_window*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for file_lane=1:ii_file_lane
        mean_odors_per_type=[];
        for ii_type=1:4
            this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file_lane(file_lane).odors);
            mean_odors_per_type=[mean_odors_per_type this_mean_odor];
            plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
        end
        plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
    end

     %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;








% x_pos=3;

% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4])
xticklabels({'Hit before','Miss before','Hit odor','Miss odor'})
xtickangle(45);

title(['Actual odor for corridor area'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 5])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridor\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridor\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.window',glm_r1.hit_miss',...
    'VariableNames',{'odor','window','hit_miss'});
mdl = fitglm(tbl,'odor~window+hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the corridor\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the corridor\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridor per file\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the corridorper file\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.window',glm_per_file.hit_miss',...
    'VariableNames',{'odor','window','hit_miss'});
mdl = fitglm(tbl,'odor~window+hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the corridor per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the corridorper file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);



%Plot the bar graph for far corner vs opposite corner actual odor merged
%lanes
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

corner_labels{2}='same corner';
corner_labels{4}='opposite corner';

lane_labels{1}='lane 1';
lane_labels{4}='lane 4';

%Plot the different R1s
all_R1s=[];

ii_time_window=2;


    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_areas=[2 4] %These are the lane and opposite far corners
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            ii_file_lane=0;
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    for ii_lane=[1 4]
                        ii_file_lane=ii_file_lane+1;
                    these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual];
                    these_odor_per_file.ii_type(ii_these_odor_per_file).file_lane(ii_file_lane).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual;
                    these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)]
                    end
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

            these_CIs(ii_these_odor_per_file).CIout=CIout;
           
            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.corner(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_areas*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' window_labels{ii_time_window}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.corner(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_areas*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[lane_labels{ii_lane} ' ' hit_miss_labels{ii_hit_miss} ' ' corner_labels{ii_time_window}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for file_lane=1:ii_file_lane 
            mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file_lane(file_lane).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)   
    end

    %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;




% x_pos=3;
% text(2,-1.5,'Lane 1')
% text(9,-1.5,'Lane 4')
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4])
xticklabels({'Hit','Miss','Hit opposite','Miss opposite'})
xtickangle(45);

title(['Actual odor for far corner'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 5])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.corner',glm_r1.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the far corner\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the far corner\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner per file\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor in the far corner per file\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.corner',glm_per_file.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the corridor per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the corridorper file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);



%Plot the bar graph comparing the corridor to far corner actual odor merging lanes 1 and 4
figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[minC_odor:0.05*(maxC_odor-minC_odor):maxC_odor];
rand_offset=0.8;

glm_r1=[];
glm_r1_ii=0;

id_r1_ii=0;
input_r1_data=[];

glm_per_file=[];
glm_per_file_ii=0;

id_per_file_ii=0;
input_per_file_data=[];



hit_miss_labels{1}='Hit';
hit_miss_labels{2}='Miss';

area_labels{1}='corridor';
area_labels{2}='far corner';


%Plot the different R1s
all_R1s=[];

ii_time_window=2;


    these_odor_per_file=[];
    these_CIs=[];
    ii_these_odor_per_file=0;
    for ii_areas=[1 2] %These are the lane and far corners
        for ii_hit_miss=1:2
            ii_these_odor_per_file=ii_these_odor_per_file+1;
            these_odors=[];
            these_mean_odors_per_file=[];
            ii_file_lane=0;
            for fileNo=1:length(handles_conc.arena_file)
                if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)
                    %load the ouptut file
                    for ii_lane=[1 4]
                        ii_file_lane=ii_file_lane+1;
                        these_odors=[these_odors odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual];
                        these_odor_per_file.ii_type(ii_these_odor_per_file).file_lane(ii_file_lane).odors=odor_decoder.lane(ii_lane).areas(ii_areas).time_window(ii_time_window).file(fileNo).hit_miss(ii_hit_miss).odor_actual;
                        these_mean_odors_per_file=[these_mean_odors_per_file mean(these_odors)];
                    end
                end
            end

            %plot bar
            switch ii_hit_miss+2*(ii_time_window-1)
                case 1
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255],'BaseValue',-10)
                case 2
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255],'BaseValue',-10)
                case 3
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255],'BaseValue',-10)
                case 4
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255],'BaseValue',-10)
                case 5
                    bar(bar_offset,mean(these_odors),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 114/255 178/255],'BaseValue',-10)
            end

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odors...
                ,edges,bar_offset,rand_offset,'k','k',1);

            these_CIs(ii_these_odor_per_file).CIout=CIout;
           
            bar_offset=bar_offset+1;

            glm_r1.data(glm_r1_ii+1:glm_r1_ii+length(these_odors))=these_odors;
            % glm_r1.lane(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_lane*ones(1,length(these_odors));
            glm_r1.corner(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_areas*ones(1,length(these_odors));
            glm_r1.hit_miss(glm_r1_ii+1:glm_r1_ii+length(these_odors))=ii_hit_miss*ones(1,length(these_odors));
            glm_r1_ii=glm_r1_ii+length(these_odors);

            id_r1_ii=id_r1_ii+1;
            input_r1_data(id_r1_ii).data=these_odors;
            input_r1_data(id_r1_ii).description=[hit_miss_labels{ii_hit_miss} ' ' area_labels{ii_areas}];

            glm_per_file.data(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=these_mean_odors_per_file;
            % glm_per_file.lane(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_lane*ones(1,length(these_mean_odors_per_file));
            glm_per_file.corner(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_areas*ones(1,length(these_mean_odors_per_file));
            glm_per_file.hit_miss(glm_per_file_ii+1:glm_per_file_ii+length(these_mean_odors_per_file))=ii_hit_miss*ones(1,length(these_mean_odors_per_file));
            glm_per_file_ii=glm_per_file_ii+length(these_mean_odors_per_file);

            id_per_file_ii=id_per_file_ii+1;
            input_per_file_data(id_per_file_ii).data=these_mean_odors_per_file;
            input_per_file_data(id_per_file_ii).description=[ hit_miss_labels{ii_hit_miss} ' ' area_labels{ii_areas}];

        end
        bar_offset=bar_offset+1;
    end

    %Plot lines between file means
    bar_offsets=[bar_offset-6 bar_offset-5 bar_offset-3 bar_offset-2];
    for file_lane=1:ii_file_lane
           mean_odors_per_type=[];
            for ii_type=1:4
                this_mean_odor=mean(these_odor_per_file.ii_type(ii_type).file_lane(file_lane).odors);
                mean_odors_per_type=[mean_odors_per_type this_mean_odor];
                plot(bar_offsets(ii_type),this_mean_odor,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',10)
            end
            plot(bar_offsets,mean_odors_per_type,'-','Color',[0.7 0.7 0.7],'LineWidth',2)
    end

    %Plot CIs again
    for ii=1:4
    plot([bar_offsets(ii) bar_offsets(ii)],these_CIs(ii).CIout,'-k','LineWidth',3)
    end

 
    bar_offset=bar_offset+1;




% x_pos=3;
% text(x_pos,0.45,'miss','Color',[240/255 228/255 66/255])
% text(x_pos,0.41,'shuffled','Color',[0/255 114/255 178/255])

xticks([0 1 3 4])
xticklabels({'Hit corridor','Miss corridor','Hit corner','Miss corner'})
xtickangle(45);

title(['Actual odor for corridor vs. far corner merging the two lanes'])
ylabel('log10(odor)')
ylim([-10 -1])
xlim([-1 5])

%Perform the glm per data point
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor for corridor vs far corner lanes merged\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor for corridor vs far corner lanes merged\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_r1.data',glm_r1.corner',glm_r1.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per data point
fprintf(1, ['\n\nRanksum or t-test p values for actual odor for corridor vs far corner lanes merged\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor for corridor vs far corner merged\n']);


[output_data] = drgMutiRanksumorTtest(input_r1_data, fileID,0);


%Perform the glm per file
fprintf(1, ['\nglm for Fig ' num2str(figNo) ' actual odor for corridor vs far corner per file lanes merged\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(figNo) ' actual odor for corridor vs far corner per file lanes merged\n']);
%
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_per_file.data',glm_per_file.corner',glm_per_file.hit_miss',...
    'VariableNames',{'odor','corner','hit_miss'});
mdl = fitglm(tbl,'odor~corner+hit_miss+corner*hit_miss'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test per file
fprintf(1, ['\n\nRanksum or t-test p values for actual odor in the far corner lanes merged per file\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values actual odor in the far corner lanes merged per file\n']);


[output_data] = drgMutiRanksumorTtest(input_per_file_data, fileID,0);


pffft=1;
