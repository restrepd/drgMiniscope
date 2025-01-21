%drgMini_analyze_batch_PathAngle
close all
clear all

is_sphgpu=0;

switch is_sphgpu
    case 0

        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';

        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        % fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end

%These are used to classify angle approach
low_angle=-130;
high_angle=-50;

figNo=0;

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

run_label{1}='0 bins before';
run_label{5}='16 bins before';

% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])

figureNo=0;
%Exclude one file with high fraction_other_angle
%We will exclude fraction_other_angle>thr_froa
thr_froa=0.2;

%Calculate fraction of horizontal angle approach for each file
fraction_other_angle=[];
all_meanAngles=[];
for fileNo=1:length(handles_conc.arena_file)
    % if sum(handles_conc.group(fileNo)==these_groups)>0
        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        trials=handles_out.trials; 

        meanAngles=[];
        for trNo=1:length(handles_out.angles.trial)
            if ~isempty(handles_out.angles.trial(trNo).mean_end_angle)
                meanAngles(trNo)=handles_out.angles.trial(trNo).mean_end_angle;
            else
                meanAngles(trNo)=0; %I need to fix this
            end
            all_meanAngles=[all_meanAngles meanAngles(trNo)];
        end

        no_hit90=0;
        no_hito=0;
        for trNo=1:trials.odor_trNo


            %Okabe_Ito colors
            switch trials.odor_trial_type(trNo)
                case 1
                    %Hit 90 degrees
                    if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                        no_hit90=no_hit90+1;
                    end

                    %Horizontal hit
                    if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                        no_hito=no_hito+1;
                    end

                case 3
                  
                       %Hit 90 degrees
                    if (meanAngles(trNo)<=high_angle)&(meanAngles(trNo)>=low_angle)
                        no_hit90=no_hit90+1;
                    end

                    %Horizontal hit
                    if (meanAngles(trNo)>high_angle)|(meanAngles(trNo)<low_angle)
                        no_hito=no_hito+1;
                    end
            end
            
        end

        fraction_other_angle=[fraction_other_angle no_hito/(no_hito+no_hit90)];
        
    % end
end

%Exclude files with p value for R1_all_trials > 0.05
thr_rho=0.05;

%Now calculate p values for conc R1 all files
ii_run=1;
ii_for_corr=0;
R1_all_trials_pre=[];
P_rho_all_trials_pre=[];
these_fileNos_pre=[];


for fileNo=1:length(handles_conc.arena_file)
    
        op_all_trials=[];
        op_decod_all_trials=[];
        

        %Load conc data
        arena_file=handles_conc.arena_file{fileNo};
        %load the ouptut file
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        odor_plume_template=handles_out.odor_plume_template;
        op_predicted=handles_out.op_predicted;

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

        last_op_predictedend=1;
        

        for trNo=1:trials.odor_trNo

            op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
            op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
            if op_predictedend>length(odor_plume_template)
                op_predictedend=length(odor_plume_template);
            end

            op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
            op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];


        end

        if ~isempty(op_all_trials)
            [R1,this_P_rho_all_trials]=corrcoef(op_all_trials,op_decod_all_trials);
            R1_all_trials_pre(fileNo)=R1(1,2);
            P_rho_all_trials_pre(fileNo)=this_P_rho_all_trials(1,2);
        else
            R1_all_trials_pre(fileNo)=NaN;
        end
end

%Now plot pseudocolor locations for last turn angles, start and end points,
%etc

%First make point maps
ii_run=1;
these_groups=[1 5];
ii_for_corr=0;

y_length=480;
x_length=500;

y_values=(y_length/20):y_length/10:y_length-(y_length/20);
x_values=(x_length/20):x_length/10:x_length-(x_length/20);

lane1_hit_turn_angle_positions=zeros(10,10);
lane1_miss_turn_angle_positions=zeros(10,10);

lane4_hit_turn_angle_positions=zeros(10,10);
lane4_miss_turn_angle_positions=zeros(10,10);

lane1_hit_start_positions=zeros(10,10);
lane1_miss_start_positions=zeros(10,10);

lane4_hit_start_positions=zeros(10,10);
lane4_miss_start_positions=zeros(10,10);

lane1_hit_end_positions=zeros(10,10);
lane1_miss_end_positions=zeros(10,10);

lane4_hit_end_positions=zeros(10,10);
lane4_miss_end_positions=zeros(10,10);


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)
       
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


         %Load conc data
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

                ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if ii_predictedend>size(XYtest,1)
                    ii_predictedend=size(XYtest,1);
                end
                ii_end=trials.odor_ii_end(trNo)
                if ii_end>size(XYtest,1)
                    ii_end=size(XYtest,1);
                end


                %Find positions for start points
                this_start_x=XYtest(trials.odor_ii_start(trNo),1);
                this_start_y=XYtest(trials.odor_ii_start(trNo),2);

                ii_start_x=ceil(this_start_x/(x_length/10));
                if ii_start_x==0
                    ii_start_x=1;
                end
                if ii_start_x>10
                    ii_start_x=10;
                end
 
                ii_start_y=ceil(this_start_y/(y_length/10));
                if ii_start_y==0
                    ii_start_y=1;
                end
                if ii_start_y>10
                    ii_start_y=10;
                end

                %Find positions for end points
                this_end_x=XYtest(ii_end,1);
                this_end_y=XYtest(ii_end,2);

                ii_end_x=ceil(this_end_x/(x_length/10));
                if ii_end_x==0
                    ii_end_x=1;
                end
                if ii_end_x>10
                    ii_end_x=10;
                end
 
                ii_end_y=ceil(this_end_y/(y_length/10));
                if ii_end_y==0
                    ii_end_y=1;
                end
                if ii_end_y>10
                    ii_end_y=10;
                end

                %Find the position for this last turn
                this_ap_ii=ii_predictedstart+ii_turns-1;

                this_ap_x=XYtest(this_ap_ii,1);
                this_ap_y=XYtest(this_ap_ii,2);

                ii_ap_x=ceil(this_ap_x/(x_length/10));
                if ii_ap_x==0
                    ii_ap_x=1;
                end
                if ii_ap_x>10
                    ii_ap_x=10;
                end
 
                ii_ap_y=ceil(this_ap_y/(y_length/10));
                if ii_ap_y==0
                    ii_ap_y=1;
                end
                if ii_ap_y>10
                    ii_ap_y=10;
                end

                %Okabe_Ito colors
                switch trials.odor_trial_type(trNo)
                    case 1
                        %Lane 1 hits vermillion
                        lane1_hit_turn_angle_positions(ii_ap_x,ii_ap_y)=lane1_hit_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane1_hit_start_positions(ii_start_x,ii_start_y)=lane1_hit_start_positions(ii_start_x,ii_start_y)+1;
                        lane1_hit_end_positions(ii_end_x,ii_end_y)=lane1_hit_end_positions(ii_end_x,ii_end_y)+1;
                        pfft=1;
                    case 2
                        %Lane 1 miss orange
                        lane1_miss_turn_angle_positions(ii_ap_x,ii_ap_y)=lane1_miss_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane1_miss_start_positions(ii_start_x,ii_start_y)=lane1_miss_start_positions(ii_start_x,ii_start_y)+1;
                        lane1_miss_end_positions(ii_end_x,ii_end_y)=lane1_miss_end_positions(ii_end_x,ii_end_y)+1;
                    case 3
                        %Lane 4 hit blue
                        lane4_hit_turn_angle_positions(ii_ap_x,ii_ap_y)=lane4_hit_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane4_hit_start_positions(ii_start_x,ii_start_y)=lane4_hit_start_positions(ii_start_x,ii_start_y)+1;
                        lane4_hit_end_positions(ii_end_x,ii_end_y)=lane4_hit_end_positions(ii_end_x,ii_end_y)+1;
                        pfft=1;
                    case 4
                        %Lane 4 miss sky blue
                        lane4_miss_turn_angle_positions(ii_ap_x,ii_ap_y)=lane4_miss_turn_angle_positions(ii_ap_x,ii_ap_y)+1;
                        lane4_miss_start_positions(ii_start_x,ii_start_y)=lane4_miss_start_positions(ii_start_x,ii_start_y)+1;
                        lane4_miss_end_positions(ii_end_x,ii_end_y)=lane4_miss_end_positions(ii_end_x,ii_end_y)+1;
                end
            end

        end
        pfft=1;
    end
end

%Normalize
lane1_hit_turn_angle_positions=lane1_hit_turn_angle_positions/sum(lane1_hit_turn_angle_positions(:));
lane1_miss_turn_angle_positions=lane1_miss_turn_angle_positions/sum(lane1_miss_turn_angle_positions(:));
lane4_hit_turn_angle_positions=lane4_hit_turn_angle_positions/sum(lane4_hit_turn_angle_positions(:));
lane4_miss_turn_angle_positions=lane4_miss_turn_angle_positions/sum(lane4_miss_turn_angle_positions(:));

lane1_hit_start_positions=lane1_hit_start_positions/sum(lane1_hit_start_positions(:));
lane1_miss_start_positions=lane1_miss_start_positions/sum(lane1_miss_start_positions(:));
lane4_hit_start_positions=lane4_hit_start_positions/sum(lane4_hit_start_positions(:));
lane4_miss_start_positions=lane4_miss_start_positions/sum(lane4_miss_start_positions(:));

lane1_hit_end_positions=lane1_hit_end_positions/sum(lane1_hit_end_positions(:));
lane1_miss_end_positions=lane1_miss_end_positions/sum(lane1_miss_end_positions(:));
lane4_hit_end_positions=lane4_hit_end_positions/sum(lane4_hit_end_positions(:));
lane4_miss_end_positions=lane4_miss_end_positions/sum(lane4_miss_end_positions(:));

%Plot turn points
%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_turn_angle_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_hit_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit last turn positions')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_turn_angle_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_miss_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss last turn positions')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_turn_angle_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_hit_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit last turn positions')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_turn_angle_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_miss_turn_angle_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss last turn positions')

%Plot start points
%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_start_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_hit_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit start positions')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_start_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_miss_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss start positions')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_start_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_hit_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit start positions')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_start_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_miss_start_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss start positions')

%Plot end points
%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_end_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_hit_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit end positions')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_end_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_miss_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss end positions')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_end_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_hit_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit end positions')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_end_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_miss_end_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss end positions')

%Now plot pseudocolor map of trajectoreis before and after turn


%First make point maps
ii_run=1;
these_groups=[1 5];
ii_for_corr=0;

y_length=480;
x_length=500;

y_values=(y_length/20):y_length/10:y_length-(y_length/20);
x_values=(x_length/20):x_length/10:x_length-(x_length/20);

lane1_hit_before_positions=zeros(10,10);
lane1_miss_before_positions=zeros(10,10);

lane4_hit_before_positions=zeros(10,10);
lane4_miss_before_positions=zeros(10,10);

lane1_hit_after_positions=zeros(10,10);
lane1_miss_after_positions=zeros(10,10);

lane4_hit_after_positions=zeros(10,10);
lane4_miss_after_positions=zeros(10,10);


for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(fraction_other_angle(fileNo)<thr_froa)&(P_rho_all_trials_pre(fileNo)<=thr_rho)

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


        %Load conc data
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

                ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                ii_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
                if ii_predictedend>size(XYtest,1)
                    ii_predictedend=size(XYtest,1);
                end
                ii_end=trials.odor_ii_end(trNo)
                if ii_end>size(XYtest,1)
                    ii_end=size(XYtest,1);
                end



                %Find the position for this last turn
                this_ap_ii=ii_predictedstart+ii_turns-1;

                %Before
                for ii=trials.odor_ii_start(trNo):this_ap_ii

                    %Find position for this point
                    this_x=XYtest(ii,1);
                    this_y=XYtest(ii,2);

                    ii_x=ceil(this_x/(x_length/10));
                    if ii_x==0
                        ii_x=1;
                    end
                    if ii_x>10
                        ii_x=10;
                    end

                    ii_y=ceil(this_y/(y_length/10));
                    if ii_y==0
                        ii_y=1;
                    end
                    if ii_y>10
                        ii_y=10;
                    end
                    %Okabe_Ito colors
                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            lane1_hit_before_positions(ii_x,ii_y)=lane1_hit_before_positions(ii_x,ii_y)+1;

                        case 2
                            %Lane 1 miss orange
                            lane1_miss_before_positions(ii_x,ii_y)=lane1_miss_before_positions(ii_x,ii_y)+1;

                        case 3
                            %Lane 4 hit blue
                            lane4_hit_before_positions(ii_x,ii_y)=lane4_hit_before_positions(ii_x,ii_y)+1;

                        case 4
                            %Lane 4 miss sky blue
                            lane4_miss_before_positions(ii_x,ii_y)=lane4_miss_before_positions(ii_x,ii_y)+1;
                    end
                end


                %After
                for ii=this_ap_ii+1:ii_end

                    %Find position for this point
                    this_x=XYtest(ii,1);
                    this_y=XYtest(ii,2);

                    ii_x=ceil(this_x/(x_length/10));
                    if ii_x==0
                        ii_x=1;
                    end
                    if ii_x>10
                        ii_x=10;
                    end

                    ii_y=ceil(this_y/(y_length/10));
                    if ii_y==0
                        ii_y=1;
                    end
                    if ii_y>10
                        ii_y=10;
                    end
                    %Okabe_Ito colors
                    switch trials.odor_trial_type(trNo)
                        case 1
                            %Lane 1 hits vermillion
                            lane1_hit_after_positions(ii_x,ii_y)=lane1_hit_after_positions(ii_x,ii_y)+1;

                        case 2
                            %Lane 1 miss orange
                            lane1_miss_after_positions(ii_x,ii_y)=lane1_miss_after_positions(ii_x,ii_y)+1;

                        case 3
                            %Lane 4 hit blue
                            lane4_hit_after_positions(ii_x,ii_y)=lane4_hit_after_positions(ii_x,ii_y)+1;

                        case 4
                            %Lane 4 miss sky blue
                            lane4_miss_after_positions(ii_x,ii_y)=lane4_miss_after_positions(ii_x,ii_y)+1;
                    end
                end

            end
            pfft=1;
        end

    end
end

%Normalize
lane1_hit_before_positions=lane1_hit_before_positions/sum(lane1_hit_before_positions(:));
lane1_miss_before_positions=lane1_miss_before_positions/sum(lane1_miss_before_positions(:));
lane4_hit_before_positions=lane4_hit_before_positions/sum(lane4_hit_before_positions(:));
lane4_miss_before_positions=lane4_miss_before_positions/sum(lane4_miss_before_positions(:));

lane1_hit_after_positions=lane1_hit_after_positions/sum(lane1_hit_after_positions(:));
lane1_miss_after_positions=lane1_miss_after_positions/sum(lane1_miss_after_positions(:));
lane4_hit_after_positions=lane4_hit_after_positions/sum(lane4_hit_after_positions(:));
lane4_miss_after_positions=lane4_miss_after_positions/sum(lane4_miss_after_positions(:));


%Plot before points
%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_before_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_hit_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit trajectories before turn')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_before_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_miss_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss trajectories before turn')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_before_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_hit_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit trajectories before turn')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_before_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_miss_before_positions(:));
caxis([minC maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss trajectories before turn')

%Plot after points
%Lane 1 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_hit_after_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_hit_after_positions(:));
caxis([minC 0.35*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 hit trejectories after turn')

%Lane 1 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane1_miss_after_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane1_miss_after_positions(:));
caxis([minC 0.35*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 1 miss trejectories after turn')

%Lane 4 hits
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_hit_after_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_hit_after_positions(:));
caxis([minC 0.35*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})

title('Lane 4 hit trejectories after turn')

%Lane 4 miss
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


drg_pcolor(repmat(x_values,10,1)',repmat(y_values,10,1),lane4_miss_after_positions)
colormap fire
shading interp

set(gca, 'YDir', 'reverse');
% xlim(x_range)
% ylim(y_range)
% Ax = gca;
% Ax.Color = 'k';
xlabel('x (mm)')
ylabel('y (mm)')

minC=0;
maxC=max(lane4_miss_after_positions(:));
caxis([minC 0.35*maxC]);

yticks([50 100 200 300 400 430])
yticklabels({'lane 4','100','200','300','400','lane 1'})


title('Lane 4 miss trejectories after turn')

% fclose(fileID);

pffft=1;