function files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
files_included=ones(1,length(handles_conc.arena_file));

is_pearson=1;

%These are used to classify angle approach
low_angle=-130;
high_angle=-50;

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

    if no_hito/(no_hito+no_hit90)>thr_froa
        files_included(fileNo)=0;
    end
    % end
end

%Exclude files with p value for R1_all_trials > 0.05
thr_rho=0.05;

%Exclude files with ROI number >300
thr_no_roi=1300;

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

    if handles_out.no_neurons>thr_no_roi
        files_included(fileNo)=0;
    end

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
        if is_pearson==1
            [R1,this_P_rho_all_trials]=corrcoef(op_all_trials,op_decod_all_trials);
            R1_all_trials_pre(fileNo)=R1(1,2);
            P_rho_all_trials_pre(fileNo)=this_P_rho_all_trials(1,2);
        else
            [R1,this_P_rho_all_trials]=corr(op_all_trials,op_decod_all_trials,'Type','Spearman');
            R1_all_trials_pre(fileNo)=R1;
            P_rho_all_trials_pre(fileNo)=this_P_rho_all_trials;
        end
    else
        R1_all_trials_pre(fileNo)=NaN;
    end
    if P_rho_all_trials_pre(fileNo)>thr_rho
        files_included(fileNo)=0;
    end
end

end