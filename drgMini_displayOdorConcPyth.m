%drgMini_displayOdorConcPyth
close all
clear all

figNo=0;
dt=0.1;

% %File 2
% save_PathConcPyth='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeBayesOdorConc05072025/';
% odorConcPythFileName='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext_pyth2_bayes.mat';

%File 3
% save_PathConcPyth='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeBayesOdorConc05072025/';
% odorConcPythFileName='20220727_FCM19_withodor_miniscope_sync_L1andL4_ncorre_ext_nonneg_pyth2_bayes.mat';


%File 20
save_PathConcPyth='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeBayesOdorConc05072025/';
odorConcPythFileName='20221117_FCM22_withodor_nearfloor_miniscope_sync_L1andL4_ncorre_fix_ext_pyth2_bayes.mat';

ii_run=1;
% 
% addpath(choiceBatchPathName)
% eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
% eval(['handles_conc_pyth=' choiceOdorConcPythFileName(1:end-2) ';'])


% arena_file=handles_conc.arena_file{fileNo}; 
% load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
% trials=handles_out.trials;
% handles_out_conc=handles_out;
% 
% dFF_file=handles_conc.dFF_file{fileNo};
load([save_PathConcPyth odorConcPythFileName])
n_shuffle=size(handles_out.op_predicted_sh,1);
handles_out_conc_pyth=handles_out;

%Prune out the trials that were not processed
%Note: In the python code Byesian processing fails with a subset of trials
odor_plume_template=[];
op_predicted=[];
op_predicted_sh=[];
trial_type=[];
ii_predictedstart=[];
ii_predictedend=[];
this_trNo=0;
this_ii=0;
for trNo=1:length(handles_out.ii_predicted_start)
    if handles_out.trial_processed(handles_out.ii_predicted_start(trNo)+1)==1
        %This trial was processed
        this_trNo=this_trNo+1;
        odor_plume_template=[odor_plume_template handles_out.op_original(handles_out.ii_predicted_start(trNo)+1:handles_out.ii_predicted_end(trNo)+1)];
        op_predicted=[op_predicted handles_out.op_predicted(handles_out.ii_predicted_start(trNo)+1:handles_out.ii_predicted_end(trNo)+1)];
        op_predicted_sh=[op_predicted_sh handles_out.op_predicted_sh(:,handles_out.ii_predicted_start(trNo)+1:handles_out.ii_predicted_end(trNo)+1)];
        trial_type=[trial_type handles_out.trial_type(handles_out.ii_predicted_start(trNo)+1:handles_out.ii_predicted_end(trNo)+1)];
        this_ii=this_ii+1;
        ii_predictedstart(this_trNo)=this_ii;
        this_ii=this_ii+(handles_out.ii_predicted_end(trNo)-handles_out.ii_predicted_start(trNo));
        ii_predictedend(this_trNo)=this_ii;
    end
end

 
%Perform convolution on the predicted output 
no_conv_points=11;
% conv_win=ones(1,no_conv_points)/no_conv_points;
conv_win_gauss = gausswin(no_conv_points);
conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

op_predicted_conv=conv(op_predicted,conv_win_gauss,'same');

%Now limit the x and y to max and min
minop=min(odor_plume_template);
maxop=max(odor_plume_template);

op_predicted_conv(op_predicted_conv<minop)=minop;
op_predicted_conv(op_predicted_conv>maxop)=maxop;


op_predicted_sh_conv=zeros(size(op_predicted_sh,1),size(op_predicted_sh,2));

for ii_sh=1:n_shuffle
    this_op_predicted_sh=zeros(1,length(op_predicted));
    this_op_predicted_sh(1,:)=op_predicted_sh(ii_sh,:);

    this_op_predicted_sh_conv=[];

    this_op_predicted_sh_conv=conv(this_op_predicted_sh,conv_win_gauss,'same');


    %Now limit the x and y to max and min
    minop=min(odor_plume_template);
    this_op_predicted_sh_conv(this_op_predicted_sh_conv<minop)=minop;
    maxop=max(odor_plume_template);
    this_op_predicted_sh_conv(this_op_predicted_sh_conv>maxop)=maxop;

    op_predicted_sh_conv(ii_sh,:)=this_op_predicted_sh_conv;
end

pfft=1;


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

op_all_hits=[];
op_all_miss=[];
op_decod_all_hits=[];
op_decod_all_miss=[];

for trNo=1:length(ii_predictedstart)

    op_predictedstart=ii_predictedstart(trNo);
    op_predictedend=ii_predictedend(trNo);


    % op_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
    % op_predictedend=trials.odor_ii_end(trNo)+handles_choices.trial_end_offset;
    % if op_predictedend>size(X_dFF,1)
    %     op_predictedend=size(X_dFF,1);
    % end


    % op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
    % op_decod_all_trials=[op_decod_all_trials; op_predicted_conv(op_predictedstart:op_predictedend)];
    %
    % op_between_trials=[op_between_trials; odor_plume_template(last_op_predictedend:op_predictedstart)'];
    % op_decod_between_trials=[op_decod_between_trials; op_predicted_conv(last_op_predictedend:op_predictedstart)];
    last_op_predictedend=op_predictedend;

    ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;



    %Okabe_Ito colors
    switch trial_type(op_predictedstart)
        case 1
            %Lane 1 hits vermillion
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[213/255 94/255 0],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_hits=[op_all_hits odor_plume_template(op_predictedstart:op_predictedend)];
            op_decod_all_hits=[op_decod_all_hits op_predicted_conv(op_predictedstart:op_predictedend)];
        case 2
            %Lane 1 miss orange
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[230/255 159/255 0],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_miss=[op_all_miss odor_plume_template(op_predictedstart:op_predictedend)];
            op_decod_all_miss=[op_decod_all_miss op_predicted_conv(op_predictedstart:op_predictedend)];

        case 3
            %Lane 4 hit blue
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[0 114/255 178/255],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_hits=[op_all_hits odor_plume_template(op_predictedstart:op_predictedend)];
            op_decod_all_hits=[op_decod_all_hits op_predicted_conv(op_predictedstart:op_predictedend)];
        case 4
            %Lane 4 miss sky blue
            plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[86/255 180/255 233/255],'LineWidth',3)
            plot(dt*[ii_start:ii_end]',op_predicted_conv(op_predictedstart:op_predictedend)','-k','LineWidth',1.5)
            plot(dt*[ii_start+10 ii_start+10],[minop maxop],'-k')
            plot(dt*[ii_end-15 ii_end-15],[minop maxop],'-k')
            op_all_miss=[op_all_miss odor_plume_template(op_predictedstart:op_predictedend)];
            op_decod_all_miss=[op_decod_all_miss op_predicted_conv(op_predictedstart:op_predictedend)];
    end
    ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;

end

ylim([min(op_all_trials)-0.1*(max(op_all_trials)-min(op_all_trials)) max(op_all_trials)+0.1*(max(op_all_trials)-min(op_all_trials))])

title('Odor concentration per trial, verm:hit1, or:mis1, b:hit4, bsky:miss4 k:predicted')