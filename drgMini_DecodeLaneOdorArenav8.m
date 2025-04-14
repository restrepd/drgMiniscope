function handles_out=drgMini_DecodeLaneOdorArenav8(handles_choices2)
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020


if exist('handles_choices2')==0
    clear all
    close all

    handles_choices2.dispay_figures=1;
    %First troubleshooting file
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    % dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';
    % pred_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm_deczdFFopt2_711103.mat';

    %File 3 for troubleshooting
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220727_FCM19/';
    dFF_file='20220727_FCM19_withodor_miniscope_sync_L1andL4_ncorre_ext_nonneg.mat';
    arena_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm.mat';
    % pred_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm_deczdFFopt2_711103.mat';

    save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
    angle_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm_aangle.mat';

    save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
    xy_arena_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mmctreexy0.mat';

    %Second troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220713_FCM6/';
    % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat';
    % pred_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm_deczdFFopt2_711103.mat';

    %No odor troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';


    handles_choices2.this_path=this_path;
    handles_choices2.dFF_file=dFF_file;
    handles_choices2.arena_file=arena_file;
    handles_choices2.save_PathAngle=save_PathAngle;
    handles_choices2.angle_file=angle_file;
    handles_choices2.save_PathXY=save_PathXY;
    handles_choices2.xy_arena_file=xy_arena_file;

    is_gpu=0;
    handles_choices2.is_gpu=is_gpu;

    


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2;
    % dt_miniscope=1/30;
    n_shuffle=10; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_choices2.dt=dt;
    % handles_choices2.dt_miniscope=dt_miniscope;
    handles_choices2.n_shuffle=n_shuffle;

    which_training_range=1;
    handles_choices2.which_training_range=which_training_range;
    % 1=all trials
    % 2=hit trials
    % 3=miss trials

    which_ml_algo=1;
    handles_choices2.which_ml_algo=which_ml_algo;
    %1 SVZ
    %2 nn
    %3 tree
    %4 bayesian
    %5 glm
    %6 linear
    %7 NN parameter optimized

    ii_cost=3;
    handles_choices2.ii_cost=ii_cost;


    neurons_included=[1:46];
    % neurons_included=[10:20];
    handles_choices2.neurons_included=neurons_included;
   

else

    % which_ROIs=4; %when this function is called the user has to specify handles_choices2.process_these_ROIs
    this_path=handles_choices2.this_path;
    dFF_file=handles_choices2.dFF_file;
    arena_file=handles_choices2.arena_file;
    save_PathAngle=handles_choices2.save_PathAngle;
    angle_file=handles_choices2.angle_file;
    save_PathXY=handles_choices2.save_PathXY;
    xy_arena_file=handles_choices2.xy_arena_file;
    dt=handles_choices2.dt;
    n_shuffle=handles_choices2.n_shuffle;
    which_training_range=handles_choices2.which_training_range;
    which_ml_algo=handles_choices2.which_ml_algo;
    ii_cost=handles_choices2.ii_cost;
    neurons_included=handles_choices2.neurons_included;
    is_gpu=handles_choices2.is_gpu;


end

algo_name{1}='SVZ';
algo_name{2}='NN';
algo_name{3}='Tree';
algo_name{4}='Bayes';
algo_name{5}='glm';
algo_name{6}='LD';
algo_name{7}='NNopt';

this_cost=[0 ii_cost;ii_cost 0];




figNo=0;

%Restart random seeds
rng('shuffle');

%Now do decoding
% tic


X_dFF=[];
pos_binned_trimmed=[];

%Get XY and dFF per trial
load([save_PathXY xy_arena_file])
trials=handles_out.trials;
no_neurons=handles_out.no_neurons;

%Now generate shuffled lane identity using additional ii_sh because in some
%cases unique(Y) is only one choice
ii_shuffled=0;
these_perms=[];
for ii_sh=1:n_shuffle*3
    got_perm=0;
    while got_perm==0
        this_perm=randperm(length(trials.lane_per_trial));
        if ii_sh==1
            got_perm=1;
        else
            got_perm=1;
            for ii_p=1:size(these_perms,1)
                this_these_perms=zeros(1,length(trials.lane_per_trial));
                this_these_perms(1,:)=these_perms(ii_p,:);
                if all(this_these_perms==this_perm)
                    got_perm=0;
                end
            end
        end
    end
    these_perms=[these_perms; this_perm];
    for ii_tr=1:length(this_perm)
        trials.shuffled(ii_sh).lane_per_trial(ii_tr)=trials.lane_per_trial(this_perm(ii_tr));
    end
end

%load the angle ouptut file
load([save_PathAngle angle_file])
angles=handles_out.angles;

%Find turn points and time vectors
no_points=5;
all_ii_turns=zeros(1,trials.odor_trNo);
all_ii_ends=zeros(1,trials.odor_trNo);
include_trial=zeros(1,trials.odor_trNo);
for trNo=1:trials.odor_trNo
    this_ii_last_turn=find(angles.trial(trNo).delta_x>=100,1,'last');
    if ~isempty(this_ii_last_turn)
        all_ii_turns(trNo)=angles.trial(trNo).ii_turns(this_ii_last_turn);
        include_trial(trNo)=1;
    else
        [maxnum, maxii]=max(angles.trial(trNo).delta_x);
        all_ii_turns(trNo)=angles.trial(trNo).ii_turns(maxii);
    end
end




%trials.odor_trial_type
%1 Hit 1
%2 Miss 1
%3 Hit 4
%4 Miss 4

% Take a subset of neurons if the user asked to do this
include_these_neurons=zeros(1,no_neurons);
for ii=1:length(neurons_included)
    include_these_neurons(neurons_included(ii))=1;
end
for trNo=1:trials.odor_trNo
    %I double checked and this is the correct location of the turn
    this_XdFFtest=[];
    this_XdFFtest=trials.trial(trNo).XdFFtest;
    trials.trial(trNo).XdFFtest=zeros(size(this_XdFFtest,1),length(neurons_included));
    trials.trial(trNo).XdFFtest(:,:)=this_XdFFtest(:,logical(include_these_neurons));
end
no_neurons=length(neurons_included);

%Now parse out the dFF
X_dFF=[];
X_dFF_trial_start=[];
X_dFF_trial_end=[];
% X_hits=[];
lane_labels=[];
trimmed_time=[-3:dt:3];
training_range_template=[];
for ii_sh=1:n_shuffle*3
    shuffled_lanes(ii_sh).lane_labels=[];
end

ii=1;
actual_trNo=0;
for trNo=1:trials.odor_trNo
    if (all_ii_turns(trNo)-(3/dt)>0)&(all_ii_turns(trNo)+(3/dt)<=size(trials.trial(trNo).XdFFtest,1))
        actual_trNo=actual_trNo+1;
        %I double checked and this is the correct location of the turn
        these_trimmed_X_dFF=zeros(length(trimmed_time),no_neurons);
        these_trimmed_X_dFF(:,:)=trials.trial(trNo).XdFFtest(all_ii_turns(trNo)-(3/dt):all_ii_turns(trNo)+(3/dt),:);

        if (trials.lane_per_trial(trNo)==1)
            these_trimmed_lanes=ones(1,length(trimmed_time));
        else
            these_trimmed_lanes=zeros(1,length(trimmed_time));
        end
        lane_labels=[lane_labels these_trimmed_lanes];

        for ii_sh=1:n_shuffle*3
            if (trials.shuffled(ii_sh).lane_per_trial(trNo)==1)
                these_trimmed_lanes=ones(1,length(trimmed_time));
            else
                these_trimmed_lanes=zeros(1,length(trimmed_time));
            end
            shuffled_lanes(ii_sh).lane_labels=[shuffled_lanes(ii_sh).lane_labels these_trimmed_lanes];
        end

        switch which_training_range
            case 1
                this_training_range_template=ones(1,length(trimmed_time));
            case 2
                %Miss trials
                if (trials.odor_trial_type(trNo)==1)||(trials.odor_trial_type(trNo)==3)
                    this_training_range_template=ones(1,length(trimmed_time));
                    % this_training_range_template(trimmed_time>0)=0;
                else
                    this_training_range_template=zeros(1,length(trimmed_time));
                end
            case 3
                %Miss trials
                if (trials.odor_trial_type(trNo)==1)||(trials.odor_trial_type(trNo)==3)
                    this_training_range_template=zeros(1,length(trimmed_time));
                    % this_training_range_template(trimmed_time>0)=0;
                else
                    this_training_range_template=ones(1,length(trimmed_time));
                end
        end

        training_range_template=[training_range_template this_training_range_template];

        X_dFF=[X_dFF; these_trimmed_X_dFF];
       
        X_dFF_trial_start(actual_trNo)=ii;
        X_dFF_trial_end(actual_trNo)=ii+length(these_trimmed_lanes)-1;
        ii=ii+length(these_trimmed_lanes);
    end
end

% training_range_template=ones(1,length(X_hits));
no_time_bins=size(X_dFF,1);

for trNo=1:actual_trNo
    % parfor trNo=1:trials.odor_trNo

    this_test_range=zeros(1,no_time_bins);


    this_test_range(X_dFF_trial_start(trNo):X_dFF_trial_end(trNo))=1;
    this_training_range=logical(training_range_template)&(~logical(this_test_range));

    XdFFtrain=X_dFF(this_training_range,:);
    XdFFtest=X_dFF(logical(this_test_range),:);


    %Decode using neural network
    tblTrn=[];
    tblTrn = array2table(XdFFtrain);
    Y=lane_labels(this_training_range);

    switch which_ml_algo
        case 1
            if is_gpu==0
                Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
            else
                %Use gpu
                tblTrn_gpu=gpuArray(XdFFtrain); %Note: gpuArray does not support tables
                Y_gpu=gpuArray(Y);
                this_cost_gpu=gpuArray(this_cost);
                Mdl_gpu = fitcsvm(tblTrn_gpu,Y_gpu,'Cost',this_cost_gpu);

                %Gather back Mdl from the gpu
                Mdl=gather(Mdl_gpu);
            end
        case 2
            if is_gpu==0
                Mdl = fitcnet(tblTrn,Y);
            else
                %Use gpu
                tblTrn_gpu=gpuArray(XdFFtrain);
                Y_gpu=gpuArray(Y);
                Mdl_gpu = fitcnet(tblTrn_gpu,Y_gpu);

                %Gather back Mdl from the gpu
                Mdl=gather(Mdl_gpu);
            end
        case 3
            if is_gpu==0
                Mdl = fitctree(tblTrn,Y);
            else

                %Use gpu
                tblTrn_gpu=gpuArray(XdFFtrain);
                Y_gpu=gpuArray(Y);
                Mdl_gpu = fitctree(tblTrn_gpu,Y_gpu);

                %Gather back Mdl from the gpu
                Mdl=gather(Mdl_gpu);

            end
        case 4
            Mdl=fitcnb(tblTrn,Y,'Cost',this_cost);
        case 5
            Mdl = fitglm(XdFFtrain,Y','Distribution','binomial');
        case 6
            Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
        case 7
            opts = struct('ShowPlots', false, ...
                'Verbose', 0, ...
                'MaxObjectiveEvaluations', 15,...
                'UseParallel',true);
            Mdl = fitcnet(tblTrn,Y,'OptimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions', opts);
    end

    label_pred(trNo).Mdl=Mdl;


    label_pred(trNo).data=predict(Mdl,XdFFtest);
    %             y_pred(trNo).data=predict(MdlY2,XdFFtest);
    pfffft=1;

end
% fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

%Parse out the parfor loop output
label_predicted=zeros(no_time_bins,1);
accuracy=zeros(no_time_bins,1);
%         y_predicted=zeros(no_time_bins,1);
for trNo=1:actual_trNo
    this_test_range=zeros(1,no_time_bins);


    this_test_range(X_dFF_trial_start(trNo):X_dFF_trial_end(trNo))=1;
    % this_training_range=logical(training_range_template)&(~logical(this_test_range));

    label_predicted(logical(this_test_range),1)=label_pred(trNo).data;
    %             y_predicted(logical(this_test_range),1)=y_pred(trNo).data;

    %Compute accuracy
    this_accuracy=zeros(size(label_pred(trNo).data,1),size(label_pred(trNo).data,2));
    if trials.lane_per_trial(trNo)==1
        these_times=logical(label_pred(trNo).data==1);
        this_accuracy(these_times)=1;
        this_accuracy(~these_times)=0;
    else
        these_times=logical(label_pred(trNo).data==1);
        this_accuracy(these_times)=0;
        this_accuracy(~these_times)=1;
    end
    accuracy(logical(this_test_range),1)=this_accuracy;
end

%Now do shuffled decoding
label_predicted_sh=zeros(no_time_bins,n_shuffle);
accuracy_sh=zeros(no_time_bins,n_shuffle);
ii_sh=0;
for ii_shuffled=1:n_shuffle
    ii_sh=ii_sh+1;
    try
        for trNo=1:actual_trNo
            y_pred(trNo).data=[];
            label_pred(trNo).data=[];
        end

      

        parfor trNo=1:trials.odor_trNo
        % for trNo=1:actual_trNo

            this_test_range=zeros(1,no_time_bins);


            this_test_range(X_dFF_trial_start(trNo):X_dFF_trial_end(trNo))=1;
            this_training_range=logical(training_range_template)&(~logical(this_test_range));

            XdFFtrain=X_dFF(this_training_range,:);
            % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
            XdFFtest=X_dFF(logical(this_test_range),:);



            %Decode using neural network
            tblTrn=[];
            tblTrn = array2table(XdFFtrain);

            Y=shuffled_lanes(ii_shuffled).lane_labels(this_training_range);

            switch which_ml_algo
                case 1
                    if is_gpu==0
                        Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                    else
                        %Use gpu
                        tblTrn_gpu=gpuArray(XdFFtrain);
                        Y_gpu=gpuArray(Y);
                        this_cost_gpu=gpuArray(this_cost);
                        Mdl_gpu = fitcsvm(tblTrn_gpu,Y_gpu,'Cost',this_cost_gpu);

                        %Gather back Mdl from the gpu
                        Mdl=gather(Mdl_gpu);
                    end
                case 2
                    if is_gpu==0
                        Mdl = fitcnet(tblTrn,Y);
                    else
                        %Use gpu
                        tblTrn_gpu=gpuArray(XdFFtrain);
                        Y_gpu=gpuArray(Y);
                        Mdl_gpu = fitcnet(tblTrn_gpu,Y_gpu);

                        %Gather back Mdl from the gpu
                        Mdl=gather(Mdl_gpu);
                    end
                case 3
                    if is_gpu==0
                        Mdl = fitctree(tblTrn,Y);
                    else
                        %Use gpu
                        tblTrn_gpu=gpuArray(XdFFtrain);
                        Y_gpu=gpuArray(Y);
                        Mdl_gpu = fitctree(tblTrn_gpu,Y_gpu);

                        %Gather back Mdl from the gpu
                        Mdl=gather(Mdl_gpu);
                    end
                case 4
                    Mdl=fitcnb(tblTrn,Y,'Cost',this_cost);
                case 5
                    Mdl = fitglm(XdFFtrain,Y','Distribution','binomial');
                case 6
                    Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                case 7
                    %Note that we do not optimize the shuffled runs
                    Mdl = fitcnet(tblTrn,Y);
            end


            label_pred(trNo).data=predict(Mdl,XdFFtest);

        end

        %Extract the output
        for trNo=1:actual_trNo
            this_test_range=zeros(1,no_time_bins);


            this_test_range(X_dFF_trial_start(trNo):X_dFF_trial_end(trNo))=1;
            % this_training_range=logical(training_range_template)&(~logical(this_test_range));

            label_predicted_sh(logical(this_test_range),ii_shuffled)=label_pred(trNo).data;
            %Compute accuracy
            this_accuracy=zeros(size(label_pred(trNo).data,1),size(label_pred(trNo).data,2));
            if trials.lane_per_trial(trNo)==1
                these_times=logical(label_pred(trNo).data==1);
                this_accuracy(these_times)=1;
                this_accuracy(~these_times)=0;
            else
                these_times=logical(label_pred(trNo).data==1);
                this_accuracy(these_times)=0;
                this_accuracy(~these_times)=1;
            end
            accuracy_sh(logical(this_test_range),ii_shuffled)=this_accuracy;
        end
    catch
        ii_sh=ii_sh+1;
    end
end


% fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])



handles_out.label_predicted_sh=label_predicted_sh;
handles_out.accuracy_sh=accuracy_sh;
handles_out.label_predicted=label_predicted;
handles_out.accuracy=accuracy;
handles_out.trimmed_time=trimmed_time;
handles_out.trials=trials;



no_conv_points=5;
conv_win_gauss = gausswin(no_conv_points);
conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

label_predicted_conv=conv(label_predicted,conv_win_gauss,'same');
accuracy_conv=conv(accuracy,conv_win_gauss,'same');


label_predicted_sh_conv=zeros(size(label_predicted_sh,1),size(label_predicted_sh,2));
accuracy_sh_conv=zeros(size(label_predicted_sh,1),size(label_predicted_sh,2));
for ii_sh=1:n_shuffle
    this_label_predicted_sh=zeros(size(label_predicted_sh,1),1);
    this_label_predicted_sh(:,1)=label_predicted_sh(:,ii_sh);
    this_label_predicted_sh_conv=[];
    this_label_predicted_sh_conv=conv(this_label_predicted_sh,conv_win_gauss,'same');
    label_predicted_sh_conv(:,ii_sh)=this_label_predicted_sh_conv;

    this_accuracy_sh=zeros(size(accuracy_sh,1),1);
    this_accuracy_sh(:,1)=accuracy_sh(:,ii_sh);
    this_accuracy_sh_conv=[];
    this_accuracy_sh_conv=conv(this_accuracy_sh,conv_win_gauss,'same');
    accuracy_sh_conv(:,ii_sh)=this_accuracy_sh_conv;
end

for trNo=1:actual_trNo
    label_predictedstart(trNo)=X_dFF_trial_start(trNo);
    label_predictedend(trNo)=X_dFF_trial_end(trNo);
end

%Show accuracy aligned to the start of the trial
align_display=0;
accuracy_all_trials=[];
accuracy_conv_all_trials=[];
accuracy_sh_conv_all_trials=[];
accuracy_sh_all_trials=[];
trimmed_trial_type=[];
% delta_odor_start_to_end=ceil(mean(trials.odor_ii_end-trials.odor_ii_start));
for trNo=1:actual_trNo
    this_display_start=label_predictedstart(trNo);
    this_display_end=label_predictedend(trNo);

    if this_display_start>=1
        trimmed_trial_type=[trimmed_trial_type trials.odor_trial_type(trNo)];
        accuracy_all_trials=[accuracy_all_trials accuracy(this_display_start:this_display_end)];
        accuracy_conv_all_trials=[accuracy_conv_all_trials accuracy_conv(this_display_start:this_display_end)];
        for ii_sh=1:n_shuffle
            accuracy_sh_all_trials=[accuracy_sh_all_trials accuracy_sh(this_display_start:this_display_end,ii_sh)];
            accuracy_sh_conv_all_trials=[accuracy_sh_conv_all_trials accuracy_sh_conv(this_display_start:this_display_end,ii_sh)];
        end
    end

end




%Plot accuracy
if handles_choices2.display_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    hold on

    CIpvsm = bootci(1000, @mean, accuracy_sh_conv_all_trials');
    meanpvsm=mean(accuracy_sh_conv_all_trials',1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(trimmed_time,mean(accuracy_sh_conv_all_trials'), CIpvsm','r');

    CIpvsm = bootci(1000, @mean, accuracy_conv_all_trials');
    meanpvsm=mean(accuracy_conv_all_trials',1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(trimmed_time,mean(accuracy_conv_all_trials'), CIpvsm','k');

    plot(trimmed_time,mean(accuracy_sh_conv_all_trials'), '-r','LineWidth',2);
    plot(trimmed_time,mean(accuracy_conv_all_trials'), '-k','LineWidth',2);

    ylim([0 1.1])

    plot([0 0],[0 1.1],'-k')
    x_pos=-2;

    text(x_pos,0.95,'All','Color',[0/255 0/255 0/255])
    % text(x_pos,0.95,'Misses','Color',[240/255 228/255 66/255])
    text(x_pos,0.9,'Shuffled','Color',[255/255 0/255 0/255])

    xlabel('Time (sec)')
    ylabel('Accuracy')
    switch which_training_range
        case 1
            title(['Accuracy, trained with all trials ' algo_name{which_ml_algo}])
        case 2
            title(['Accuracy, trained with hits ' algo_name{which_ml_algo}])
        case 3
            title(['Accuracy, trained with misses ' algo_name{which_ml_algo}])
    end


    %Plot accuracy for hits
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    hold on

    CIpvsm = bootci(1000, @mean, accuracy_sh_conv_all_trials');
    meanpvsm=mean(accuracy_sh_conv_all_trials',1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(trimmed_time,mean(accuracy_sh_conv_all_trials'), CIpvsm','cmap',[255/255 0/255 0/255]);

    no_hits=sum((trimmed_trial_type==1)|(trimmed_trial_type==3));
    accuracy_conv_all_hits=zeros(no_hits,size(accuracy_conv_all_trials,2));
    accuracy_conv_all_hits=accuracy_conv_all_trials(:,(trimmed_trial_type==1)|(trimmed_trial_type==3));
    CIpvsm = bootci(1000, @mean, accuracy_conv_all_hits');
    meanpvsm=mean(accuracy_conv_all_hits',1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(trimmed_time,mean(accuracy_conv_all_hits'), CIpvsm','cmap',[86/255 180/255 233/255]);


    plot(trimmed_time,mean(accuracy_sh_conv_all_trials'), 'Color',[255/255 0/255 0/255],'LineWidth',2);
    % plot(trimmed_time,mean(accuracy_conv_all_misses'), 'Color', [240/255 228/255 66/255]);
    plot(trimmed_time,mean(accuracy_conv_all_hits'), 'Color',[86/255 180/255 233/255],'LineWidth',2);

    x_pos=-2;
    text(x_pos,0.95,'Hits','Color',[86/255 180/255 233/255])
    % text(x_pos,0.95,'Misses','Color',[240/255 228/255 66/255])
    text(x_pos,0.9,'Shuffled','Color',[255/255 0/255 0/255])

    ylim([0 1.1])

    plot([0 0],[0 1.1],'-k')


    xlabel('Time (sec)')
    ylabel('Accuracy')
    switch which_training_range
        case 1
            title(['Accuracy for Hits, trained with all trials ' algo_name{which_ml_algo}])
        case 2
            title(['Accuracy for Hits, trained with hits ' algo_name{which_ml_algo}])
        case 3
            title(['Accuracy for Hits, trained with misses ' algo_name{which_ml_algo}])
    end



    %Plot accuracy for miss
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    hold on

    CIpvsm = bootci(1000, @mean, accuracy_sh_conv_all_trials');
    meanpvsm=mean(accuracy_sh_conv_all_trials',1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(trimmed_time,mean(accuracy_sh_conv_all_trials'), CIpvsm','cmap',[255/255 0/255 0/255]);

    no_misses=sum((trimmed_trial_type==2)|(trimmed_trial_type==4));
    accuracy_conv_all_misses=zeros(no_misses,size(accuracy_conv_all_trials,2));
    accuracy_conv_all_misses=accuracy_conv_all_trials(:,(trimmed_trial_type==2)|(trimmed_trial_type==4));
    CIpvsm = bootci(1000, @mean, accuracy_conv_all_misses');
    meanpvsm=mean(accuracy_conv_all_misses',1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(trimmed_time,mean(accuracy_conv_all_misses'), CIpvsm','cmap',[0/255 158/255 115/255]);


    plot(trimmed_time,mean(accuracy_sh_conv_all_trials'), 'Color',[255/255 0/255 0/255],'LineWidth',2);
    plot(trimmed_time,mean(accuracy_conv_all_misses'), 'Color', [0/255 158/255 115/255],'LineWidth',2);
    % plot(trimmed_time,mean(accuracy_conv_all_hits'), 'Color',[0/255 158/255 115/255]);

    x_pos=-2;
    % text(x_pos,1,'Hits','Color',[0/255 158/255 115/255])
    text(x_pos,0.95,'Misses','Color',[0/255 158/255 115/255])
    text(x_pos,0.9,'Shuffled','Color',[255/255 0/255 0/255])

    ylim([0 1.1])

    plot([0 0],[0 1.1],'-k')


    xlabel('Time (sec)')
    ylabel('Accuracy')
    switch which_training_range
        case 1
            title(['Accuracy for Miss, trained with all trials ' algo_name{which_ml_algo}])
        case 2
            title(['Accuracy for Miss, trained with hits ' algo_name{which_ml_algo}])
        case 3
            title(['Accuracy for Miss, trained with misses ' algo_name{which_ml_algo}])
    end


    %Plot accuracy separately for each event
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    hold on

    %hit1 red
    plot(trimmed_time,mean(accuracy_conv_all_trials(:,trimmed_trial_type==1),2), '-r','LineWidth',2);

    %miss1 cyan
    plot(trimmed_time,mean(accuracy_conv_all_trials(:,trimmed_trial_type==2),2), '-c','LineWidth',2);

    %hit4 blue
    plot(trimmed_time,mean(accuracy_conv_all_trials(:,trimmed_trial_type==3),2), '-b','LineWidth',2);

    %miss4 magenta
    plot(trimmed_time,mean(accuracy_conv_all_trials(:,trimmed_trial_type==4),2), '-m','LineWidth',2);

    plot([0 0],[0 1.1],'-k')

    x_pos=-2;
    text(x_pos,1,'hit1','Color','r')
    text(x_pos,0.95,'miss1','Color','c')
    text(x_pos,0.9,'hit4','Color','b')
    text(x_pos,0.85,'miss4','Color','m')

    ylim([0 1.1])
    xlabel('Time (sec)')
    ylabel('Accuracy')

    switch which_training_range
        case 1
            title(['Smoothed accuracy, trained with all trials ' algo_name{which_ml_algo}])
        case 2
            title(['Smoothed accuracy, trained with hits ' algo_name{which_ml_algo}])
        case 3
            title(['Smoothed accuracy, trained with misses ' algo_name{which_ml_algo}])
    end
end

handles_out.accuracy_conv_all_trials=accuracy_conv_all_trials;
handles_out.accuracy_sh_conv_all_trials=accuracy_sh_conv_all_trials;
handles_out.trimmed_trial_type=trimmed_trial_type;
handles_out.accuracy_sh_all_trials=accuracy_sh_all_trials;
handles_out.accuracy_all_trials=accuracy_all_trials;
handles_out.mean_accuracy_after=mean(mean(accuracy_all_trials(trimmed_time>=0,:),2));
handles_out.mean_sh_accuracy_after=mean(mean(accuracy_sh_all_trials(trimmed_time>=0,:),2));
handles_out.mean_accuracy_miss_after=mean(mean(accuracy_all_trials(trimmed_time>=0,(trimmed_trial_type==2)|(trimmed_trial_type==4))));
handles_out.mean_accuracy_hit_after=mean(mean(accuracy_all_trials(trimmed_time>=0,(trimmed_trial_type==1)|(trimmed_trial_type==3))));
handles_out.mean_accuracy_before=mean(mean(accuracy_all_trials(trimmed_time<0,:),2));
handles_out.mean_sh_accuracy_before=mean(mean(accuracy_sh_all_trials(trimmed_time<0,:),2));
handles_out.mean_accuracy_miss_before=mean(mean(accuracy_all_trials(trimmed_time<0,(trimmed_trial_type==2)|(trimmed_trial_type==4))));
handles_out.mean_accuracy_hit_before=mean(mean(accuracy_all_trials(trimmed_time<0,(trimmed_trial_type==1)|(trimmed_trial_type==3))));


pffft=1;
