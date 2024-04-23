%drgDecodeOdorArenaKording.m
%does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020
close all
clear all

%Let's do only some of the trajectory
last_point=5000;

figNo=0;

% isKording=0;

%Note: The data brought into the Kording lab jupyter notebbok seems to be
%binned in 200 msec bins
dt=0.2;

%Define the different ranges (training, valid and testing)
training_range=[0, 0.5];
valid_range=[0.5,0.65];
testing_range=[0.65, 0.8];

%The user can define what time period to use spikes from (with respect to the output).
bins_before=10; %How many bins of neural data prior to the output are used for decoding, 4
bins_current=1; %Whether to use concurrent time bin of neural data, 1
bins_after=10; %How many bins of neural data after the output are used for decoding, 0


%Load Kording's data
this_path='/Users/restrepd/Documents/Projects/Neural_Decoding/';
example_data_hc='hc_data_raw.mat';
load([this_path example_data_hc])


%Plot the spike times
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on

for ii_spikes=1:20
    plot([spike_times{1}(ii_spikes) spike_times{1}(ii_spikes)],[0 1],'-b','LineWidth',2)
end

ylim([0 1.2])

%bin the data to dt
no_neurons=length(spike_times);
no_time_bins=ceil(pos_times(end)/dt);
neural_data=zeros(no_time_bins,no_neurons);

for ii_neuron=1:no_neurons
    for ii_spike=1:length(spike_times{ii_neuron})
        this_bin=ceil(spike_times{ii_neuron}(ii_spike)/dt);
        neural_data(this_bin,ii_neuron)=neural_data(this_bin,ii_neuron)+1;
    end
end

pos_binned=zeros(no_time_bins,2);
time_binned=[1:no_time_bins]*dt-dt/2;

for ii_bins=1:no_time_bins
    time_from=time_binned(ii_bins)-dt/2;
    time_to=time_binned(ii_bins)+dt/2;
    pos_binned(ii_bins,:)=mean(pos((pos_times>=time_from)&(pos_times<time_to),:),1);
end
% else
%
%     % this_path='/Users/restrepd/Documents/Projects/SFTP/Miniscope/20220713_sync/';
%     %
%     % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_8bits_moco_dfof.csv';
%     %
%     % arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';
%
%     this_path='/Users/restrepd/Documents/Projects/SFTP/Miniscope/20220518_sync/';
%
%     dFF_file='20220518_FCM6_WithOdor_miniscope_Results.csv';
%
%     arena_file='20220518_FCM6_WithOdor_odorarena.mat';
%
%
%
%     dFF=readmatrix([this_path dFF_file]);
%
%     load([this_path arena_file])
%
%     pos=[];
%     pos(:,1)=arena.xsync;
%     pos(:,2)=arena.ysync;
%     no_time_points=size(pos,1);
%
%     dt_miniscope=1/30;
%     dFF_times=[1:no_time_points]*dt_miniscope;
%
%     no_neurons=size(dFF,2)-1;
%     no_time_bins=ceil(dFF_times(end)/dt);
%     time_binned=[1:no_time_bins]*dt-dt/2;
%     neural_data=zeros(no_time_bins,no_neurons);
%     pos_binned=zeros(no_time_bins,2);
%
%     for ii_time_bin=1:no_time_bins
%         time_from=time_binned(ii_time_bin)-dt/2;
%         time_to=time_binned(ii_time_bin)+dt/2;
%         pos_binned(ii_time_bin,:)=mean(pos((dFF_times>=time_from)&(dFF_times<time_to),:),1);
%         for ii_neuron=1:no_neurons
%             neural_data(ii_time_bin,ii_neuron)=mean(dFF((dFF_times>=time_from)&(dFF_times<time_to),ii_neuron+1),1);
%         end
%     end
% end

no_time_points=size(pos,1);

%Plot the trajectory of the mouse
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

plot(pos(:,1),pos(:,2),'-b')
title('Mouse trajectory')


%Plot the binned trajectory of the mouse (this should look like Kording's
%plots
figNo=figNo+1;
try
    close(figNo)
catch
end
 
hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

plot(pos_binned(2000:5000,1),pos_binned(2000:5000,2),'-b','LineWidth',2)
title('Mouse trajectory binned')

%Plot the binned neural activity, this should look like Kording's
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

plot(neural_data(1:2000,1),'-b','LineWidth',2)
title('Spikes for neuron 1 binned')


% if isKording==1
    %Remove neurons with too few spikes in HC dataset
    neuron_mask=ones(1,no_neurons);
    for ii_neuron=1:no_neurons
        this_no_spikes=length(spike_times{ii_neuron}); %Total number of spikes of each neuron
        if this_no_spikes<100
            neuron_mask(ii_neuron)=0;
        end
    end

    %Some positions are NaN
    nan_mask_1=zeros(1,no_time_bins);
    nan_mask_1(1,:)=isnan(pos_binned(:,1));
    nan_mask_2=zeros(1,no_time_bins);
    nan_mask_2(1,:)=isnan(pos_binned(:,2));
    nan_mask_all=logical(nan_mask_1)|logical(nan_mask_2);

    neural_data_trimmed=zeros(sum(~nan_mask_all),sum(neuron_mask));
    neural_data_trimmed(:,:)=neural_data(logical(~nan_mask_all),logical(neuron_mask));
    no_time_bins=size(neural_data_trimmed,1);
    no_neurons=size(neural_data_trimmed,2);

    pos_binned_trimmed=zeros(no_time_bins,2);
    pos_binned_trimmed(:,:)=pos_binned(logical(~nan_mask_all),:);
% else
%     nan_mask=logical(ones(1,no_time_bins));
%     for ii_neuron=1:no_neurons
%         this_nan_mask=ones(1,no_time_bins);
%         this_nan_mask(1,:)=~isnan(neural_data(:,ii_neuron));
%         nan_mask=nan_mask&this_nan_mask;
%     end
%     neural_data_trimmed=neural_data(nan_mask,:);
%     pos_binned_trimmed=pos_binned(nan_mask,:);
%     no_time_bins=sum(nan_mask);
% end

%Do z scores
mean_neural_data_trimmed_col=mean(neural_data_trimmed,1);
mean_neural_data_trimmed=repmat(mean_neural_data_trimmed_col,no_time_bins,1);

std_neural_data_trimmed_col=std(neural_data_trimmed,1);
std_neural_data_trimmed=repmat(std_neural_data_trimmed_col,no_time_bins,1);

neural_data_trimmed=(neural_data_trimmed-mean_neural_data_trimmed)./std_neural_data_trimmed;

pffft=1;

% Format for Wiener Filter, Wiener Cascade, XGBoost, and Dense Neural Network
% Put in "flat" format, so each "neuron / time" is a single feature
% i.e. each time point in the before and after window becomes a different
% "neuron"
all_bins_per_window=bins_before+bins_after+bins_current;
no_x_flat_neurons=no_neurons*all_bins_per_window;
X_flat=zeros(no_time_bins,no_x_flat_neurons);

for ii_t=bins_before+1:no_time_bins-bins_after
    ii_n=0;
    for no_win=1:all_bins_per_window
        ii_this_t=ii_t-bins_before+no_win-1;
        X_flat(ii_t,ii_n+1:ii_n+no_neurons)=neural_data_trimmed(ii_this_t,:);
        ii_n=ii_n+no_neurons;
    end
end

%Set what part of data should be part of the training/testing/validation sets
%Note that there was a long period of no movement after about 80% of recording, so I did not use this data.

ii_train_range=ceil(training_range*no_time_bins);
if ii_train_range(1)<1
    ii_train_range(1)=1;
end

ii_valid_range=ceil(valid_range*no_time_bins);
ii_test_range=ceil(valid_range*no_time_bins);

Xtrain=X_flat(ii_train_range(1):ii_train_range(2),:);
Xvalid=X_flat(ii_valid_range(1):ii_valid_range(2),:);
Xtest=X_flat(ii_test_range(1):ii_test_range(2),:);

Ytrain=pos_binned_trimmed(ii_train_range(1):ii_train_range(2),:);
Yvalid=pos_binned_trimmed(ii_valid_range(1):ii_valid_range(2),:);
Ytest=pos_binned_trimmed(ii_test_range(1):ii_test_range(2),:);

%Decode using the generatlized linear model
MdlY1 = fitglm(Xtrain,Ytrain(:,1));
MdlY2 = fitglm(Xtrain,Ytrain(:,2));

Y1_valid=predict(MdlY1,Xvalid);
Y2_valid=predict(MdlY2,Xvalid);

fprintf(1, 'R2 glm: %d %d\n',drgGetR2(Yvalid(:,1),Y1_valid),drgGetR2(Yvalid(:,2),Y2_valid));
R1=corrcoef(Yvalid(:,1),Y1_valid);
R2=corrcoef(Yvalid(:,2),Y2_valid);
fprintf(1, 'Correlation coefficient glm: %d %d\n\n',R1(1,2),R2(1,2));

% if isKording==1
y1v_start=2000;
y1v_end=4000;
% else
%     y1v_start=1;
%     y1v_end=length(Y1_valid(:,1));
% end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on
plot(Y1_valid(y1v_start:y1v_end,1),Y2_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('xy for glm')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y1_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),'-b','LineWidth',1.5)
title('x for glm')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y2_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('y for glm')

%Decode using suport vector machine
MdlY1 = fitrsvm(Xtrain,Ytrain(:,1));
MdlY2 = fitrsvm(Xtrain,Ytrain(:,2));

Y1_valid=predict(MdlY1,Xvalid);
Y2_valid=predict(MdlY2,Xvalid);

fprintf(1, 'R2 svm: %d %d\n',drgGetR2(Yvalid(:,1),Y1_valid),drgGetR2(Yvalid(:,2),Y2_valid));
R1=corrcoef(Yvalid(:,1),Y1_valid);
R2=corrcoef(Yvalid(:,2),Y2_valid);
fprintf(1, 'Correlation coefficient svm: %d %d\n\n',R1(1,2),R2(1,2));

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on
plot(Y1_valid(y1v_start:y1v_end,1),Y2_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('xy for svm')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y1_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),'-b','LineWidth',1.5)
title('x for svm')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y2_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('y for svm')

%Decode using neural network
MdlY1 = fitrnet(Xtrain,Ytrain(:,1));
MdlY2 = fitrnet(Xtrain,Ytrain(:,2));

Y1_valid=predict(MdlY1,Xvalid);
Y2_valid=predict(MdlY2,Xvalid);

fprintf(1, 'R2 nn: %d %d\n',drgGetR2(Yvalid(:,1),Y1_valid),drgGetR2(Yvalid(:,2),Y2_valid));
R1=corrcoef(Yvalid(:,1),Y1_valid);
R2=corrcoef(Yvalid(:,2),Y2_valid);
fprintf(1, 'Correlation coefficient nn: %d %d\n\n',R1(1,2),R2(1,2));

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y1_valid(y1v_start:y1v_end,1),Y2_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('xy for neural network')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y1_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),'-b','LineWidth',1.5)
title('x for neural network')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y2_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('y for neural network')



no_conv_points=11;
% conv_win=ones(1,no_conv_points)/no_conv_points;
conv_win_gauss = gausswin(no_conv_points);
conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

Y1_valid_conv=conv(Y1_valid,conv_win_gauss,'same');
Y2_valid_conv=conv(Y2_valid,conv_win_gauss,'same');

%Now limit the x and y to max and min
minY1=min(Yvalid(:,1));
Y1_valid_conv(Y1_valid_conv<minY1)=minY1;
maxY1=max(Yvalid(:,1));
Y1_valid_conv(Y1_valid_conv>maxY1)=maxY1;

minY2=min(Yvalid(:,2));
Y2_valid_conv(Y2_valid_conv<minY2)=minY2;
maxY2=max(Yvalid(:,2));
Y2_valid_conv(Y2_valid_conv>maxY2)=maxY2;


fprintf(1, 'R2 nn conv: %d %d\n',drgGetR2(Yvalid(:,1),Y1_valid_conv),drgGetR2(Yvalid(:,2),Y2_valid_conv));
R1=corrcoef(Yvalid(:,1),Y1_valid_conv);
R2=corrcoef(Yvalid(:,2),Y2_valid_conv);
fprintf(1, 'Correlation coefficient nn conv %d %d\n\n',R1(1,2),R2(1,2));


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y1_valid_conv(y1v_start:y1v_end,1),Y2_valid_conv(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('xy for neural network convolved')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y1_valid_conv(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,1),'-b','LineWidth',1.5)
title('x for neural network convolved')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(Y2_valid_conv(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
plot(Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
title('y for neural network convolved')

pffft=1;

%
% %Decode using nonlinear
% modelfun = @(b,x)b(1) + b(2)*x(:,1) + ...
%     b(3)*x(:,2).^2;
% beta0 = [1 1 1];
% MdlY1 = fitnlm(Xtrain,Ytrain(:,1),modelfun,beta0);
% MdlY2 = fitnlm(Xtrain,Ytrain(:,2),modelfun,beta0);
%
% Y1_valid=predict(MdlY1,Xvalid);
% Y2_valid=predict(MdlY2,Xvalid);
%
% fprintf(1, 'R2 nonlinear: %d %d\n',drgGetR2(Yvalid(:,1),Y1_valid),drgGetR2(Yvalid(:,2),Y2_valid));
% R1=corrcoef(Yvalid(:,1),Y1_valid);
% R2=corrcoef(Yvalid(:,2),Y2_valid);
% fprintf(1, 'Correlation coefficient nonlinear: %d %d\n\n',R1(1,2),R2(1,2));
%
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
%
%
% hold on
% plot(Y1_valid(y1v_start:y1v_end,1),Y2_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
% plot(Yvalid(y1v_start:y1v_end,1),Yvalid(y1v_start:y1v_end,2),'-b','LineWidth',1.5)
% title('xy for neural network')
%
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
%
% hFig = figure(figNo);
%
% set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
%
%
% hold on
% plot(Y1_valid(y1v_start:y1v_end,1),'-r','LineWidth',1.5)
% plot(Yvalid(y1v_start:y1v_end,1),'-b','LineWidth',1.5)
% title('x for neural network')

pfft=1;
%
% %Leave one out
% %This is very slow!
% Y1=[];
% Y2=[];
%
%
% tic
% parfor ii=2000:5000
%     disp(ii)
%     time_mask=ones(1,no_time_bins);
%     time_mask(ii)=0;
%     training_X=zeros(no_time_bins-1,no_x_flat_neurons);
%     training_X(:,:)=X_flat(logical(time_mask),:);
%     training_Y1=zeros(1,no_time_bins-1);
%     training_Y1(1,:)=pos_binned_trimmed(logical(time_mask),1);
%     training_Y2=zeros(1,no_time_bins-1);
%     training_Y2(1,:)=pos_binned_trimmed(logical(time_mask),2);
%
%     MdlY1 = fitglm(training_X,training_Y1);
%     MdlY2 = fitglm(training_X,training_Y2);
%
% %     Mdlx = fitcecoc(these_training_dFFs,x(logical(time_mask)));
% %     Mdly = fitcecoc(these_training_dFFs,y(logical(time_mask)));
%
%     this_X=zeros(1,no_x_flat_neurons);
%     this_X(1,:)=X_flat(ii,:);
%
%
%     Y1(ii)=predict(MdlY1,this_X);
%     Y2(ii)=predict(MdlY2,this_X);
% %     Y_pred(ii-2000+1,:)=[Y1 Y2];
% end
% toc
%
% figure(3)
%
% plot(arena.xsync(1:last_point),arena.ysync(1:last_point),'-b')
% title('Mouse trajectory')
% hold on
% plot(xy_pred(:,1),xy_pred(:,2),'-r')

pffft=1;
