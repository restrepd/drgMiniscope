function drgDisplayDecodeOdorArenav2
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020
close all
clear all


[FileName,PathName] = uigetfile({'*.mat'},'Select the dec .mat file for analysis');

load([PathName FileName])


fprintf(1,['Decoding analysis for ' PathName FileName])
x_predicted_sh=handles_out.x_predicted_sh;
y_predicted_sh=handles_out.y_predicted_sh;
x_predicted=handles_out.x_predicted;
y_predicted=handles_out.y_predicted;
XYtest=handles_out.XYtest;
trials=handles_out.trials;

which_training_algorithm=handles_choices.which_training_algorithm;
dt=handles_choices.dt;
dt_miniscope=handles_choices.dt_miniscope;
n_shuffle=handles_choices.n_shuffle;

x_predictedstart=1;
x_predictedend=length(x_predicted(:,1));

switch which_training_algorithm
    case 1
        fprintf(1,['\nTrained with data for the entire session\n\n'])
    case 2
        fprintf(1,['\nTrained with within trial data\n\n'])
    case 3
        fprintf(1,['\nTrained with between trial data\n\n'])
end

figNo=0;


no_conv_points=11;
% conv_win=ones(1,no_conv_points)/no_conv_points;
conv_win_gauss = gausswin(no_conv_points);
conv_win_gauss=conv_win_gauss/sum(conv_win_gauss);

x_predicted_conv=conv(x_predicted,conv_win_gauss,'same');
y_predicted_conv=conv(y_predicted,conv_win_gauss,'same');

%Now limit the x and y to max and min
minY1=min(XYtest(:,1));
x_predicted_conv(x_predicted_conv<minY1)=minY1;
maxY1=max(XYtest(:,1));
x_predicted_conv(x_predicted_conv>maxY1)=maxY1;

minY2=min(XYtest(:,2));
y_predicted_conv(y_predicted_conv<minY2)=minY2;
maxY2=max(XYtest(:,2));
y_predicted_conv(y_predicted_conv>maxY2)=maxY2;

x_predicted_sh_conv=zeros(size(x_predicted_sh,1),size(x_predicted_sh,2));
y_predicted_sh_conv=zeros(size(y_predicted_sh,1),size(y_predicted_sh,2));
for ii_sh=1:n_shuffle
    this_x_predicted_sh=zeros(size(x_predicted_sh,1),1);
    this_x_predicted_sh(:,1)=x_predicted_sh(:,ii_sh);
    this_y_predicted_sh=zeros(size(y_predicted_sh,1),1);
    this_y_predicted_sh(:,1)=y_predicted_sh(:,ii_sh);
    this_x_predicted_sh_conv=[];
    this_y_predicted_sh_conv=[];
    this_x_predicted_sh_conv=conv(this_x_predicted_sh,conv_win_gauss,'same');
    this_y_predicted_sh_conv=conv(this_y_predicted_sh,conv_win_gauss,'same');

    %Now limit the x and y to max and min
    minY1=min(XYtest(:,1));
    this_x_predicted_sh_conv(this_x_predicted_sh_conv<minY1)=minY1;
    maxY1=max(XYtest(:,1));
    this_x_predicted_sh_conv(this_x_predicted_sh_conv>maxY1)=maxY1;

    minY2=min(XYtest(:,2));
    this_y_predicted_sh_conv(this_y_predicted_sh_conv<minY2)=minY2;
    maxY2=max(XYtest(:,2));
    this_y_predicted_sh_conv(this_y_predicted_sh_conv>maxY2)=maxY2;

    x_predicted_sh_conv(:,ii_sh)=this_x_predicted_sh_conv;
    y_predicted_sh_conv(:,ii_sh)=this_y_predicted_sh_conv;
end




% fprintf(1, 'R2 nn conv x, y entire run: %d %d\n',drgGetR2(XYtest(:,1),x_predicted_conv),drgGetR2(XYtest(:,2),y_predicted_conv));
R1=corrcoef(XYtest(:,1),x_predicted_conv);
R2=corrcoef(XYtest(:,2),y_predicted_conv);
fprintf(1, 'Correlation coefficient nn conv x, y entire run %d %d\n\n',R1(1,2),R2(1,2));


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(x_predicted_conv(x_predictedstart:x_predictedend,1),y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1.5)
plot(XYtest(x_predictedstart:x_predictedend,1),XYtest(x_predictedstart:x_predictedend,2),'-b','LineWidth',1.5)
switch which_training_algorithm
    case 1
        title('xy for neural network convolved (trained per session)')
    case 2
        title('xy for neural network convolved (trained per trial)')
    case 3
        title('xy for neural network convolved (trained between)')
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on

x_predictedstart=1;
x_predictedend=length(x_predicted(:,1));
plot(XYtest(x_predictedstart:x_predictedend,1),'-b','LineWidth',3)
plot(x_predicted_conv(x_predictedstart:x_predictedend,1),'-r','LineWidth',1)

for trNo=1:trials.odor_trNo
    this_x_predictedstart=trials.odor_ii_start(trNo)-10;
    this_x_predictedend=trials.odor_ii_end(trNo)+15;
    plot([this_x_predictedstart:this_x_predictedend],75*ones(this_x_predictedend-this_x_predictedstart+1,1),'-k','LineWidth',2)
end


switch which_training_algorithm
    case 1
        title('x for nn, b:original, r:predicted (trained per session)')
    case 2
        title('x for nn, b:original, r:predicted (trained per trial)')
    case 3
        title('x for nn, b:original, r:predicted (trained between)')
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(XYtest(x_predictedstart:x_predictedend,1),x_predicted_conv(x_predictedstart:x_predictedend,1),'.b')

xlabel('Actual x')
ylabel('Decoded x')

switch which_training_algorithm
    case 1
        title('x for nn (trained per session)')
    case 2
        title('x for nn (trained per trial)')
    case 3
        title('x for nn (trained between)')
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
plot(XYtest(x_predictedstart:x_predictedend,2),'-b','LineWidth',3)
plot(y_predicted_conv(x_predictedstart:x_predictedend,1),'-r','LineWidth',1)

for trNo=1:trials.odor_trNo
    this_x_predictedstart=trials.odor_ii_start(trNo)-10;
    this_x_predictedend=trials.odor_ii_end(trNo)+15;
    plot([this_x_predictedstart:this_x_predictedend],75*ones(this_x_predictedend-this_x_predictedstart+1,1),'-k','LineWidth',2)
end


switch which_training_algorithm
    case 1
        title('y for nn, b:original, r:predicted (trained per session)')
    case 2
        title('y for nn, b:original, r:predicted (trained per trial)')
    case 3
        title('y for nn, b:original, r:predicted (trained between)')
end


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on
plot(XYtest(x_predictedstart:x_predictedend,2),y_predicted_conv(x_predictedstart:x_predictedend,1),'.b')
xlabel('Actual y')
ylabel('Decoded y')

switch which_training_algorithm
    case 1
        title('y for nn(trained per session)')
    case 2
        title('y for nn (trained per trial)')
    case 3
        title('y for nn (trained between)')
end


%Keep track of the per trial decoding
x_all_trials=[];
y_all_trials=[];
x_decod_all_trials=[];
y_decod_all_trials=[];
x_all_trials_sh=[];
y_all_trials_sh=[];
x_decod_all_trials_sh=[];
y_decod_all_trials_sh=[];

x_between_trials=[];
y_between_trials=[];
x_decod_between_trials=[];
y_decod_between_trials=[];
x_between_trials_sh=[];
y_between_trials_sh=[];
x_decod_between_trials_sh=[];
y_decod_between_trials_sh=[];

%Plot the per trial results for y for permuted input
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;
last_x_predictedend=1;
for trNo=1:trials.odor_trNo
    
    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    y_all_trials_sh=[y_all_trials_sh; XYtest(x_predictedstart:x_predictedend,2)];
    y_decod_all_trials_sh=[y_decod_all_trials_sh; y_predicted_sh_conv(x_predictedstart:x_predictedend,1)];

    y_between_trials_sh=[y_between_trials_sh; XYtest(last_x_predictedend:x_predictedstart,2)];
    y_decod_between_trials_sh=[y_decod_between_trials_sh; y_predicted_sh_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;

    %Plot accuracy per trial
    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    %Plot accuracy per trial for permuted training control
    CIsp = bootci(1000, @mean, y_predicted_sh_conv(x_predictedstart:x_predictedend,:)');
    meansp=mean(y_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;

    [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(y_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1)', CIsp', '-k');

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-r','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-c','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-b','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 4
           %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-m','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;


end

switch which_training_algorithm
    case 1
        title('y for nn, permuted b:original, r:predicted (trained per session)')
    case 2
        title('y for nn, permuted b:original, r:predicted (trained per trial)')
    case 3
        title('y for nn, permuted b:original, r:predicted (trained between)')
end


%Plot the per trial results for y
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;
last_x_predictedend=1;
for trNo=1:trials.odor_trNo
    
    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    y_all_trials=[y_all_trials; XYtest(x_predictedstart:x_predictedend,2)];
    y_decod_all_trials=[y_decod_all_trials; y_predicted_conv(x_predictedstart:x_predictedend,1)];

    y_between_trials=[y_between_trials; XYtest(last_x_predictedend:x_predictedstart,2)];
    y_decod_between_trials=[y_decod_between_trials; y_predicted_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;

    %Plot accuracy per trial
    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-r','LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-c','LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-b','LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 4
           %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,2),'-m','LineWidth',3)
            plot(dt*[ii_start:ii_end]',y_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;


end

switch which_training_algorithm
    case 1
        title('y for nn per trial, b:original, r:predicted (trained per session)')
    case 2
        title('y for nn per trial, b:original, r:predicted (trained per trial)')
    case 3
        title('y for nn per trial, b:original, r:predicted (trained between)')
end

%Plot the per trial results for x
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .7 .3])


hold on
ii_start=0;

last_x_predictedend=1;
for trNo=1:trials.odor_trNo
    
    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    x_all_trials=[x_all_trials; XYtest(x_predictedstart:x_predictedend,1)];
    x_decod_all_trials=[x_decod_all_trials; x_predicted_conv(x_predictedstart:x_predictedend,1)];

    x_between_trials=[x_between_trials; XYtest(last_x_predictedend:x_predictedstart,1)];
    x_decod_between_trials=[x_decod_between_trials; x_predicted_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;
    
    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-r','LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-c','LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-b','LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 4
           %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-m','LineWidth',3)
            plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;
end

switch which_training_algorithm
    case 1
        title('x for nn per trial, r: hit1, c: mis1, b: hit4, m: miss4 (trained per session)')
    case 2
        title('x for nn per trial, r: hit1, c: mis1, b: hit4, m: miss4 (trained per trial)')
    case 3
        title('x for nn per trial, r: hit1, c: mis1, b: hit4, m: miss4 (trained between)')
end

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
last_x_predictedend=1;
for trNo=1:trials.odor_trNo
    
    x_predictedstart=trials.odor_ii_start(trNo)-10;
    x_predictedend=trials.odor_ii_end(trNo)+15;

    x_all_trials_sh=[x_all_trials_sh; XYtest(x_predictedstart:x_predictedend,1)];
    x_decod_all_trials_sh=[x_decod_all_trials_sh; x_predicted_sh_conv(x_predictedstart:x_predictedend,1)];

    x_between_trials_sh=[x_between_trials_sh; XYtest(last_x_predictedend:x_predictedstart,1)];
    x_decod_between_trials_sh=[x_decod_between_trials_sh; x_predicted_sh_conv(last_x_predictedend:x_predictedstart,1)];
    last_x_predictedend=x_predictedend;
    
    ii_end=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))-1;

    %Plot accuracy per trial for permuted training control
    CIsp = bootci(1000, @mean, x_predicted_sh_conv(x_predictedstart:x_predictedend,:)');
    meansp=mean(x_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;

    [hlsp, hpsp] = boundedline(dt*[ii_start:ii_end]',mean(x_predicted_sh_conv(x_predictedstart:x_predictedend,:)',1)', CIsp', '-k');


    switch trials.odor_trial_type(trNo)
        case 1
            %Lane 1 hits
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-r','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 2
            %Lane 1 miss
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-c','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 3
            %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-b','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
        case 4
           %Lane 4 hit
            plot(dt*[ii_start:ii_end]',XYtest(x_predictedstart:x_predictedend,1),'-m','LineWidth',3)
%             plot(dt*[ii_start:ii_end]',x_predicted_conv(x_predictedstart:x_predictedend,1),'-k','LineWidth',1)
            plot(dt*[ii_start+10 ii_start+10],[0 250],'-k')
            plot(dt*[ii_end-15 ii_end-15],[0 250],'-k')
    end
    ii_start=ii_start+length(XYtest(x_predictedstart:x_predictedend,2))+20;
end

switch which_training_algorithm
    case 1
        title('x for nn per trial permuted, r: hit1, c: mis1, b: hit4, m: miss4 (trained per session)')
    case 2
        title('x for nn per trial permuted, r: hit1, c: mis1, b: hit4, m: miss4 (trained per trial)')
    case 3
        title('x for nn per trial permuted, r: hit1, c: mis1, b: hit4, m: miss4 (trained between)')
end


% fprintf(1, 'R2 nn conv x, y per trial run: %d %d\n',drgGetR2(x_all_trials,x_decod_all_trials),drgGetR2(y_all_trials,y_decod_all_trials));
R1=corrcoef(x_all_trials,x_decod_all_trials);
R2=corrcoef(y_all_trials,y_decod_all_trials);
fprintf(1, 'Correlation coefficient nn conv x, y per trial %d %d\n',R1(1,2),R2(1,2));

R1=corrcoef(x_between_trials,x_decod_between_trials);
R2=corrcoef(y_between_trials,y_decod_between_trials);
fprintf(1, 'Correlation coefficient nn conv x, y between trials %d %d\n\n',R1(1,2),R2(1,2));

% fprintf(1, 'R2 nn permuted x, y per trial run: %d %d\n',drgGetR2(x_all_trials_sh,x_decod_all_trials_sh),drgGetR2(y_all_trials_sh,y_decod_all_trials_sh));
R1=corrcoef(x_all_trials_sh,x_decod_all_trials_sh);
R2=corrcoef(y_all_trials_sh,y_decod_all_trials_sh);
fprintf(1, 'Correlation coefficient nn permuted x, y per trial %d %d\n',R1(1,2),R2(1,2));

R1=corrcoef(x_between_trials_sh,x_decod_between_trials_sh);
R2=corrcoef(y_between_trials_sh,y_decod_between_trials_sh);
fprintf(1, 'Correlation coefficient nn permuted x, y between trial %d %d\n',R1(1,2),R2(1,2));

%Plot x vs decoded
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

plot(x_all_trials,x_decod_all_trials,'.b')
xlabel('x')
ylabel('x decoded')
switch which_training_algorithm
    case 1
        title('x for nn per trial(trained per session)')
    case 2
        title('x for nn per trial (trained per trial)')
    case 3
        title('x for nn per trial (trained between)')
end

%Plot y vs decoded
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

plot(y_all_trials,y_decod_all_trials,'.b')
xlabel('y')
ylabel('y decoded')


switch which_training_algorithm
    case 1
        title('y for nn per trial (trained per session)')
    case 2
        title('y for nn per trial (trained per trial)')
    case 3
        title('y for nn per trial (trained between)')
end

pffft=1;
