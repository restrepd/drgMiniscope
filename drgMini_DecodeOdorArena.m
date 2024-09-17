close all
clear all

%Let's do only some of the trajectory
last_point=5000;

this_path='/Users/restrepd/Documents/Projects/SFTP/Miniscope/20220713_sync/';

dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_8bits_moco_dfof.csv';

arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';

% this_path='/Users/restrepd/Documents/Projects/SFTP/Miniscope/20220518_sync/';
% 
% dFF_file='20220518_FCM6_WithOdor_miniscope_Results.csv';
% 
% arena_file='20220518_FCM6_WithOdor_odorarena.mat';



dFF=readmatrix([this_path dFF_file]);

load([this_path arena_file])

%Kording data

% this_path='/Users/restrepd/Documents/Projects/Neural_Decoding/';
% 
% example_data_hc='hc_data_raw.mat';

% 
% load([this_path example_data_hc])
% 
% arena.xsync=pos(:,1);
% arena.ysync=pos(:,2);

pos=[arena.xsync arena.ysync];


figure(3)

plot(pos(:,1),pos(:,2),'-b')
title('Mouse trajectory')

figure(2)
plot(dFF(1:last_point,3),'-b')

xy_pred=zeros(last_point,2);
x=arena.xsync(1:last_point);
y=arena.ysync(1:last_point);
trimmed_dFF=zeros(last_point,size(dFF,2));
trimmed_dFF(:,:)=dFF(1:last_point,:);

%Do leave one out

for ii=1:last_point
    time_mask=ones(1,last_point);
    time_mask(ii)=0;
    these_training_dFFs=zeros(last_point-1,size(dFF,2)-1);
    these_training_dFFs(:,:)=trimmed_dFF(logical(time_mask),2:end);
    this_xy=zeros(last_point-1,2);
    this_xy(:,1)=x(logical(time_mask));
    this_xy(:,2)=y(logical(time_mask)); 
    
    Mdlx = fitglm(these_training_dFFs,x(logical(time_mask)));
    Mdly = fitglm(these_training_dFFs,y(logical(time_mask)));

%     Mdlx = fitcecoc(these_training_dFFs,x(logical(time_mask)));
%     Mdly = fitcecoc(these_training_dFFs,y(logical(time_mask)));

    this_dFF=zeros(1,size(dFF,2)-1);
    this_dFF(1,:)=trimmed_dFF(ii,2:end);
   
    this_dFF(1,:)=trimmed_dFF(ii,2:end);
    xy_pred(ii,1)=predict(Mdlx,this_dFF);
    xy_pred(ii,2)=predict(Mdly,this_dFF);
end

figure(3)

plot(arena.xsync(1:last_point),arena.ysync(1:last_point),'-b')
title('Mouse trajectory')
hold on
plot(xy_pred(:,1),xy_pred(:,2),'-r')

pffft=1;
