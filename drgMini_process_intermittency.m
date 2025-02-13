%drgMini_process_intermittency.m
% data = 
% 
%   struct with fields:
% 
%      xGrid: [xy_span×xy_span×151 double]
%      yGrid: [xy_span×xy_span×151 double]
%      zGrid: [xy_span×xy_span×151 double]
%     c1Grid: [xy_span×xy_span×151 double]
%     c2Grid: [xy_span×xy_span×151 double]
%     c3Grid: [xy_span×xy_span×151 double]
%     c4Grid: [xy_span×xy_span×151 double]

% sources 1 & 2: x1 = x2 = 0.02 m, y1 = y2 = -0.2 m (5 cm from wall at y = -0.25 m)
% sources 3 & 4: x3 = x4 = 0.02 m, y3 = y4 = 0.2 m (5 cm from wall at y = 0.25 m)

% z1 = z3 = -0.115 m (1 cm from floor at z = - 0.115 m)
% z2 = z4 = -0.105 m (2 cm from floor at z = - 0.105 m)


close all
clear all


    handles_conc.group=1;

    %Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4
    %Group 2 is rewarded, with odor lane 4, no odor in lane 1
    %Group 3 is rewarded, with odor lane 1, no odor in lane 4
    %Group 4 is rewarded, with no odor in lane 1 and lane 4
    %No odor troubleshooting files
    %     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    %     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    %     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';

    handles_conc.is_sphgpu=0;
    is_sphgpu=handles_conc.is_sphgpu;

    % handles_conc.this_path=this_path;
    % handles_conc.dFF_file=dFF_file;
    % handles_conc.arena_file=arena_file;

   

    %     isKording=0;

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    %     dt=0.2;

    %Define the different ranges (training, valid and testing)
    % training_fraction=0.9;
    % handles_conc.training_fraction=training_fraction;

    %     training_range=[0, 0.5];
    %     valid_range=[0.5,0.65];
    %     test_range=[0.5, 1];

    %The user can define what time period to use spikes from (with respect to the output).
    bins_before=0; %How many bins of neural data prior to the output are used for decoding, 10
    bins_current=1; %Whether to use concurrent time bin of neural data, 1
    bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    handles_conc.bins_before=bins_before;
    handles_conc.bins_current=bins_current;
    handles_conc.bins_after=bins_after;


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2; %Time bins for decoding
    dt_miniscope=1/30;
    n_shuffle=5; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_conc.dt=dt;
    handles_conc.dt_miniscope=dt_miniscope;
    handles_conc.n_shuffle=n_shuffle;

    which_training_algorithm=3;
    handles_conc.which_training_algorithm=which_training_algorithm;
    %1=fitrnet
    %2=fitrgp
    %3=fitrtree
    %4=fitglm
    %5=fitrsvm
    %6=fitrtree gpu

    handles_conc.cm_from_floor=2;

    % handles_conc.weber_fechner=1;
    % handles_conc.alpha=1;
    % handles_conc.multiplier=1;
    % handles_conc.lowest_conc=-20;
    % %0 is Stevens Law, R proportional to C^alpha
    % %1 is Weber-Flechner law R proportional to multiplier*log(C)
    % %See Copelli et al DOI: 10.1103/PhysRevE.65.060901


   %Weber-Frechner or Stevens
    handles_conc.weber_fechner=1; 
    %0 is Stevens Law, R proportional to C^alpha
    %1 is Weber-Flechner law R proportional to multiplier*log(C)
    %See Copelli et al DOI: 10.1103/PhysRevE.65.060901

    handles_conc.alpha=1;
    handles_conc.multiplier=1;
    handles_conc.lowest_conc=-200;

    %Hill transform
    handles_conc.hill=0; %0=no Hill transform, 1=Hill transform
    handles_conc.k_half=10^-8; %Hill equation K1/2
    % handles_conc.actual_maxC=(handles_conc.k_half^handles_conc.n_hill); %Measured from the simulated data 0.0043
    
    handles_conc.maxC=10^-6.5;%maxC=0.0043 maximum of simulated odor plume
    handles_conc.n_hill=2; %Hill coefficient


    handles_conc.displayFigures=1;

    handles_conc.trial_start_offset=-15; %This was -10
    handles_conc.trial_end_offset=15;

    % handles_conc.save_tag='OdorConc';


odor_plumes=[];

load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumes.mat')
xy_span=size(data.xGrid,1);
z_span=size(data.xGrid,3);
first_xGrid=zeros(xy_span,xy_span);
first_xGrid(:,:)=data.xGrid(:,:,1);
first_x=zeros(1,xy_span);
first_x(1,:)=first_xGrid(1,:);


first_yGrid=zeros(xy_span,xy_span);
first_yGrid(:,:)=data.yGrid(:,:,1);
first_y=zeros(1,xy_span);
first_y(1,:)=first_yGrid(:,1);
iiy_souce1=find(first_y==-0.2);
iiy_souce4=find(first_y==0.2);



first_z=zeros(1,z_span);
first_z(1,:)=data.zGrid(1,1,:);
iiz_1cm=find(first_z==first_z(1)+0.01);


% %plot the average concentration vs x
% c1Grid_z_average=zeros(xy_span,xy_span);
% c1Grid_z_average(:,:)=mean(data.c1Grid(:,:,1:iiz_1cm),3);
% 
% c1_vs_x=zeros(1,xy_span);
% c1_vs_x(1,:)=c1Grid_z_average(iiy_souce1,:);
% figure(1)
% plot(first_x,c1_vs_x,'-r')
% xlabel('x (mt)')
% ylabel('Concentration (A.U.)')
% title('source 1')



%Do pseudocolor for each source

%Calculate average plumes
c1Grid_z_average=zeros(xy_span,xy_span);
c1Grid_z_average(:,:)=mean(data.c1Grid(:,:,1:iiz_1cm),3);

c2Grid_z_average=zeros(xy_span,xy_span);
c2Grid_z_average(:,:)=mean(data.c2Grid(:,:,1:iiz_1cm),3);

c3Grid_z_average=zeros(xy_span,xy_span);
c3Grid_z_average(:,:)=mean(data.c3Grid(:,:,1:iiz_1cm),3);

c4Grid_z_average=zeros(xy_span,xy_span);
c4Grid_z_average(:,:)=mean(data.c4Grid(:,:,1:iiz_1cm),3);

maxC=max([max(c1Grid_z_average(:)) max(c2Grid_z_average(:)) max(c3Grid_z_average(:)) max(c4Grid_z_average(:))]);
minC=min([min(c1Grid_z_average(:)) min(c2Grid_z_average(:)) min(c3Grid_z_average(:)) min(c4Grid_z_average(:))]);
 
figNo=0;

%Source 1 
%5 cm from wall at y = -0.25 m, 1 cm from floor at z=-0.115
figNo=figNo+1;
try
    close figNo
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.2 .2 .4 .4])


x=100*repmat(first_x,xy_span,1);
y=100*repmat(first_y',1,xy_span);
drg_pcolor(x,y,c1Grid_z_average)


colormap hot
shading interp
% caxis([minC maxC]);
xlabel('x(cm)')
ylabel('y(cm)');
title('Average odor plume source 1, 1 cm from floor')


source=1;
odor_plumes.source(source).mean_plume=c1Grid_z_average;
odor_plumes.source(source).x=zeros(1,size(x,1));
odor_plumes.source(source).x(1,:)=x(1,:);
odor_plumes.source(source).y=zeros(1,size(x,1));
odor_plumes.source(source).y(1,:)=y(:,1);
odor_plumes.source(source).lane=4;
odor_plumes.source(source).cm_from_floor=1;

%Source 2 
%5 cm from wall at y = -0.25 m, 2 cm from floor
figNo=figNo+1;
try
    close figNo
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.2 .2 .4 .4])


x=100*repmat(first_x,xy_span,1);
y=100*repmat(first_y',1,xy_span);
drg_pcolor(x,y,c2Grid_z_average)


colormap hot
shading interp
% caxis([minC maxC]);
xlabel('x(cm)')
ylabel('y(cm)');
title('Average odor plume source 2, 2 cm from floor')

source=2;
odor_plumes.source(source).mean_plume=c2Grid_z_average;
odor_plumes.source(source).x=zeros(1,size(x,1));
odor_plumes.source(source).x(1,:)=x(1,:);
odor_plumes.source(source).y=zeros(1,size(x,1));
odor_plumes.source(source).y(1,:)=y(:,1);
odor_plumes.source(source).lane=4;
odor_plumes.source(source).cm_from_floor=2;

%Source 3 
%5 cm from wall at y = 0.25 m, 1 cm from floor at z=-0.115
figNo=figNo+1;
try
    close figNo
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.2 .2 .4 .4])


x=100*repmat(first_x,xy_span,1);
y=100*repmat(first_y',1,xy_span);
drg_pcolor(x,y,c3Grid_z_average)


colormap hot
shading interp
% caxis([minC maxC]);
xlabel('x(cm)')
ylabel('y(cm)');
title('Average odor plume source 3, 1 cm from floor')

source=3;
odor_plumes.source(source).mean_plume=c3Grid_z_average;
odor_plumes.source(source).x=zeros(1,size(x,1));
odor_plumes.source(source).x(1,:)=x(1,:);
odor_plumes.source(source).y=zeros(1,size(x,1));
odor_plumes.source(source).y(1,:)=y(:,1);
odor_plumes.source(source).lane=1;
odor_plumes.source(source).cm_from_floor=1;

%Source 4 
%5 cm from wall at y = 0.25 m, 2 cm from floor
figNo=figNo+1;
try
    close figNo
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.2 .2 .4 .4])


x=100*repmat(first_x,xy_span,1);
y=100*repmat(first_y',1,xy_span);
drg_pcolor(x,y,c4Grid_z_average)


colormap hot
shading interp
% caxis([minC maxC]);
xlabel('x(cm)')
ylabel('y(cm)');
title('Average odor plume source 4, 2 cm from floor')

source=4;
odor_plumes.source(source).mean_plume=c4Grid_z_average;
odor_plumes.source(source).x=zeros(1,size(x,1));
odor_plumes.source(source).x(1,:)=x(1,:);
odor_plumes.source(source).y=zeros(1,size(x,1));
odor_plumes.source(source).y(1,:)=y(:,1);
odor_plumes.source(source).lane=1;
odor_plumes.source(source).cm_from_floor=2;

% save('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumesDR2.mat','odor_plumes')
   
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
  
pfft=1