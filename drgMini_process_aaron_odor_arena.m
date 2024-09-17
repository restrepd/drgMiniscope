%drgMini_process_aaron_odor_arena.m
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

save('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumesDR.mat','odor_plumes')
 
pffft=1;