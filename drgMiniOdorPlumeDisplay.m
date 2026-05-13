close all
clear all
load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumesDR.mat')

figNo=0;

%Use Aaron's simulation
cm_from_floor=2;

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
 


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l4')
    colormap fire
    shading interp
    % caxis([minC maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')

    title('Mean odor plume lane 4 before shift')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l1')
    colormap fire
    shading interp
    % caxis([minC maxC]);
    set(gca, 'YDir', 'reverse');
    % xlim(x_range)
    % ylim(y_range)
    % Ax = gca;
    % Ax.Color = 'k';
    xlabel('x (mm)')
    ylabel('y (mm)')

    title('Mean odor plume lane 1 before shift')
