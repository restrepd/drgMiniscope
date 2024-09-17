function handles_out=drgMini_dFF_class_prediction(handles_choices)
%Does decoding using a classificaiton decoder where dFF> threshold = 1, 0 otherwise
%following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020

warning('off')

if exist('handles_choices')==0
    clear all
    close all

    handles_choices.distances_in_mm=1; %If the distance is already in mm there will not be morphing of space

    %Files for dFF data
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    % dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % 
    % if handles_choices.distances_in_mm==0
    %     arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';
    % else
    %     arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';
    % end


%     %Second troubleshooting files
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';

    if handles_choices.distances_in_mm==0
        arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';
    else
        arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat'; %Distances converted to mm
    end

    %No odor troubleshooting files
%     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
%     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
%     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';



    handles_choices.which_odor_plume=2;   %1=use odor plume from the publication, 2=use odor plume simulated by Aaron
    handles_choices.cm_from_floor=2;

    handles_choices.weber_fechner=1; 
    handles_choices.alpha=1;
    handles_choices.multiplier=1;
    %0 is Stevens Law, R proportional to C^alpha
    %1 is Weber-Flechner law R proportional to log(C)
    %See Copelli et al DOI: 10.1103/PhysRevE.65.060901

    handles_choices.algo=1; 
    %1 ann, fitcnet
    %2 glm, fitglm
    %3 tree, fittctree
    %4 svm with radial basis, fitcsvm
  

    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

    handles_choices.group=1; %1=odor plume, 2=spatial
    
  


    
    handles_choices.repeats=1;
    handles_choices.sh_repeats=handles_choices.repeats;
    
    handles_choices.max_overlap=3; %Maximum overlap of the shuffled segments
    handles_choices.displayFigures=1;
   


    %     isKording=0;

    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    %     dt=0.2;

    %Define the different ranges (training, valid and testing)
    training_fraction=0.9;
    handles_choices.training_fraction=training_fraction;

%     training_range=[0, 0.5];
    %     valid_range=[0.5,0.65];
%     test_range=[0.5, 1];

    %The user can define what time period to use spikes from (with respect to the output).
    bins_before=10; %How many bins of neural data prior to the output are used for decoding, 4
    bins_current=1; %Whether to use concurrent time bin of neural data, 1
    bins_after=0; %How many bins of neural data after the output are used for decoding, 10
    handles_choices.bins_before=bins_before;
    handles_choices.bins_current=bins_current;
    handles_choices.bins_after=bins_after;


    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2;
    dt_miniscope=1/30;
    n_shuffle=5; %Note that n_shuffle is changed to a maximum of ii_n_training

    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    handles_choices.n_shuffle=n_shuffle;

    % which_training_algorithm=2; 
    % handles_choices.which_training_algorithm=which_training_algorithm;
    %1=entire session is used for training
    %2=only per trial is used for training
    %3=trained with data between trials

else
    
    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;
    training_fraction=handles_choices.training_fraction;
    bins_before=handles_choices.bins_before;
    bins_current=handles_choices.bins_current;
    bins_after=handles_choices.bins_after;
    dt=handles_choices.dt;
    dt_miniscope=handles_choices.dt_miniscope;
    n_shuffle=handles_choices.n_shuffle;
    
    % which_training_algorithm=handles_choices.which_training_algorithm;

    op_path=handles_choices.op_path;
    op_file=handles_choices.op_file;

end

 

% switch which_training_algorithm
%     case 1
%         fprintf(1,['\nTrained with data for the entire session\n\n'])
%     case 2
%         fprintf(1,['\nTrained with within trial data\n\n'])
%     case 3
%         fprintf(1,['\nTrained with between trial data\n\n'])
% end

figNo=0;

%Load odor plume and calculate the odor plume for lane 1 and lane 4
switch handles_choices.which_odor_plume
    case 1
        
        %File for odor plume
        %11282017_20cms_unbounded__README
        %         Quantification of airborne odor plumes using planar laser‐induced fluorescence
        % E. G. Connor, M. K. McHugh and J. P. Crimaldi
        % Experiments in Fluids 2018 Vol. 59 Pages 137-
        % DOI: https://doi.org/10.1007/s00348-018-2591-3
        op_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Data_Plume/';
        op_file='dataset.mat';

        %Use the odor plume in the publication
        %Plot odor plume
        load([op_path op_file])
        x_range=[0 500]; %in mm
        y_range=[0 480];
        y_lane1=410;    %Note that lane 1 is at 410
        y_lane4=50;

        mean_plume=mean_plume-min(mean_plume(:));
        mean_plume(mean_plume==0)=min(mean_plume(mean_plume~=0));
        if handles_choices.weber_fechner==1
            %Weber-Frechner
            mean_plume=log10(mean_plume);
        end


        if handles_choices.displayFigures==1
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the timecourse
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
            

            drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),mean_plume')
            colormap fire
            shading interp
            % xlim(x_range)
            % ylim(y_range)
            % Ax = gca;
            % Ax.Color = 'k';
            xlabel('x (mm)')
            ylabel('y (mm)')
            title('Mean odor plume')
        end

        %Now make mean_plume matrices for lanes 1 and 4
        delta_x=x(2)-x(1); %please note this is the same as delta_y
        x_for_plume=0:delta_x:500;
        y_for_plume=0:delta_x:480;


        [minabsy,minabsy_ii]=min(abs(y));

        %lane 4
        % mean_plume_lane4=prctile(mean_plume(:),1)*ones(length(y_for_plume),length(x_for_plume));
        if handles_choices.weber_fechner==0
            %Stevens
            mean_plume_lane4=zeros(length(y_for_plume),length(x_for_plume));
        else
            %Weber-Frechner
            mean_plume_lane4=prctile(mean_plume(:),1)*ones(length(y_for_plume),length(x_for_plume));
        end


        [minabsy4,minabsy4_ii]=min(abs(y_for_plume-y_lane4));

        iiy4_to=minabsy4_ii+(length(y)-minabsy_ii);
        delta_iiy4_to=(length(y)-minabsy_ii);
        if iiy4_to>size(mean_plume_lane4,1)
            delta_iiy4_to=size(mean_plume_lane4,1)-minabsy4_ii;
            iiy4_to=size(mean_plume_lane4,1);
        end

        iiy4_from=minabsy4_ii-(length(y)-minabsy_ii);
        delta_iiy4_from=-(length(y)-minabsy_ii);
        if iiy4_from<1
            iiy4_from=1;
            delta_iiy4_from=-(minabsy4_ii-1);
        end

        if handles_choices.weber_fechner==0
            %Stevens
            mean_plume_lane4(iiy4_from:iiy4_to,1:size(mean_plume,2))=(mean_plume(minabsy_ii+delta_iiy4_from:minabsy_ii+delta_iiy4_to,:)).^handles_choices.alpha;
        else
            %Weber-Frechner
            mean_plume_lane4(iiy4_from:iiy4_to,1:size(mean_plume,2))=(mean_plume(minabsy_ii+delta_iiy4_from:minabsy_ii+delta_iiy4_to,:));
        end

        %lane 1
        if handles_choices.weber_fechner==0
            %Stevens
            mean_plume_lane1=zeros(length(y_for_plume),length(x_for_plume));
        else
            %Weber-Frechner
            mean_plume_lane1=prctile(mean_plume(:),1)*ones(length(y_for_plume),length(x_for_plume));
        end


        [minabsy1,minabsy1_ii]=min(abs(y_for_plume-y_lane1));

        iiy1_to=minabsy1_ii+(length(y)-minabsy_ii);
        delta_iiy1_to=(length(y)-minabsy_ii);
        if iiy1_to>size(mean_plume_lane1,1)
            delta_iiy1_to=size(mean_plume_lane1,1)-minabsy1_ii;
            iiy1_to=size(mean_plume_lane1,1);
        end

        iiy1_from=minabsy1_ii-(length(y)-minabsy_ii);
        delta_iiy1_from=-(length(y)-minabsy_ii);
        if iiy1_from<1
            iiy1_from=1;
            delta_iiy1_from=minabsy1_ii-1;
        end

        if handles_choices.weber_fechner==0
            %Stevens
            mean_plume_lane1(iiy1_from:iiy1_to,1:size(mean_plume,2))=(mean_plume(minabsy_ii+delta_iiy1_from:minabsy_ii+delta_iiy1_to,:)).^handles_choices.alpha;
        else
            %Weber-Frechner
            mean_plume_lane1(iiy1_from:iiy1_to,1:size(mean_plume,2))=(mean_plume(minabsy_ii+delta_iiy1_from:minabsy_ii+delta_iiy1_to,:));
        end

        minC=min(mean_plume(:));
        maxC=max(mean_plume(:));

        if handles_choices.displayFigures==1
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_lane4')
            colormap fire
            shading interp
            caxis([minC maxC]);
            % xlim(x_range)
            % ylim(y_range)
            % Ax = gca;
            % Ax.Color = 'k';
            xlabel('x (mm)')
            ylabel('y (mm)')

            title('Mean odor plume shifted to lane 4')

            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_lane1')
            colormap fire
            shading interp
            caxis([minC maxC]);
            % xlim(x_range)
            % ylim(y_range)
            % Ax = gca;
            % Ax.Color = 'k';
            xlabel('x (mm)')
            ylabel('y (mm)')

            title('Mean odor plume shifted to lane 1')
        end
        pffft=1;
    case 2
        %Use Aaron's simulation
        %Plot odor plume
        load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumesDR.mat')
 
        %Now make mean_plume matrices for lanes 1 and 4

        %lane 1
        this_lane=1;
        for ii_source=1:length(odor_plumes.source)
            if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==handles_choices.cm_from_floor)
                this_source=ii_source;
            end
        end

        x_for_plume=10*odor_plumes.source(this_source).x;
        y_for_plume=10*(odor_plumes.source(this_source).y-min(odor_plumes.source(this_source).y));
        mean_plume_lane1=odor_plumes.source(this_source).mean_plume;
        mean_plume_lane1=mean_plume_lane1-min(mean_plume_lane1(:));
        mean_plume_lane1(mean_plume_lane1==0)=min(mean_plume_lane1(mean_plume_lane1~=0));

        %Shift the plume to 7 cm (20 mm)
       

        if handles_choices.weber_fechner==0
             mean_plume_lane1=handles_choices.multiplier*mean_plume_lane1.^handles_choices.alpha;
        else
            %Weber-Frechner
            mean_plume_lane1=handles_choices.multiplier*log10(mean_plume_lane1);
        end

        %lane 4
        this_lane=4;
        for ii_source=1:length(odor_plumes.source)
            if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==handles_choices.cm_from_floor)
                this_source=ii_source;
            end
        end

        % x_for_plume=10*odor_plumes.source(this_source).x;
        % y_for_plume=10*(odor_plumes.source(this_source).y-min(odor_plumes.source(this_source).y));
        mean_plume_lane4=odor_plumes.source(this_source).mean_plume;
        mean_plume_lane4=mean_plume_lane4-min(mean_plume_lane4(:));
        mean_plume_lane4(mean_plume_lane4==0)=min(mean_plume_lane4(mean_plume_lane4~=0));

        if handles_choices.weber_fechner==0
             mean_plume_lane4=handles_choices.multiplier*mean_plume_lane4.^handles_choices.alpha;
        else
            %Weber-Frechner
            mean_plume_lane4=handles_choices.multiplier*log10(mean_plume_lane4);
        end

        minC=min([min(mean_plume_lane1(:)) min(mean_plume_lane4(:))]);
        maxC=min([max(mean_plume_lane1(:)) max(mean_plume_lane4(:))]);

        %Now shift to the actual dimensions of the odorant arena
        new_mean_plume_lane4=mean_plume_lane4(y_for_plume<=480,:);
        mean_plume_lane4=[];
        mean_plume_lane4=new_mean_plume_lane4;
        new_mean_plume_lane1=mean_plume_lane1(y_for_plume>=20,:);

        %Now shift the center to 7 cm from end
        dy=abs(y_for_plume(2)-y_for_plume(1));
        ii_dy_shift=20/dy;
        mean_plume_lane1=zeros(size(new_mean_plume_lane1,1),size(new_mean_plume_lane1,2));
        mean_plume_lane1(ii_dy_shift+1:end,:)=new_mean_plume_lane1(1:size(new_mean_plume_lane1,1)-ii_dy_shift,:);
        new_y_for_plume=y_for_plume(:,y_for_plume<=480);
        y_for_plume=[];
        y_for_plume=new_y_for_plume;
        ii_y_center=find(y_for_plume==410);
        mean_plume_lane1(1:ii_y_center-1,:)=mean_plume_lane1(ii_y_center+ii_y_center-1:-1:ii_y_center+1,:);

        if handles_choices.distances_in_mm==1
            mean_plume_lane1(:,:)=mean_plume_lane1(end:-1:1,:);
            mean_plume_lane4(:,:)=mean_plume_lane4(end:-1:1,:);
        end

        if handles_choices.displayFigures==1
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_lane4')
            colormap fire
            shading interp
            caxis([minC maxC]);
            set(gca, 'YDir', 'reverse');
            % xlim(x_range)
            % ylim(y_range)
            % Ax = gca;
            % Ax.Color = 'k';
            xlabel('x (mm)')
            ylabel('y (mm)')

            title('Mean odor plume shifted to lane 4')

            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_lane1')
            colormap fire
            shading interp
            caxis([minC maxC]);
            set(gca, 'YDir', 'reverse');
            % xlim(x_range)
            % ylim(y_range)
            % Ax = gca;
            % Ax.Color = 'k';
            xlabel('x (mm)')
            ylabel('y (mm)')

            title('Mean odor plume shifted to lane 1')
        end

        pffft=1;
end

%Restart random seeds
rng('shuffle');


dFF=[];
if strcmp(dFF_file(end-3:end),'.mat')
    %This reads the extract file
    load([this_path dFF_file])
    dFF=zeros(size(output.temporal_weights,1),size(output.temporal_weights,2));
    for traceNo=1:size(output.temporal_weights,2)
        dFF(:,traceNo)=output.temporal_weights(:,traceNo);
    end
else
    if strcmp(dFF_file(end-4:end),'.hdf5')
        %This is an hdf5 generated by CaImAn
        dFF=h5read([this_path dFF_file],'/estimates/F_dff' );
    else
        %This is a csv file created from ImageJ
        dFF=readmatrix([this_path dFF_file]);
    end
end





load([this_path arena_file])

%Note that Ryan divided by 4
if handles_choices.distances_in_mm~=1
    arena.xsync=4*arena.xsync;
    arena.ysync=4*arena.ysync;
end

%Extract trials
trials=[];


%Extract odor on using the camera sync
at_end=0;
ii=0;
jj=0;
jj_l1=0;
jj_l4=0;
while at_end==0
    next_ii=find(arena.odorsync(ii+1:end)==1,1,'first');
    if ~isempty(next_ii)
        jj=jj+1;
        trials.ii_odor(jj)=ii+next_ii;
        trials.x_odor(jj)=arena.xsync(ii+next_ii);
        trials.y_odor(jj)=arena.ysync(ii+next_ii);
      
        ii=ii+next_ii;
        ii_mini=arena.index_flirsynctominiscope(ii);

        if sum(arena.laneodor1(ii_mini-5:ii_mini+5)==1)>0
            %Note: laneodor4 is 1 only for one time point
            jj_l1=jj_l1+1;
            trials.ii_laneodor1(jj_l1)=ii;
            trials.x_laneodor1(jj_l1)=trials.x_odor(jj);
            trials.y_laneodor1(jj_l1)=trials.y_odor(jj);
        end

        if sum(arena.laneodor4(ii_mini-5:ii_mini+5)==1)>0
            %Note: laneodor4 is 1 only for one time point
            jj_l4=jj_l4+1;
            trials.ii_laneodor4(jj_l4)=ii;
            trials.x_laneodor4(jj_l4)=trials.x_odor(jj);
            trials.y_laneodor4(jj_l4)=trials.y_odor(jj);
        end

        next_ii=find(arena.odorsync(ii+1:end)==0,1,'first');
        if ~isempty(next_ii)
            ii=ii+next_ii;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end


%Extract lanewater1
at_end=0;
ii_flir=0;
jj=0;
while at_end==0
    next_ii_flir=find(arena.lanewater1(ii_flir+1:end)==1,1,'first');
    if ~isempty(next_ii_flir)
        next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
        jj=jj+1;
        trials.ii_lanewater1(jj)=next_ii;
        trials.x_lanewater1(jj)=arena.xsync(next_ii);
        trials.y_lanewater1(jj)=arena.ysync(next_ii);
        next_ii_flir2=find(arena.lanewater1(ii_flir+next_ii_flir+1:end)==0,1,'first');
%         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
        if ~isempty(next_ii_flir2)
            ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
        else
            at_end=1;
        end
    else
         at_end=1;
    end
end


%Extract lanewater4
at_end=0;
ii_flir=0;
jj=0;
while at_end==0
    next_ii_flir=find(arena.lanewater4(ii_flir+1:end)==1,1,'first');
    if ~isempty(next_ii_flir)
        next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir+ii_flir,1,'last');
        jj=jj+1;
        trials.ii_lanewater4(jj)=next_ii;
        trials.x_lanewater4(jj)=arena.xsync(next_ii);
        trials.y_lanewater4(jj)=arena.ysync(next_ii);
        next_ii_flir2=find(arena.lanewater4(ii_flir+next_ii_flir+1:end)==0,1,'first');
%         next_ii=find(arena.index_flirsynctominiscope<=next_ii_flir,1,'last');
        if ~isempty(next_ii_flir2)
            ii_flir=ii_flir+next_ii_flir+next_ii_flir2;
        else
            at_end=1;
        end
    else
         at_end=1;
    end
end


if handles_choices.displayFigures==1
    %Display the location of trial start and water delivery
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on
    plot(trials.x_laneodor1,trials.y_laneodor1,'or')
    plot(trials.x_laneodor4,trials.y_laneodor4,'ob')
    % plot([95 95],[200 250],'-r')
    % plot([95 95],[0 50],'-b')

    plot(trials.x_lanewater1,trials.y_lanewater1,'xr')
    plot(trials.x_lanewater4,trials.y_lanewater4,'xb')
    set(gca, 'YDir', 'reverse');
    xlabel('x (pixels)')
    ylabel('y (pixels)')
    title('Location of trial start (o) and water delivery (x)')
    

    %Plot mouse pathway
        figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on
    these_no_points=150;
    plot(arena.xsync(end),arena.ysync(end),'ok')
    % plot(arena.xsync(end-these_no_points:end),arena.ysync(end-these_no_points:end),'-b')
    plot(arena.xsync,arena.ysync,'-b')
    set(gca, 'YDir', 'reverse');
    % title(['End of trajectory for ' arena_file])
    title(['Trajectory for ' arena_file])
    ylabel('y')
    xlabel('x')
    legend('end of trajectory','trajectory')
 
end

if handles_choices.distances_in_mm==1
    %The distances are already in mm

    %Bin positions into dt time bins
    pos=[];
    pos(:,1)=arena.xsync;
    pos(:,2)=arena.ysync;
    no_time_points=size(pos,1);

else
    %The physical odor arena is 50 cm along the air flow and 48 cm wide,
    %with lanes 1 and 4, 5 cm from the wall
    %We need to morph the video camera locations to the physical location
    %For now I will do this isometrically
    y_mean_video_lane1=mean(trials.y_lanewater1);
    y_mean_video_lane4=mean(trials.y_lanewater4);
    y_mean_video_lanes14=mean([y_mean_video_lane1 y_mean_video_lane4]);
    delta_y_physical_lane14=380;
    y_half_physical=480/2;


    %Bin positions into dt time bins
    pos=[];
    pos(:,1)=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(arena.xsync-95);
    pos(:,2)=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(arena.ysync-y_mean_video_lanes14)+y_half_physical;
    no_time_points=size(pos,1);

    %Plot morphed end locations
    if handles_choices.displayFigures==1
        %Display the location of trial start and water delivery
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

        hold on

        morphed_x1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_laneodor1-95);
        morphed_y1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_laneodor1-y_mean_video_lanes14)+y_half_physical;
        plot(morphed_x1,morphed_y1,'or')

        morphed_x4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_laneodor4-95);
        morphed_y4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_laneodor4-y_mean_video_lanes14)+y_half_physical;
        plot(morphed_x4,morphed_y4,'ob')

        % plot([95 95],[200 250],'-r')
        % plot([95 95],[0 50],'-b')

        morphed_x1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_lanewater1-95);
        morphed_y1=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_lanewater1-y_mean_video_lanes14)+y_half_physical;
        plot(morphed_x1,morphed_y1,'xr')

        morphed_x4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.x_lanewater4-95);
        morphed_y4=(delta_y_physical_lane14/(y_mean_video_lane1-y_mean_video_lane4))*(trials.y_lanewater4-y_mean_video_lanes14)+y_half_physical;
        plot(morphed_x4,morphed_y4,'xb')
        set(gca, 'YDir', 'reverse');
        xlabel('x (mm)')
        ylabel('y (mm)')
        title('Location of trial start (o) and water delivery (x), x and y morphed')
        %Note: This is clearly not correct, the mice move to 25 mm and 225 mm
        %This is half as large as the physical odor arena!!


    end

end

dFF_times=[1:no_time_points]*dt_miniscope;

no_neurons=size(dFF,2)-1;
no_time_bins=ceil(dFF_times(end)/dt);
time_binned=[1:no_time_bins]*dt-dt/2;
neural_data=zeros(no_time_bins,no_neurons);
pos_binned=zeros(no_time_bins,2);

for ii_time_bin=1:no_time_bins
    time_from=time_binned(ii_time_bin)-dt/2;
    time_to=time_binned(ii_time_bin)+dt/2;
    pos_binned(ii_time_bin,:)=mean(pos((dFF_times>=time_from)&(dFF_times<time_to),:),1);
    for ii_neuron=1:no_neurons
        neural_data(ii_time_bin,ii_neuron)=mean(dFF((dFF_times>=time_from)&(dFF_times<time_to),ii_neuron+1),1);
    end
end

%Convert neural_data to binary

%Calculate the threshold for binary conversion
neural_data_threshold=prctile(neural_data(:),5);
neural_data_binary=zeros(size(neural_data,1),size(neural_data,2));

logical_threshold=neural_data > neural_data_threshold;
neural_data_binary(logical_threshold)=1;

logical_threshold=neural_data <= neural_data_threshold;
neural_data_binary(logical_threshold)=0;

% dFF=readmatrix([this_path dFF_file]); %Timepoints x ROIs

%Now calculate the behavioral performance
trials.hit1=0;
trials.miss1=0;
trials.hit4=0;
trials.miss4=0;
trials.odor_trNo=0;

trim_factor=no_time_bins/no_time_points;

dii_trial=[];
for trNo=1:length(trials.ii_odor)

    this_ii_laneodor1=find(abs(trials.ii_odor(trNo)-trials.ii_laneodor1)<3);
    if ~isempty(this_ii_laneodor1)
        if trNo==length(trials.ii_odor)
            ii_next=length(arena.xsync);
        else
            ii_next=trials.ii_odor(trNo+1);
        end
        this_water=find((trials.ii_lanewater1>trials.ii_odor(trNo))&(trials.ii_lanewater1<ii_next));
        if ~isempty(this_water)
            trials.hit1=trials.hit1+1;
            trials.hit1_ii_start(trials.hit1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.hit1_ii_end(trials.hit1)=ceil(trim_factor*trials.ii_lanewater1(this_water));
            dii_trial=[dii_trial trials.hit1_ii_end(trials.hit1)-trials.hit1_ii_start(trials.hit1)];
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.hit1_ii_start(trials.hit1);
            trials.odor_ii_end(trials.odor_trNo)=trials.hit1_ii_end(trials.hit1);
            trials.odor_trial_type(trials.odor_trNo)=1;
            trials.odor_lane(trials.odor_trNo)=1;
        else
            trials.miss1=trials.miss1+1;
            trials.miss1_ii_start(trials.miss1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss1_ii_start(trials.miss1);
            trials.odor_trial_type(trials.odor_trNo)=2;
            trials.odor_lane(trials.odor_trNo)=1;
        end
    end

    this_ii_laneodor4=find(abs(trials.ii_odor(trNo)-trials.ii_laneodor4)<3);
    if ~isempty(this_ii_laneodor4)
        if trNo==length(trials.ii_odor)
            ii_next=length(arena.xsync);
        else
            ii_next=trials.ii_odor(trNo+1);
        end
        this_water=find((trials.ii_lanewater4>trials.ii_odor(trNo))&(trials.ii_lanewater4<ii_next));
        if ~isempty(this_water)
            trials.hit4=trials.hit4+1;
            trials.hit4_ii_start(trials.hit4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.hit4_ii_end(trials.hit4)=ceil(trim_factor*trials.ii_lanewater4(this_water));
            dii_trial=[dii_trial trials.hit4_ii_end(trials.hit4)-trials.hit4_ii_start(trials.hit4)];
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.hit4_ii_start(trials.hit4);
            trials.odor_ii_end(trials.odor_trNo)=trials.hit4_ii_end(trials.hit4);
            trials.odor_trial_type(trials.odor_trNo)=3;
            trials.odor_lane(trials.odor_trNo)=4;
        else
            trials.miss4=trials.miss4+1;
            trials.miss4_ii_start(trials.miss4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.miss4_ii_start(trials.miss4);
            trials.odor_trial_type(trials.odor_trNo)=4;
            trials.odor_lane(trials.odor_trNo)=4;
        end
    end

end

for ii_miss=1:trials.miss1
    trials.miss1_ii_end(ii_miss)=trials.miss1_ii_start(ii_miss)+ceil(mean(dii_trial));
end

for ii_miss=1:trials.miss4
    trials.miss4_ii_end(ii_miss)=trials.miss4_ii_start(ii_miss)+ceil(mean(dii_trial));
end

for ii_odor=1:trials.odor_trNo
    if trials.odor_trial_type(ii_odor)==4
        trials.odor_ii_end(ii_odor)=trials.odor_ii_start(ii_odor)+ceil(mean(dii_trial));
    end
    if trials.odor_trial_type(ii_odor)==2
        trials.odor_ii_end(ii_odor)=trials.odor_ii_start(ii_odor)+ceil(mean(dii_trial));
    end
end

percent_correct=100*(trials.hit4+trials.hit1)/(trials.hit4+trials.hit1+trials.miss4+trials.miss1);
percent_correct1=100*(trials.hit1)/(trials.hit1+trials.miss1);
percent_correct4=100*(trials.hit4)/(trials.hit4+trials.miss4);
fprintf(1,['\nPercent correct ' num2str(percent_correct) ' percent correct lane 1 ' num2str(percent_correct1) ' percent correct lane 4 ' num2str(percent_correct4) '\n\n'])

lane1_trials=trials.hit1+trials.miss1;
lane4_trials=trials.hit4+trials.miss4;

[p_hit_miss, Q]= chi2test([trials.hit1, trials.miss1; trials.hit4, trials.miss4]);
fprintf(1,['Chi squared testing difference in hit vs miss between lane 1 and 4 ' num2str(p_hit_miss) '\n\n'])

[p_hit_miss_lane1, Q]= chi2test([trials.hit1, trials.miss1; (trials.hit1+trials.miss1)/2, (trials.hit1+trials.miss1)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for lane 1 ' num2str(p_hit_miss_lane1) '\n\n'])

[p_hit_miss_lane4, Q]= chi2test([trials.hit4, trials.miss4; (trials.hit4+trials.miss4)/2, (trials.hit4+trials.miss4)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for lane 4 ' num2str(p_hit_miss_lane4) '\n\n'])

[p_hit_miss_both_lanes, Q]= chi2test([trials.hit1+trials.hit4, trials.miss1+trials.miss4; (trials.hit4+trials.hit1+trials.miss4+trials.miss1)/2, (trials.hit4+trials.hit1+trials.miss4+trials.miss1)/2]);
fprintf(1,['Chi squared testing hit/miss vs random for both lanes ' num2str(p_hit_miss_both_lanes) '\n\n'])

% nan_mask=logical(ones(1,no_time_bins));
% for ii_neuron=1:no_neurons
%     this_nan_mask=ones(1,no_time_bins);
%     this_nan_mask(1,:)=~isnan(neural_data(:,ii_neuron));
%     nan_mask=nan_mask&this_nan_mask;
% end
% neural_data=neural_data(nan_mask,:);
% pos_binned=pos_binned(nan_mask,:);
% no_time_bins=sum(nan_mask);


no_time_bins=size(neural_data,1);

%Do z scores
mean_neural_data_col=mean(neural_data,1);
mean_neural_data=repmat(mean_neural_data_col,no_time_bins,1);

std_neural_data_col=std(neural_data,1);
std_neural_data=repmat(std_neural_data_col,no_time_bins,1);

neural_data=(neural_data-mean_neural_data)./std_neural_data;

pffft=1;

% % Format for Wiener Filter, Wiener Cascade, XGBoost, and Dense Neural Network
% % Put in "flat" format, so each "neuron / time" is a single feature
% % i.e. each time point in the before and after window becomes a different
% % "neuron"
% all_bins_per_window=bins_before+bins_after+bins_current;
% no_X_dFF_neurons=no_neurons*all_bins_per_window;
% X_dFF=zeros(no_time_bins,no_X_dFF_neurons);
% 
% for ii_t=bins_before+1:no_time_bins-bins_after
%     ii_n=0;
%     for no_win=1:all_bins_per_window
%         ii_this_t=ii_t-bins_before+no_win-1;
%         X_dFF(ii_t,ii_n+1:ii_n+no_neurons)=neural_data(ii_this_t,:);
%         ii_n=ii_n+no_neurons;
%     end
% end

%Now do the neural networks
tic
  

%Train with per trial data using a leave one out approach
per_ROI=[];
for ii_repeat=1:handles_choices.repeats
    for this_ROI=1:no_neurons
        %training_range_template has all the trials
        training_range_template=zeros(1,no_time_bins);
        lane_template=zeros(1,no_time_bins);
        odor_plume_template=zeros(1,no_time_bins);
        for trNo=1:trials.odor_trNo
            dFF_pred(trNo).dFFpred_xy=[];
            % dFF_pred(trNo).dFFpred_xyl=[];
            % dFF_pred(trNo).dFFpred_xyop=[];
            dFF_pred(trNo).dFFpred_op=[];
            dFF_pred(trNo).dFF=[];
            dFF_pred(trNo).dFFpred_xyop=[];
            
            % x_pred(trNo).data=[];

            x_predictedstart=trials.odor_ii_start(trNo)-10;
            x_predictedend=trials.odor_ii_end(trNo)+15;
            training_range_template(x_predictedstart:x_predictedend)=1;
            if trials.odor_lane(trNo)==1
                lane_template(x_predictedstart:x_predictedend)=1;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned(x_ii,1);
                    this_y=pos_binned(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane1(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            else
                lane_template(x_predictedstart:x_predictedend)=4;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned(x_ii,1);
                    this_y=pos_binned(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane4(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            end
        end

        % parfor trNo=1:trials.odor_trNo
        start_toc=toc;
        gcp;
        parfor trNo=1:trials.odor_trNo
        % for trNo=1:trials.odor_trNo

            this_test_range=zeros(1,no_time_bins);
            if trNo==1
                ii_test_range_start=1;
            else
                ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
            end

            if trNo==trials.odor_trNo
                ii_test_range_end=no_time_bins;
            else
                ii_test_range_end=trials.odor_ii_end(trNo)+15;
            end

            this_test_range(ii_test_range_start:ii_test_range_end)=1;
            this_training_range=logical(training_range_template)&(~logical(this_test_range));

            this_bin_dFFtrain=zeros(sum(this_training_range),1);
            this_bin_dFFtrain=neural_data_binary(this_training_range,this_ROI);
            
            % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
            this_dFFtest=zeros(sum(this_test_range),1);
            this_dFFtest(:,1)=neural_data(logical(this_test_range),this_ROI);
            dFF_pred(trNo).dFF=this_dFFtest;

            this_bin_dFFtest=zeros(sum(this_test_range),1);
            this_bin_dFFtest(:,1)=neural_data_binary(logical(this_test_range),this_ROI);
            dFF_pred(trNo).bindFF=this_bin_dFFtest;
            

            XYtrain=pos_binned(this_training_range,:);
            % Yvalid=pos_binned(ii_valid_range(1):ii_valid_range(2),:);
            XYtest=pos_binned(logical(this_test_range),:);

            lane_train=lane_template(this_training_range);
            lane_test=lane_template(logical(this_test_range));

            odor_plume_train=odor_plume_template(this_training_range);
            odor_plume_test=odor_plume_template(logical(this_test_range));

            %Decode using neural network

            %Decode using x and y only
            switch handles_choices.algo
                case 1
                    Mdl_dFFxy = fitcnet(XYtrain,this_bin_dFFtrain);
                case 2
                    Mdl_dFFxy = fitglm(XYtrain,this_bin_dFFtrain);
                case 3
                    Mdl_dFFxy = fitctree(XYtrain,this_bin_dFFtrain);
                case 4
                    Mdl_dFFxy = fitcsvm(XYtrain,this_bin_dFFtrain, 'KernelFunction','rbf','Standardize', true);
                % case 5
                %     Mdl_dFFxy = fitrgp(XYtrain,this_bin_dFFtrain);
                % case 6
                %     Mdl_dFFxy = fitrgp(XYtrain,this_bin_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
            end
            dFF_pred(trNo).bindFFpred_xy=predict(Mdl_dFFxy,XYtest);


            %Decode using odor plume only
            switch handles_choices.algo
                case 1
                    Mdl_dFFop = fitcnet([odor_plume_train'],this_bin_dFFtrain);
                case 2
                    Mdl_dFFop = fitglm([odor_plume_train'],this_bin_dFFtrain);
                case 3
                    Mdl_dFFop = fitctree([odor_plume_train'],this_bin_dFFtrain);
                case 4
                    Mdl_dFFop = fitcsvm([odor_plume_train'],this_bin_dFFtrain, 'KernelFunction','rbf','Standardize', true);
                % case 5
                %     Mdl_dFFop = fitrgp([odor_plume_train'],this_dFFtrain);
                % case 6
                %     Mdl_dFFop = fitrgp([odor_plume_train'],this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
            end
            dFF_pred(trNo).bindFFpred_op=predict(Mdl_dFFop,[odor_plume_test']);
 
            %Decode using both xy and odor plume only
            switch handles_choices.algo
                case 1
                    Mdl_dFFxyop = fitcnet([XYtrain odor_plume_train'],this_bin_dFFtrain);
                case 2
                    Mdl_dFFxyop = fitglm([XYtrain odor_plume_train'],this_bin_dFFtrain);
                case 3
                    Mdl_dFFxyop = fitctree([XYtrain odor_plume_train'],this_bin_dFFtrain);
                case 4
                    Mdl_dFFxyop = fitcsvm([XYtrain odor_plume_train'],this_bin_dFFtrain, 'KernelFunction','rbf','Standardize', true);
                % case 5
                %     Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train'],this_dFFtrain);
                % case 6
                %     Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train'],this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
            end
            dFF_pred(trNo).bindFFpred_xyop=predict(Mdl_dFFxyop,[XYtest odor_plume_test']);

            pffft=1;
            % y_pred(trNo).data=predict(MdlY2,XdFFtest);
        end
        fprintf(1,['Elapsed time ' num2str((toc-start_toc)/(60)) ' mins\n\n'])

        %Parse out the parfor loop output
        all_bindFFpred_xy=[];
        all_bindFFpred_xyop=[];
        all_bindFFpred_op=[];

        all_bindFFpredl1_xy=[];
        all_bindFFpredl1_xyop=[];
        all_bindFFpredl1_op=[];

        all_bindFFpredl4_xy=[];
        all_bindFFpredl4_xyop=[];
        all_bindFFpredl4_op=[];

        all_dFF=[];
        all_dFFl1=[];
        all_dFFl4=[];

        all_bindFF=[];
        all_bindFFl1=[];
        all_bindFFl4=[];

        for trNo=1:trials.odor_trNo
            per_ROI(this_ROI).repeats(ii_repeat).trial(trNo).bindFFpred_xy=dFF_pred(trNo).bindFFpred_xy;
            per_ROI(this_ROI).repeats(ii_repeat).trial(trNo).bindFFpred_xyop=dFF_pred(trNo).bindFFpred_xyop;
            per_ROI(this_ROI).repeats(ii_repeat).trial(trNo).bindFFpred_op=dFF_pred(trNo).bindFFpred_op;
            per_ROI(this_ROI).repeats(ii_repeat).trial(trNo).dFF=dFF_pred(trNo).dFF;

            all_bindFFpred_xy=[all_bindFFpred_xy dFF_pred(trNo).bindFFpred_xy'];
            all_bindFFpred_xyop=[all_bindFFpred_xyop dFF_pred(trNo).bindFFpred_xyop'];
            all_bindFFpred_op=[all_bindFFpred_op dFF_pred(trNo).bindFFpred_op'];
            all_dFF=[all_dFF dFF_pred(trNo).dFF'];
            all_bindFF=[all_bindFF dFF_pred(trNo).bindFF'];

            if trials.odor_lane(trNo)==1
                all_bindFFpredl1_xy=[all_bindFFpredl1_xy dFF_pred(trNo).bindFFpred_xy'];
                all_bindFFpredl1_xyop=[all_bindFFpredl1_xyop dFF_pred(trNo).bindFFpred_xyop'];
                all_bindFFpredl1_op=[all_bindFFpredl1_op dFF_pred(trNo).bindFFpred_op'];
                all_dFFl1=[all_dFFl1 dFF_pred(trNo).dFF'];
                all_bindFFl1=[all_bindFFl1 dFF_pred(trNo).bindFF'];
            else
                all_bindFFpredl4_xy=[all_bindFFpredl4_xy dFF_pred(trNo).bindFFpred_xy'];
                all_bindFFpredl4_xyop=[all_bindFFpredl4_xyop dFF_pred(trNo).bindFFpred_xyop'];
                all_bindFFpredl4_op=[all_bindFFpredl4_op dFF_pred(trNo).bindFFpred_op'];
                all_dFFl4=[all_dFFl4 dFF_pred(trNo).dFF'];
                all_bindFFl4=[all_bindFFl4 dFF_pred(trNo).bindFF'];
            end
        end

        per_ROI(this_ROI).repeats(ii_repeat).all_dFF=all_dFF;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFF=all_bindFF;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpred_xy=all_bindFFpred_xy;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpred_xyop=all_bindFFpred_xyop;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpred_op=all_bindFFpred_op;

        per_ROI(this_ROI).repeats(ii_repeat).all_dFFl1=all_dFFl1;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFl1=all_bindFFl1;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl1_xy=all_bindFFpredl1_xy;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl1_xyop=all_bindFFpredl1_xyop;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl1_op=all_bindFFpredl1_op;

        per_ROI(this_ROI).repeats(ii_repeat).all_dFFl4=all_dFFl4;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFl4=all_bindFFl4;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl4_xy=all_bindFFpredl4_xy;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl4_xyop=all_bindFFpredl4_xyop;
        per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl4_op=all_bindFFpredl4_op;

        thisRxy=corrcoef([all_bindFF' all_bindFFpred_xy']);
        per_ROI(this_ROI).repeats(ii_repeat).Rxy=thisRxy(1,2);
        % thisRxyl=corrcoef([all_dFF' all_dFFpred_xyl']);
        % per_ROI(this_ROI).repeats(ii_repeat).Rxyl=thisRxyl(1,2);
        thisRxyop=corrcoef([all_bindFF' all_bindFFpred_xyop']);
        per_ROI(this_ROI).repeats(ii_repeat).Rxyop=thisRxyop(1,2);
        thisRop=corrcoef([all_bindFF' all_bindFFpred_op']);
        per_ROI(this_ROI).repeats(ii_repeat).Rop=thisRop(1,2);
        % fprintf(1, 'Repeat %d ROI No %d, rho for prediction dFFxy, dFFxyl, dFFxyop, dFFop %d %d %d %d\n\n'...
        %     ,ii_repeat,this_ROI, per_ROI(this_ROI).repeats(ii_repeat).Rxy,per_ROI(this_ROI).repeats(ii_repeat).Rxyl,per_ROI(this_ROI).repeats(ii_repeat).Rxyop,per_ROI(this_ROI).repeats(ii_repeat).Rop);
        fprintf(1, 'Repeat %d ROI No %d, rho for prediction dFFxy, dFFop dFFxyop%d %d\n\n'...
            ,ii_repeat,this_ROI, per_ROI(this_ROI).repeats(ii_repeat).Rxy,per_ROI(this_ROI).repeats(ii_repeat).Rop...
            ,per_ROI(this_ROI).repeats(ii_repeat).Rxyop);


    end
end

%Train with per trial with shuffled data using a leave one out approach
per_ROI_sh=[];
for ii_repeat=1:handles_choices.sh_repeats
    for this_ROI=1:no_neurons
        %training_range_template has all the trials
        training_range_template=zeros(1,no_time_bins);
        lane_template=zeros(1,no_time_bins);
        odor_plume_template=zeros(1,no_time_bins);
        for trNo=1:trials.odor_trNo
            dFF_pred(trNo).dFFpred_xy=[];
            % dFF_pred(trNo).dFFpred_xyl=[];
            dFF_pred(trNo).dFF=[];
            % x_pred(trNo).data=[];

            x_predictedstart=trials.odor_ii_start(trNo)-10;
            x_predictedend=trials.odor_ii_end(trNo)+15;
            training_range_template(x_predictedstart:x_predictedend)=1;
            if trials.odor_lane(trNo)==1
                lane_template(x_predictedstart:x_predictedend)=1;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned(x_ii,1);
                    this_y=pos_binned(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane1(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            else
                lane_template(x_predictedstart:x_predictedend)=4;
                for x_ii=x_predictedstart:x_predictedend
                    this_x=pos_binned(x_ii,1);
                    this_y=pos_binned(x_ii,1);
                    [minabx,ii_minax]=min(abs(x_for_plume-this_x));
                    [minaby,ii_minay]=min(abs(y_for_plume-this_y));
                    this_ca=mean_plume_lane4(ii_minay,ii_minax);
                    odor_plume_template(x_ii)=this_ca;
                end
            end
        end

        % parfor trNo=1:trials.odor_trNo
        start_toc=toc;
        max_overlap=handles_choices.max_overlap;
        gcp;
        parfor trNo=1:trials.odor_trNo

            this_test_range=zeros(1,no_time_bins);
            if trNo==1
                ii_test_range_start=1;
            else
                ii_test_range_start=trials.odor_ii_end(trNo-1)+15;
            end

            if trNo==trials.odor_trNo
                ii_test_range_end=no_time_bins;
            else
                ii_test_range_end=trials.odor_ii_end(trNo)+15;
            end

            this_test_range(ii_test_range_start:ii_test_range_end)=1;
            this_training_range=logical(training_range_template)&(~logical(this_test_range));

            % % ii_valid_range=ceil(valid_range*no_time_bins);
            % %     ii_test_range=ceil(test_range*no_time_bins);
            % this_dFFtrain=zeros(sum(this_training_range),1);
            % this_dFFtrain=neural_data(this_training_range,this_ROI);
            % % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
            % this_dFFtest=zeros(sum(this_test_range),1);
            % this_dFFtest(:,1)=neural_data(logical(this_test_range),this_ROI);
            % dFF_pred(trNo).dFF=this_dFFtest;

            this_bin_dFFtrain=zeros(sum(this_training_range),1);
            this_bin_dFFtrain=neural_data_binary(this_training_range,this_ROI);
            
            % Xvalid=X_dFF(ii_valid_range(1):ii_valid_range(2),:);
            this_dFFtest=zeros(sum(this_test_range),1);
            this_dFFtest(:,1)=neural_data(logical(this_test_range),this_ROI);
            dFF_pred(trNo).dFF=this_dFFtest;

            this_bin_dFFtest=zeros(sum(this_test_range),1);
            this_bin_dFFtest(:,1)=neural_data_binary(logical(this_test_range),this_ROI);
            dFF_pred(trNo).bindFF=this_bin_dFFtest;

            %Now shuffle
            XYtrain_pre=pos_binned(this_training_range,:);
            shuffle_length=floor(size(XYtrain_pre,1)/(trials.odor_trNo*2));

            this_overlap=10;
            while this_overlap>max_overlap
                perm_indices=randperm(trials.odor_trNo*2);
                this_overlap=0;
                for ii_p=1:trials.odor_trNo*2
                    if perm_indices(ii_p)==ii_p
                        this_overlap=this_overlap+1;
                    end
                end
            end

            % Yvalid=pos_binned(ii_valid_range(1):ii_valid_range(2),:);
            XYtest=pos_binned(logical(this_test_range),:);

            lane_train_pre=lane_template(this_training_range);
            lane_train_pre=lane_train_pre';
            lane_test=lane_template(logical(this_test_range));

            odor_plume_train_pre=odor_plume_template(this_training_range);
            odor_plume_train_pre=odor_plume_train_pre';
            odor_plume_test=odor_plume_template(logical(this_test_range));

            XYtrain=zeros(size(XYtrain_pre,1),size(XYtrain_pre,2));
            lane_train=zeros(size(XYtrain_pre,1),1);
            odor_plume_train=zeros(size(XYtrain_pre,1),1);
            for sh_ii=1:trials.odor_trNo*2
                XYtrain((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
                    XYtrain_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
                lane_train((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
                    lane_train_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
                odor_plume_train((sh_ii-1)*shuffle_length+1:sh_ii*shuffle_length,:)=...
                    odor_plume_train_pre((perm_indices(sh_ii)-1)*shuffle_length+1:perm_indices(sh_ii)*shuffle_length,:);
            end

            %Make the last few equal to the first few
            if trials.odor_trNo*2*shuffle_length<size(XYtrain,1)
                delta_ii=size(XYtrain,1)-trials.odor_trNo*2*shuffle_length;
                XYtrain(trials.odor_trNo*2*shuffle_length+1:end,:)=...
                    XYtrain_pre(1:delta_ii,:);
                lane_train(trials.odor_trNo*2*shuffle_length+1:end,:)=...
                    lane_train_pre(1:delta_ii,:);
                odor_plume_train(trials.odor_trNo*2*shuffle_length+1:end,:)=...
                    odor_plume_train_pre(1:delta_ii,:);
            end


            %Decode using neural network

            %Decode using x and y only
            switch handles_choices.algo
                case 1
                    Mdl_dFFxy = fitcnet(XYtrain,this_bin_dFFtrain');
                case 2
                    Mdl_dFFxy = fitglm(XYtrain,this_bin_dFFtrain');
                case 3
                    Mdl_dFFxy = fitctree(XYtrain,this_bin_dFFtrain');
                case 4
                    Mdl_dFFxy = fitcsvm(XYtrain,this_bin_dFFtrain', 'KernelFunction','rbf','Standardize', true);
                % case 5
                %     Mdl_dFFxy = fitrgp(XYtrain,this_dFFtrain');
                % case 6
                %     Mdl_dFFxy = fitrgp(XYtrain,this_dFFtrain','KernelFunction', 'squaredexponential', 'Standardize', true);
            end
            dFF_pred(trNo).bindFFpred_xy=predict(Mdl_dFFxy,XYtest);



            %Decode using odor plume
            switch handles_choices.algo
                case 1
                    Mdl_dFFop = fitcnet([odor_plume_train],this_bin_dFFtrain);
                case 2
                    Mdl_dFFop = fitglm([odor_plume_train],this_bin_dFFtrain);
                case 3
                    Mdl_dFFop = fitctree([odor_plume_train],this_bin_dFFtrain);
                case 4
                    Mdl_dFFop = fitcsvm([odor_plume_train],this_bin_dFFtrain, 'KernelFunction','rbf','Standardize', true);
                % case 5
                %     Mdl_dFFop = fitrgp([odor_plume_train],this_dFFtrain);
                % case 6
                %     Mdl_dFFop = fitrgp([odor_plume_train],this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
            end
            dFF_pred(trNo).bindFFpred_op=predict(Mdl_dFFop,[odor_plume_test']);

            %Decode using xy and odor plume
            switch handles_choices.algo
                case 1
                    Mdl_dFFxyop = fitcnet([XYtrain odor_plume_train],this_bin_dFFtrain);
                case 2
                    Mdl_dFFxyop = fitglm([XYtrain odor_plume_train],this_bin_dFFtrain);
                case 3
                    Mdl_dFFxyop = fitctree([XYtrain odor_plume_train],this_bin_dFFtrain);
                case 4
                    Mdl_dFFxyop = fitcsvm([XYtrain odor_plume_train],this_bin_dFFtrain, 'KernelFunction','rbf','Standardize', true);
                % case 5
                %     Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train],this_dFFtrain);
                % case 6
                %     Mdl_dFFxyop = fitrgp([XYtrain odor_plume_train],this_dFFtrain,'KernelFunction', 'squaredexponential', 'Standardize', true);
            end
            dFF_pred(trNo).bindFFpred_xyop=predict(Mdl_dFFxyop,[XYtest odor_plume_test']);



            pffft=1;
            % y_pred(trNo).data=predict(MdlY2,XdFFtest);
        end
        fprintf(1,['Elapsed time ' num2str((toc-start_toc)/(60)) ' mins\n\n'])

        %Parse out the parfor loop output
        all_bindFFpred_xy=[];
        all_bindFFpred_xyl=[];
        all_bindFFpred_xyop=[];
        all_bindFFpred_op=[];
        all_bindFF=[];
        for trNo=1:trials.odor_trNo
            per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).bindFFpred_xy=dFF_pred(trNo).bindFFpred_xy;
            all_bindFFpred_xy=[all_bindFFpred_xy dFF_pred(trNo).bindFFpred_xy'];
            % per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).dFFpred_xyl=dFF_pred(trNo).dFFpred_xyl;
            % all_dFFpred_xyl=[all_dFFpred_xyl dFF_pred(trNo).dFFpred_xyl'];
            all_bindFFpred_xyop=[all_bindFFpred_xyop dFF_pred(trNo).bindFFpred_xyop'];
            all_bindFFpred_op=[all_bindFFpred_op dFF_pred(trNo).bindFFpred_op'];
            per_ROI_sh(this_ROI).repeats(ii_repeat).trial(trNo).bindFF=dFF_pred(trNo).bindFF;
            all_bindFF=[all_bindFF dFF_pred(trNo).bindFF'];
        end

        per_ROI_sh(this_ROI).repeats(ii_repeat).all_bindFF=all_bindFF;
        per_ROI_sh(this_ROI).repeats(ii_repeat).all_bindFFpred_xy=all_bindFFpred_xy;
        % per_ROI_sh(this_ROI).repeats(ii_repeat).all_dFFpred_xyl=all_dFFpred_xyl;
        per_ROI_sh(this_ROI).repeats(ii_repeat).all_bindFFpred_xyop=all_bindFFpred_xyop;
        per_ROI_sh(this_ROI).repeats(ii_repeat).all_bindFFpred_op=all_bindFFpred_op;

        thisRxy=corrcoef([all_bindFF' all_bindFFpred_xy']);
        per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy=thisRxy(1,2);
        % thisRxyl=corrcoef([all_dFF' all_dFFpred_xyl']);
        % per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyl=thisRxyl(1,2);
        thisRxyop=corrcoef([all_bindFF' all_bindFFpred_xyop']);
        per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop=thisRxyop(1,2);
        thisRop=corrcoef([all_bindFF' all_bindFFpred_op']);
        per_ROI_sh(this_ROI).repeats(ii_repeat).Rop=thisRop(1,2);
        % fprintf(1, 'Repeat %d ROI No %d, rho for shuffled prediction dFFxy, dFFxyl dFFxyop, dFFop %d %d %d %d\n\n',ii_repeat,...
        %     this_ROI,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyl,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop,per_ROI_sh(this_ROI).repeats(ii_repeat).Rop);
                fprintf(1, 'Repeat %d ROI No %d, rho for shuffled prediction dFFxy, dFFop, dFFxyop %d %d \n\n',ii_repeat,...
            this_ROI,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxy,per_ROI_sh(this_ROI).repeats(ii_repeat).Rop...
            ,per_ROI_sh(this_ROI).repeats(ii_repeat).Rxyop);

 
    end
end

fprintf(1,['Elapsed time ' num2str(toc/(60*60)) ' hrs\n\n'])

handles_out.per_ROI=per_ROI;
handles_out.per_ROI_sh=per_ROI_sh;
handles_out.trials=trials;

if handles_choices.displayFigures==1
    % figNo=figNo+1;
    % try
    %     close(figNo)
    % catch
    % en
    %Calculate 95th percentile
    % figNo=0;
    noROIs=length(handles_out.per_ROI);
    if isfield(handles_out.per_ROI(1),'sh_repeats')
        noshRepeats=length(handles_out.per_ROI(1).sh_repeats);
    else
        noshRepeats=length(handles_out.per_ROI(1).repeats);
    end

    shRops=[];
    shRxys=[];
    shRxyops=[];
    for ii_ROI=1:noROIs
        for ii_sh_repeats=1:noshRepeats
            shRops=[shRops handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rop];
            shRxys=[shRxys handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rxy];
            shRxyops=[shRxyops handles_out.per_ROI_sh(ii_ROI).repeats(ii_sh_repeats).Rxyop];
        end
    end

    % Bootstrap resampling op
    bootstrap_samples_op = bootstrp(1000, @prctile, shRops, 95);

    % Compute the 95th percentile of the bootstrap samples
    Ropsh_bootstrap_95th_percentile = prctile(bootstrap_samples_op, 95);

    % Bootstrap resampling xy
    bootstrap_samples_xy = bootstrp(1000, @prctile, shRxys, 95);

    % Compute the 95th percentile of the bootstrap samples
    Rxysh_bootstrap_95th_percentile = prctile(bootstrap_samples_xy, 95);

        % Bootstrap resampling xyop
    bootstrap_samples_xyop = bootstrp(1000, @prctile, shRxyops, 95);

    % Compute the 95th percentile of the bootstrap samples
    Rxyopsh_bootstrap_95th_percentile = prctile(bootstrap_samples_xyop, 95);


    %Find the Rops that are above 95th percentile
    noRepeats=length(handles_out.per_ROI(1).repeats);

    Rops=[];
    Rxys=[];
    Rxyops=[];
    iiROI_Rops=[];
    iiROI_Rxys=[];
    iiROI_Rxyops=[];
    fileNo_Rops=[];
    fileNo_Rxys=[];
    fileNo_Rxyops=[];
    Rops_above_95=[];
    Rxys_above_95=[];
    Rxyops_above_95=[];
    Rops_for_Rxys_above_95=[];
    Rxys_for_Rops_above_95=[];
    Rxys_for_Rxyops_above_95=[];
    Rxyops_for_Rxys_above_95=[];
    for ii_ROI=1:noROIs
        these_Rops=[];
        these_Rxys=[];
        these_Rxyops=[];
        for ii_repeats=1:noRepeats
            these_Rops=[these_Rops handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rop];
            these_Rxys=[these_Rxys handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rxy];
            these_Rxyops=[these_Rxyops handles_out.per_ROI(ii_ROI).repeats(ii_repeats).Rxyop];
        end
        Rops=[Rops mean(these_Rops)];
        Rxys=[Rxys mean(these_Rxys)];
        Rxyops=[Rxyops mean(these_Rxyops)];
        if mean(these_Rops)>Ropsh_bootstrap_95th_percentile
            iiROI_Rops=[iiROI_Rops ii_ROI];
            % fileNo_Rops=[fileNo_Rops fileNo];
            Rops_above_95=[Rops_above_95 mean(these_Rops)];
            Rxys_for_Rops_above_95=[Rxys_for_Rops_above_95 mean(these_Rxys)];
        end
        if mean(these_Rxys)>Rxysh_bootstrap_95th_percentile
            iiROI_Rxys=[iiROI_Rxys ii_ROI];
            % fileNo_Rxys=[fileNo_Rxys fileNo];
            Rxys_above_95=[Rxys_above_95 mean(these_Rxys)];
            Rops_for_Rxys_above_95=[Rops_for_Rxys_above_95 mean(these_Rops)];
            Rxyops_for_Rxys_above_95=[Rxyops_for_Rxys_above_95 mean(these_Rxyops)];
        end
        if mean(these_Rxyops)>Rxyopsh_bootstrap_95th_percentile
            iiROI_Rxyops=[iiROI_Rxyops ii_ROI];
            % fileNo_Rxys=[fileNo_Rxys fileNo];
            Rxyops_above_95=[Rxyops_above_95 mean(these_Rxyops)];
            Rxys_for_Rxyops_above_95=[Rxys_for_Rxyops_above_95 mean(these_Rxys)];
        end

    end


    %Display Rho op vs Rho xy figure for this ii_wf and fileNo
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    hold on

    %Plot Rops above 95
    for ii_rho1=1:length(Rops_above_95)
        match_found=0;
        for ii_rho2=1:length(Rxys_above_95)
            if (iiROI_Rxys(ii_rho2)==iiROI_Rops(ii_rho1))
                match_found=1;
            end
        end
        if match_found==1
            plot(Rxys_for_Rops_above_95(ii_rho1),Rops_above_95(ii_rho1),'ok')
        else
            plot(Rxys_for_Rops_above_95(ii_rho1),Rops_above_95(ii_rho1),'ob')
        end
    end
 
    %Plot Rxys above 95
    for ii_rho1=1:length(Rxys_above_95)
        match_found=0;
        for ii_rho2=1:length(Rops_above_95)
            if (iiROI_Rops(ii_rho2)==iiROI_Rxys(ii_rho1))
                match_found=1;
            end
        end
        if match_found==1
            plot(Rxys_above_95(ii_rho1),Rops_for_Rxys_above_95(ii_rho1),'ok')
        else
            plot(Rxys_above_95(ii_rho1),Rops_for_Rxys_above_95(ii_rho1),'or')
        end
    end

    this_xlim=xlim;
    this_ylim=ylim;
    text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.15*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rop>95%','Color',[0 0 1])
    text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.10*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy>95%','Color',[1 0 0])
    text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.05*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy and Rop>95%','Color',[0 0 0])


    plot([this_xlim(1) this_xlim(2)],[this_xlim(1) this_xlim(2)],'-k')

    xlabel('Rho xy')
    ylabel('Rho op')
    title(['Rho op vs. Rho xy log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha) ' op= ' num2str(handles_choices.which_odor_plume)])

  %Display Rho xyop vs Rho xy figure for this ii_wf and fileNo
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


    hold on
 
    %Plot Rxyops above 95
    for ii_rho1=1:length(Rxyops_above_95)
        match_found=0;
        for ii_rho2=1:length(Rxys_above_95)
            if (iiROI_Rxys(ii_rho2)==iiROI_Rxyops(ii_rho1))
                match_found=1;
            end
        end
        if match_found==1
            plot(Rxys_for_Rxyops_above_95(ii_rho1),Rxyops_above_95(ii_rho1),'ok')
        else
            plot(Rxys_for_Rxyops_above_95(ii_rho1),Rxyops_above_95(ii_rho1),'ob')
        end
    end
 
    %Plot Rxys above 95
    for ii_rho1=1:length(Rxys_above_95)
        match_found=0;
        for ii_rho2=1:length(Rops_above_95)
            if (iiROI_Rxyops(ii_rho2)==iiROI_Rxys(ii_rho1))
                match_found=1;
            end
        end
        if match_found==1
            plot(Rxys_above_95(ii_rho1),Rxyops_for_Rxys_above_95(ii_rho1),'ok')
        else
            plot(Rxys_above_95(ii_rho1),Rxyops_for_Rxys_above_95(ii_rho1),'or')
        end
    end

    this_xlim=xlim;
    this_ylim=ylim;
    text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.15*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxyop>95%','Color',[0 0 1])
    text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.10*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy>95%','Color',[1 0 0])
    text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.05*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy and Rxyop>95%','Color',[0 0 0])
  
    plot([this_xlim(1) this_xlim(2)],[this_xlim(1) this_xlim(2)],'-k')

    xlabel('Rho xy')
    ylabel('Rho xyop')
    title(['Rho xyop vs. Rho xy log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha) ' op= ' num2str(handles_choices.which_odor_plume)])


end
    

save([this_path arena_file(1:end-4) 'decdFFa' num2str(handles_choices.algo) ...
    'wf' num2str(handles_choices.weber_fechner) 'g' num2str(handles_choices.group) ...
    'a' num2str(handles_choices.alpha) 'op' num2str(handles_choices.which_odor_plume)...
    '.mat'],'handles_out','handles_choices','-v7.3')
 
%Plot the traces and predicted traces
%Rxys above 95 percentile
to_sort=[Rxys_above_95' iiROI_Rxys'];
sorted_rows=sortrows(to_sort,1);
iiROI_Rxys_to_plot=zeros(1,length(iiROI_Rxys));
iiROI_Rxys_to_plot(1,:)=sorted_rows(:,2);
% 
% per_ROI(this_ROI).repeats(ii_repeat).all_dFF=all_dFF;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpred_xy=all_dFFpred_xy;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpred_xyop=all_dFFpred_xyop;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpred_op=all_dFFpred_op;
% 
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFl1=all_dFFl1;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpredl1_xy=all_dFFpredl1_xy;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpredl1_xyop=all_dFFpredl1_xyop;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpredl1_op=all_dFFpredl1_op;
% 
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFl4=all_dFFl4;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpredl4_xy=all_dFFpredl4_xy;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpredl4_xyop=all_dFFpredl4_xyop;
% per_ROI(this_ROI).repeats(ii_repeat).all_dFFpredl4_op=all_dFFpredl4_op;

%Note: I will only show the results for repeat 1
these_alldFF=[];
ii_repeat=1;
for ii=1:length(iiROI_Rxys)
    these_alldFF=[these_alldFF per_ROI(iiROI_Rxys(ii)).repeats(ii_repeat).all_dFF];
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on

% Determine the y spacing of the traces
y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));


for iiROI=1:length(iiROI_Rxys)
    this_ROI=iiROI_Rxys_to_plot(iiROI);
    no_time_bins1=length(per_ROI(this_ROI).repeats(ii_repeat).all_dFFl1);
    no_time_bins4=length(per_ROI(this_ROI).repeats(ii_repeat).all_dFFl4);

    plot([1:no_time_bins1], per_ROI(this_ROI).repeats(ii_repeat).all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).repeats(ii_repeat).all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)

    plot([1:no_time_bins1], per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl1_xy*0.8*y_shift+y_shift*iiROI-0.05*y_shift,'-r','LineWidth',1)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl4_xy*0.8*y_shift+y_shift*iiROI-0.05*y_shift,'-r','LineWidth',1)
end

%Show the last few
ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])

xlabel('time(sec)')
title(['All dFF timecourses for ROIs with xy predictions above 95 percentile'])

to_sort=[Rxyops_above_95' iiROI_Rxyops'];
sorted_rows=sortrows(to_sort,1);
iiROI_Rxyops_to_plot=zeros(1,length(iiROI_Rxyops));
iiROI_Rxyops_to_plot(1,:)=sorted_rows(:,2);

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on

% Determine the y spacing of the traces
y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));


for iiROI=1:length(iiROI_Rxyops_to_plot)
    this_ROI=iiROI_Rxyops_to_plot(iiROI);
    no_time_bins1=length(per_ROI(this_ROI).repeats(ii_repeat).all_dFFl1);
    no_time_bins4=length(per_ROI(this_ROI).repeats(ii_repeat).all_dFFl4);

    plot([1:no_time_bins1], per_ROI(this_ROI).repeats(ii_repeat).all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).repeats(ii_repeat).all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)

    plot([1:no_time_bins1], per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl1_xyop*0.8*y_shift+y_shift*iiROI-0.05*y_shift,'-r','LineWidth',1)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).repeats(ii_repeat).all_bindFFpredl4_xyop*0.8*y_shift+y_shift*iiROI-0.05*y_shift,'-r','LineWidth',1)
end

%Show the last few
ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])

xlabel('time(sec)')
title(['All dFF timecourses for ROIs with xyop predictions above 95 percentile'])

%Now plot pseudocolor maps of activity


pffft=1;