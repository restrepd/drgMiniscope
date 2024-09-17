function handles_out=drgMini_dFFcrosscorr(handles_choices)
%Does decoding following Glaser et al, 2020 https://doi.org/10.1523/ENEURO.0506-19.2020

warning('off')

if exist('handles_choices')==0
    clear all
    close all

     %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=0.2;
    handles_choices.dt=dt;
    handles_choices.cross_shift_dt=10;
    dt_miniscope=1/30;
    handles_choices.dt_miniscope=dt_miniscope;
    %Note that n_shuffle is changed to a maximum of ii_n_training
    n_shuffle=5;
    handles_choices.n_shuffle=n_shuffle;
    handles_choices.normalize=0; %Normalize the calcium transients before doing the analysis
    handles_choices.no_to_display=20;
    no_to_display=handles_choices.no_to_display;
    handles_choices.trim_neural_data=0;

    handles_choices.distances_in_mm=1; %If the distance is already in mm there will not be morphing of space

    %Files for dFF data
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    dFF_file='20220804_FCM22_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';

    if handles_choices.distances_in_mm==0
        arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync.mat';
    else
        arena_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm.mat';
    end


%     %Second troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
    % dFF_file='20220713_FCM6_withodor_miniscope_sync_L1andL4_ncorre_ext.mat';
    % 
    % if handles_choices.distances_in_mm==0
    %     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';
    % else
    %     arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm.mat'; %Distances converted to mm
    % end
    
    %No odor troubleshooting files
%     this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/';
%     dFF_file='20220824_FCM6_withoutodor_miniscope_sync_L1andL4_ncorre_ext.mat';
%     arena_file='20220824_FCM6withoutodor_odorarena_L1andL4_sync.mat';


    
    handles_choices.this_path=this_path;
    handles_choices.dFF_file=dFF_file;
    handles_choices.arena_file=arena_file;

    handles_choices.which_odor_plume=2;   %1=use odor plume from the publication, 2=use odor plume simulated by Aaron
    handles_choices.cm_from_floor=2;

    handles_choices.weber_fechner=0; 
    handles_choices.alpha=1;
    handles_choices.multiplier=1;
    %0 is Stevens Law, R proportional to C^alpha
    %1 is Weber-Flechner law R proportional to log(C)
    %See Copelli et al DOI: 10.1103/PhysRevE.65.060901

    handles_choices.algo=1; 
    %1 ann
    %2 glm
    %3 tree
    %4 svm with radial basis
    %5 Gaussian process regression
    %6 GPR with 'KernelFunction', 'squaredexponential', 'Standardize', true
    %7 ann optimize hyper parameters, this is VERY slow!


    handles_choices.group=1; %1=odor plume, 2=spatial
    
  
    handles_choices.save_tag='decdFF'; %This will be used in the name of the save file

    
    handles_choices.repeats=1;
    handles_choices.sh_repeats=handles_choices.repeats;
    handles_choices.best_repeat=1; %If multiple runs are performed choose the best result
    
    handles_choices.max_overlap=3; %Maximum overlap of the shuffled segments
    handles_choices.displayFigures=1;

   


    handles_choices.dt_decoding_op=2; %All bins from -dt_decoding to 0 will be included in the decoder
    %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
    fprintf(1,['Number of time bins for op decoding ' num2str(no_dec_time_bins_op) '\n\n'])


    handles_choices.dt_decoding_xy=0.6;
    no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
    fprintf(1,['Number of time bins for xy decoding ' num2str(no_dec_time_bins_xy) '\n\n'])


    handles_choices.dt=dt;
    handles_choices.dt_miniscope=dt_miniscope;
    handles_choices.n_shuffle=n_shuffle;

  

else
    
    %Note: The data brought into the Kording lab jupyter notebbok seems to be
    %binned in 200 msec bins
    dt=handles_choices.dt;
    dt_miniscope=handles_choices.dt_miniscope;
    n_shuffle=handles_choices.n_shuffle;

    distances_in_mm=handles_choices.distances_in_mm;

    this_path=handles_choices.this_path;
    dFF_file=handles_choices.dFF_file;
    arena_file=handles_choices.arena_file;


    which_odor_plume=handles_choices.which_odor_plume;   %1=use odor plume from the publication, 2=use odor plume simulated by Aaron
    cm_from_floor=handles_choices.cm_from_floor;

    weber_fechner=handles_choices.weber_fechner; 
    alpha=handles_choices.alpha;
    multiplier=handles_choices.multiplier;
    %0 is Stevens Law, R proportional to C^alpha
    %1 is Weber-Flechner law R proportional to log(C)
    %See Copelli et al DOI: 10.1103/PhysRevE.65.060901

    algo=handles_choices.algo; 
    %1 ann
    %2 glm
    %3 tree
    %4 svm with radial basis
    %5 Gaussian process regression
    %6 GPR with 'KernelFunction', 'squaredexponential', 'Standardize', true
    %7 ann optimize hyper parameters, this is VERY slow!


    group=handles_choices.group; %1=odor plume, 2=spatial
    
  
    save_tag=handles_choices.save_tag; %This will be used in the name of the save file

    
    repeats=handles_choices.repeats;
    handles_choices.sh_repeats=handles_choices.repeats;
    best_repeat=handles_choices.best_repeat; %If multiple runs are performed choose the best result
    
    max_overlap=handles_choices.max_overlap; %Maximum overlap of the shuffled segments
    displayFigures=handles_choices.displayFigures;

   


    dt_decoding_op=handles_choices.dt_decoding_op; %All bins from -dt_decoding to 0 will be included in the decoder
    %i.e. given that the dt bin is 0.2 sec dt_decoding=0.1 will use only the
    %current bin, dt_decoding=2 will use the current bin and the last 9 bins
    no_dec_time_bins_op=ceil(handles_choices.dt_decoding_op/dt);
    fprintf(1,['Number of time bins for op decoding ' num2str(no_dec_time_bins_op) '\n\n'])


    dt_decoding_xy=handles_choices.dt_decoding_xy;
    no_dec_time_bins_xy=ceil(handles_choices.dt_decoding_xy/dt);
    fprintf(1,['Number of time bins for xy decoding ' num2str(no_dec_time_bins_xy) '\n\n'])


end

 


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


% dFF=readmatrix([this_path dFF_file]); %Timepoints x ROIs

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

%Now calculate the behavioral performance
trials.hit1=0;
trials.miss1=0;
trials.lane1=0;
trials.lane4=0;
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
            trials.lane1=trials.lane1+1;
            trials.lane1_ii_start(trials.lane1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.lane1_ii_end(trials.lane1)=ceil(trim_factor*trials.ii_lanewater1(this_water));
            dii_trial=[dii_trial trials.hit1_ii_end(trials.hit1)-trials.hit1_ii_start(trials.hit1)];
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.hit1_ii_start(trials.hit1);
            trials.odor_ii_end(trials.odor_trNo)=trials.hit1_ii_end(trials.hit1);
            trials.odor_trial_type(trials.odor_trNo)=1;
            trials.odor_lane(trials.odor_trNo)=1;
        else
            trials.miss1=trials.miss1+1;
            trials.miss1_ii_start(trials.miss1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.lane1=trials.lane1+1;
            trials.lane1_ii_start(trials.lane1)=floor(trim_factor*trials.ii_odor(trNo));
            trials.lane1_ii_end(trials.lane1)=-1;
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
            trials.lane4=trials.lane4+1;
            trials.lane4_ii_start(trials.lane4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.lane4_ii_end(trials.lane4)=ceil(trim_factor*trials.ii_lanewater4(this_water));
            dii_trial=[dii_trial trials.hit4_ii_end(trials.hit4)-trials.hit4_ii_start(trials.hit4)];
            trials.odor_trNo=trials.odor_trNo+1;
            trials.odor_ii_start(trials.odor_trNo)=trials.hit4_ii_start(trials.hit4);
            trials.odor_ii_end(trials.odor_trNo)=trials.hit4_ii_end(trials.hit4);
            trials.odor_trial_type(trials.odor_trNo)=3;
            trials.odor_lane(trials.odor_trNo)=4;
        else
            trials.miss4=trials.miss4+1;
            trials.miss4_ii_start(trials.miss4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.lane4=trials.lane4+1;
            trials.lane4_ii_start(trials.lane4)=floor(trim_factor*trials.ii_odor(trNo));
            trials.lane4_ii_end(trials.lane4)=-1;
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

for iil1=1:length(trials.lane1_ii_end)
    if trials.lane1_ii_end(iil1)==-1
        trials.lane1_ii_end(iil1)=trials.lane1_ii_start(iil1)+ceil(mean(dii_trial));
    end
end

for iil1=1:length(trials.lane4_ii_end)
    if trials.lane4_ii_end(iil1)==-1
        trials.lane4_ii_end(iil1)=trials.lane4_ii_start(iil1)+ceil(mean(dii_trial));
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


if handles_choices.displayFigures==1
    
    %Display the lane 1 trials
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on

    %These are hits
    for iil1=1:trials.hit1
        plot(pos_binned(trials.hit1_ii_start(iil1):trials.hit1_ii_end(iil1),1),pos_binned(trials.hit1_ii_start(iil1):trials.hit1_ii_end(iil1),2),'-r')
        plot(pos_binned(trials.hit1_ii_start(iil1),1),pos_binned(trials.hit1_ii_start(iil1),2),'or')
        plot(pos_binned(trials.hit1_ii_end(iil1),1),pos_binned(trials.hit1_ii_end(iil1),2),'xr')
    end

     %These are miss
    for iil1=1:trials.miss1
        plot(pos_binned(trials.miss1_ii_start(iil1):trials.miss1_ii_end(iil1),1),pos_binned(trials.miss1_ii_start(iil1):trials.miss1_ii_end(iil1),2),'-b')
        plot(pos_binned(trials.miss1_ii_start(iil1),1),pos_binned(trials.miss1_ii_start(iil1),2),'ob')
        plot(pos_binned(trials.miss1_ii_end(iil1),1),pos_binned(trials.miss1_ii_end(iil1),2),'xb')
    end
    
    set(gca, 'YDir', 'reverse');
    xlabel('x (pixels)')
    ylabel('y (pixels)')
    title('Lane 1 trajectories hit (red) and miss (blue)')
    

   %Display the lane 4 trials
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on

    %These are hits
    for iil4=1:trials.hit4
        plot(pos_binned(trials.hit4_ii_start(iil4):trials.hit4_ii_end(iil4),1),pos_binned(trials.hit4_ii_start(iil4):trials.hit4_ii_end(iil4),2),'-r')
        plot(pos_binned(trials.hit4_ii_start(iil4),1),pos_binned(trials.hit4_ii_start(iil4),2),'or')
        plot(pos_binned(trials.hit4_ii_end(iil4),1),pos_binned(trials.hit4_ii_end(iil4),2),'xr')
    end

    %These are miss
    for iil4=1:trials.miss4
        plot(pos_binned(trials.miss4_ii_start(iil4):trials.miss4_ii_end(iil4),1),pos_binned(trials.miss4_ii_start(iil4):trials.miss4_ii_end(iil4),2),'-b')
        plot(pos_binned(trials.miss4_ii_start(iil4),1),pos_binned(trials.miss4_ii_start(iil4),2),'ob')
        plot(pos_binned(trials.miss4_ii_end(iil4),1),pos_binned(trials.miss1_ii_end(iil4),2),'xb')
    end

    set(gca, 'YDir', 'reverse');
    xlabel('x (mm)')
    ylabel('y (mm)')
    title('Lane 4 trajectories hit (red) and miss (blue)')

end

% nan_mask=logical(ones(1,no_time_bins));
% for ii_neuron=1:no_neurons
%     this_nan_mask=ones(1,no_time_bins);
%     this_nan_mask(1,:)=~isnan(neural_data(:,ii_neuron));
%     nan_mask=nan_mask&this_nan_mask;
% end
% neural_data_trimmed=neural_data(nan_mask,:);
% pos_binned_trimmed=pos_binned(nan_mask,:);
% no_time_bins=sum(nan_mask);

if handles_choices.trim_neural_data==1
    training_range_template=zeros(1,no_time_bins);
    for trNo=1:trials.odor_trNo
        x_predictedstart=trials.odor_ii_start(trNo);
        x_predictedend=trials.odor_ii_end(trNo);
        training_range_template(x_predictedstart:x_predictedend)=1;
    end
else
    training_range_template=ones(1,no_time_bins);
end

% no_time_bins=size(neural_data,1);

no_ROIs=size(neural_data,2);

if handles_choices.normalize==0
    
    trimmed_neural_data=zeros(sum(training_range_template),no_ROIs);
    trimmed_neural_data(:,:)=neural_data(logical(training_range_template),:);

    no_time_bins_trimmed=sum(training_range_template);
    % dt=0.2;
    %   handles_choices.dt=dt;
    %   handles_choices.cross_shift_dt=5;
    pair_rhos=[];
    pair_rhos.ii_pair_rho=0;

    no_bins_shift=handles_choices.cross_shift_dt/dt;
    max_rho=[];
    max_rho_ii=[];
    for ii_ROI1=1:no_ROIs
        for ii_ROI2=ii_ROI1+1:no_ROIs
            ii_bin=0;
            pair_rhos.ii_pair_rho=pair_rhos.ii_pair_rho+1;
            pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI1=ii_ROI1;
            pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI2=ii_ROI2;
            for bin_shift=-no_bins_shift:no_bins_shift
                this_dFF1=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                this_dFF1(1,:)=trimmed_neural_data(no_bins_shift+1:no_time_bins_trimmed-no_bins_shift,ii_ROI1)';
                this_dFF2=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                this_dFF2(1,:)=trimmed_neural_data(no_bins_shift+1+bin_shift:no_time_bins_trimmed-no_bins_shift+bin_shift,ii_ROI2)';
                thisRho=corrcoef([this_dFF1' this_dFF2']);
                ii_bin=ii_bin+1;
                pair_rhos.pair(pair_rhos.ii_pair_rho).rho(ii_bin)=thisRho(1,2);
            end
            [M,I]=max(pair_rhos.pair(pair_rhos.ii_pair_rho).rho);
            max_rho=[max_rho M];
            max_rho_ii=[max_rho_ii I];
        end
    end

    %Now do the same calculation for circularly permuted dFF
    max_rho_sh=[];
    for ii_shuffle=1:n_shuffle
        ii_pair_rho=0;
        shift_amount = randi([1, no_time_bins_trimmed]);
        no_bins_shift=handles_choices.cross_shift_dt/dt;
        for ii_ROI1=1:no_ROIs
            for ii_ROI2=ii_ROI1+1:no_ROIs
                % Perform the circular shift
                shifted_trimmed_neural_dataROI2=zeros(size(trimmed_neural_data,1),1);
                this_dFF_ROI2=zeros(size(trimmed_neural_data,1),1);
                this_dFF_ROI2(:,1)=trimmed_neural_data(:,ii_ROI2);
                shifted_trimmed_neural_dataROI2(:,1) = circshift(this_dFF_ROI2, shift_amount);
                ii_bin=0;
                ii_pair_rho=ii_pair_rho+1;
                pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI1=ii_ROI1;
                pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI2=ii_ROI2;
                for bin_shift=-no_bins_shift:no_bins_shift
                    this_dFF1=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                    this_dFF1(1,:)=trimmed_neural_data(no_bins_shift+1:no_time_bins_trimmed-no_bins_shift,ii_ROI1)';
                    this_dFF2=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                    this_dFF2(1,:)=shifted_trimmed_neural_dataROI2(no_bins_shift+1+bin_shift:no_time_bins_trimmed-no_bins_shift+bin_shift,1)';
                    thisRho=corrcoef([this_dFF1' this_dFF2']);
                    ii_bin=ii_bin+1;
                    pair_rhos.shuffle(ii_shuffle).pair(ii_pair_rho).rho(ii_bin)=thisRho(1,2);
                end
                max_rho_sh=[max_rho_sh max(pair_rhos.shuffle(ii_shuffle).pair(ii_pair_rho).rho)];
            end
        end
    end

    %Plot the histograms
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    histogram(max_rho_sh, 'BinWidth', 0.01, 'Normalization','probability')
    hold on
    histogram(max_rho, 'BinWidth', 0.01, 'Normalization','probability')
    title('Histograms for max rho for pair dF/F')

    p_thr=prctile(max_rho_sh,99);



    ii_ROI_sig_rho=zeros(1,no_ROIs);
    sig_rhos=[];
    sig_rho_dts=[];
    dts=dt*([1:2*no_bins_shift+1]-(no_bins_shift+1));
    no_sig_pairs=0;
    dt_thr=1;
    ii_pair_sig_rhos=[];
    for ii_pair=1:length(pair_rhos.pair)
        if (max(pair_rhos.pair(ii_pair).rho)>p_thr)&(dts(max_rho_ii(ii_pair))>dt_thr)
            ii_ROI_sig_rho(pair_rhos.pair(ii_pair).ii_ROI1)=1;
            ii_ROI_sig_rho(pair_rhos.pair(ii_pair).ii_ROI2)=1;
            sig_rhos=[sig_rhos max_rho(ii_pair)];
            sig_rho_dts=[sig_rho_dts dts(max_rho_ii(ii_pair))];
            no_sig_pairs=no_sig_pairs+1;
            ii_pair_sig_rhos=[ii_pair_sig_rhos ii_pair];
        end
    end

    fprintf(1,['Number of ROIs with significant pair crosscorr ' num2str(sum(ii_ROI_sig_rho)), ' out of ' num2str(no_ROIs) '\n'])
    fprintf(1,['Number of ROI pairs with significant pair crosscorr ' num2str(no_sig_pairs) '\n\n'])


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    histogram(sig_rhos,'BinWidth', 0.05)
    xlabel('rhos')
    title('Histograms for significant max rhos for pair dF/F')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    histogram(abs(sig_rho_dts),'BinWidth', 1)
    xlabel('Time shift (sec)')
    title('Histograms for significant max rhos dts (sec)')

    %For the top few pairs plot the dF/Fs and the cross correlograms
    to_sort=[sig_rhos' ii_pair_sig_rhos' sig_rho_dts'];
    sorted_out=sortrows(to_sort,'descend');
    ii_pair_sig_rhos_sorted=zeros(1,length(sig_rhos));
    ii_pair_sig_rhos_sorted(1,:)=sorted_out(:,2);

    this_time=dt*[1:size(trimmed_neural_data,1)];
    for ii=1:no_to_display
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

        hold on

        this_ii_ROI1=pair_rhos.pair(ii_pair_sig_rhos_sorted(ii)).ii_ROI1;
        this_ii_ROI2=pair_rhos.pair(ii_pair_sig_rhos_sorted(ii)).ii_ROI2;

        plot(this_time,trimmed_neural_data(:,this_ii_ROI1),'-b')
        plot(this_time,trimmed_neural_data(:,this_ii_ROI2)+max(trimmed_neural_data(:,this_ii_ROI1)),'-r')

        xlabel('Time (sec)')
        title(['dF/F for Rois ' num2str(this_ii_ROI1) ' and ' num2str(this_ii_ROI2)])

    end

    pffft=1;

else
    %The dF/F per transient is quite uneven leading to focusing on the large
    %transients. Let's make all transients above a threshold the same size

    %Calculate the threshold for binary conversion
    neural_data_threshold=prctile(neural_data(:),5);
    neural_data_binary=zeros(size(neural_data,1),size(neural_data,2));

    logical_threshold=neural_data > neural_data_threshold;
    neural_data_binary(logical_threshold)=1;

    logical_threshold=neural_data <= neural_data_threshold;
    neural_data_binary(logical_threshold)=0;

    normalized_neural_data=zeros(size(neural_data,1),size(neural_data,2));

    for ii_ROI=1:no_ROIs
        at_end=0;
        ii=1;
        this_neural_data_binary=zeros(size(neural_data,1),1);
        this_neural_data_binary(:,1)=neural_data_binary(:,ii_ROI);
        this_neural_data=zeros(size(neural_data,1),1);
        this_neural_data(:,1)=neural_data(:,ii_ROI);
        while at_end==0
            if ~isempty(find(this_neural_data_binary(ii:end)==1,1,'first'))
                %Find the start of this transient
                delta_ii_next=find(this_neural_data_binary(ii:end)==1,1,'first');
                %Find the end of this transient
                if ~isempty(find(this_neural_data_binary(ii-1+delta_ii_next:end)==0,1,'first'))
                    delta_ii_next_end=find(this_neural_data_binary(ii-1+delta_ii_next:end)==0,1,'first');
                    this_max_dFF=max(this_neural_data(ii-1+delta_ii_next:ii-1+delta_ii_next+delta_ii_next_end-1));
                    normalized_neural_data(ii-1+delta_ii_next:ii-1+delta_ii_next+delta_ii_next_end-1,ii_ROI)=...
                        this_neural_data(ii-1+delta_ii_next:ii-1+delta_ii_next+delta_ii_next_end-1)/this_max_dFF;
                else
                    at_end=1;
                end
                ii=ii-1+delta_ii_next+delta_ii_next_end-1;
            else
                at_end=1;
            end
        end
    end

    trimmed_normalized_neural_data=zeros(sum(training_range_template),no_ROIs);
    trimmed_normalized_neural_data(:,:)=normalized_neural_data(logical(training_range_template),:);
    no_time_bins_trimmed=sum(training_range_template);
    % dt=0.2;
    %   handles_choices.dt=dt;
    %   handles_choices.cross_shift_dt=5;
    pair_rhos=[];
    pair_rhos.ii_pair_rho=0;

    no_bins_shift=handles_choices.cross_shift_dt/dt;
    max_rho=[];
    max_rho_ii=[];
    for ii_ROI1=1:no_ROIs
        for ii_ROI2=ii_ROI1+1:no_ROIs
            ii_bin=0;
            pair_rhos.ii_pair_rho=pair_rhos.ii_pair_rho+1;
            pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI1=ii_ROI1;
            pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI2=ii_ROI2;
            for bin_shift=-no_bins_shift:no_bins_shift
                this_dFF1=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                this_dFF1(1,:)=trimmed_normalized_neural_data(no_bins_shift+1:no_time_bins_trimmed-no_bins_shift,ii_ROI1)';
                this_dFF2=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                this_dFF2(1,:)=trimmed_normalized_neural_data(no_bins_shift+1+bin_shift:no_time_bins_trimmed-no_bins_shift+bin_shift,ii_ROI2)';
                thisRho=corrcoef([this_dFF1' this_dFF2']);
                ii_bin=ii_bin+1;
                pair_rhos.pair(pair_rhos.ii_pair_rho).rho(ii_bin)=thisRho(1,2);
            end
            [M,I]=max(pair_rhos.pair(pair_rhos.ii_pair_rho).rho);
            max_rho=[max_rho M];
            max_rho_ii=[max_rho_ii I];
        end
    end

    %Now do the same calculation for circularly permuted dFF
    max_rho_sh=[];
    for ii_shuffle=1:n_shuffle
        ii_pair_rho=0;
        shift_amount = randi([1, no_time_bins_trimmed]);
        no_bins_shift=handles_choices.cross_shift_dt/dt;
        for ii_ROI1=1:no_ROIs
            for ii_ROI2=ii_ROI1+1:no_ROIs
                % Perform the circular shift
                shifted_trimmed_normalized_neural_dataROI2=zeros(size(trimmed_normalized_neural_data,1),1);
                this_dFF_ROI2=zeros(size(trimmed_normalized_neural_data,1),1);
                this_dFF_ROI2(:,1)=trimmed_normalized_neural_data(:,ii_ROI2);
                shifted_trimmed_normalized_neural_dataROI2(:,1) = circshift(this_dFF_ROI2, shift_amount);
                ii_bin=0;
                ii_pair_rho=ii_pair_rho+1;
                pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI1=ii_ROI1;
                pair_rhos.pair(pair_rhos.ii_pair_rho).ii_ROI2=ii_ROI2;
                for bin_shift=-no_bins_shift:no_bins_shift
                    this_dFF1=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                    this_dFF1(1,:)=trimmed_normalized_neural_data(no_bins_shift+1:no_time_bins_trimmed-no_bins_shift,ii_ROI1)';
                    this_dFF2=zeros(1,no_time_bins_trimmed-2*no_bins_shift);
                    this_dFF2(1,:)=shifted_trimmed_normalized_neural_dataROI2(no_bins_shift+1+bin_shift:no_time_bins_trimmed-no_bins_shift+bin_shift,1)';
                    thisRho=corrcoef([this_dFF1' this_dFF2']);
                    ii_bin=ii_bin+1;
                    pair_rhos.shuffle(ii_shuffle).pair(ii_pair_rho).rho(ii_bin)=thisRho(1,2);
                end
                max_rho_sh=[max_rho_sh max(pair_rhos.shuffle(ii_shuffle).pair(ii_pair_rho).rho)];
            end
        end
    end

    %Plot the histograms
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    histogram(max_rho_sh, 'BinWidth', 0.01, 'Normalization','probability')
    hold on
    histogram(max_rho, 'BinWidth', 0.01, 'Normalization','probability')
    title('Histograms for max rho for pair dF/F')

    p_thr=prctile(max_rho_sh,99);



    ii_ROI_sig_rho=zeros(1,no_ROIs);
    sig_rhos=[];
    sig_rho_dts=[];
    dts=dt*([1:2*no_bins_shift+1]-(no_bins_shift+1));
    no_sig_pairs=0;
    dt_thr=1;
    ii_pair_sig_rhos=[];
    for ii_pair=1:length(pair_rhos.pair)
        if (max(pair_rhos.pair(ii_pair).rho)>p_thr)&(dts(max_rho_ii(ii_pair))>dt_thr)
            ii_ROI_sig_rho(pair_rhos.pair(ii_pair).ii_ROI1)=1;
            ii_ROI_sig_rho(pair_rhos.pair(ii_pair).ii_ROI2)=1;
            sig_rhos=[sig_rhos max_rho(ii_pair)];
            sig_rho_dts=[sig_rho_dts dts(max_rho_ii(ii_pair))];
            no_sig_pairs=no_sig_pairs+1;
            ii_pair_sig_rhos=[ii_pair_sig_rhos ii_pair];
        end
    end

    fprintf(1,['Number of ROIs with significant pair crosscorr ' num2str(sum(ii_ROI_sig_rho)), ' out of ' num2str(no_ROIs) '\n'])
    fprintf(1,['Number of ROI pairs with significant pair crosscorr ' num2str(no_sig_pairs) '\n\n'])


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    histogram(sig_rhos,'BinWidth', 0.05)
    xlabel('rhos')
    title('Histograms for significant max rhos for pair dF/F')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    histogram(abs(sig_rho_dts),'BinWidth', 1)
    xlabel('Time shift (sec)')
    title('Histograms for significant max rhos dts (sec)')

    %For the top few pairs plot the dF/Fs and the cross correlograms
    to_sort=[sig_rhos' ii_pair_sig_rhos' sig_rho_dts'];
    sorted_out=sortrows(to_sort,'descend');
    ii_pair_sig_rhos_sorted=zeros(1,length(sig_rhos));
    ii_pair_sig_rhos_sorted(1,:)=sorted_out(:,2);

    this_time=dt*[1:size(trimmed_normalized_neural_data,1)];
    for ii=1:no_to_display
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

        hold on

        this_ii_ROI1=pair_rhos.pair(ii_pair_sig_rhos_sorted(ii)).ii_ROI1;
        this_ii_ROI2=pair_rhos.pair(ii_pair_sig_rhos_sorted(ii)).ii_ROI2;

        plot(this_time,trimmed_normalized_neural_data(:,this_ii_ROI1),'-b')
        plot(this_time,trimmed_normalized_neural_data(:,this_ii_ROI2)+max(trimmed_normalized_neural_data(:,this_ii_ROI1)),'-r')

        xlabel('Time (sec)')
        title(['dF/F for Rois ' num2str(this_ii_ROI1) ' and ' num2str(this_ii_ROI2)])

    end

    pffft=1;
end
 
pffft=1;