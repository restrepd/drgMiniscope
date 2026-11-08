%drgMini_information_content_hit_vs_miss_odor_concv2
close all
clear all

is_sphgpu=1;

switch is_sphgpu
    case 0



        %Trained with all trials and mult=1 taking on account when mouse detects the odor
        save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc_1_07162026/';
        choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_1_07162026.m';
    
        %This one has the odor encounter
        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle05152025/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_05102025.m';
        
        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
        fileID = fopen([choiceBatchPathName 'odor_conc_per_trial_stats.txt'],'w');

    case 1
     %Trained with all trials and mult=1 taking on account when mouse detects the odor
        save_PathConc='/data2/SFTP/DecodeDynOdorConc_1_07162026/';
        choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_bin_all_1_07162026.m';
    
        %This one has the odor encounter
        save_PathAngle='/data2/SFTP/Angle05152025/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_05102025.m';
        
        choiceBatchPathName='/data2/SFTP/PreProcessed/';
        fileID = fopen([choiceBatchPathName 'odor_conc_per_trial_stats.txt'],'w');

        fileID = fopen('/data2/SFTP/PreProcessed/odor_conc_per_trial_stats.txt','w');
        addpath('/data2/DRMatlab/drgMiniscope')
        addpath('/data2/DRMatlab/m new/Chi Squared')
        addpath('/data2/DRMatlab/drgMaster')
        addpath(genpath('/data2/DRMatlab/m new/kakearney-boundedline-pkg-32f2a1f'))
end

%Speed for background air flow in mm/sec
air_flow_speed=50; %5 cm/sec = 50 mm/sec


percentile_stringency=95; %mean_imps at this percentile stringency are included as cells of high prediction importance

%Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4, 2 cm from floor
%Group 2 is rewarded, with odor lane 4, no odor in lane 1
%Group 3 is rewarded, with odor lane 1, no odor in lane 4
%Group 4 is rewarded, with no odor in lane 1 and lane 4
%Group 5 is rewarded, with ISO1 in both lane 1 and lane 4, 1 cm from floor

group_label{1}='Odor both lanes 2 cm';
group_label{2}='Odor lane 4';
group_label{3}='Odor lane 1';
group_label{4}='No odor';
group_label{5}='Odor both lanes 1 cm';


% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
% eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])
% eval(['handles_Moser=' choiceMoserFileName(1:end-2) ';'])

figureNo=0;

try
    colormap fire
    this_cmap=colormap;
catch
    %This is here because of the lack of definition of fire in Matlab 2025b
    % Load the saved colormap matrix from the .mat file
    load('/Users/restrepd/Documents/GitHub/drgMiniscope/fireColormap2023b.mat');
    this_cmap=fireMap;
end


% this_cmap(1,:)=[0.7 0.7 0.7];

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

these_groups=[1 5];
ii_run=1;

n_shuffle_SI=100;
%
% %Load the odor plumes
% load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odor_plume_patternsDR.mat')
%
% lowest_conc=-200;
% op_threshold=-3.5;
% odor_c_bounds=[-10:8/10:-2];
% for cm_from_floor=1:2
%     mean_plume_l1=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1';
%     mean_plume_l4=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4';
%     x_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume;
%     y_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume;
%
%     %Generate binary plumes
%     binary_plumel1=zeros(10,10);
%     binary_plumel4=zeros(10,10);
%
%     for ii_x=1:10
%         for ii_y=1:10
%             these_opl1=[];
%             these_opl1=mean_plume_l1( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
%             this_mean_opl1=mean(these_opl1(:));
%             if this_mean_opl1<op_threshold
%                 binary_plumel1(ii_x,ii_y)=0;
%             else
%                 binary_plumel1(ii_x,ii_y)=1;
%             end
%
%             these_opl4=[];
%             these_opl4=mean_plume_l4( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
%             this_mean_opl4=mean(these_opl4(:));
%             if this_mean_opl4<op_threshold
%                 binary_plumel4(ii_x,ii_y)=0;
%             else
%                 binary_plumel4(ii_x,ii_y)=1;
%             end
%         end
%     end
%
%     odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1=binary_plumel1;
%     odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4=binary_plumel4;
%
%     minC=prctile([mean_plume_l4(:); mean_plume_l1(:)],0.5);
%     if minC<lowest_conc
%         minC=lowest_conc;
%     end
%     maxC=max([max(mean_plume_l4(:)) max(mean_plume_l1(:))]);
%     if minC==maxC
%         maxC=minC+0.1;
%     end
%     odor_plume_patterns.cm_from_floor(cm_from_floor).minC=minC;
% end


%Now calculate odor concentration time courses per trial
figNo=0;
x=25:50:475;
y=24:48:456;

all_trial_OC=[];
all_trial_x=[];
all_trial_y=[];
all_trial_ii=0;
all_trial_file=[];
all_trial_mouse=[];
all_trial_hit_vs_miss=[];

%Let's do the glm and ZETA stats from -3 to 3 sec in 1 sec bins
time_range=[-3 3];
ii_for_time_bins=(time_range(2)-time_range(1))/handles_conc.dt;
odor_aligned_time_bins=[0:ii_for_time_bins]*handles_conc.dt+time_range(1);

handles_outic=[];

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)

        %Get XY and dFF per trial
        arena_file=handles_conc.arena_file{fileNo};
        load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        handles_out_conc=handles_out;
        dt=handles_conc.dt;

        %
        %
        % %OC will be calculated on a 0.5 sec bins
        % trimmed_dt=0.5;
        % trimmed_time_bins=[time_range(1)+(trimmed_dt/2):trimmed_dt:time_range(2)-(trimmed_dt/2)];

        angle_file=handles_Angle.arena_file{fileNo};
        %load the ouptut file
        load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        
        handles_out_angles=handles_out;
        angles=handles_out.angles;

        %Load odor plume and calculate the odor plume for lane 1 and lane 4

        %Use Aaron's simulation
        %Plot odor plume
        if is_sphgpu==1
            load('/data/SFTP/PreProcessedDR/Odor Arena Plumes/odorArenaPlumesDR.mat')
        else
            load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odorArenaPlumesDR.mat')
        end

        %Now make mean_plume matrices for lanes 1 and 4
        %Determine which lanes have odorants
        %Group 1 is rewarded, odor ISO1 in both lane 1 and lane 4
        %Group 2 is rewarded, with odor lane 4, no odor in lane 1
        %Group 3 is rewarded, with odor lane 1, no odor in lane 4
        %Group 4 is rewarded, with no odor in lane 1 and lane 4

        switch handles_conc.group(fileNo)
            case {1,5}
                %Groups 1 and 5 are rewarded, odor ISO1 in both lane 1 and lane 4
                %They differe in cm from floor
                lane1_odor_on=1;
                lane4_odor_on=1;
            case 2
                %Group 2 is rewarded, with odor lane 4, no odor in lane 1
                lane1_odor_on=0;
                lane4_odor_on=1;
            case 3
                %Group 3 is rewarded, with odor lane 1, no odor in lane 4
                lane1_odor_on=1;
                lane4_odor_on=0;
            case 4
                %Group 4 is rewarded, with no odor in lane 1 and lane 4
                %Note that we are testing whether we can "predict" the odor when
                %there is no odor applied
                lane1_odor_on=1;
                lane4_odor_on=1;
        end

        %Calculate mean for lane 4
        this_lane=4; %Note that I had mistakenly made lane 1 close to y=0, this is =1, not 4 even though I am calcualting lane 4
        for ii_source=1:length(odor_plumes.source)
            if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==handles_conc.cm_from_floor(fileNo))
                this_source=ii_source;
            end
        end
        x_for_plume=10*odor_plumes.source(this_source).x;
        y_for_plume=10*(odor_plumes.source(this_source).y-min(odor_plumes.source(this_source).y));
        mean_plume_l4=odor_plumes.source(this_source).mean_plume;

        %Calculate mean lane 1
        this_lane=1; %Note that I had mistakenly made lane 4 close to y=480, this is =4, not 1 even though I am calcualting lane 4
        for ii_source=1:length(odor_plumes.source)
            if (odor_plumes.source(ii_source).lane==this_lane)&(odor_plumes.source(ii_source).cm_from_floor==handles_conc.cm_from_floor(fileNo))
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

        minC=min([min(mean_plume_l4(:)) min(mean_plume_l1(:))]);
        % minC=prctile([mean_plume_l4(:); mean_plume_l1(:)],0.5);
        % minC=-9.6332;
        if minC<handles_conc.lowest_conc
            minC=handles_conc.lowest_conc;
        end
        maxC=max([max(mean_plume_l4(:)) max(mean_plume_l1(:))]);
        if minC==maxC
            maxC=minC+0.1;
        end

        % if handles_conc.displayFigures==1
        %     figNo=figNo+1;
        %     try
        %         close(figNo)
        %     catch
        %     end
        % 
        %     %Plot the shifted odor plume
        %     hFig = figure(figNo);
        %     set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
        % 
        %     drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l4')
        %     colormap fire
        %     shading interp
        %     caxis([minC maxC]);
        %     set(gca, 'YDir', 'reverse');
        %     % xlim(x_range)
        %     % ylim(y_range)
        %     % Ax = gca;
        %     % Ax.Color = 'k';
        %     xlabel('x (mm)')
        %     ylabel('y (mm)')
        % 
        %     title('Mean odor plume lane 4 before shift')
        % 
        %     figNo=figNo+1;
        %     try
        %         close(figNo)
        %     catch
        %     end
        % 
        %     %Plot the shifted odor plume
        %     hFig = figure(figNo);
        %     set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
        % 
        %     drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l1')
        %     colormap fire
        %     shading interp
        %     caxis([minC maxC]);
        %     set(gca, 'YDir', 'reverse');
        %     % xlim(x_range)
        %     % ylim(y_range)
        %     % Ax = gca;
        %     % Ax.Color = 'k';
        %     xlabel('x (mm)')
        %     ylabel('y (mm)')
        % 
        %     title('Mean odor plume lane 1 before shift')
        % end

        %Now shift to the actual dimensions of the odorant arena

        %Old code
        % new_mean_plume_l1=mean_plume_l1(y_for_plume>=20,:);
        % mean_plume_l1=[];
        % mean_plume_l1=new_mean_plume_l1;
        % new_mean_plume_l4=mean_plume_l4(y_for_plume<=480,:);
        % mean_plume_l4=new_mean_plume_l4;
        %
        % new_y_for_plume=y_for_plume(:,y_for_plume<=480);
        % y_for_plume=[];
        % y_for_plume=new_y_for_plume;
        %
        %Now shift the center of lane 4 to 70 mm
        % new_mean_plume_l4(y_for_plume>=20,:)=mean_plume_l4(y_for_plume<=460,:);
        % from_y_ii=find(y_for_plume<70,1,'first');
        % new_mean_plume_l4(length(y_for_plume):-1:from_y_ii,:)=mean_plume_l4((y_for_plume>50)&(y_for_plume<=120),:);
        % mean_plume_l4=new_mean_plume_l4;

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

        % %Note that I am replacing minC with prctile because there are a small
        % %number of points with very low min
        % minC=prctile([mean14_plume_l4(:); mean14_plume_l1(:)],0.5);
        % maxC=min([max(mean14_plume_l4(:)) max(mean14_plume_l1(:))]);
        % if maxC==minC
        %     maxC=minC+0.1;
        % end

        mean_plume_l4(mean_plume_l4<minC)=minC;
        mean_plume_l1(mean_plume_l1<minC)=minC;

        % if handles_conc.displayFigures==1
        %Setup the figures for the paths
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .7 .25 .25])


            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l4')
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

            title(['Lane 4 hits, file ' num2str(fileNo)])

             figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .35 .25 .25])


            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l4')
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

            title(['Lane 4 miss, file ' num2str(fileNo)])

            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.35 .7 .25 .25])

            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l1')
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

            title(['Lane 1 hit, file ' num2str(fileNo)])

            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.35 .35 .25 .25])

            drg_pcolor(repmat(x_for_plume,length(y_for_plume),1)',repmat(y_for_plume,length(x_for_plume),1),mean_plume_l1')
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

            title(['Lane 1 miss, file ' num2str(fileNo)])

            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.4 .25 .05 .3])

            prain=[minC:(maxC-minC)/99:maxC];
            drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
            colormap fire
            shading interp
            ax=gca;
            set(ax,'XTickLabel','')

            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            %Plot the shifted odor plume
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

            hold on
            histogram ([mean_plume_l1(:) mean_plume_l4(:)])

        % end
 

        odor_plume_template=handles_out_conc.odor_plume_template;
        pos_binned=handles_out_conc.pos_binned;
        odor_encounter_ii=trials.odor_encounter_ii;



        %Now find the calcium concentrations and the odor encounter
        lane1_hits=0;
        lane4_hits=0;
        lane1_miss=0;
        lane4_miss=0;
        for trNo=1:trials.odor_trNo
            op_predictedstart=odor_encounter_ii(trNo)-(ii_for_time_bins/2);
            op_predictedend=odor_encounter_ii(trNo)+(ii_for_time_bins/2);
            if (op_predictedstart>=1)&(op_predictedend<=length(odor_plume_template))

            

                these_x=[];
                these_y=[];
                these_ca=[];
                these_x=pos_binned(op_predictedstart:op_predictedend,1);
                these_y=pos_binned(op_predictedstart:op_predictedend,2);
                these_ca=odor_plume_template(op_predictedstart:op_predictedend);


                

                %Plot path taken by the mouse
                path_before_oc_x=these_x(1:1+(ii_for_time_bins/2));
                path_before_oc_y=these_y(1:1+(ii_for_time_bins/2));
                path_after_oc_x=these_x(1+(ii_for_time_bins/2):end);
                path_after_oc_y=these_y(1+(ii_for_time_bins/2):end);

                if trials.hit1(trNo)==1
                    figure(3)
                    lane1_hits=lane1_hits+1;
                end
                 if trials.hit4(trNo)==1
                    figure(1)
                    lane4_hits=lane4_hits+1;
                 end
                 if trials.miss1(trNo)==1
                    figure(4)
                    lane1_miss=lane1_miss+1;
                end
                 if trials.miss4(trNo)==1
                    figure(2)
                    lane4_miss=lane4_miss+1;
                end
                

                hold on
                plot(path_before_oc_x,path_before_oc_y,'-','Color',[0.7 0.7 0.7],'LineWidth',1)
                plot(path_after_oc_x,path_after_oc_y,'-k','LineWidth',1)

                %Assign values
                all_trial_ii=all_trial_ii+1;
                all_trial_OC(all_trial_ii,1:length(these_ca))=these_ca;
                all_trial_x(all_trial_ii,1:length(these_ca))=these_x;
                all_trial_y(all_trial_ii,1:length(these_ca))=these_y;
                all_trial_file(all_trial_ii)=fileNo;
                all_trial_mouse(all_trial_ii)=handles_conc.mouse(fileNo);
                if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                    all_trial_hit_vs_miss(all_trial_ii)=1;
                else
                    all_trial_hit_vs_miss(all_trial_ii)=0;
                end

            end
        end
    end

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .5 .2])


    hold on
    ii_start=0;

    last_op_predictedend=1;
    wall_threshold=50;
    op_all_trials=[];
    op_all_hits=[];
    op_all_miss=[];
    ii_start=1;

    for trNo=1:trials.odor_trNo

        op_predictedstart=odor_encounter_ii(trNo)-(ii_for_time_bins/2);
        op_predictedend=odor_encounter_ii(trNo)+(ii_for_time_bins/2);
        if (op_predictedstart>=1)&(op_predictedend<=length(odor_plume_template))
            op_all_trials=[op_all_trials; odor_plume_template(op_predictedstart:op_predictedend)'];
            ii_end=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))-1;


            %Okabe_Ito colors
            switch trials.odor_trial_type(trNo)
                case 1
                    %Lane 1 hits vermillion
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[230/255 159/255 0/255],'LineWidth',2)
                    op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];

                case 2
                    %Lane 1 miss orange
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend),'Color',[86/255 180/255 233/255],'LineWidth',2)
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];

                case 3
                    %Lane 4 hit blue
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[230/255 159/255 0/255],'LineWidth',2)
                    op_all_hits=[op_all_hits; odor_plume_template(op_predictedstart:op_predictedend)'];
                case 4
                    %Lane 4 miss sky blue
                    plot(dt*[ii_start:ii_end]',odor_plume_template(op_predictedstart:op_predictedend)','Color',[86/255 180/255 233/255],'LineWidth',2)
                    op_all_miss=[op_all_miss; odor_plume_template(op_predictedstart:op_predictedend)'];
            end
         
            plot(dt*[ii_start+(ii_for_time_bins/2)+1 ii_start+(ii_for_time_bins/2)+1],ylim,'-k')
            ii_start=ii_start+length(odor_plume_template(op_predictedstart:op_predictedend))+20;
        end
    end

    title(['Odor concentration per trial, file ' num2str(fileNo) ' vermillion: HIT, sky blue: MISS'])

    pfft=1;
    figNo=0;
    if fileNo==20
        pffft=1;
    end
end

all_trial_hits=zeros(sum(all_trial_hit_vs_miss==1),ii_for_time_bins+1);
all_trial_hits(:,:)=all_trial_OC(logical(all_trial_hit_vs_miss==1),:);
ii_hits=sum(all_trial_hit_vs_miss==1);

all_trial_miss=zeros(sum(all_trial_hit_vs_miss==0),ii_for_time_bins+1);
all_trial_miss(:,:)=all_trial_OC(logical(all_trial_hit_vs_miss==0),:);
ii_miss=sum(all_trial_hit_vs_miss==0);


%Plot all hit time courses sorted by the logC integral
figNo=7;
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .2 .7])

sum_all_trial_hits=sum(all_trial_hits,2);
to_sort_aths=[all_trial_hits sum_all_trial_hits];
sorted_aths=sortrows(to_sort_aths,ii_for_time_bins+2);
sorted_all_trial_hits=zeros(sum(all_trial_hit_vs_miss==1),ii_for_time_bins+1);
sorted_all_trial_hits(:,:)=sorted_aths(:,1:end-1);

%Pseudocolor
time_span_mat=repmat(odor_aligned_time_bins,ii_hits,1);
trial_mat=repmat(1:ii_hits,length(odor_aligned_time_bins),1)';

pcolor(time_span_mat,trial_mat,sorted_all_trial_hits)
% pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_hit1)
colormap(this_cmap)
shading flat

hold on
plot([0 0],[0 ii_hits],'-r')


xlim([-3 3])
ylim([1 ii_hits])
title('logC for hits, all mice')
xlabel('Time (sec)')
ylabel('Trial number')


%Plot all miss time courses sorted by the logC integral
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .2 .7])

sum_all_trial_miss=sum(all_trial_miss,2);
to_sort_aths=[all_trial_miss sum_all_trial_miss];
sorted_aths=sortrows(to_sort_aths,ii_for_time_bins+2);
sorted_all_trial_miss=zeros(sum(all_trial_hit_vs_miss==0),ii_for_time_bins+1);
sorted_all_trial_miss(:,:)=sorted_aths(:,1:end-1);

%Pseudocolor
time_span_mat=repmat(odor_aligned_time_bins,ii_miss,1);
trial_mat=repmat(1:ii_miss,length(odor_aligned_time_bins),1)';


pcolor(time_span_mat,trial_mat,sorted_all_trial_miss)
% pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_mis1)
colormap(this_cmap)
shading flat

hold on
plot([0 0],[0 ii_miss],'-r')


xlim([-3 3])
ylim([1 ii_miss])
title('logC time for miss, all mice')
xlabel('Time (sec)')
ylabel('Trial number')

%Now plot the time courses per mouse
glm_oc=[];
glm_oc_ii=0;

id_oc_ii=0;
input_oc_data=[];

edges=[0:10:250];
rand_offset=0.5;
bar_offset=0;

bar_figNo=figNo+7;
try
    close(bar_figNo)
catch
end

hFig = figure(bar_figNo);

set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
hold on

for ii_mouse=1:max(all_trial_mouse)

    all_trial_hits=zeros(sum((all_trial_hit_vs_miss==1)&(all_trial_mouse==ii_mouse)),ii_for_time_bins+1);
    all_trial_hits(:,:)=all_trial_OC(logical((all_trial_hit_vs_miss==1)&(all_trial_mouse==ii_mouse)),:);
    ii_hits=sum((all_trial_hit_vs_miss==1)&(all_trial_mouse==ii_mouse));
    cum_all_trial_hits=zeros(1,ii_hits);
    cum_all_trial_hits(1,:)=(sum(odor_aligned_time_bins>=0))*((sum(all_trial_hits(:,odor_aligned_time_bins>=0),2)/sum(odor_aligned_time_bins>=0))...
        -(sum(all_trial_hits(:,odor_aligned_time_bins<0),2)/sum(odor_aligned_time_bins<0)));

    %Hit bar
    hFig = figure(bar_figNo);
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(cum_all_trial_hits),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
    [mean_out, CIout]=drgViolinPoint(cum_all_trial_hits...
        ,edges,bar_offset,rand_offset,'k','k',3);

    glm_oc.data(glm_oc_ii+1:glm_oc_ii+ii_hits)=cum_all_trial_hits;
    glm_oc.mouse(glm_oc_ii+1:glm_oc_ii+ii_hits)=ii_mouse*ones(1,ii_hits);
    glm_oc.hit_miss(glm_oc_ii+1:glm_oc_ii+ii_hits)=1*ones(1,ii_hits);
    glm_oc_ii=glm_oc_ii+ii_hits;

    id_oc_ii=id_oc_ii+1;
    input_oc_data(id_oc_ii).data=cum_all_trial_hits;
    input_oc_data(id_oc_ii).description=['Hit, mouse no ' num2str(ii_mouse)];

    all_trial_miss=zeros(sum((all_trial_hit_vs_miss==0)&(all_trial_mouse==ii_mouse)),ii_for_time_bins+1);
    all_trial_miss(:,:)=all_trial_OC(logical((all_trial_hit_vs_miss==0)&(all_trial_mouse==ii_mouse)),:);
    ii_miss=sum((all_trial_hit_vs_miss==0)&(all_trial_mouse==ii_mouse));
    cum_all_trial_miss=zeros(1,ii_miss);
    cum_all_trial_miss(1,:)=(sum(odor_aligned_time_bins>=0))*((sum(all_trial_miss(:,odor_aligned_time_bins>=0),2)/sum(odor_aligned_time_bins>=0))...
        -(sum(all_trial_miss(:,odor_aligned_time_bins<0),2)/sum(odor_aligned_time_bins<0)));

        %Miss bar
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(cum_all_trial_miss),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
    [mean_out, CIout]=drgViolinPoint(cum_all_trial_miss...
        ,edges,bar_offset,rand_offset,'k','k',3);

    glm_oc.data(glm_oc_ii+1:glm_oc_ii+ii_miss)=cum_all_trial_miss;
    glm_oc.mouse(glm_oc_ii+1:glm_oc_ii+ii_miss)=ii_mouse*ones(1,ii_miss);
    glm_oc.hit_miss(glm_oc_ii+1:glm_oc_ii+ii_miss)=0*ones(1,ii_miss);
    glm_oc_ii=glm_oc_ii+ii_miss;

    id_oc_ii=id_oc_ii+1;
    input_oc_data(id_oc_ii).data=cum_all_trial_miss;
    input_oc_data(id_oc_ii).description=['Miss, mouse no ' num2str(ii_mouse)];

    bar_offset=bar_offset+1;

    %Plot all hit time courses sorted by the logC integral
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .2 .7])

    sum_all_trial_hits=sum(all_trial_hits,2);
    to_sort_aths=[all_trial_hits sum_all_trial_hits];
    sorted_aths=sortrows(to_sort_aths,ii_for_time_bins+2);
    sorted_all_trial_hits=zeros(sum((all_trial_hit_vs_miss==1)&(all_trial_mouse==ii_mouse)),ii_for_time_bins+1);
    sorted_all_trial_hits(:,:)=sorted_aths(:,1:end-1);

    %Pseudocolor
    time_span_mat=repmat(odor_aligned_time_bins,ii_hits,1);
    trial_mat=repmat(1:ii_hits,length(odor_aligned_time_bins),1)';

    pcolor(time_span_mat,trial_mat,sorted_all_trial_hits)
    % pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_hit1)
    colormap(this_cmap)
    shading flat

    hold on
    plot([0 0],[0 ii_hits],'-r')


    xlim([-3 3])
    ylim([1 ii_hits])
    title(['logC for hits, mouse ' num2str(ii_mouse)])
    xlabel('Time (sec)')
    ylabel('Trial number')


    %Plot all miss time courses sorted by the logC integral
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .2 .7])

    sum_all_trial_miss=sum(all_trial_miss,2);
    to_sort_aths=[all_trial_miss sum_all_trial_miss];
    sorted_aths=sortrows(to_sort_aths,ii_for_time_bins+2);
    sorted_all_trial_miss=zeros(sum((all_trial_hit_vs_miss==0)&(all_trial_mouse==ii_mouse)),ii_for_time_bins+1);
    sorted_all_trial_miss(:,:)=sorted_aths(:,1:end-1);

    %Pseudocolor
    time_span_mat=repmat(odor_aligned_time_bins,ii_miss,1);
    trial_mat=repmat(1:ii_miss,length(odor_aligned_time_bins),1)';


    pcolor(time_span_mat,trial_mat,sorted_all_trial_miss)
    % pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_mis1)
    colormap(this_cmap)
    shading flat

    hold on
    plot([0 0],[0 ii_miss],'-r')


    xlim([-3 3])
    ylim([1 ii_miss])
    title(['logC for miss, mouse ' num2str(ii_mouse)])
    xlabel('Time (sec)')
    ylabel('Trial number')

end

hFig = figure(bar_figNo);

xticks([1.5 4.5 7.5])
xticklabels({'Mouse 1','Mouse 2','Mouse 3'})


title(['logC integral per mouse'])
ylabel('logC integral')
xlim([0 9])
ylim([0 300])
text(4,280,'Reached OP','Color',[230/255 159/255 0/255])
text(4,260,'Missed OP','Color',[86/255 180/255 233/255])
 
%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(bar_figNo) ' logC integral per mouse\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(bar_figNo) ' logC integral per mouse\n']);


% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_oc.data',glm_oc.mouse',glm_oc.hit_miss',...
    'VariableNames',{'logC','mouse','reached'});
mdl = fitglm(tbl,'logC~mouse+reached'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for logC integral vs. mouse\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for logC integral vs. mouse\n']);


[output_data] = drgMutiRanksumorTtest(input_oc_data, fileID,0);

%Plot the speed of the mouse for each session
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .8 .8])
hold on
y_offset=0;
delta_y=200;
delta_y_mouse=200;

all_turn_speed=[];
all_turn_ii=0;
all_turn_file=[];
all_turn_mouse=[];
all_turn_hit_vs_miss=[];

%These are used to calculate the significance of the increase in speed
all_turn_peaks=[];
all_turn_before_peaks=[];
all_turn_peaks_rnd=[];
all_turn_before_peaks_rnd=[];


for ii_mouse=1:max(all_trial_mouse)
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(handles_conc.mouse(fileNo)==ii_mouse)


            %Get XY and dFF per trial
            arena_file=handles_conc.arena_file{fileNo};
            load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            trials=handles_out.trials;
            handles_out_conc=handles_out;
            dt=handles_conc.dt;

            pos_binned=handles_out_conc.pos_binned;
            odor_encounter_ii=trials.odor_encounter_ii;

            %Calculate speed
            speed=zeros(1,size(pos_binned,1));
            these_x=zeros(1,size(pos_binned,1));
            these_y=zeros(1,size(pos_binned,1));
            these_x(1,:)=pos_binned(:,1);
            these_y(1,:)=pos_binned(:,2);

            speed(1,2:end)=sqrt((these_x(1,2:end)-these_x(1,1:end-1)).^2 + (these_y(1,2:end)-these_y(1,1:end-1)).^2)/dt;
            speed(1)=speed(2);
            speed=speed/10; %cm/sec
            % max_speed=prctile(speed,99);
            plot([1:length(speed)]*dt,speed+y_offset,'-k')

            %Plot the turns
            angle_file=handles_Angle.arena_file{fileNo};
            %load the ouptut file
            load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
            handles_out_angles=handles_out;
            for trNo=1:trials.odor_trNo

                %Find the last turn and save the per turn speed
                this_ii_turn=find(handles_out_angles.angles.trial(trNo).delta_x>100,1,'last');
                if ~isempty(this_ii_turn)

                    ii_last_turn=handles_out_angles.angles.trial(trNo).ii_turns(this_ii_turn);
                    ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                    ii_last_turn=ii_last_turn+ii_predictedstart-1;

                    if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                        plot([ii_last_turn ii_last_turn]*dt,[y_offset y_offset+delta_y],'-','Color',[230/255 159/255 0/255])
                    else
                        plot([ii_last_turn ii_last_turn]*dt,[y_offset y_offset+delta_y],'-','Color',[86/255 180/255 233/255])
                    end

                    %Save the per turn speed
                    op_predictedstart=ii_last_turn-(ii_for_time_bins/2);
                    op_predictedend=ii_last_turn+(ii_for_time_bins/2);
                    if (op_predictedstart>=1)&(op_predictedend<=length(speed))
                        these_speeds=zeros(1,ii_for_time_bins+1);
                        these_speeds(1,:)=speed(op_predictedstart:op_predictedend);
                        all_turn_ii=all_turn_ii+1;
                        all_turn_speed(all_turn_ii,:)=these_speeds;
                        all_turn_file(all_turn_ii)=fileNo;
                        all_turn_mouse(all_turn_ii)=ii_mouse;
                        if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                            all_turn_hit_vs_miss(all_turn_ii)=1;
                        else
                            all_turn_hit_vs_miss(all_turn_ii)=0;
                        end
                        all_turn_peaks=[all_turn_peaks these_speeds((ii_for_time_bins/2)+2)];
                        all_turn_before_peaks=[all_turn_before_peaks these_speeds((ii_for_time_bins/2))];
                        rnd_ii=randperm(length(speed),1);
                        all_turn_peaks_rnd=[all_turn_peaks_rnd speed(rnd_ii)];
                        all_turn_before_peaks_rnd=[all_turn_before_peaks_rnd speed(rnd_ii)];
                    end
                end
            end

            y_offset=y_offset+delta_y;
        end
    end
    y_offset=y_offset+ delta_y_mouse+500;
end

xlabel('Time (sec)')
ylabel('Speed (cm/sec)')
title('Mouse speed and turn times')

[h, p_peak, ci,stats]=ttest2(all_turn_peaks,all_turn_peaks_rnd,'Vartype','unequal');
[h, p_before, ci,stats]=ttest2(all_turn_before_peaks,all_turn_before_peaks_rnd,'Vartype','unequal');

fprintf(1, ['\np value for t test for peak speed vs random speeds ' num2str(p_peak) ', n= ' num2str(length(all_turn_peaks)) '\n'])
fprintf(fileID, ['\np value for t test for peak speed vs random speeds ' num2str(p_peak) ', n= ' num2str(length(all_turn_peaks)) '\n']);

fprintf(1, ['\np value for t test for peak speed vs random speeds ' num2str(p_before) ', n= ' num2str(length(all_turn_peaks)) '\n'])
fprintf(fileID,['\np value for t test for peak speed vs random speeds ' num2str(p_before) ', n= ' num2str(length(all_turn_peaks)) '\n']);

x = all_turn_peaks;
y = all_turn_peaks_rnd;

bootfun = @(xbs, ybs) mean(xbs) - mean(ybs);
bootdiff = bootstrp(10000, bootfun, x, y);

ci = prctile(bootdiff, [2.5 97.5]);
obsdiff = mean(x) - mean(y);

p_peak_boot = 2 * min(mean(bootdiff <= 0), mean(bootdiff >= 0));

fprintf(1, ['\np value for bootstrap for difference in means for peak speed vs random speeds ' num2str(p_peak_boot) ', n= ' num2str(length(all_turn_peaks)) '\n'])
fprintf(fileID, ['\np value for bootstrap for difference in means for peak speed vs random speeds ' num2str(p_peak_boot) ', n= ' num2str(length(all_turn_peaks)) '\n']);

%Show the time course for speed aligned to turn
delta_y=120;
all_turn_speed_hits=zeros(sum((all_turn_hit_vs_miss==1)),ii_for_time_bins+1);
all_turn_speed_hits(:,:)=all_turn_speed(logical((all_turn_hit_vs_miss==1)),:)+delta_y;
ii_turn_hits=sum(all_turn_hit_vs_miss==1);

all_turn_speed_misses=zeros(sum((all_turn_hit_vs_miss==0)),ii_for_time_bins+1);
all_turn_speed_misses(:,:)=all_turn_speed(logical((all_turn_hit_vs_miss==0)),:);
ii_turn_misses=sum(all_turn_hit_vs_miss==0);

figNo = figNo + 1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
%Miss
CIpvsm = bootci(1000, @mean, all_turn_speed_misses);
meanpvsm=mean(all_turn_speed_misses,1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(all_turn_speed_misses), CIpvsm','cmap',[86/255 180/255 233/255]);

%Hits
CIpvsm = bootci(1000, @mean, all_turn_speed_hits);
meanpvsm=mean(all_turn_speed_hits,1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

[hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(all_turn_speed_hits), CIpvsm','cmap',[230/255 159/255 0/255]);

plot([0 0],[ylim],'-k')

text(-1,80,'Hits','Color',[230/255 159/255 0/255])
text(-1,70,'Misses','Color',[86/255 180/255 233/255])

ylim([-20 240])
xlabel('Time (sec)')
ylabel('Speed (cm/sec)')
title('Mouse speed aligned to turn')


%Show per mouse the time course for speed aligned to turn
delta_y=120;

for ii_mouse=1:max(all_trial_mouse)
    
    these_turns=(all_turn_hit_vs_miss==1)&(all_turn_mouse==ii_mouse);
    all_turn_speed_hits=zeros(sum(these_turns),ii_for_time_bins+1);
    all_turn_speed_hits(:,:)=all_turn_speed(logical(these_turns),:)+delta_y;
    ii_turn_hits=sum(these_turns);

    these_turns=(all_turn_hit_vs_miss==0)&(all_turn_mouse==ii_mouse);
    all_turn_speed_misses=zeros(sum(these_turns),ii_for_time_bins+1);
    all_turn_speed_misses(:,:)=all_turn_speed(logical(these_turns),:);
    ii_turn_misses=sum(these_turns);

    figNo = figNo + 1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    %Miss
    CIpvsm = bootci(1000, @mean, all_turn_speed_misses);
    meanpvsm=mean(all_turn_speed_misses,1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(all_turn_speed_misses), CIpvsm','cmap',[86/255 180/255 233/255]);

    %Hits
    CIpvsm = bootci(1000, @mean, all_turn_speed_hits);
    meanpvsm=mean(all_turn_speed_hits,1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(all_turn_speed_hits), CIpvsm','cmap',[230/255 159/255 0/255]);

    plot([0 0],[ylim],'-k')

    text(-1,80,'Hits','Color',[230/255 159/255 0/255])
    text(-1,70,'Misses','Color',[86/255 180/255 233/255])

    ylim([-20 240])
    xlabel('Time (sec)')
    ylabel('Speed (cm/sec)')
    title(['Mouse speed for mouse ' num2str(ii_mouse)])
end


%Now plot a bar graph of speed at t=0.1 sec per mouse
glm_sp=[];
glm_sp_ii=0;

id_sp_ii=0;
input_sp_data=[];

edges=[0:10:250];
rand_offset=0.5;
bar_offset=0;

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
hold on

for ii_mouse=1:max(all_trial_mouse)

    %Hits
    these_turns=(all_turn_hit_vs_miss==1)&(all_turn_mouse==ii_mouse);
    all_turn_speed_hits=zeros(sum(these_turns),ii_for_time_bins+1);
    all_turn_speed_hits(:,:)=all_turn_speed(logical(these_turns),:);
    peak_speed=zeros(1,size(all_turn_speed_hits,1));
    peak_speed(1,:)=all_turn_speed_hits(:, (ii_for_time_bins/2)+2);
    ii_turn_hits=sum(these_turns);


    %Hit bar
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(peak_speed),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
    [mean_out, CIout]=drgViolinPoint(peak_speed...
        ,edges,bar_offset,rand_offset,'k','k',3);

    glm_sp.data(glm_sp_ii+1:glm_sp_ii+ii_turn_hits)=peak_speed;
    glm_sp.mouse(glm_sp_ii+1:glm_sp_ii+ii_turn_hits)=ii_mouse*ones(1,ii_turn_hits);
    glm_sp.hit_miss(glm_sp_ii+1:glm_sp_ii+ii_turn_hits)=1*ones(1,ii_turn_hits);
    glm_sp_ii=glm_sp_ii+ii_turn_hits;

    id_sp_ii=id_sp_ii+1;
    input_sp_data(id_sp_ii).data=peak_speed;
    input_sp_data(id_sp_ii).description=['Reached OP, mouse no ' num2str(ii_mouse)];



    %Miss
    these_turns=(all_turn_hit_vs_miss==0)&(all_turn_mouse==ii_mouse);
    all_turn_speed_miss=zeros(sum(these_turns),ii_for_time_bins+1);
    all_turn_speed_miss(:,:)=all_turn_speed(logical(these_turns),:);
    peak_speed=zeros(1,size(all_turn_speed_miss,1));
    peak_speed(1,:)=all_turn_speed_miss(:, (ii_for_time_bins/2)+2);
    ii_turn_miss=sum(these_turns);


    %Hit bar
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(peak_speed),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
    [mean_out, CIout]=drgViolinPoint(peak_speed...
        ,edges,bar_offset,rand_offset,'k','k',3);

    glm_sp.data(glm_sp_ii+1:glm_sp_ii+ii_turn_miss)=peak_speed;
    glm_sp.mouse(glm_sp_ii+1:glm_sp_ii+ii_turn_miss)=ii_mouse*ones(1,ii_turn_miss);
    glm_sp.hit_miss(glm_sp_ii+1:glm_sp_ii+ii_turn_miss)=0*ones(1,ii_turn_miss);
    glm_sp_ii=glm_sp_ii+ii_turn_miss;

    id_sp_ii=id_sp_ii+1;
    input_sp_data(id_sp_ii).data=peak_speed;
    input_sp_data(id_sp_ii).description=['Missed OP, mouse no ' num2str(ii_mouse)];

    bar_offset=bar_offset+1;

end

xticks([1.5 4.5 7.5])
xticklabels({'Mouse 1','Mouse 2','Mouse 3'})


title(['Speed per mouse'])
ylabel('cm/sec')
xlim([0 9])
ylim([0 500])
text(4,420,'Reached OP','Color',[230/255 159/255 0/255])
text(4,460,'Missed OP','Color',[86/255 180/255 233/255])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(bar_figNo) ' speed per mouse\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(bar_figNo) ' speed per mouse\n']);


tbl = table(glm_sp.data',glm_sp.mouse',glm_sp.hit_miss',...
    'VariableNames',{'speed','mouse','reached'});
mdl = fitglm(tbl,'speed~mouse+reached'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for speed vs. mouse\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for speed vs. mouse\n']);


[output_data] = drgMutiRanksumorTtest(input_sp_data, fileID,0);
  
%Plot the turn times per session
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .8 .8])
hold on
y_offset=0;
delta_y=200;
delta_y_mouse=200;

all_turn_speed=[];
all_turn_x=[];
all_turn_y=[];
all_turn_ii=0;
all_turn_times=[];
all_turn_file=[];
all_turn_mouse=[];
all_turn_hit_vs_miss=[];
all_turn_lane=[];

for ii_mouse=1:max(all_trial_mouse)
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(handles_conc.mouse(fileNo)==ii_mouse)


            %Get XY and dFF per trial
            arena_file=handles_conc.arena_file{fileNo};
            load([save_PathConc arena_file(1:end-4) handles_conc.save_tag{ii_run} '.mat'])
            trials=handles_out.trials;
            handles_out_conc=handles_out;
            dt=handles_conc.dt;

            pos_binned=handles_out_conc.pos_binned;
            odor_encounter_ii=trials.odor_encounter_ii;

            %Calculate speed
            speed=zeros(1,size(pos_binned,1));
            these_x=zeros(1,size(pos_binned,1));
            these_y=zeros(1,size(pos_binned,1));
            these_x(1,:)=pos_binned(:,1);
            these_y(1,:)=pos_binned(:,2);

            speed(1,2:end)=sqrt((these_x(1,2:end)-these_x(1,1:end-1)).^2 + (these_y(1,2:end)-these_y(1,1:end-1)).^2)/dt;
            speed(1)=speed(2);
            speed=speed/10; %cm/sec
            % max_speed=prctile(speed,99);
            plot([1 length(speed)]*dt,[0 0]+y_offset,'-k')

            %Plot the turns
            angle_file=handles_Angle.arena_file{fileNo};
            %load the ouptut file
            load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
            handles_out_angles=handles_out;
            for trNo=1:trials.odor_trNo

                %Find the last turn and save the per turn speed
                this_ii_turn=find(handles_out_angles.angles.trial(trNo).delta_x>100,1,'last');
                if ~isempty(this_ii_turn)

                    ii_last_turn=handles_out_angles.angles.trial(trNo).ii_turns(this_ii_turn);
                    ii_predictedstart=trials.odor_ii_start(trNo)+handles_choices.trial_start_offset;
                    ii_last_turn=ii_last_turn+ii_predictedstart-1;

                    if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                        plot([ii_last_turn ii_last_turn]*dt,[y_offset y_offset+delta_y],'-','Color',[230/255 159/255 0/255])
                    else
                        plot([ii_last_turn ii_last_turn]*dt,[y_offset y_offset+delta_y],'-','Color',[86/255 180/255 233/255])
                    end

                    %Save the per turn speed
                    op_predictedstart=ii_last_turn-(ii_for_time_bins/2);
                    op_predictedend=ii_last_turn+(ii_for_time_bins/2);
                    if (op_predictedstart>=1)&(op_predictedend<=length(speed))
                        these_speeds=zeros(1,ii_for_time_bins+1);
                        these_speeds(1,:)=speed(op_predictedstart:op_predictedend);
                        these_x=zeros(1,ii_for_time_bins+1);
                        these_x(1,:)=pos_binned(op_predictedstart:op_predictedend,1);
                        these_y=zeros(1,ii_for_time_bins+1);
                        these_y(1,:)=pos_binned(op_predictedstart:op_predictedend,2);
                        all_turn_ii=all_turn_ii+1;
                        all_turn_speed(all_turn_ii,:)=these_speeds;
                        all_turn_x(all_turn_ii,:)=these_x;
                        all_turn_y(all_turn_ii,:)=these_y;
                        all_turn_times(all_turn_ii)=ii_last_turn*dt;
                        all_turn_file(all_turn_ii)=fileNo;
                        all_turn_mouse(all_turn_ii)=ii_mouse;
                        if (trials.hit1(trNo)==1)||(trials.hit4(trNo)==1)
                            all_turn_hit_vs_miss(all_turn_ii)=1;
                        else
                            all_turn_hit_vs_miss(all_turn_ii)=0;
                        end
                         if (trials.hit1(trNo)==1)||(trials.miss1(trNo)==1)
                            all_turn_lane(all_turn_ii)=1;
                        else
                            all_turn_lane(all_turn_ii)=4;
                        end
                    end
                end
            end

            y_offset=y_offset+delta_y;
        end
    end
    y_offset=y_offset+ delta_y_mouse+500;
end
 
%Trial rate
start_date=datetime('01/01/2022','InputFormat','MM/dd/yyyy');
days_per_session=[];
rates_per_session=[];
mouse_per_session=[];
ii_session=0;

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .8 .8])
hold on
for ii_mouse=1:max(all_trial_mouse)
    for fileNo=1:length(handles_conc.arena_file)
        if (sum(handles_conc.group(fileNo)==these_groups)>0)&(handles_conc.mouse(fileNo)==ii_mouse)
            these_turn_times=all_turn_times(all_turn_file==fileNo);
            these_rates=[];
            for ii_tt=2:length(these_turn_times)
                this_trial_rate=60/(all_turn_times(ii_tt)-all_turn_times(ii_tt-1));
                if ii_tt==2
                    last_trial_rate=this_trial_rate;
                end
                switch ii_mouse
                    case 1
                        plot([these_turn_times(ii_tt-1) these_turn_times(ii_tt)]/60, [last_trial_rate this_trial_rate]*60,'o-','LineWidth', 2 ...
                            ,'Color',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'MarkerFaceColor',[0/255 158/255 115/255])
                    case 2
                         plot([these_turn_times(ii_tt-1) these_turn_times(ii_tt)]/60, [last_trial_rate this_trial_rate]*60,'o-','LineWidth', 2 ...
                            ,'Color',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'MarkerFaceColor',[230/255 159/255 0/255])
                    case 3
                         plot([these_turn_times(ii_tt-1) these_turn_times(ii_tt)]/60, [last_trial_rate this_trial_rate]*60,'o-','LineWidth', 2 ...
                            ,'Color',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'MarkerFaceColor',[86/255 180/255 233/255])
                end
                last_trial_rate=this_trial_rate;
                these_rates=[these_rates this_trial_rate];
            end
            ii_session=ii_session+1;
            days_per_session(ii_session)=days(handles_conc.date{fileNo}-start_date);
            rates_per_session(ii_session)=mean(these_rates);
            mouse_per_session(ii_session)=handles_conc.mouse(fileNo);
        end
    end
    pffft=1;
end

xlabel('Time (min)')
ylabel('Trials per minute')
title('Trial rate per session')

%PLot trial rates across sessions
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

all_rates=[];
all_days=[];
 
for ii_mouse=1:max(all_trial_mouse)
    these_days=days_per_session(mouse_per_session==ii_mouse);
    these_days=these_days-min(these_days);
    these_rates=rates_per_session(mouse_per_session==ii_mouse);
    all_rates=[all_rates these_rates];
    all_days=[all_days these_days];
    switch ii_mouse
        case 1
            plot(these_days, these_rates,'o','LineWidth', 2 ...
                ,'Color',[0/255 158/255 115/255],'MarkerEdgeColor',[0/255 158/255 115/255],'MarkerFaceColor',[0/255 158/255 115/255])
        case 2
            plot(these_days, these_rates,'o','LineWidth', 2 ...
                ,'Color',[230/255 159/255 0/255],'MarkerEdgeColor',[230/255 159/255 0/255],'MarkerFaceColor',[230/255 159/255 0/255])
        case 3
            plot(these_days, these_rates,'o','LineWidth', 2 ...
                ,'Color',[86/255 180/255 233/255],'MarkerEdgeColor',[86/255 180/255 233/255],'MarkerFaceColor',[86/255 180/255 233/255])
    end
end

p=polyfit(all_days,all_rates,1);
slope=p(1);
intercept=p(2);
yfit=polyval(p,[-10 150]);
plot([-10 150],yfit,'-k')
ylim([0 4])
xlim([-10 150])
xlabel('Days')
ylabel('Trials per minute')
title('Trial rate')

[rho,pval]=corr(all_days',all_rates');
fprintf(1, ['\nRho for trial rate vs days ' num2str(rho) ', p value ' num2str(pval) '\n'])
fprintf(fileID, ['\nRho for trial rate vs days ' num2str(rho) ', p value ' num2str(pval) '\n']);


%Plot path mean and CI or SD per mouse
use_CI=0; %If use_CI==1 use CI, otherwise, use SD
for ii_mouse=1:max(all_trial_mouse)
    for ii_lane=[1 4]
        for ii_hit_miss=[0 1]
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);

            set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
            hold on

            set(gca, 'YDir', 'reverse');

            % Assume you already have:
            % paths: [nPoints x nTrials] matrix of x,y pairs
            % Here: nPoints = 800, nTrials = 61, sampled every 0.1 s, centered at turn

            % PARAMETERS
            nBoot = 1000;       % number of bootstrap resamples
            arenaSize = [500 480];     % mm, for axis scaling (square 50x48 arena)

            % [nPoints, nCols] = size(paths);

            % If paths is organized as [x;y] in rows or in some other way, adjust here.
            % Below assumes: each column is one trial, containing interleaved [x,y] pairs
            % or sequentially x then y. Since you said paths is 800x61, I'll assume:
            %   paths = [x1 y1 x2 y2 ...] is NOT the case.
            % More common is that each column is one trial and rows alternate x,y or
            % hold x and y separately. If your actual layout is:
            %   pathsX: [nPoints x nTrials]
            %   pathsY: [nPoints x nTrials]
            % then skip the parsing and set pathsX, pathsY directly.

            % --- Example layout: rows are time, columns alternate x and y ---
            % If your matrix is actually 800 x 122 and columns 1,3,5,... are x and
            % 2,4,6,... are y (61 trials), use this instead:
            %   xMat = paths(:,1:2:end);  % [nPoints x nTrials]
            %   yMat = paths(:,2:2:end);
            %
            % For your stated size (800x61), I'll assume:
            %   pathsX = pathsX_raw; % 800x61, x-coordinates
            %   pathsY = pathsY_raw; % 800x61, y-coordinates
            % Replace next two lines with your actual x,y matrices.
            
            pathsX=zeros(sum((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss)),size(all_turn_x,2));
            pathsX(:,:)=all_turn_x((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss),:);
            pathsY=zeros(sum((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss)),size(all_turn_x,2));
            pathsY(:,:)=all_turn_y((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss),:);

            pathsX=pathsX';
            pathsY=pathsY';
            
            [nPoints, nTrials] = size(pathsX);

            % Time vector (centered at turn)
            t  = (-nPoints/2 : nPoints/2-1) * dt;      % length 800, centered at 0

            % MEAN PATH
            meanX = mean(pathsX, 2);   % [nPoints x 1]
            meanY = mean(pathsY, 2);   % [nPoints x 1]

            if use_CI==1
                % BOOTSTRAP 95% CI
                bootMeanX = zeros(nPoints, nBoot);
                bootMeanY = zeros(nPoints, nBoot);

                for b = 1:nBoot
                    % resample trials with replacement
                    idx = randi(nTrials, [1 nTrials]);
                    sampleX = pathsX(:, idx);  % [nPoints x nTrials]
                    sampleY = pathsY(:, idx);

                    bootMeanX(:, b) = mean(sampleX, 2);
                    bootMeanY(:, b) = mean(sampleY, 2);
                end

                % 95% CI (2.5th and 97.5th percentiles) at each time point
                ciLowX  = prctile(bootMeanX, 2.5, 2);   % [nPoints x 1]
                ciHighX = prctile(bootMeanX, 97.5, 2);
                ciLowY  = prctile(bootMeanY, 2.5, 2);
                ciHighY = prctile(bootMeanY, 97.5, 2);

            else
                %Use SD
                % MEAN PATH
                sdX = std(pathsX,0,2);   % [nPoints x 1]
                sdY = std(pathsY,0,2);   % [nPoints x 1]

                ciLowX  = meanX -sdX;   % [nPoints x 1]
                ciHighX = meanX +sdX;
                ciLowY  = meanY-sdY;
                ciHighY = meanY+sdY;
            end

            % PLOT IN ARENA COORDINATES
            

            % Shaded confidence band around mean path.
            % We can shade in 2D by plotting the polygon formed by upper and lower CI.
            % Here we shade separately for x and y around the mean path.

            % Construct polygon for X CI vs Y mean (one way to visualize longitudinal uncertainty)
            px = [ciLowX; flipud(ciHighX)];
            py = [meanY;  flipud(meanY)];
            fill(px, py, [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

            % Construct polygon for Y CI vs X mean (transverse uncertainty)
            px2 = [meanX; flipud(meanX)];
            py2 = [ciLowY; flipud(ciHighY)];
            fill(px2, py2, [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

            % Plot mean path on top
            plot(meanX, meanY, 'k', 'LineWidth', 2);

            % Arena boundaries (optional)
            plot([0 arenaSize(1) arenaSize(1) 0 0], [0 0 arenaSize(2) arenaSize(2) 0], 'k--');

            % axis equal;
            xlim([0 arenaSize(1)]);
            ylim([0 arenaSize(2)]);
            xlabel('X (mm)');
            ylabel('Y (mm)');
            if ii_hit_miss==1
                title(['Mean path with CI, mouse ' num2str(ii_mouse) ' lane ' num2str(ii_lane) ' hit']);
            else
                title(['Mean path with CI, mouse ' num2str(ii_mouse) ' lane ' num2str(ii_lane) ' miss']);
            end
            % set(gca, 'Box', 'on');

        end
    end
end


%Plot path mean and CI for all mice

for ii_lane=[1 4]
    for ii_hit_miss=[0 1]
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
        hold on

        set(gca, 'YDir', 'reverse');

        % Assume you already have:
        % paths: [nPoints x nTrials] matrix of x,y pairs
        % Here: nPoints = 800, nTrials = 61, sampled every 0.1 s, centered at turn

        % PARAMETERS
        nBoot = 1000;       % number of bootstrap resamples
        arenaSize = [500 480];     % mm, for axis scaling (square 50x48 arena)

        % [nPoints, nCols] = size(paths);

        % If paths is organized as [x;y] in rows or in some other way, adjust here.
        % Below assumes: each column is one trial, containing interleaved [x,y] pairs
        % or sequentially x then y. Since you said paths is 800x61, I'll assume:
        %   paths = [x1 y1 x2 y2 ...] is NOT the case.
        % More common is that each column is one trial and rows alternate x,y or
        % hold x and y separately. If your actual layout is:
        %   pathsX: [nPoints x nTrials]
        %   pathsY: [nPoints x nTrials]
        % then skip the parsing and set pathsX, pathsY directly.

        % --- Example layout: rows are time, columns alternate x and y ---
        % If your matrix is actually 800 x 122 and columns 1,3,5,... are x and
        % 2,4,6,... are y (61 trials), use this instead:
        %   xMat = paths(:,1:2:end);  % [nPoints x nTrials]
        %   yMat = paths(:,2:2:end);
        %
        % For your stated size (800x61), I'll assume:
        %   pathsX = pathsX_raw; % 800x61, x-coordinates
        %   pathsY = pathsY_raw; % 800x61, y-coordinates
        % Replace next two lines with your actual x,y matrices.

        pathsX=zeros(sum((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss)),size(all_turn_x,2));
        pathsX(:,:)=all_turn_x((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss),:);
        pathsY=zeros(sum((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss)),size(all_turn_x,2));
        pathsY(:,:)=all_turn_y((all_turn_lane==ii_lane)&(all_turn_hit_vs_miss==ii_hit_miss),:);

        pathsX=pathsX';
        pathsY=pathsY';

        [nPoints, nTrials] = size(pathsX);

        % Time vector (centered at turn)
        t  = (-nPoints/2 : nPoints/2-1) * dt;      % length 800, centered at 0

        % MEAN PATH
        meanX = mean(pathsX, 2);   % [nPoints x 1]
        meanY = mean(pathsY, 2);   % [nPoints x 1]


        if use_CI==1
            % BOOTSTRAP 95% CI
            bootMeanX = zeros(nPoints, nBoot);
            bootMeanY = zeros(nPoints, nBoot);

            for b = 1:nBoot
                % resample trials with replacement
                idx = randi(nTrials, [1 nTrials]);
                sampleX = pathsX(:, idx);  % [nPoints x nTrials]
                sampleY = pathsY(:, idx);

                bootMeanX(:, b) = mean(sampleX, 2);
                bootMeanY(:, b) = mean(sampleY, 2);
            end

            % 95% CI (2.5th and 97.5th percentiles) at each time point
            ciLowX  = prctile(bootMeanX, 2.5, 2);   % [nPoints x 1]
            ciHighX = prctile(bootMeanX, 97.5, 2);
            ciLowY  = prctile(bootMeanY, 2.5, 2);
            ciHighY = prctile(bootMeanY, 97.5, 2);

        else
            %Use SD
            % MEAN PATH
            sdX = std(pathsX,0,2);   % [nPoints x 1]
            sdY = std(pathsY,0,2);   % [nPoints x 1]

            ciLowX  = meanX -sdX;   % [nPoints x 1]
            ciHighX = meanX +sdX;
            ciLowY  = meanY-sdY;
            ciHighY = meanY+sdY;
        end

        % PLOT IN ARENA COORDINATES


        % Shaded confidence band around mean path.
        % We can shade in 2D by plotting the polygon formed by upper and lower CI.
        % Here we shade separately for x and y around the mean path.

        % Construct polygon for X CI vs Y mean (one way to visualize longitudinal uncertainty)
        for ii=1:length(ciLowX)-1
            px = [ciLowX(ii) ciLowX(ii+1) ciHighX(ii+1) ciHighX(ii)];
            py = [meanY(ii) meanY(ii+1) meanY(ii+1) meanY(ii)];
            fill(px, py, [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        end

        % Construct polygon for Y CI vs X mean (transverse uncertainty)
        for ii=1:length(ciLowX)-1
            px2 = [meanX(ii) meanX(ii+1) meanX(ii+1) meanX(ii)];
            py2 = [ciLowY(ii) ciLowY(ii+1) ciHighY(ii+1) ciHighY(ii)];
            fill(px2, py2, [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        end

        % Plot mean path on top
        plot(meanX, meanY, 'k', 'LineWidth', 2);

        % Arena boundaries (optional)
        plot([0 arenaSize(1) arenaSize(1) 0 0], [0 0 arenaSize(2) arenaSize(2) 0], 'k--');

        % axis equal;
        xlim([0 arenaSize(1)]);
        ylim([0 arenaSize(2)]);
        xlabel('X (mm)');
        ylabel('Y (mm)');
        if ii_hit_miss==1
            title(['Mean path with CI, all mice, lane ' num2str(ii_lane) ' hit']);
        else
            title(['Mean path with CI, all mice, lane ' num2str(ii_lane) ' miss']);
        end
        % set(gca, 'Box', 'on');

    end
end
%Calculate measures of straightness


% If turn point is essentially the same across trials you can also
% fix it as e.g. turnX0, turnY0; here we allow small variability.
straightness=[];
meanAngleDev=[];
min_stepLen=20; %When the mouse stays in one place do not use the dx,dy for angle calculation
for i = 1:size(all_turn_x,1)

    pathsX=zeros(size(all_turn_x,1),size(all_turn_x,2));
    pathsX(:,:)=all_turn_x;
    pathsY=zeros(size(all_turn_x,1),size(all_turn_x,2));
    pathsY(:,:)=all_turn_y;

    pathsX=pathsX';
    pathsY=pathsY';

    x = pathsX(:, i);
    y = pathsY(:, i);

    % --- Define start (turn) and end (spout) ---

    ii_turn=1+((size(all_turn_x,2)-1)/2);
    turnPt  = [x(ii_turn), y(ii_turn)];
    if all_turn_lane(i)==1
        spoutPt = [1, 430];  % assumed fixed
    else
        spoutPt = [1, 70];  % assumed fixed
    end

    x_from_turn=x(ii_turn:end);
    y_from_turn=y(ii_turn:end);

    % Path length (sum of step-wise distances)
    dx = diff(x_from_turn);     % [nPoints-1 x 1]
    dy = diff(y_from_turn);
    stepLen = hypot(dx, dy);
    L_path = sum(stepLen);

    % Net displacement from turn to spout
    D_net = hypot(spoutPt(1) - turnPt(1), spoutPt(2) - turnPt(2));

  

    % Straightness index
    straightness(i) = D_net/ L_path; %Mouse is judged to reach the spout 50 mm away from the spout



    % --- Velocity and spout direction angles ---

    % Instantaneous velocity vectors
    vx = dx;
    vy = dy;

    % Heading angle of velocity (radians)
    heading = atan2(vy, vx);  % [-pi, pi]
    heading_valid=heading(stepLen>=min_stepLen);

    % Direction from current position to spout at each step
    spVecX = spoutPt(1) - x_from_turn(1:end-1);  % same length as vx
    spVecY = spoutPt(2) - y_from_turn(1:end-1);
    spAng  = atan2(spVecY, spVecX);    % [-pi, pi]
    spAng_valid=spAng(stepLen>=min_stepLen);

    % Angular deviation between velocity and spout direction
    % angDiff = heading - spAng;

    % Wrap to [-pi, pi]
    angDiff = atan2(sin(heading_valid-spAng_valid),cos(heading_valid-spAng_valid));
    % Convert to degrees for reporting
    angDiffDeg = rad2deg(angDiff);
   

    meanAngleDev(i) = mean(abs(angDiffDeg), 'omitnan');
    % % --- Heading circular variance ---
    %    % Use standard circular statistics: circVar = 1 - R
    %    % where R is mean resultant length of unit vectors e^(i*heading)
    %    cosH = cos(heading);
    %    sinH = sin(heading);
    %    R = hypot(mean(cosH), mean(sinH));   % mean resultant length
    %    headingCircVar(i) = 1 - R;
    %
    %    % --- Longitudinal vs lateral speed components ---
    %    % Define the straight reference direction from turn to spout
    %    refVec = spoutPt - turnPt;
    %    refLen = hypot(refVec(1), refVec(2));
    %    if refLen > 0
    %        refUnit = refVec / refLen;
    %    else
    %        refUnit = [NaN, NaN];
    %    end
    %
    %    % Project each velocity step onto refUnit
    %    vMat   = [vx, vy];      % [nPoints-1 x 2]
    %    vLong  = vMat * refUnit';   % component along ref direction
    %    vSpeed = stepLen;           % total speed per step
    %
    %    % Fraction of speed along spout direction (use absolute longitudinal component)
    %    longitudinalFrac(i) = sum(abs(vLong)) / sum(vSpeed);
end

%Plot meanAngleDev
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

bar_offset=0;

glm_st=[];
glm_st_ii=0;

id_st_ii=0;
input_st_data=[];

%Hit bar 
these_meanAngleDev=meanAngleDev((all_turn_hit_vs_miss==1)&(~isnan(meanAngleDev)));
these_mice=all_turn_mouse((all_turn_hit_vs_miss==1)&(~isnan(meanAngleDev)));
ii_hit=length(these_meanAngleDev);
bar_offset=bar_offset+1;
bar(bar_offset,mean(these_meanAngleDev),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
[mean_out, CIout]=drgViolinPoint(these_meanAngleDev...
    ,edges,bar_offset,rand_offset,'k','k',3);

glm_st.data(glm_st_ii+1:glm_st_ii+ii_hit)=these_meanAngleDev;
glm_st.mouse(glm_st_ii+1:glm_st_ii+ii_hit)=these_mice;
glm_st.hit_miss(glm_st_ii+1:glm_st_ii+ii_hit)=1*ones(1,ii_hit);
glm_st_ii=glm_st_ii+ii_hit;

% id_st_ii=id_st_ii+1;
% input_st_data(id_st_ii).data=cum_all_trial_miss;
% input_st_data(id_st_ii).description=['Miss, mouse no ' num2str(ii_mouse)];




%Miss bar
these_meanAngleDev=meanAngleDev((all_turn_hit_vs_miss==0)&(~isnan(meanAngleDev)));
these_mice=all_turn_mouse((all_turn_hit_vs_miss==0)&(~isnan(meanAngleDev)));
ii_miss=length(these_meanAngleDev);
bar_offset=bar_offset+1;
bar(bar_offset,mean(these_meanAngleDev),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
[mean_out, CIout]=drgViolinPoint(these_meanAngleDev...
    ,edges,bar_offset,rand_offset,'k','k',3);

glm_st.data(glm_st_ii+1:glm_st_ii+ii_miss)=these_meanAngleDev;
glm_st.mouse(glm_st_ii+1:glm_st_ii+ii_miss)=these_mice;
glm_st.hit_miss(glm_st_ii+1:glm_st_ii+ii_miss)=0*ones(1,ii_miss);
glm_st_ii=glm_st_ii+ii_miss;

% id_st_ii=id_st_ii+1;
% input_st_data(id_st_ii).data=cum_all_trial_miss;
% input_st_data(id_st_ii).description=['Miss, mouse no ' num2str(ii_mouse)];


xticks([1 2])
xticklabels({'present','absent'})
title('Heading angle')
ylabel('Degrees')

tbl = table(glm_st.data',glm_st.mouse',glm_st.hit_miss',...
    'VariableNames',{'logC','mouse','reached'});

mdl = fitglm(tbl,'logC~mouse+reached'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);

%Plot the heading per mouse


%Now plot a bar graph of heading per mouse
glm_hd=[];
glm_hd_ii=0;

id_hd_ii=0;
input_hd_data=[];

edges=[0:10:250];
rand_offset=0.5;
bar_offset=0;

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.3 .3 .2 .3])
hold on

for ii_mouse=1:max(all_trial_mouse)

    %Hits
    these_turns=(all_turn_hit_vs_miss==1)&(all_turn_mouse==ii_mouse)&(~isnan(meanAngleDev));
    all_turn_heading_hits=zeros(1,sum(these_turns));
    all_turn_heading_hits(1,:)=meanAngleDev(1,logical(these_turns));
    ii_turn_hits=sum(these_turns);

    %Hit bar
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(all_turn_heading_hits),'LineWidth', 3,'EdgeColor','none','FaceColor',[230/255 159/255 0/255])
    [mean_out, CIout]=drgViolinPoint(all_turn_heading_hits...
        ,edges,bar_offset,rand_offset,'k','k',2);

    glm_hd.data(glm_hd_ii+1:glm_hd_ii+ii_turn_hits)=all_turn_heading_hits;
    glm_hd.mouse(glm_hd_ii+1:glm_hd_ii+ii_turn_hits)=ii_mouse*ones(1,ii_turn_hits);
    glm_hd.hit_miss(glm_hd_ii+1:glm_hd_ii+ii_turn_hits)=1*ones(1,ii_turn_hits);
    glm_hd_ii=glm_hd_ii+ii_turn_hits;

    id_hd_ii=id_hd_ii+1;
    input_hd_data(id_hd_ii).data=all_turn_heading_hits;
    input_hd_data(id_hd_ii).description=['Reached OP, mouse no ' num2str(ii_mouse)];



    %Miss
    these_turns=(all_turn_hit_vs_miss==0)&(all_turn_mouse==ii_mouse)&(~isnan(meanAngleDev));
    all_turn_heading_miss=zeros(1,sum(these_turns));
    all_turn_heading_miss(1,:)=meanAngleDev(1,logical(these_turns));
    ii_turn_miss=sum(these_turns);


    %Miss bar
    bar_offset=bar_offset+1;
    bar(bar_offset,mean(all_turn_heading_miss),'LineWidth', 3,'EdgeColor','none','FaceColor',[86/255 180/255 233/255])
    [mean_out, CIout]=drgViolinPoint(all_turn_heading_miss...
        ,edges,bar_offset,rand_offset,'k','k',2);

    glm_hd.data(glm_hd_ii+1:glm_hd_ii+ii_turn_miss)=all_turn_heading_miss;
    glm_hd.mouse(glm_hd_ii+1:glm_hd_ii+ii_turn_miss)=ii_mouse*ones(1,ii_turn_miss);
    glm_hd.hit_miss(glm_hd_ii+1:glm_hd_ii+ii_turn_miss)=0*ones(1,ii_turn_miss);
    glm_hd_ii=glm_hd_ii+ii_turn_miss;

    id_hd_ii=id_hd_ii+1;
    input_hd_data(id_hd_ii).data=all_turn_heading_miss;
    input_hd_data(id_hd_ii).description=['Missed OP, mouse no ' num2str(ii_mouse)];

    bar_offset=bar_offset+1;

end

xticks([1.5 4.5 7.5])
xticklabels({'Mouse 1','Mouse 2','Mouse 3'})


title(['Heading per mouse'])
ylabel('degrees')
xlim([0 9])
ylim([0 180])
text(4,120,'Reached OP','Color',[230/255 159/255 0/255])
text(4,130,'Missed OP','Color',[86/255 180/255 233/255])

%Perform the glm  for errors
fprintf(1, ['\nglm for Fig ' num2str(bar_figNo) ' heading per mouse\n'])
fprintf(fileID, ['\nglm for Fig ' num2str(bar_figNo) ' heading per mouse\n']);


tbl = table(glm_hd.data',glm_hd.mouse',glm_hd.hit_miss',...
    'VariableNames',{'heading','mouse','reached'});
mdl = fitglm(tbl,'heading~mouse+reached'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for heading vs. mouse\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for heading vs. mouse\n']);


[output_data] = drgMutiRanksumorTtest(input_hd_data, fileID,0);

% 
% %Do the ranksum/t-test
% fprintf(1, ['\n\nRanksum or t-test p values for logC va mouse\n'])
% fprintf(fileID, ['\n\nRanksum or t-test p values for logC va mouse\n']);



%Plot all miss time courses
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .25 .7])

all_trial_miss=zeros(sum(all_trial_hit_vs_miss==0),ii_for_time_bins+1);
all_trial_miss(:,:)=all_trial_OC(logical(all_trial_hit_vs_miss==0),:);
ii_miss=sum(all_trial_hit_vs_miss==0);

%Plot hits

time_span_mat=repmat(odor_aligned_time_bins,ii_miss,1);
trial_mat=repmat(1:ii_miss,length(odor_aligned_time_bins),1)';

pcolor(time_span_mat,trial_mat,all_trial_miss)
% pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_hit1)
colormap(this_cmap)
shading flat

hold on
plot([0 0],[0 ii_miss],'-r')


xlim([-3 3])
ylim([1 ii_miss])
title('logC time courses per trial for miss')
xlabel('Time (sec)')
ylabel('Trial number')

all_trial_miss=zeros(sum(all_trial_hit_vs_miss==0),ii_for_time_bins+1);
all_trial_miss(:,:)=all_trial_OC(logical(all_trial_hit_vs_miss==0),:);


% croscorr_traces=corrcoef(all_trial_hits);
croscorr_traces=corrcoef(all_trial_OC);

Z = linkage(croscorr_traces,'complete','correlation');

no_clusters=3;
handles_out2.clusters = cluster(Z,'Maxclust',no_clusters);
figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

%Do cutoff for no_clusters
cutoff = median([Z(end-(no_clusters-1),no_clusters) Z(end-(no_clusters-2),no_clusters)]);
[H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
set(H,'LineWidth',2)
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

%re-sort the matrix
for ii=1:size(croscorr_traces,1)
    for jj=1:size(croscorr_traces,1)
        perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
    end
end



figureNo=figureNo+1;
try
    close(figureNo)
catch
end

hFig = figure(figureNo);

set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
hold on
pcolor(perm_croscorr_traces)
colormap(this_cmap)
shading flat

% caxis([-1  1])
clim([-1 1])
% eval(['title([''Cross correlations for all ROIs ' trial_type_comp_labels{ii_comp}  '''])'])
title('Cross correlations for OC time course for hits')

ROIs_included=size(all_div_to_sort,1);
xlim([1 ROIs_included])
ylim([1 ROIs_included])

for ii_type=1:2
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
    hold on

    eval(['sorted_handles_out.' these_adiv_names{ii_type} '=[];'])
    eval(['sorted_handles_out.' these_adiv_names{ii_type} '= ' these_adiv_names{ii_type} '(outperm,:);'])

    % sorted_handles_out.all_div_hit1=[];
    % sorted_handles_out.all_div_hit1=all_div_hit1(outperm,:);

    eval(['ii_included=size(sorted_handles_out.' these_adiv_names{ii_type} ',1);'])

    % ii_included=size(sorted_handles_out.all_div_hit1,1);

    time_span_mat=repmat(odor_aligned_time_bins,ii_included,1);
    ROI_mat=repmat(1:ii_included,length(odor_aligned_time_bins),1)';

    eval(['pcolor(time_span_mat,ROI_mat,sorted_handles_out.' these_adiv_names{ii_type} ')'])
    % pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_hit1)
    colormap(this_cmap)
    shading flat

    % caxis([prctile(sorted_handles_out.all_div_hit1(:),1) prctile(sorted_handles_out.all_div_hit1(:),99.9)])

    caxis([prctile(all_div(:),1) prctile(all_div(:),99.9)])

    plot([0 0],[0 ii_included],'-r')


    xlim([-3 3])
    ylim([1 ii_included])
    eval(['title([''' trial_type_comp_labels{ii_comp} ': ' trial_type_labels{ii_type} '''])'])
    xlabel('Time (sec)')
    ylabel('ROI number')
end




%Plot the average timecourses per cluster
do_std=0;
these_ylim=[];
for clus=1:no_clusters
    %Report the number of ROIs per cluster
    fprintf(1, ['\nFor '  trial_type_labels{trial_type_comp_ii_type1(ii_comp)}...
        ' vs ' trial_type_labels{trial_type_comp_ii_type2(ii_comp)} ' cluster No '...
        num2str(clus) ' has ' num2str(sum(handles_out2.clusters==clus)) ' ROIs\n'])

    clusters=handles_out2.clusters;

    %Plot the average of this cluster for the type of trials compared
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .2 .18 .18])


    hold on

    %get the dF/F

    % these_adiv_names{1}='all_div_hit1';
    % these_adiv_names{2}='all_div_miss1';
    % these_adiv_names{3}='all_div_hit4';
    % these_adiv_names{4}='all_div_miss1';

    % this_adiv_name1=these_adiv_names{trial_type_comp_ii_type1(ii_comp)};
    % this_adiv_name2=these_adiv_names{trial_type_comp_ii_type2(ii_comp)};


    %Miss
    this_cluster_miss=zeros(sum(clusters==clus),size(all_div_miss,2));
    ii_included=0;
    no_ROIs_this_clus= size(all_div_miss,1);
    for ii=1:no_ROIs_this_clus
        if clusters(ii)==clus
            ii_included=ii_included+1;
            this_cluster_miss(ii_included,:) = all_div_miss(ii,:);
        end
    end

    % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_miss);
        meanpv=mean(this_cluster_miss,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_miss);
        meanpv=mean(this_cluster_miss,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end


    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins,mean(this_cluster_miss), CIpv','cmap',[158/255 31/255 99/255]);

    %Hits
    this_cluster_hit=zeros(sum(clusters==clus),size(all_div_miss,2));

    ii_included=0;
    no_ROIs_this_clus= size(all_div_hit,1);
    for ii=1:no_ROIs_this_clus
        if clusters(ii)==clus
            ii_included=ii_included+1;
            this_cluster_hit(ii_included,:) = all_div_hit(ii,:);
        end
    end

    % dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

    if do_std==0
        CIpv = bootci(1000, @mean, this_cluster_hit);
        meanpv=mean(this_cluster_hit,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
    else
        STDpv = std(this_cluster_hit);
        meanpv=mean(this_cluster_hit,1);
        CIpv(1,:)=STDpv;
        CIpv(2,:)=STDpv;
    end



    [hlpvl, hppvl] = boundedline(odor_aligned_time_bins, mean(this_cluster_hit), CIpv','cmap',[0 114/255 178/255]);


    plot(odor_aligned_time_bins',mean(this_cluster_miss)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
    plot(odor_aligned_time_bins',mean(this_cluster_hit)','Color',[0 114/255 178/255],'LineWidth',1.5);


    % ylim([-0.5 1.5])
    xlim([-3 3])
    these_ylim=[these_ylim; ylim];




    xlabel('Time(sec)')
    ylabel('dFF')
    title(['dFF mean for hit vs miss cluster No ' num2str(clus)])
    % eval(['title([''' trial_type_comp_labels{ii_comp} ''' '' significantly diff, cluster no ''  num2str(clus) ])'])

end


fig_minus=0;
for clus=1:no_clusters

    figNo=figureNo+fig_minus;
    figure(figNo)
    ylim([min(these_ylim(:)) max(these_ylim(:))])
    this_ylim=ylim;

    %Odor on markers
    plot([0 0],this_ylim,'-k')

    text(-2.5,this_ylim(1)+0.9*(this_ylim(2)-this_ylim(1)),'Hit','Color',[158/255 31/255 99/255])
    text(-2.5,this_ylim(1)+0.8*(this_ylim(2)-this_ylim(1)),'Miss','Color',[0 114/255 178/255])


    %Odor on markers
    plot([0 0],this_ylim,'-k')

    fig_minus=fig_minus-1;

end


fclose(fileID);

pffft=1;

