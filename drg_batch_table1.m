%drgMini_batch_table1
close all
clear all

is_sphgpu=0;

switch is_sphgpu
    case 0

        % %There was a bug and the shuffled runs were the same for all shuffled runs
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc12192024/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_12192024.m';
        %
        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput12192024/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_12192024.m';

        %For the files below the shuffled runs should be different
        %Trained with all trials
        % save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeOdorConc01062025/';
        % choiceOdorConcFileName='drgOdorConcChoices_Fabio_Good_01062025.m';

         %Trained with hits only
         save_PathConc='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/DecodeDynOdorConc04192024/';
         choiceOdorConcFileName='drgDynamicOdorConcChoices_Fabio_Good_04192024.m';

        % save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01062925/';
        % choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01062025.m';

        %This one has the dFF per trial
        save_PathXY='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/OdorArenaOutput01122925/';
        choiceXYFileName='drgOdorArenaChoices_Fabio_Good_01122025.m';

        %Angle file
        % save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle12212024/';
        % choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_12212024.m';

        %This one has the odor encounter
        save_PathAngle='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Angle05152025/';
        choiceAngleFileName='drgMiniAngleChoices_Fabio_Good_05102025.m';

        % save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser12212024/';
        % choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_12192024.m';

        %This is not used here, all the Moser infor is input through PredImp
        save_PathMoser='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Moser02032025/';
        choiceMoserFileName='drgMiniMoserChoices_Fabio_Good_02032025.m';


        choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/CurrentChoices/';
        fileID = fopen([choiceBatchPathName 'decode_XYandconc_stats.txt'],'w');

        %The imps file with predictive importance values is be saved here
        save_PathPredImp='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FilePredImp='outputPredictionImportanceHitMiss.mat';

        %The output of drgMini_information_contentv2 is saved here
        save_PathIC='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FileIC='outputPerROIInformationContentHitMiss.mat';

        save_PathTable='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
        save_FileTable='Table1.xlsx';

    case 1
        fileID = fopen('/data2/SFTP/PreProcessed/decoder_odor_conc_hitmiss_stats.txt','w');
        addpath('/home/restrepd/Documents/MATLAB/drgMiniscope')
        addpath('/home/restrepd/Documents/MATLAB/m new/Chi Squared')
        addpath('/home/restrepd/Documents/MATLAB/drgMaster')
        addpath(genpath('/home/restrepd/Documents/MATLAB/m new/kakearney-boundedline-pkg-32f2a1f'))
end


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


%Mouse 1 FCM6
%Mouse 2 FCM19
%Mouse 3 FCM22


% handles.bins_before=[0 0 0 0 0 0 0 0 1 2 4];

addpath(choiceBatchPathName)
eval(['handles_conc=' choiceOdorConcFileName(1:end-2) ';'])
eval(['handles_XY=' choiceXYFileName(1:end-2) ';'])
eval(['handles_Angle=' choiceAngleFileName(1:end-2) ';'])
eval(['handles_Moser=' choiceMoserFileName(1:end-2) ';'])

figureNo=0;

colormap fire
this_cmap=colormap;
this_cmap(1,:)=[0.3 0.3 0.3];

%Find which files are included in the analysis
files_included = drgMini_included_files(handles_Angle,save_PathAngle, handles_conc, save_PathConc);

these_groups=[1 5];
ii_run=1;

n_shuffle_SI=100;

%Load the odor plumes
load('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Odor Arena Plumes/odor_plume_patternsDR.mat')

op_threshold=-3.5;
for cm_from_floor=1:2
    mean_plume_l1=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l1';
    mean_plume_l4=odor_plume_patterns.cm_from_floor(cm_from_floor).mean_plume_l4';
    x_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).x_for_plume;
    y_for_plume=odor_plume_patterns.cm_from_floor(cm_from_floor).y_for_plume;

    %Generate binary plumes
    binary_plumel1=zeros(10,10);
    binary_plumel4=zeros(10,10);

    for ii_x=1:10
        for ii_y=1:10
            these_opl1=[];
            these_opl1=mean_plume_l1( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
            this_mean_opl1=mean(these_opl1(:));
            if this_mean_opl1<op_threshold
                binary_plumel1(ii_x,ii_y)=0;
            else
                binary_plumel1(ii_x,ii_y)=1;
            end

             these_opl4=[];
            these_opl4=mean_plume_l4( (x_for_plume>=((ii_x-1)*50))&(x_for_plume<(ii_x*50)),(y_for_plume>=((ii_y-1)*48))&(y_for_plume<(ii_y*48)) );
            this_mean_opl4=mean(these_opl4(:));
            if this_mean_opl4<op_threshold
                binary_plumel4(ii_x,ii_y)=0;
            else
                binary_plumel4(ii_x,ii_y)=1;
            end
        end
    end

    odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel1=binary_plumel1;
    odor_plume_patterns.cm_from_floor(cm_from_floor).binary_plumel4=binary_plumel4;

end

%Now calculate mutual information
figNo=figureNo+1;
x=25:50:475;
y=24:48:456;

all_info_ii=0;
all_info_fileNo=[];
all_info_ii_ROI=[];
all_info_hit=[];
all_info_miss=[];
all_info_mutual_info14=[];

all_info_hit_sh=[];
all_info_miss_sh=[];
all_info_mutual_info14_sh=[];

all_info_op_bin=[];
all_info_mutual_info_dFFbin_xy_op_bin=[];
all_info_mutual_xy_op_bin=[];
handles_outic=[];

no_trials=[];
no_neurons=[];
mouse_no=[];
cm_from_floor=[];
ii_files=0;

for fileNo=1:length(handles_conc.arena_file)
    if (sum(handles_conc.group(fileNo)==these_groups)>0)&(files_included(fileNo)==1)

        ii_files=ii_files+1;

        %Get XY and dFF per trial
        arena_file=handles_XY.arena_file{fileNo};
        load([save_PathXY arena_file(1:end-4) handles_XY.save_tag{ii_run} '.mat'])
        trials=handles_out.trials;
        no_neurons(ii_files)=handles_out.no_neurons;
        no_trials(ii_files)=trials.odor_trNo;
        mouse_no(ii_files)=handles_conc.mouse(fileNo);

        if handles_conc.group(fileNo)==1
            cm_from_floor(ii_files)=2;
        end

        if handles_conc.group(fileNo)==5
            cm_from_floor(ii_files)=1;
        end

        % angle_file=handles_Angle.arena_file{fileNo};
        % %load the ouptut file
        % load([save_PathAngle angle_file(1:end-4) handles_Angle.save_tag '.mat'])
        % angles=handles_out.angles;

        pfft=1;
    end
end

% Create table
T = table([1:ii_files]',mouse_no', no_trials', no_neurons', cm_from_floor','VariableNames', {'File #','mouse #','# of trials','# of neurons','cm from floor'});

% Save so it opens in Excel
writetable(T, [save_PathTable save_FileTable]);

pffft=1;
