%drgMini_FindTrialsFromArenaFile
close all
clear all

% %20221117_FCM22_lanes_1_4 for first figure in manuscript
% this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20221117_FCM22_lanes_1_4/';
% % dFF_file='20221117_FCM22_withodor_nearfloor_miniscope_sync_L1andL4_ncorre_fix_ext.mat';
% arena_file='20221117_FCM22withodor_nearfloor_odorarena_L1andL4_fix_sync_mm.mat';


this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Fabio New Data/20220929_FCM22_withodor_farfromfloor/';
arena_file='20220929_FCM22withtodor_farfromfloor_odorarena_L1andL4_sync.mat';

load([this_path arena_file])

%Extract trials
trials=[];

%Extract lanes using FLIR data
at_end=0;
ii=0;
trNo=0;
trNo_l1=0;
trNo_l4=0;
while at_end==0
    next_ii=find(arena.odor(ii+1:end)==1,1,'first');
    if ~isempty(next_ii)
        trNo=trNo+1;
        % trials.odor_ii(trNo)=ii+next_ii;
        % trials.x_odor(trNo)=arena.xsync(ii+next_ii);
        % trials.y_odor(trNo)=arena.ysync(ii+next_ii);

        ii=ii+next_ii;
        % ii_mini=arena.index_flirsynctominiscope(ii);

        if sum(arena.laneodor1(ii-3:ii+3)==1)>0
            %Note: laneodor1 is 1 only for one time point
            trials.lane_per_trial(trNo)=1;
        end

        if sum(arena.laneodor4(ii-3:ii+3)==1)>0
            %Note: laneodor4 is 1 only for one time point
            trials.lane_per_trial(trNo)=4;
        end

        next_ii=find(arena.odor(ii+1:end)==0,1,'first');
        if ~isempty(next_ii)
            ii=ii+next_ii;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end

fprintf(1,['\nThis file has ' num2str(trNo) ' trials\n\n'])

