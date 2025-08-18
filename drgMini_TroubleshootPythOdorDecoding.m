%drgMiniTroubleshootPythOdorDecoding
%Used to troubleshoot decoding using Kording's python code

close all
clear all

%Load the data
%This is already working in the python code
this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220727_FCM19/';
this_file='20220727_FCM19_withodor_miniscope_sync_L1andL4_ncorre_ext_nonneg_pyth.mat';

%We load handles_outp with the following fields
%xdFF is the z scored neuronal activity (time points x number of neurons)
%op is the odor concetration in the odor plume that will be decoded from
%xdFF
%trials is a structure that has all the information on which trials are hits vs miss, etc
%trial_start_ii and trial_end_ii are the time points when each trial starts
%and ends in the time course 
%In python print(xdFF.shape) yields (46, 17888), which is transposed from
%Matlab's 1788, 46 for size(xdFF)
load([this_path this_file])

%These are already defined in the python code
bins_before=4; %How many bins of neural data prior to the output are used for decoding
bins_current=1; %Whether to use concurrent time bin of neural data
bins_after=0; %How many bins of neural data after the output are used for decoding

%We need to translate the code below to python
%We decode each trial and we train with the other hit trials
trial_start_ii=handles_outp.trial_start_ii;
trial_end_ii=handles_outp.trial_end_ii;
trials=handles_outp.trials;
xdFF=handles_outp.xdFF';
no_neurons=size(xdFF,1);
no_bins=size(xdFF,2);
op=handles_outp.op';
for trNo=1:length(trial_start_ii)

    %Get testing data
    X_test=xdFF(:,trial_start_ii(trNo)+bins_before:trial_end_ii(trNo)-bins_after);
    y_test=op(1,trial_start_ii(trNo)+bins_before:trial_end_ii(trNo)-bins_after);

    %Get training data
    ii=0;
    X_train=zeros(no_neurons,no_bins);
    y_train=zeros(1,no_bins);
    for this_trNo=1:length(trial_start_ii)
        if this_trNo~=trNo %Leave one out
            if (trials.hit1(this_trNo)==1)||(trials.hit4(this_trNo)==1) %Must be a hit             
                this_xdFF=xdFF(:,trial_start_ii(this_trNo):trial_end_ii(this_trNo));
                this_op=op(1,trial_start_ii(this_trNo):trial_end_ii(this_trNo));
                X_train(:,ii+1:ii+length(this_op))=this_xdFF;
                y_train(1,ii+1:ii+length(this_op))=this_op;
                ii=ii+length(this_op);
            end
        end
    end
    X_train=X_train(:,1:ii);
    y_train=y_train(1,1:ii);
    
    %Decoding code will go here

end

pffft=1;