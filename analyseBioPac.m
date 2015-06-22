clear all 
close all

NRUN = 5;
NTR = 350;

% load regressors - valence and arousal data
load('/triton/becs/scratch/braindata/shared/2pnhyperEMO/BioPacTesting/val_aro.mat'); % ARO, VAL and tempreg

% load biopac data preprocessed with DRIFTER
load('/triton/becs/scratch/braindata/shared/2pnhyperEMO/BioPacTesting/biopac_data.mat'); % breath and heartrate

subjects = {
'Emo_listening_4';
'Emo_listening_6';
'Emo_listening_7';
'Emo_listening_8';
'Emo_listening_9';
'Emo_listening_10';
'Emo_listening_12';
'Emo_listening_13';
'Emo_listening_15';
'Emo_listening_16';
'Emo_listening_17';
'Emo_listening_18';
'Emo_listening_19';
'Emo_listening_5_2';
'Emo_listening_10_2';
'Emo_listening_11_2';
'Emo_listening_12_2';
'Emo_listening_14_2';
'Emo_listening_15_2';
'Emo_listening_16_2';
'Emo_listening_17_2';
};


% drop last two subjects, as we don't have behavioral data for them
% and listening_14, as they have less data (maybe biopac turned off in the
% course of scanning)
breath = breath([1:8 10:22],:);
heartrate = heartrate([1:8 10:22],:);

%% Subsample
br = zeros(length(subjects),1015);
hr = zeros(length(subjects),1015);
aro = zeros(length(subjects),1015);
val = zeros(length(subjects),1015);
for s = 1:length(subjects)
    temp_br = [];
    temp_hr = [];
    
    %% prepro arousal
    TR=1.7;
    sample=0.2; %200ms, 5 times per sec
    TRms=TR*10; % Tenth of second
    samplems=sample*10;
    % 1) zscore
    taro = zscore(ARO{s});
    % 2) downsample
    taro = resample(taro,samplems,TRms);
    % 3) zscore again
    taro = zscore(taro);
    aro(s,:) = taro(find(tempreg));
    
    %% prepro valence
    % 1) zscore
    tval = zscore(VAL{s});
    % 2) downsample
    tval = resample(tval,samplems,TRms);
    % 3) zscore again
    tval = zscore(tval);
    val(s,:) = tval(find(tempreg));

    
    for r = 1:NRUN
        % breath
        tempfreq = breath{s,r}.frequency;
        dtempfreq = resample(tempfreq,1,170);
        dtempfreq = dtempfreq(1:NTR);
        temp_br = cat(2,temp_br,dtempfreq);
        % heartrate
        tempfreq = heartrate{s,r}.frequency;
        dtempfreq = resample(tempfreq,1,170);
        dtempfreq = dtempfreq(1:NTR);
        temp_hr = cat(2,temp_hr,dtempfreq);
    end    
    br(s,:) = zscore(temp_br(find(tempreg)),[],2);
    hr(s,:) = zscore(temp_hr(find(tempreg)),[],2);
    
end

%% Run correlation
addpath(genpath('/triton/becs/scratch/braindata/shared/toolboxes/bramila/bramila/'));
% run correlation
for s = 1:length(subjects)
    % Breath and arousal
    % Get rho and pval
    [rho_br_aro(s,1),pval_br_aro] = corr(br(s,:)',aro(s,:)','Type','Spearman');
    % Get adjusted df
    newdf_br_aro = bramila_autocorr(br(s,:)',aro(s,:)');
    % Get adjusted pval
    corr_pval_br_aro(s,1) = pvalPearson('b',rho_br_aro(s,1),newdf_br_aro);
    
    % Breath and valence
    [rho_br_val(s,1),pval_br_val] = corr(br(s,:)',val(s,:)','Type','Spearman');
    newdf_br_val = bramila_autocorr(br(s,:)',val(s,:)');
    corr_pval_br_val(s,1) = pvalPearson('b',rho_br_val(s,1),newdf_br_val);
    
    % Heartrate and arousal
    [rho_hr_aro(s,1),pval_hr_aro] = corr(hr(s,:)',aro(s,:)','Type','Spearman');
    newdf_hr_aro = bramila_autocorr(hr(s,:)',aro(s,:)');
    corr_pval_hr_aro(s,1) = pvalPearson('b',rho_hr_aro(s,1),newdf_hr_aro);
    
    % Heartrate and valence
    [rho_hr_val(s,1),pval_hr_val] = corr(hr(s,:)',val(s,:)','Type','Spearman');
    newdf_hr_val = bramila_autocorr(hr(s,:)',val(s,:)');
    corr_pval_hr_val(s,1) = pvalPearson('b',rho_hr_val(s,1),newdf_hr_val);
end
%% Get group p value
% PERMUTATION
data = rho_hr_val;
permnum = 5000;
[corrsum,permdist,permmap] = sum_sign_permutator(data,permnum);
hist(permdist);
thresh = prctile(permdist,[5 95]);
mean(permdist > corrsum);

% 
function [corrsum,permdist,permmap] = sum_sign_permutator(data,permnum)
% Parameters:
% data - vector of values to sum
% permnum - number of permutations
% Output:
% corrsum - sum of unpermuted data
% permdist - permutation sums
% permmap - map of permutations (to reproduce)

% Get real correlation sum
corrsum = sum(data);

permlen = length(data);
permmap = zeros(permlen,permnum);
for p = 1:permnum
    % store the permutation
    tperm = ones(permlen,1);
    iperm = randi([0 1],permlen,1) * -2;
    tperm = tperm + iperm;
    permmap(:,p) = tperm;
    
    % get correlation sum
    permcorr = data.*permmap(:,p);
    permdist(p) = sum(permcorr);
end


