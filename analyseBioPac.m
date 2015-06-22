clear all 
close all

NRUN = 5;
NTR = 350;

% load regressors - valence and arousal data
load('/triton/becs/scratch/braindata/shared/2pnhyperEMO/isps_regressors.mat'); % finARO, finVAL and tempreg

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
%     try
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
%     catch
%         disp(subjects{s})
%     end
    br(s,:) = temp_br(find(tempreg));
    hr(s,:) = temp_hr(find(tempreg));
    
    % unpack behavior
    val(s,:) = finVAL{s}(find(tempreg));
    aro(s,:) = finARO{s}(find(tempreg));
end

%% Run correlation
load('biopac_data.mat')
load('data_processed.mat')
br = zscore(br,[],2);
hr = zscore(hr,[],2);
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
data = corr_pval_hr_val;
permnum = 5000;

function permdist = correlation_sign_permutator(data,permnum)
% Get real correlation sum
corrsum = sum(data);
% Auxilary variables
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
hist(permdist);
thresh = prctile(permdist,[5 95]);
mean(permdist > corrsum);




