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
% run correlation, look how it is on average, what are the maximums and minimums
plot(br(1,:),aro(1,:),'.')
[rho,pval] = corr(br(1,:)',aro(1,:)','Type','Spearman');


%% Probably some sort of autocorrelation model makes more sense, as we want to check if frequency increases every time the signal in aro or val starts to change



