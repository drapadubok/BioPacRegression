clear all
close all

addpath(genpath('/triton/becs/scratch/braindata/shared/toolboxes/bramila/bramila'));
addpath('/triton/becs/scratch/braindata/shared/toolboxes/NIFTI');
addpath('/triton/becs/scratch/braindata/shared/toolboxes/bramila/bramila/external/DRIFTER')

rmvframes = 3; % How many volumes to remove in the beginning (sync trial)
nrun = 5;
dataroot = '/triton/becs/scratch/braindata/shared/2pnhyperEMO';

subjects = {
'Emo_listening_4';
'Emo_listening_6';
'Emo_listening_7';
'Emo_listening_8';
'Emo_listening_9';
'Emo_listening_10';
'Emo_listening_12';
'Emo_listening_13';
'Emo_listening_14';
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
'Emo_listening_18_2';
'Emo_listening_19_2'
};


%% Biopac settings, these are setup according to the template that Dima and Heini use
% Which channels were active 1 - active, 0 - not active
biopacfile.CH1=1; %usually breath data from belt
biopacfile.CH2=0; %usually ECG
biopacfile.CH3=0; %usually GSR
biopacfile.CH4=1; %usually heart rate from pulse oxymeter_left
biopacfile.CH5=0; %usually heart rate from pulse oxymeter_right
biopacfile.CH35=1; %usually MRI scan off/on information
biopacfile.dtMRI = 1.7;
biopacfile.dt=0.001; % in seconds/sample
biopacfile.breath=0.01;%new sampling interval for breath it should not be higher than 2.4( for 25 inhales/min)
biopacfile.HR=0.01; %new sampling interval for heart rate, should not be higher than 0.75(for 80 bpm)
% set the range for frequencies in signals in bpm (try to keep those as narrow as posible)
biopacfile.freqBreath=10:25; % in breaths per min
biopacfile.freqHR=40:90; % in beats per minutes
biopacfile.controlplot=0;
       
%%
breath = cell(length(subjects),nrun);
heartrate = breath;
cleanbreath = breath;
clearnheartrate = breath;
for s = 1:length(subjects)
    for r = 1:nrun
        clear refdata
        % get filename
        dirpath = sprintf('%s/%s/run%i/',dataroot,subjects{s},r);
        
        % specify filename for this biopac file
        if ~isempty(dir([dirpath '/*.acq'])) % if .acq file
            tempfn = dir([dirpath '/*.acq']);
            biopacfile.name = fullfile(dirpath,tempfn.name);
            % Loads data and drops unused channels
            refdata=acq2mat(biopacfile);
        elseif ~isempty(dir([dirpath '/*.mat'])) % if .mat file
            tempfn = dir([dirpath '/*.mat' ]);
            name=ones(1,length(tempfn));
            if length(tempfn)>1
                for i=1:length(tempfn)
                    if tempfn(i).bytes < 200000 % euristica idiotica
                        name(i)=0;
                    end
                end
            end
            refdata = load(fullfile(dirpath,tempfn(name==1).name));
            refdata = refdata.refdata;
        end
        
        % Cut out the beginning corresponding to rmvtrials, add default params
        samptoremove = biopacfile.dtMRI*(refdata{1}.dt^-1)*rmvframes;
        refdata{1}.data = refdata{1}.data(samptoremove+1:end,1);
        refdata{2}.data = refdata{2}.data(samptoremove+1:end,1);
        
        % Fill in necessary parameters
        refdata{1}.downdt=biopacfile.breath;
        refdata{2}.downdt=biopacfile.HR;
        refdata{1}.freqlist=biopacfile.freqBreath;
        refdata{2}.freqlist=biopacfile.freqHR;
        refdata{1}.filter = 1;
        refdata{2}.filter = 1;
        refdata{1}.name = 'breath';
        refdata{2}.name = 'heartrate';
        
        % Prepare random data
        T = refdata{1}.dt*(0:size(refdata{1}.data,1)-2);
        data.data = randn(1,numel(T));
        data.dt = refdata{1}.dt;

        % Run DRIFTER, I use my own modified function because plotting is
        % crashing it and I don't need it.
        [~,refdata] = drifter2(data,refdata);

        % save raw downsampled
        breath{s,r} = refdata{1};
        heartrate{s,r} = refdata{2};
        
    end
end
save('biopac_data.mat','breath','heartrate','-v7.3');


