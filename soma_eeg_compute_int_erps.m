%% Compute ERPs

clear 
close all

EXP = 1;        % 1 = DRT | 2 = MT
mydir = '...';  % Directory containing project folder

%% Directories etc.

if EXP == 1
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' };
    data_dir = [mydir '\SomA_EEG\DRT\data'];   
    prefix = 'bfraeTMdfsoma_DRT_';
elseif EXP == 2
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S27' 'S28' };
    data_dir = [mydir '\SomA_EEG\MT\data'];  
    prefix = 'bfraeTMdfsoma_MT_';
end

condition = 'int';
n = 10;
leg = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'};

trg_dir = fullfile(data_dir,'2nd level','ERP');

% Load EEG data dummy
D = spm_eeg_load(fullfile(data_dir, SJs{1},'preprocessed',[prefix SJs{1} '.mat']));

channels = D.indchantype('EEG');
nChannels = length(channels);
nSamples = D.nsamples;

%% Get ERPs
ERP_file = fullfile(trg_dir,['ERPs_' condition '.mat']);
disp(['Computing grand mean ERPs: ' condition])

ERP = nan(nChannels, nSamples, n, numel(SJs));

for s = 1:numel(SJs)
    
    disp(SJs{s})
    sj_data_dir = fullfile(data_dir, SJs{s},'preprocessed');
    
    % Get EEG data
    D = spm_eeg_load(fullfile(sj_data_dir,[prefix SJs{s} '.mat']));
    
    % Get trial definitions
    trlog = fullfile(data_dir, SJs{s}, 'logs', [SJs{s}, '_trial_log.mat']);
    load(trlog)
    
    % make sure trial numbers match
    if D.ntrials ~= size(trial_log.data,2)
        error(['Trial numbers do not match!!! ' SJs{s}])
    end
    
    % Get trial conditions
    [~,idx] = ismember(condition,trial_log.labels);
    conds = unique(trial_log.data(idx,:));
    
    for cc = conds
        filt = find(trial_log.data(idx,:) == cc);
        cond_data = D(channels,:,filt);
        ERP(:,:,cc,s) = mean(cond_data,3);
    end
end
save(ERP_file,'ERP');

