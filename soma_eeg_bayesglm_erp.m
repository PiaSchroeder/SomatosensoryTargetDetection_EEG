%% Run Bayesian GLM estimation on electrode data

clear 
close all

disp('Initialising...')
spm('defaults','eeg')
spm_jobman('initcfg');
spm_get_defaults('cmdline',true)

EXP = 1;        % 1 = DRT | 2 = MT
mydir = '...';  % Directory containing project folder

t_start = tic;
info_str = ['BayesGLM on all SJs\n' ...
           'Models: null, int, det, pf, unc, rep, cue \n'];
       
%% Directories etc.

if EXP == 1
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S21' 'S22' 'S23' };
    data_dir = [mydir '\SomA_EEG\DRT\data'];   
    prefix = 'bfraeTMdfsoma_DRT_';
    log_f = 'trial_log_DRT_';
elseif EXP == 2
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S09' 'S10' 'S11' 'S12' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S24' 'S25' 'S26' 'S27'};
    data_dir = [mydir '\SomA_EEG\MT\data'];  
    prefix = 'bfraeTMdfsoma_MT_';
    log_f = 'trial_log_MT_';
end

models = {'null', 'int', 'det', 'pf', 'unc', 'rep', 'cue'}; % [ null | int | det | pf | unc | rep ]  (must be label in trial_log or null)

chantype = 'EEG';

src_folder  = 'preprocessed';
log_folder  = 'logs';
trg_folder  = 'bayesglm';

%% Loop through subjects and models
for s = 1:numel(SJs)
    
    sj_eeg_file = fullfile(data_dir, SJs{s}, src_folder, [prefix,SJs{s},'.mat']);
    sj_trial_log = fullfile(data_dir, SJs{s}, log_folder, [log_f SJs{s} '.mat']);
    
    D = spm_eeg_load(sj_eeg_file);
    load(sj_trial_log);
    
    for m = 1:numel(models)
        
        fprintf('%s %s ', SJs{s}, models{m});
        
        % Get regressor
        if strcmp(models{m},'null')
            reg = [];
        else
            model_idx = strcmp(trial_log.labels, models{m});
            reg = zscore(trial_log.data(model_idx,:)');
        end
        
        S = [];
        S.D = D;
        S.chantype = chantype;
        S.reg = reg;
        S.alpha = 200;
        
        results = bayesglm_sensors(S);
        
        results.sj = SJs{s};
        results.model = models{m};
        results.sensors = D.chanlabels(D.indchantype('EEG'));
        results.time = D.time;
    
        sj_trg_dir = fullfile(data_dir, SJs{s}, trg_folder);
        if ~exist(sj_trg_dir,'dir')
            mkdir(sj_trg_dir)
        end
        save(fullfile(sj_trg_dir, ['bayesglm_results_' SJs{s} '_' models{m} '.mat']),'results')
        
    end
    
end

%%
t_elapsed = toc(t_start)/60;
fprintf('\n---------------------------- DONE! -----------------------------\n')
fprintf(info_str)
fprintf('Elapsed time: %.1f min\n',t_elapsed)
