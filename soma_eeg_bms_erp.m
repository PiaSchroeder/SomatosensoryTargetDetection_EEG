%% Sensor level BMS
% On sensor level log evidence (from bayesglm_sensorlevel.m)

clear 
close all

EXP = 1;        % 1 = DRT | 2 = MT
mydir = '...';  % Directory containing project folder

%% Directories etc.

if EXP == 1
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S21' 'S22' 'S23' };
    data_dir = [mydir '\SomA_EEG\DRT\data'];   
    
elseif EXP == 2
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S09' 'S10' 'S11' 'S12' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S24' 'S25' 'S26' 'S27'};
    data_dir = [mydir '\SomA_EEG\MT\data'];  
end

models = {'null', 'int', 'det', 'pf', 'unc', 'rep', 'cue'}; % [ null | int | det | pf | unc | rep | cue ]  (must be label in trial_log or null)
partition = [1 2 2 2 3 4 5];

mod_names = strcat(models, '-');
mod_names = strcat(mod_names{:});
mod_names = mod_names(1:end-1);

glm_dir = 'bayesglm';

trg_dir = fullfile(data_dir, '2nd level', 'BMS', mod_names);
if ~exist(trg_dir, 'dir')
    mkdir(trg_dir)
end

%% Get dummy EEG file to retrieve channels and timing
results_file = fullfile(data_dir, SJs{1}, glm_dir, ['bayesglm_results_' SJs{1} '_' models{1}] );
load(results_file);
channels = results.sensors;
time = results.time;
nTimepoints = length(time);
nSensors = length(channels);

%% Assemble data

disp('----------------------- BMS -----------------------')
disp(mod_names)

disp('Retrieving data...')
lme = nan(numel(SJs),numel(models),nTimepoints,nSensors);

for s = 1:numel(SJs)
    for m = 1:numel(models)
        
        load(fullfile(data_dir, SJs{s},glm_dir,['bayesglm_results_' SJs{s} '_' models{m}]));
        
        lme(s,m,:,:) = results.LogEv;
        
        clear results
    end
end
fprintf('Done!\n')

%% Run BMS per sensor and time point

fprintf('Running BMS\n')

exp_r       = nan(numel(models), nTimepoints, nSensors);
xp          = nan(numel(models), nTimepoints, nSensors);
exp_r_fam   = nan(numel(unique(partition)), nTimepoints, nSensors);
xp_fam      = nan(numel(unique(partition)), nTimepoints, nSensors);

% Set up BMS structure
BMS.models      = models;
BMS.SJs         = SJs;
BMS.channels    = channels;
BMS.time        = time;
BMS.partition   = partition;
BMS.exp_r       = exp_r;
BMS.xp          = xp;
BMS.exp_r_fam   = exp_r_fam;
BMS.xp_fam      = xp_fam;

for sensor = 1:nSensors
    
    fprintf('Sensor %d/%d ',sensor,nSensors)
    
    % Compute posterior model probabilities
    nbytes = fprintf('Time point 0/%d ',nTimepoints);
    for t = 1:nTimepoints
        while nbytes > 0
            fprintf('\b');
            nbytes = nbytes - 1;
        end
        nbytes = fprintf('Time point %d/%d',t,nTimepoints);
        
        [exp_r(:,t,sensor),exp_r_fam(:,t,sensor),xp(:,t,sensor),xp_fam(:,t,sensor)] = spm_BMS_family(lme(:,:,t,sensor),[],partition);
        
    end
    
    % Update BMS structure and save
    BMS.exp_r       = exp_r;
    BMS.xp          = xp;
    BMS.exp_r_fam   = exp_r_fam;
    BMS.xp_fam      = xp_fam;
    
    save(fullfile(trg_dir, 'BMS_FamXPs.mat'),'BMS');
    fprintf('\n')
    
    fprintf('Sensor %d/%d done!\n',sensor,nSensors)
end



