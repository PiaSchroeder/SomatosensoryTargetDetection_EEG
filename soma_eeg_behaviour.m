%% SomA EEG behaviour
% Detection rates
% Reaction times
% Response associations: Detection - Match reports

clear 
close all

EXP = 1;        % 1 = DRT | 2 = MT
mydir = '...';  % Directory containing project folder

%% Directories etc

if EXP == 1
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S21' 'S22' 'S23' };
    data_dir = [mydir '\SomA_EEG\DRT\data'];   
    log_f = 'trial_log_DRT_';
elseif EXP == 2
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S09' 'S10' 'S11' 'S12' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S24' 'S25' 'S26' 'S27'};
    data_dir = [mydir '\SomA_EEG\MT\data'];  
    log_f = 'trial_log_MT_';
end

trg_dir = fullfile(data_dir,'2nd level','Behaviour');
log_dir = 'logs';

%% Detection rates
det_rates = nan(numel(SJs),1);

for s = 1:numel(SJs)
    
    sj_trl_log = fullfile(data_dir, SJs{s}, log_dir, [log_f SJs{s} '.mat']);
    load(sj_trl_log)
    
    [~,det_idx] = ismember('det',trial_log.labels);
    
    det_rates(s) = sum(trial_log.data(det_idx,:))/size(trial_log.data,2)*100;
    
end

mean_det_rates = mean(det_rates);
sd_det_rates = std(det_rates);

DetRates.all = det_rates;
DetRates.mean = mean_det_rates;
DetRates.sd = sd_det_rates;

save(fullfile(trg_dir, 'DetRates.mat'),'DetRates')

%% Reaction times
RTs.all = [];
RTs.sj_mean = nan(1,numel(SJs));
RTs.group_mean = [];
RTs.group_sd = [];

for s = 1:numel(SJs)
    
    sj_trl_log = fullfile(data_dir, SJs{s}, log_dir, [log_f SJs{s} '.mat']);
    load(sj_trl_log)
    
    [~,rt_idx] = ismember('rts',trial_log.labels);
    rt_data = trial_log.data(rt_idx,:)*1000;
    
    RTs.all = [RTs.all rt_data];
    RTs.sj_mean(s) = mean(rt_data);
    
    [~,cidx] = ismember('det',trial_log.labels);
    c_data = trial_log.data(cidx,:);
    cs = unique(c_data);
    
    for cc = 1:length(cs)
        RTs.det.sj_mean(s,cc) = mean(rt_data(c_data==cs(cc)));
    end
end

RTs.group_mean = mean(RTs.sj_mean);
RTs.group_sd = std(RTs.sj_mean);

RTs.det.labels = {'miss' 'hit'};
RTs.det.group_mean = mean(RTs.det.sj_mean);
RTs.det.group_sd = std(RTs.det.sj_mean);
RTs.det.BF10 = bf_ttest(RTs.det.sj_mean(:,1),RTs.det.sj_mean(:,2)); % Uses the ttest.m function from the bayesFactor toolbox (Bart Krekelberg (2021).
                                                                    % bayesFactor (https://github.com/klabhub/bayesFactor)), renamed to bf_ttest.m
                                                                    % to avoid naming conflicts  

save(fullfile(trg_dir,'RTs.mat'),'RTs')

%% Response associations: Detection (det) - Match reports (rep)
BF10_det_rep = nan(numel(SJs),1);

for s = 1:numel(SJs)
    
    sj_trl_log = fullfile(data_dir, SJs{s}, log_dir, [log_f SJs{s} '.mat']);
    load(sj_trl_log)
    
    [~,det_idx] = ismember('det',trial_log.labels);
    [~,rep_idx] = ismember('rep',trial_log.labels);
    
    dets = trial_log.data(det_idx,:);
    reps = trial_log.data(rep_idx,:);

    det_rep(1,1) = sum(dets==0 & reps==0);
    det_rep(1,2) = sum(dets==0 & reps==1);
    det_rep(2,1) = sum(dets==1 & reps==0);
    det_rep(2,2) = sum(dets==1 & reps==1);
    BF10_det_rep(s) = c_table(det_rep);     % Uses the c_table.m function from the companion software to Johnson, V. E., & Albert, J. H. 
                                            % (2006). Ordinal data modeling. Springer Science & Business Media. Available at: 
                                            % https://de.mathworks.com/matlabcentral/fileexchange/2264-ordinal-data-modeling 
    
end

BF01_det_rep = 1./BF10_det_rep;

Det_Rep.BF10 = BF10_det_rep;
Det_Rep.BF01 = BF01_det_rep;

save(fullfile(trg_dir, 'ResponseAssociations.mat'),'Det_Rep')


