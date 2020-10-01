%% Subsample ERPs to get intensity matched detection ERP 

clear 
close all

rng('shuffle')

EXP = 1;        % 1 = DRT | 2 = MT
mydir = '...';  % Directory containing project folder

useBMS = 1;     % Mark detection effects from BMS

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

condition = 'det_intmatched';
n = 2;
leg = {'Miss' 'Hit'};
  
fig_pos = [1,2,5.6,3.5];
chansel = {'CP4' 'C6' 'CPz'};

xlims = [-50 600];

bms_dir = fullfile(data_dir,'2nd level','BMS','null-int-det-pf-unc-rep-cue');
erp_dir = fullfile(data_dir,'2nd level','ERP');
trg_dir = fullfile(erp_dir,'Plots','det_intmatched');
if ~exist(trg_dir, 'dir')
    mkdir(trg_dir)
end

% Load EEG data dummy
D = spm_eeg_load(fullfile(data_dir, SJs{1},'preprocessed',[prefix SJs{1} '.mat']));

% Load BMS and beta test
if useBMS
    load(fullfile(bms_dir,'BMS_FamXPs.mat'))
    load(fullfile(bms_dir,'beta_test.mat'))
end

channels = D.indchantype('EEG');
nChannels = length(channels);
nSamples = D.nsamples;

%% Get ERPs
ERP_file = fullfile(erp_dir,['ERPs_' condition '.mat']);

if ~exist(ERP_file,'file')
    
    disp(['Computing grand mean ERPs: ' condition])
    
    ERP = nan(nChannels, nSamples, n, numel(SJs));
    trial_samples = cell(2,10,numel(SJs));
    ntrials = nan(numel(SJs),1);
    
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
        
        % Get detection and intensity conditions
        [~,int_idx] = ismember('int',trial_log.labels);
        [~,det_idx] = ismember('det',trial_log.labels);
        int_conds = unique(trial_log.data(int_idx,:));
        det_conds = unique(trial_log.data(det_idx,:));
        
        % Get int x det trial numbers
        n_intdet = nan(length(det_conds),length(int_conds));
        for d = 1:length(det_conds)
            for i = 1:length(int_conds)
                n_intdet(d,i) = sum(trial_log.data(det_idx,:) == det_conds(d) & trial_log.data(int_idx,:) == int_conds(i));
            end
        end
        min_trials = min(n_intdet);
        ntrials(s) = sum(min_trials);
        
        % Sample trials
        trials = nan(nChannels,nSamples,n,ntrials(s));
        idx = 1;
        for i = 1:length(int_conds)
            idx1 = idx;
            idx2 = idx + min_trials(i) - 1;
            for d = 1:length(det_conds)
                filt = find(trial_log.data(det_idx,:) == det_conds(d) & trial_log.data(int_idx,:) == int_conds(i));
                trial_sample = sort(randsample(n_intdet(d,i),min_trials(i)));
                trial_idx = filt(trial_sample);
                trials(:,:,d,idx1:idx2) = D(channels,:,trial_idx);
                trial_samples{d,i,s} = trial_idx;
            end
            idx = idx + min_trials(i);
        end
        
        % Compute ERP
        for d = 1:length(det_conds)
            ERP(:,:,d,s) = mean(trials(:,:,d,:),4);
        end
    end
    disp('Saving...')
    save(ERP_file,'ERP','trial_samples','ntrials');
    disp('Done!')
else
    disp('Loading ERP file...')
    load(ERP_file)
    disp('Done!')
end

%% Get significant time points: detection
if useBMS
    xp_thr = .99;
    beta_thr = 3;
    model_idx = 2;      % detection model
    fam_idx = 2;        % +family
    fam_mods = [2,3,4];
    
    % +fam wins
    fam_filt = squeeze(BMS.xp_fam(fam_idx,:,:) >= xp_thr);
    
    % det wins
    mod_xps = BMS.xp(fam_mods,:,:);
    [~,max_mod] = max(mod_xps);
    mod_filt = squeeze(max_mod == model_idx);
    
    % model beta significant
    mod_betas = beta_test.bf10(fam_mods(model_idx),:,:);
    beta_filt = squeeze(mod_betas >= beta_thr);
    
    % overall significance
    sign_filt = fam_filt & mod_filt & beta_filt;
end

%% Plot 

mean_ERP = squeeze(mean(ERP,4));
sd_ERP = squeeze(std(ERP,[],4));
se_ERP = sd_ERP./sqrt(numel(SJs));

erp_cols = [ .9 .8  0
              0 .5  .9 ];
erp_lines = {':','-'};

offset = .1;
width = .8;
height = .88;

cidx = D.indchannel(chansel);
for c = cidx'
   
    f = figure('units','centimeters','position',fig_pos);
    subplot('position',[offset,offset,width,height])
    hold on
    set(gca,'FontName','Calibri','FontSize',10)
    
    h = nan(n,1);
    min_y = 0;
    max_y = 0;
    for d = 1:n
        x = (D.time*1000)';
        y = mean_ERP(c,:,d)';
        se = se_ERP(c,:,d)';
        fill([x;flipud(x)],[y-se;flipud(y+se)],erp_cols(d,:),'linestyle','none','facealpha',.2);
        h(d) = line(x,y,'Color',erp_cols(d,:),'Linestyle',erp_lines{d},'LineWidth',2);
        max_y = max([max_y max(y(x<=xlims(2))+se(x<=xlims(2)))]);
        min_y = min([min_y min(y(x<=xlims(2))-se(x<=xlims(2)))]);
    end
    
    ylim([min_y max_y])    
    ylims = ylim;
    
    % mark significant time points
    if useBMS
        chan_sign = sign_filt(:,c);
        diff_chan_sign = [0; diff(chan_sign)];
        start_sign = find(diff_chan_sign==1);
        end_sign = find(diff_chan_sign==-1)-1;
        if length(end_sign) < length(start_sign)
            end_sign = [end_sign; length(chan_sign)];
        end
        mark_sign = reshape([start_sign end_sign]',1,length(start_sign)*2);
        if ~isempty(mark_sign)
            sign_time = D.time(mark_sign)*1000;
            vfill(sign_time,[.8 .8 .8],'edgecolor','none','bottom');
        end
    end
    
    xlim(xlims)
    
    set(gca,...
        'xtick',0:100:600,...
        'ticklength',[.015 .025],...
        'tickdir','both',...
        'xminortick','on')  
    
    legend(h,leg);
    name = [D.chanlabels{c} '_det_intmatched.tiff'];
    
    line([0 0],ylims,'color',[0 0 0],'linewidth',1)
    line(xlims,[0 0],'color',[0 0 0],'linewidth',1)
    line(xlims,[ylims(1) ylims(1)],'color',[0 0 0],'linewidth',1)
    line([xlims(1) xlims(1)],[ylims(1) ylims(2)],'color',[0 0 0],'linewidth',1)
    ylim(ylims);
    
    saveas(f,fullfile(trg_dir,name))
end