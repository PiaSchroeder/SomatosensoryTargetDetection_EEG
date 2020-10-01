%% Plot BMS results as topoplot

clear 
close all

EXP = 1;        % 1 = DRT | 2 = MT
mydir = '...';  % Directory containing project folder

%% Directories etc.

if EXP == 1
    data_dir = [mydir '\SomA_EEG\DRT\data']; 
    prefix = 'bfraeTMdfsoma_DRT_';
elseif EXP == 2
    data_dir = [mydir '\SomA_EEG\MT\data'];  
    prefix = 'bfraeTMdfsoma_MT_';
end

topo_template = fullfile(data_dir,'templates','topo.set');  % Use template EEGLAB data set to get channel topography

bms_dir = fullfile(data_dir,'2nd level','BMS','null-int-det-pf-unc-rep-cue');
trg_dir = fullfile(bms_dir,'topos');
if ~exist(trg_dir, 'dir')
    mkdir(trg_dir)
end

models      = {'null' 'int', 'det', 'pf', 'unc', 'rep', 'cue'};
isfam       = [0 1 1 1 0 0 0];
plus_fam    = 2;

name    = 'BMS_topo';
THRESH  = 1;

% name    = 'BMS_topo_unthresholded';
% THRESH  = 0;

% BMS file and beta_test
bms_file = fullfile(bms_dir, 'BMS_FamXPs.mat');
load(bms_file)
beta_test_file = fullfile(bms_dir, 'beta_test.mat');
load(beta_test_file)

% dummy data sets
D = spm_eeg_load(fullfile(data_dir,BMS.SJs{1},'preprocessed',[prefix BMS.SJs{1} '.mat']));
EEG = pop_loadset(topo_template);

xp_theshold = .99;
beta_threshold = 3;

%% Prepare plots      

fam_cols = [.5 .5 .5   % null - grey
            .8 .8 .8   % +fam - light grey
             0  1  1   % unc - cyan
             1  0  1   % rep - magenta
             1  1  0   % cue - yellow
             ];

% Family info
partition   = BMS.partition;
nFam        = length(unique(partition)); 
fam_idx     = cell(1,nFam);
fam_size    = nan(1,nFam);
for i = 1:nFam
    fam_idx{i} = find(partition == i);
    fam_size(i) = length(fam_idx{i});
end

fam_xps = BMS.xp_fam;
xps = BMS.xp(isfam==1,:,:);
xps = xps.*repmat(fam_xps(plus_fam,:,:),fam_size(plus_fam),1,1);

%% Plot
toi = [0 0.05 0.09 0.12 0.15 0.25 0.28 0.35 0.5];
plotsamples = D.indsample(toi);

for s = 1:length(plotsamples)
    
    fig = figure;
    
    ps = plotsamples(s);
    
    bf10 = nan(64,1);
    tp_chan_cols = cell(64,1);
    for c = 1:64
        [~,max_fam] = max(fam_xps(:,ps,c),[],1);  % Determine winning family
        
        if max_fam == plus_fam
            [~,max_mod] = max(xps(:,ps,c),[],1);
            bf10(c) = beta_test.bf10(fam_idx{max_fam}(max_mod),ps,c);
            if fam_size(fam_size == 3)
                tp_chan_cols{c} = xps([2,1,3],ps,c)';
            else
                tp_chan_cols{c} = [xps(2,ps,c) xps(1,ps,c) 0];
            end
        else
            bf10(c) = beta_test.bf10(fam_idx{max_fam},ps,c);
            tp_chan_cols{c} = fam_cols(max_fam,:)*fam_xps(max_fam,ps,c);
        end
    end
    
    if THRESH
        data = any(squeeze(fam_xps(:,ps,:))' >= xp_theshold,2) & bf10 >= beta_threshold;
    else
        data = ones(64,1);
    end
    
    topoplot(data,EEG.chanlocs,'style','blank','plotdisk','on','efontsize',20,'emarkercolors',tp_chan_cols);
%     title(sprintf('%.0f ms',toi(s)*1000),'position',[-0.5 0.5],'FontSize',18','FontWeight','b')
    set(gcf,'Units','centimeters');
    set(gcf,'Units','centimeters','position',[0,1,20,16])
    
    fig_name = sprintf('%s_%.0fms.tiff',name,toi(s)*1000);
    saveas(fig,fullfile(trg_dir, fig_name)); 
end
