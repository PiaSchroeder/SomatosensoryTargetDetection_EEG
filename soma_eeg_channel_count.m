%% Plot temporal evolution of model fits across the scalp
% For each model: significant electrodes/all electrodes

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

bms_dir = fullfile(data_dir,'2nd level','BMS','null-int-det-pf-unc-rep-cue');
trg_dir = fullfile(bms_dir,'channel count');
if ~exist(trg_dir, 'dir')
    mkdir(trg_dir)
end

load(fullfile(bms_dir,'BMS_FamXPs.mat'))
load(fullfile(bms_dir,'beta_test.mat'))

nChannels = length(BMS.channels);
nSamples = length(BMS.time);

xp_threshold = .99;
beta_threshold = 3;

uselims = [-50 600];
xticks = 0:100:600;
yticks = 0:10:100;
plot_mods = [1 5 6 7 2 4 3]; % plot models in this order

%% Family info
partition   = BMS.partition;
nFam        = length(unique(partition)); 
fam_idx     = cell(1,nFam);
fam_size    = nan(1,nFam);
for i = 1:nFam
    fam_idx{i} = find(partition == i);
    fam_size(i) = length(fam_idx{i});
end

%% Count channels
XPs = BMS.xp(fam_idx{fam_size > 1},:,:);
FamXPs = BMS.xp_fam;
model_count = nan(7,nSamples);
[~,max_fams] = max(FamXPs,[],1);

idx = 1;
for f = 1:nFam
    
    if fam_size(f) > 1
        [~,max_mods] = max(XPs,[],1);
        beta = nan(size(max_mods));
        for s = 1:nSamples
            for c = 1:nChannels
                beta(1,s,c) = beta_test.bf10(fam_idx{f}(max_mods(1,s,c)),s,c);
            end
        end
        %     winning family & XP threshold exceedance      & beta threshold exceedance 
        filt = max_fams == f & FamXPs(f,:,:) >= xp_threshold & beta >= beta_threshold;    
        max_mods(~filt) = 0;
        for ff = 1:fam_size(f)
            model_count(idx,:) = sum(max_mods == ff,3)./nChannels;
            idx = idx+1;
        end
    else
        beta = beta_test.bf10(fam_idx{f},:,:);
        filt = max_fams == f & FamXPs(f,:,:) >= xp_threshold & beta >= beta_threshold;
        model_count(idx,:) = sum(filt,3)./nChannels;
        idx = idx+1;
    end
end
 
%% Plot

cols = [.7 .7 .7
         0  1  0
         1  0  0
         0  0  1
         0  1  1
         1  0  1
         1  1  0];
  
fig_pos = [1 2 16 4];
line_width = 2;
subplot_pos = [.05 .1 .9 .8];
fig = figure('units','centimeters','position',fig_pos);
subplot('Position',subplot_pos)
hold on
set(gca,'FontName','Calibri','FontSize',12,'ytick',yticks,'xtick',xticks,'xminortick','on','yminortick','on')

for f = plot_mods
    plot(BMS.time*1000,model_count(f,:)*100,'color',cols(f,:),'linewidth',line_width)
end

xlabel('time (ms)')
ylabel('% of channels')
xlim(uselims)
ylim([0 45])

line([0 0],ylim,'color',[0 0 0],'linewidth',1)
line([-50 -50],ylim,'color',[0 0 0],'linewidth',1)
line(xlim,[0 0],'color',[0 0 0],'linewidth',1)
    
saveas(fig,fullfile(trg_dir, 'BMS_channel_count.tiff'))    
    