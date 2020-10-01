%% Plot PFs 
% Prior to any data analysis, mean psychometric functions were plotted and
% data exclusion criteria applied (min P(det) > 10% || max P(det) < 90%)

clear 
close all

EXP = 1;        % 1 = DRT | 2 = MT
mydir = '...';  % Directory containing project folder

%% 1. Load data
% =========================================================================
if EXP == 1
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' };
    data_dir = [mydir '\SomA_EEG\DRT\data'];   
    
elseif EXP == 2
    SJs  = { 'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S27' 'S28' };
    data_dir = [mydir '\SomA_EEG\MT\data'];  
end

runs = repmat({1:6},1,numel(SJs));

ana_dir  = [mydir '\SomA_EEG\analysis\behaviour'];
trg_dir  = fullfile(data_dir,'2nd level','Behaviour');

nSub = numel(SJs);
nRuns = 8;
nInts = 10;
nTrials = 200;

cd(ana_dir)
[Data, PFs] = load_logs(SJs,nRuns,runs,data_dir,ana_dir);

%% 2. Get psychometric functions
% =========================================================================

logistic = @(c,x) (1./(1+exp(-c(2)*(x-c(1)))));

% Normalise PFs to intensity range 1-10
iniT = 5.5;
ints = cell(nSub,nRuns);
resp = cell(nSub,nRuns);
norm_PFs = cell(nSub,nRuns);
norm_slopes = nan(nSub,nRuns);
norm_threshs = nan(nSub,nRuns);

for s = 1:nSub
    disp(SJs{s})
    for r = runs{s}
        [~,ints{s,r}] = ismember(Data{s,r}.behaviour.PF.Intensities,Data{s,r}.Exp.Intensities); % Get intensity levels
        resp{s,r} = Data{s,r}.behaviour.PF.Responses;                                           % Get responses 
        norm_PFs{s,r} = fit_logistic(iniT,ints{s,r},resp{s,r},SJs{s},0);                        % Fit logistic function --> normalised PF
        norm_slopes(s,r) = norm_PFs{s,r}.fit_logistic(1,2);                                     % Get normalised slope
        norm_threshs(s,r) = norm_PFs{s,r}.fit_logistic(1,1);                                    % Get normalised T50
    end
end

mean_norm_slopes = nanmean(norm_slopes,2);
mean_norm_threshs = nanmean(norm_threshs,2);

normPFs.norm_PFs = norm_PFs;
normPFs.norm_slopes = norm_slopes;
normPFs.norm_threshs = norm_threshs;
normPFs.mean_norm_slopes = mean_norm_slopes;
normPFs.mean_norm_threshs = mean_norm_threshs;

% Fit mean PF for every participant
plot_range = 1:nInts;
for s = 1:nSub
    mean_fit(s,:) = logistic([mean_norm_threshs(s), mean_norm_slopes(s)],plot_range);
end

normPFs.mean_fit = mean_fit;

% Apply exclusion criterion: mean fitted detection probability min>10% or max<90%
excl = round(mean_fit(:,1)*100) > 10 | round(mean_fit(:,end)*100) < 90;

%% 3. Plot PFs
fig_pos = [1,2,4.5,4];
plot_range = 0:0.01:nInts;
f = figure('Units','centimeters','position',fig_pos); 
hold on
set(gca,'FontName','calibri','FontSize',12)

area([0 10],[1 1],'FaceColor',[1 1 1], 'LineStyle', 'none')
area([0 10],[.9 .9],'FaceColor',[.7 .7 .7], 'LineStyle', 'none')
area([0 10],[.1 .1],'FaceColor',[1 1 1], 'LineStyle', 'none')

for s = 1:nSub
    sub_pf = logistic([mean_norm_threshs(s), mean_norm_slopes(s)],plot_range);
    if excl(s)
        fprintf('%s min: %f max: %f\n',SJs{s}, mean_fit(s,1), mean_fit(s,10));
        line(plot_range,sub_pf,'Color',[.8 0 .1],'LineWidth',1,'LineStyle','--');
    elseif EXP == 2 && strcmp(SJs{s},'S23') % Bad EEG data quality
        line(plot_range,sub_pf,'Color',[.8 0 .1],'LineWidth',1,'LineStyle',':');
    else
        line(plot_range,sub_pf,'Color',[0 0 0],'LineWidth',1);
    end
end
axis([1 10 0 1]);
line([1 1],[0 1],'color',[0 0 0])
line([1 10],[0 0],'color',[0 0 0])
ylabel('Detection probability');
xlabel('Intensity level');
set(gca,...
    'xtick',1:10,...
    'ytick',[0,.5,1],...
    'ticklength',[.015 .025],...
    'tickdir','both',...
    'yminortick','on')

saveas(f,fullfile(trg_dir,'PFs.tiff'))

