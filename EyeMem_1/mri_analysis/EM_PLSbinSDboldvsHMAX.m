% function accuracy = EM_PLSbinSDboldvsHMAX(cfg)
% Corr SDBOLD and HMAX across trials for each subject, then
% 1) Adaptation to feature richness: task PLS OA vs YA: positive link YA, weaker OA?
% 2) FR adaptation vs memory: behav PLS: positive link YA, weaker OA? 

mri = ft_read_mri('/Users/kloosterman/projectdata/eyemem/Standards/MNI152_T1_3mm_brain.nii.gz');
cfg=[]; ft_sourceplot(cfg, mri); %, mri_mask

mri_mask = ft_read_mri('/Users/kloosterman/projectdata/eyemem/Standards/MNI152_T1_3mm_brain_mask.nii.gz');
% mri_mask.anatomy = round(mri_mask.anatomy); 
mri_mask.anatomy(mri_mask.anatomy < 1) = 0; % binarize
% cfg=[]; ft_sourceplot(cfg, mri_mask); %, mri_mask

%% get ROIs
region_of_interest = 'Hippocampus' % 'Hippocampus' 'Visual cortex'
region_of_interest = 'Visual cortex' % 'Hippocampus' 'Visual cortex'
% region_of_interest = 'V1' % 'Hippocampus' 'Visual cortex'
% read Juelich atlas, select ROIs
atlas = ft_read_atlas('/Users/kloosterman/fsl/data/atlases/Juelich.xml');
cfg=[];
cfg.funparameter = 'tissue';
cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, atlas, mri);
% roiInds = find(contains(atlas.tissuelabel, 'Visual cortex')); % get roi Indices
roiInds = find(contains(atlas.tissuelabel, region_of_interest));
atlas.tissuelabel(roiInds)
mri_mask = atlas;
mri_mask.tissue = ismember( mri_mask.tissue, roiInds);
% resample mask to 3 mm
cfg=[];
cfg.parameter = 'tissue';
mri_mask3mm = ft_sourceinterpolate(cfg, mri_mask, mri);

cfg=[];
cfg.funparameter = 'tissue';
cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, mri_mask3mm, mri);

%%
%%
%% load all subjects, compute SD over time, correlate to HMAX, append
rmpath('/Users/kloosterman/Documents/GitHub/plscmd')
datapath = '/Users/kloosterman/projectdata/eyemem/variability2/5TRspertrial/ftsource';
datapath = '/Users/kloosterman/projectdata/eyemem/variability2/1TRspertrial/ftsource';

nbins = 3;
sources = [];
sources.young = cell(1,nbins); % 
sources.old = cell(1,nbins); % nbins
cd(datapath)
subjects = dir(fullfile(datapath, 'source_sub*'));
labels = [];
load participantinfo.mat % TODO make this reliable

stat_mse = {}; corrlist = []; inside = []; hmax_dat = [];
for isub = 1:length(subjects) % TODO leave out subject under investigation
  subj = subjects(isub).name(8:end-4);
  subjinfo = Participants(Participants.participant_id == subj, :);     % give different outfolder for OA and YA
  disp(subj)
  % get eye data for labels
  eyedat = load(['eye' subjects(isub).name(7:end)]);
  look_region_dat = eyedat.trialinfo.HMAX_fix_lookregion_mean; % gives
  % Bin2 > 1
 
  % look_region_dat = eyedat.trialinfo.HMAX_fix; % HMAX_fix: HMAX at fixation location
  % look_region_dat = eyedat.trialinfo.HMAX_C1; % HMAX_C1: mean HMAX over whole picture

  source = load(subjects(isub).name);

  % source.pow = zscore(source.pow,0,2); %normalize across trials within subject!
  source.pow = normalize(source.pow, 2, "range"); % [0 1] normalization

  % source.inside = logical(mri_mask.anatomy(:)); 
  inside(:,end+1) = source.inside;
  source.freq = 1;

  % Define bin edges using quantiles (equal-sized bins)
  edges = quantile(look_region_dat, 1/nbins:1/nbins:(nbins-1)/nbins);
  edges = [-Inf, edges, Inf]; % Extend edges to cover all values
  bin_idx = discretize(look_region_dat, edges); % Assign bin indices based on edges
  SDbold = [];
  for ibin=1:nbins
    trialdat = source.pow(:,bin_idx == ibin,:);
    source_bin = source;
    source_bin.pow = squeeze(std(trialdat(:,:), 0, 2));% compute SD over time;
    source_bin.powdimord = 'pos';
    sources.(string(subjinfo.group)){ibin}{end+1} = source_bin; %TODO change order of dims
    hmax_dat(isub,ibin) = mean(look_region_dat(bin_idx == ibin,1));
  end

  plotit=0;
  if plotit
    tmp2 = source;
    % tmp2.pow = tmp2.pow(:,itrial);
    cfg=[];
    cfg.funparameter = 'pow';
    cfg.method = 'ortho'; % slice ortho glassbrain vertex
    load colormap_jetlightgray.mat
    cfg.funcolormap = cmap;
    cfg.funcolorlim = 'maxabs';
    cfg.location = [25 11 30];
    cfg.locationcoordinates = 'voxel';
    ft_sourceplot(cfg, tmp2, mri_mask); %, mri_mask
  end
end
% figure; histogram(corrlist,10)
figure; boxplot(hmax_dat); ylabel('Saliency')

cfg=[];
cfg.keepindividual = 'yes';
agegroups = fieldnames(sources);
for iage = 1:length(agegroups)
  for ibin = 1:nbins
    C = sources.(agegroups{iage}){ibin};
    C = cellfun(@(s) setfield(s, 'inside', all(inside,2)), C, 'UniformOutput', false);
    sources.(agegroups{iage}){ibin} = ft_sourcegrandaverage(cfg, C{:});
  end
end
% remove standard brain parts not measured
mri.anatomy(not(all(inside,2))) = 0;

plotit=0;
if plotit
  allsourceplot = sources.(agegroups{iage}){ibin};
  allsourceplot.pow = mean(allsourceplot.pow)';
  allsourceplot.dimord = 'pos';
  cfg=[];
  cfg.funparameter = 'pow';
  cfg.method = 'ortho'; % slice ortho glassbrain vertex
  load colormap_jetlightgray.mat
  cfg.funcolormap = cmap;
  cfg.funcolorlim = 'maxabs';
  % cfg.funcolorlim = 'zeromax';
  cfg.location = [25 11 30];
  cfg.locationcoordinates = 'voxel';
  ft_sourceplot(cfg, allsourceplot); %, mri_mask
end
addpath('/Users/kloosterman/Documents/GitHub/plscmd')

%%
%% Task PLS 2 group, cond (nbins)
cfg = [];
cfg.parameter = 'pow';
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 10;                         % Number of permutation
cfg.num_boot = 100;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 1;                         % 1 is taskPLS; 3 is behavPLS
cfg.num_cond = nbins;                       % Number of conditions
cfg.num_subj_lst = [size(sources.young{1}.pow,1) size(sources.old{1}.pow,1)];                 % Number of subjects per condition
cfg.design = zeros(sum(cfg.num_subj_lst)*nbins,1); % Placeholder TODO omit design for Task PLS
stat_taskPLS = ft_sourcestatistics(cfg, sources.young{:}, sources.old{:});
stat_taskPLS.results.perm_result.sprob

%% Task PLS 2 group, cond (nbins): plot Brainscores, bootstrapratios
% bardat = -reshape(stat_taskPLS.brainscores, [], 6);
% bardat = 
% figure; bar(mean(bardat));
% figure; boxplot(bardat);
figure; tiledlayout(1,2)
for i=1:2
  nexttile;  boxplot(stat_taskPLS.brainscores{i})
  % nexttile;  bar(mean(stat_taskPLS.brainscores{i}))
end

stat_taskPLS.bootstrapratios = zeros(size(stat_taskPLS.inside));
stat_taskPLS.bootstrapratios(stat_taskPLS.inside,1) = stat_taskPLS.results.boot_result.compare_u(:,1);
stat_taskPLS.mask = (stat_taskPLS.bootstrapratios < -2 | stat_taskPLS.bootstrapratios > 2);
% stat_taskPLS.bootstrapratios_mask(stat_taskPLS.bootstrapratios > -2 & stat_taskPLS.bootstrapratios < 2) = 0;

stat_taskPLS.bootstrapratios_masked = stat_taskPLS.bootstrapratios;
stat_taskPLS.bootstrapratios_masked(~(stat_taskPLS.mask)) = 0;

cfg=[];
cfg.funparameter = 'bootstrapratios_masked';
cfg.maskparameter = 'mask';
cfg.anaparameter = 'anatomy';
cfg.method = 'ortho'; % slice ortho glassbrain vertex
load colormap_jetlightgray.mat
cfg.funcolormap = cmap;
cfg.funcolorlim = 'maxabs';
cfg.funcolorlim = [-6 6];
% cfg.funcolorlim = 'zeromax';
cfg.location = [25 11 30];
cfg.locationcoordinates = 'voxel';
ft_sourceplot(cfg, stat_taskPLS, mri); %, mri_mask , mri

%% Task PLS 1 group , nbins cond (nbins)
%%
addpath('/Users/kloosterman/Documents/GitHub/plscmd')
cfg = [];
cfg.parameter = 'pow';
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 10;                         % Number of permutation
cfg.num_boot = 10;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 1;                         % 1 is taskPLS; 3 is behavPLS
cfg.num_cond = nbins;                       % Number of conditions
cfg.num_subj_lst = 39;                 % Number of subjects per condition
cfg.design = zeros(sum(cfg.num_subj_lst)*nbins,1); % Placeholder TODO omit design for Task PLS
stat_taskPLS_YA = ft_sourcestatistics(cfg, sources.young{:});
stat_taskPLS_YA.results.perm_result.sprob

%% plot Task PLS 1 group , nbins cond
bardat = -reshape(stat_taskPLS_YA.brainscores, [], nbins);  

figure; tiledlayout(1,2); 
nexttile; boxplot(bardat);
% nexttile; bar(mean(bardat));
% title('Overall brainscore')

% get ROI-specific brain scores
stat_taskPLS_YA.brainweights = zeros(size(stat_taskPLS_YA.inside));
stat_taskPLS_YA.brainweights(stat_taskPLS_YA.inside,1) = -stat_taskPLS_YA.results.u(:,1);
HCvoxels = logical(mri_mask3mm.tissue(:));
stat_taskPLS_YA.brainscores_HC = [];
for ibin=1:nbins
  stat_taskPLS_YA.brainscores_HC(:,ibin) = sum(stat_taskPLS_YA.brainweights(HCvoxels)' .* sources.young{ibin}.pow(:,HCvoxels),2);
end
nexttile; boxplot(stat_taskPLS_YA.brainscores_HC); 
title(region_of_interest)

%% plot bootstrapratios
stat_taskPLS_YA.bootstrapratios = zeros(size(stat_taskPLS_YA.inside));
stat_taskPLS_YA.bootstrapratios(stat_taskPLS_YA.inside,1) = stat_taskPLS_YA.results.boot_result.compare_u(:,1);
stat_taskPLS_YA.bootstrapratios(stat_taskPLS_YA.bootstrapratios > -2 & stat_taskPLS_YA.bootstrapratios < 2) = 0;

stat_taskPLS_YA.mask = double(stat_taskPLS_YA.bootstrapratios < -2 | stat_taskPLS_YA.bootstrapratios > 2);
mri.mask = stat_taskPLS_YA.mask;

stat_taskPLS_YA.bootstrapratios_masked = stat_taskPLS_YA.bootstrapratios;
stat_taskPLS_YA.bootstrapratios_masked(~(stat_taskPLS_YA.mask)) = NaN;

cfg=[];
cfg.funparameter = 'bootstrapratios';
cfg.maskparameter = 'mask';
cfg.anaparameter = 'anatomy';
cfg.method = 'ortho'; % slice ortho glassbrain vertex
load colormap_jetlightgray.mat
cfg.funcolormap = cmap;
cfg.funcolorlim = 'maxabs';
cfg.funcolorlim = [-4 4];
cfg.location = [25 11 30];
cfg.locationcoordinates = 'voxel';
ft_sourceplot(cfg, stat_taskPLS_YA, mri); %, mri_mask , mri










%% CBPM corr vs 0
Nsubj = size(allsource.pow,1);
design = ones(1, Nsubj); % All subjects in one group

cfg = [];
cfg.dim   = allsource.dim;
cfg.parameter   = 'pow';  % Specify the parameter for testing
cfg.method      = 'montecarlo'; % Use Monte Carlo permutation method
cfg.statistic   = 'ft_statfun_depsamplesT'; % INDependent samples t-test (against zero)
cfg.correctm    = 'cluster'; % Use cluster-based correction
% cfg.correctm    = 'no'; % Use cluster-based correction
cfg.clusteralpha = 0.05; % Cluster significance threshold
cfg.numrandomization = 1000; % Number of permutations

% Cluster definition settings
cfg.clusterstatistic = 'maxsum';
% cfg.clustercritval   = 0.05;
% cfg.minnbchan        = 2; % Minimum neighboring sources to form a cluster
cfg.tail             = 0; % Two-tailed test
cfg.alpha            = 0.025; % Significance level (two-tailed)
% cfg.correcttail      = 'alpha';

% Define design matrix
cfg.design(1,:) = [1:Nsubj 1:Nsubj];
cfg.design(2,:) = [ones(1,Nsubj)*1 ones(1,Nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

% Define neighborhood structure
% cfg.connectivity = 'nearest'; % Use a nearest-neighbor connectivity criterion
% cfg.connectivity = 1;

% allsource.inside = true(size(allsource.inside ))

allsourcezeros = allsource;
allsourcezeros.pow = zeros(size(allsource.pow));
% Run cluster permutation test
stat = ft_sourcestatistics(cfg, allsource, allsourcezeros);

% stat.mask = double(stat.prob < 0.05);
stat.mask = double(stat.posclusterslabelmat > 21);
cfg = [];
cfg.funparameter = 'stat';
cfg.maskparameter = 'mask';
ft_sourceplot(cfg, stat) % mri


%% predict saliency of the presented image from trialwise SDbold maps
% does SDbold represent image saliency?
% need for training format: Ntrials x Nmaps, ca: 78*150 = 11700 trials (exemplars), 259200 voxels (features) 

param=[];
param.svm_type= 3;
param.kernel= 'linear'; % linear rbf
param.degree= 3;
param.gamma= []; % 
% param.gamma= 1/size(Xdat,1); % 1/nfeatures, only for rbf kernel
param.coef0= 0;
param.cost= 0.1; % 1
param.nu= 0.5000;
param.epsilon= 1;  % 0.1
param.cachesize= 4000; % 100
param.eps= 1e-03; % 1e-03
param.shrinking= 0;
param.probability_estimates= 0;
param.weight= 1;
param.cv= [];
param.quiet= 0;
param.kernel_type= 0;

X = permute(allsource.pow(:,allsource.inside,:), [1 3 2]);
ypred_subj = []; corr_predvsactual = [];
nsub = size(X,1);
for isub = 1:3 %nsub
  subjbool = true([nsub 1]);
  subjbool(isub) = false; % remove isub from training set
  Xdat = reshape(X(subjbool,:,:), [], size(X,3)); % samples X features
  labeldat =  labels(:,subjbool);
  
  tic
  cf = train_libsvm(param, Xdat, labeldat); % call train_libsvm and test_libsvm directly
  toc
  % for testing: 150 exemplars * N features
  testdat = double(squeeze(X(isub,:,:)));
  [ypred, dval] = test_libsvm(cf, testdat);

  testdatlabels = labels(:, isub);
  corr1 = corr(ypred, testdatlabels, 'Type','Spearman');
  figure; scatter(ypred, testdatlabels); lsline; title(corr1)
  figure; histogram(ypred); hold on;  histogram(testdatlabels)

  ypred_subj(isub,:) = ypred;
  yactual_subj(isub,:) = testdatlabels;
  corr_predvsactual(isub,:) = corr1;
end
figure; histogram(corr_predvsactual,10);
mean(corr_predvsactual)
[h,p] = ttest(corr_predvsactual)



%% OLD
% nvoxels = 3000;
% Xdat = Xdat(:,1:nvoxels);
% testdat = testdat(:,1:nvoxels);
% cfg = [];
% % cfg.classifier = 'multiclass_lda';
% cfg.classifier = 'libsvm';
% cfg.hyperparameter.svm_type = 3;
% cfg.hyperparameter.kernel = 'linear';
% % cfg.metric = 'accuracy'; %auc accuracy
% cfg.metric = 'none'; %auc accuracy
% % cfg.cv = 'leaveout'; % leaveout kfold
% cfg.sample_dimension = 1;
% cfg.feature_dimension = 2;
% % cfg.k = 10;
% % cfg.repeat = 1;
% tic
% accuracy = mv_classify(cfg, Xdat, labeldat, testdat, testdatlabels);
% toc
% figure; scatter(accuracy{1}, testdatlabels); lsline; title(corr(accuracy{1}, testdatlabels, 'Type','Spearman'))
% figure; histogram(accuracy{1})
% acc = mean(accuracy{1} == testdatlabels)*100
% chance = 1/150*100

%% OLDER
% %% decode image complexity from fMRI entropy
% % To give a concrete code example, consider the ?faces vs. houses? experiment.
% % For each trial, the BOLD response has been recorded for all voxels. This yields a
% % [samples x voxels] data matrix for one subject, where the samples correspond to
% % trials and the voxels serve as features. The matrix is denoted as X. Each trial
% % corresponds to either a ?face? or a ?house? stimulus. This is encoded in a vector
% % of class labels, denoted as clabel, that contains 1?s and 2?s (?face? = 1, ?house? = 2).
% % Then the following piece of code performs 10-fold cross-validation with 2 repetitions.
% % LDA is used as classifier and area under the ROC curve (AUC) is calculated as a classification metric.
% 
% % get data into trials X voxels
% source = src.source(2);
% accuracy = [];
% for itim = 1:length(source.time)
%   tmp = source;
%   % tmp.pow = squeeze(source.pow(:,source.inside,3));
%   % tmp.pow = squeeze(mean(source.pow(:,source.inside,3:end),3));
%   tmp.pow = source.pow(:,source.inside,itim);
%   tmp.powdimord = 'rpt_pos';
%   X = tmp.pow;
% 
%   % get HMAX bins and run decoding
%   disp 'make bins of trials based on hmax'
%   [sortHMAX, sortinds] = sort(source.trialinfo(:,10));  %hmax in 10, ascending, trial inds
%   nbins = 2;
%   [~,~,clabel] = histcounts(sortinds, nbins);
%   cfg = [];
%   cfg.model = 'lda';
%   cfg.metric = 'accuracy'; %auc accuracy
%   cfg.cv = 'kfold'; % leaveout kfold
%   cfg.k = 10;
%   cfg.repeat = 20;
%   accuracy(1,itim) = mv_classify(cfg, X, clabel)
% end
% disp(accuracy)
% %%
% disp(outfile)
% save(outfile, 'accuracy')