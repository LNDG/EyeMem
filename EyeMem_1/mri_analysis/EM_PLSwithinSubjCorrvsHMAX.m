function accuracy = EM_PLSwithinSubjCorrvsHMAX(cfg)
% Corr SDBOLD and HMAX across trials for each subject, then
% 1) Adaptation to feature richness: task PLS OA vs YA: positive link YA, weaker OA?
% 2) FR adaptation vs memory: behav PLS: positive link YA, weaker OA? 

mri = ft_read_mri('/Users/kloosterman/projectdata/eyemem/Standards/MNI152_T1_3mm_brain.nii.gz');

mri_mask = ft_read_mri('/Users/kloosterman/projectdata/eyemem/Standards/MNI152_T1_3mm_brain_mask.nii.gz');
% mri_mask.anatomy = round(mri_mask.anatomy); 
mri_mask.anatomy(mri_mask.anatomy < 1) = 0; % binarize
cfg=[];
ft_sourceplot(cfg, mri_mask); %, mri_mask
%% format PLS toolbox multiple conds:
% Condition	Subject	Feature 1	Feature 2	...	Feature N
% C1	S1	X_11	X_12	...	X_1N
% C1	S2	X_21	X_22	...	X_2N
% C2	S3	X_31	X_32	...	X_3N
% C2	S4	X_41	X_42	...	X_4N
% C3	S5	X_51	X_52	...	X_5N
% C3	S6	X_61	X_62	...	X_6N

%% load all subjects, compute SD over time, correlate to HMAX, append
datapath = '/Users/kloosterman/projectdata/eyemem/variability2/5TRspertrial/ftsource';

sources = cell(2,1);
cd(datapath)
subjects = dir(fullfile(datapath, 'source_sub*'));
labels = [];
load participantinfo.mat % TODO make this reliable

stat_mse = {}; corrlist = [];
for isub = 1:length(subjects) % TODO leave out subject under investigation
  subj = subjects(isub).name(8:end-4);
  subjinfo = Participants(Participants.participant_id == subj, :);     % give different outfolder for OA and YA
  % if subjinfo.group == 'old';    continue;   end
  if subjinfo.group == 'young';
    iage = 1;
  else
    iage =2;
  end % TODO get both

  disp(isub)
  % get eye data for labels
  eyedat = load(['eye' subjects(isub).name(7:end)]);
  % look_region_dat = eyedat.trialinfo.HMAX_fix_lookregion_mean; % _lookregion_std
  look_region_dat = eyedat.trialinfo.HMAX_fix; % _lookregion_std

  source = load(subjects(isub).name);
  source.inside = logical(mri_mask.anatomy(:)); 
  source.freq = 1;

  SDbold = squeeze(std(source.pow, 0, 3));% compute SD over time
  corrdat = NaN(size(source.inside));
  corrdat(source.inside,1) = corr(SDbold(source.inside,:)', look_region_dat, 'Type', 'Spearman');
  source.pow = corrdat;
  source.powdimord = 'pos';

  sources{iage}{end+1} = source;

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

cfg=[];
cfg.keepindividual = 'yes';
allsource_YA = ft_sourcegrandaverage(cfg, sources{1}{:});
allsource_OA = ft_sourcegrandaverage(cfg, sources{2}{:});
allsource_all = [sources{:}];
allsource_all = ft_sourcegrandaverage(cfg, allsource_all{:});
% remove Nan voxels where correlation went wrong
allsource_all.inside = all(not(isnan(allsource_all.pow)))';
allsource_YA.inside = all(not(isnan(allsource_all.pow)))';
allsource_OA.inside = all(not(isnan(allsource_all.pow)))';

plotit=0;
if plotit
  allsourceplot = allsource;
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


%% Task PLS 2 group, 
cfg = [];
cfg.parameter = 'pow';
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 100;                         % Number of permutation
cfg.num_boot = 100;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 1;                         % 1 is taskPLS; 3 is behavPLS
cfg.num_cond = 1;                           % Number of conditions
% cfg.num_subj_lst = size(allsource.pow,1);   % Number of subjects per condition
cfg.num_subj_lst = [39 39];   % Number of subjects per condition
cfg.design = zeros(size(allsource_all.pow,1),1); % Placeholder TODO omit design for Task PLS
stat_taskPLS = ft_sourcestatistics(cfg, allsource_OA, allsource_YA);

stat_taskPLS.results.perm_result.sprob
bardat = [stat_taskPLS.brainscores(1:39,1) stat_taskPLS.brainscores(40:end,1)];
figure; bar(mean(bardat));
figure; boxplot(bardat);

%%
stat_taskPLS.bootstrapratios = zeros(size(stat_taskPLS.inside));
stat_taskPLS.bootstrapratios(stat_taskPLS.inside,1) = stat_taskPLS.results.boot_result.compare_u(:,1);

cfg=[];
cfg.funparameter = 'bootstrapratios';
cfg.method = 'ortho'; % slice ortho glassbrain vertex
load colormap_jetlightgray.mat
cfg.funcolormap = cmap;
cfg.funcolorlim = 'maxabs';
% cfg.funcolorlim = 'zeromax';
cfg.location = [25 11 30];
cfg.locationcoordinates = 'voxel';
ft_sourceplot(cfg, stat_taskPLS); %, mri_mask

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

%
% stat.mask = double(stat.prob < 0.05);
stat.mask = double(stat.posclusterslabelmat > 21);
cfg = [];
cfg.funparameter = 'stat';
cfg.maskparameter = 'mask';
ft_sourceplot(cfg, stat) % mri




%% Task PLS 1 cond, does not work?
% cfg = [];
% cfg.parameter = 'pow';
% cfg.statistic = 'ft_statfun_pls';           % PLS statistics
% cfg.num_perm = 100;                         % Number of permutation
% cfg.num_boot = 10;
% cfg.method = 'analytic';                    % analytic method for statistics
% cfg.pls_method = 1;                         % 1 is taskPLS; 3 is behavPLS
% cfg.num_cond = 1;                           % Number of conditions
% cfg.num_subj_lst = size(allsource.pow,1);   % Number of subjects per condition
% cfg.design = zeros(size(allsource.pow,1),1); % Placeholder TODO omit design for Task PLS
% stat_mse{end+1} = ft_sourcestatistics(cfg, allsource);



if ismac
  cfg=[];
  cfg.funparameter = 'pow';
  cfg.method = 'ortho'; % slice ortho glassbrain vertex
  load colormap_jetlightgray.mat
  %     cfg.funcolormap = cmap(129:end,:);
  cfg.funcolormap = cmap;
  %         cfg.funcolorlim = 'zeromax';
  cfg.funcolorlim = 'maxabs';
  cfg.location = [25 11 30];
  cfg.locationcoordinates = 'voxel';
  %     cfg.colorbar = 'yes'
  %     ft_sourceplot(cfg, tmp, anat)
  allsourceplot = allsource;

  allsourceplot.pow = mean(allsourceplot.pow,3);
  allsourceplot.pow = mean(allsourceplot.pow)';
  allsourceplot.dimord = 'pos';
  ft_sourceplot(cfg, allsourceplot);
  figure; histogram(allsourceplot.pow); max(allsourceplot.pow)
  figure; histogram(allsourceplot.pow(allsourceplot.inside)); numel(allsourceplot.pow(allsourceplot.inside))% ca 55000 voxels, 3500 around zero
end

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