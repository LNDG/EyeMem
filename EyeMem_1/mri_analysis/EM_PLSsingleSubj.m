function accuracy = EM_PLSsingleSubj(cfg)

% infile = cfg.infile;
% outfile = cfg.outfile;
% 
% disp(infile)
% source = load(infile);

mri_mask = ft_read_mri('/Users/kloosterman/projectdata/eyemem/Standards/MNI152_T1_3mm_brain_mask.nii.gz');
mri_mask.anatomy = round(mri_mask.anatomy); 

%% load all subjects, compute SD over time, append
datapath = '/Users/kloosterman/projectdata/eyemem/variability2/5TRspertrial/ftsource';

sources = {};
cd(datapath)
subjects = dir(fullfile(datapath, 'source_sub*'));
labels = [];
load participantinfo.mat % TODO make this reliable

stat_mse = {}; corrlist = [];
for isub = 1:length(subjects) % TODO leave out subject under investigation
  subj = subjects(isub).name(8:end-4);
  subjinfo = Participants(Participants.participant_id == subj, :);     % give different outfolder for OA and YA
  % if subjinfo.group == 'old';    continue;   end
  if subjinfo.group == 'young';    continue;   end

  disp(isub)
  % get eye data for labels
  eyedat = load(['eye' subjects(isub).name(7:end)]);
  label_subj = eyedat.trialinfo.HMAX_fix_lookregion_mean;
  labels = [labels label_subj];

  source = load(subjects(isub).name);
  source.pow = squeeze(std(source.pow, 0, 3));% compute SD over time

  % source.pow = zscore(source.pow,0,2); % normalize each feaure separately across subject's trials
  
  source.powdimord = 'pos_rpt';
  source.inside = logical(mri_mask.anatomy(:)); 
  % source.inside = logical(mri_mask3mm.tissue(:)); % use visual cortex mask  
  source.freq = 1;
  % run single subject PLS

  cfg = [];
  cfg.parameter = 'pow';
  cfg.statistic = 'ft_statfun_pls';           % PLS statistics
  cfg.num_perm = 0;                         % Number of permutation
  cfg.num_boot = 0;
  cfg.method = 'analytic';                    % analytic method for statistics
  cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
  cfg.cormode = 8;                            % 0 is Pearson corr, 8 is Spearman
  cfg.num_cond = 1;                           % Number of conditions
  cfg.design = label_subj;
  cfg.num_subj_lst = size(source.pow,2); % Number of subjects per condition


  % Step 3: Compute statistics
  stat_mse{end+1} = ft_sourcestatistics(cfg, source);
  % corrlist = [corrlist stat_mse{end}.results.lvlvcorr];
  corrlist = [corrlist mean(stat_mse{end}.results.usc)];




  % plotit=0;
  % if plotit
  %   for itrial=2:5:20
  %     tmp2 = source; 
  %     tmp2.pow = tmp2.pow(:,itrial);
  %     cfg=[];
  %     cfg.funparameter = 'pow';
  %     cfg.method = 'ortho'; % slice ortho glassbrain vertex
  %     load colormap_jetlightgray.mat
  %     cfg.funcolormap = cmap;
  %     cfg.funcolorlim = 'maxabs';
  %     cfg.location = [25 11 30];
  %     cfg.locationcoordinates = 'voxel';
  %     ft_sourceplot(cfg, tmp2, mri); %, mri_mask
  %   end
  % end
  % sources{end+1} = tmp;
end
figure; histogram(corrlist,10)

cfg=[];
cfg.keepindividual = 'yes';
allsource = ft_sourcegrandaverage(cfg, sources{:});

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