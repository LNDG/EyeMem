function accuracy = EM_decodeSaliency(cfg)

%  TODO only nans: load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/LOTO/ERPremoved/YA/source_sub-50.mat
% if nargin == 0
%   load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/LOTO/ERPremoved/YA/source_sub-65.mat')
% end

infile = cfg.infile;
outfile = cfg.outfile;

disp(infile)
source = load(infile);

% make inside field uniform across subjects
mri = ft_read_mri('/Users/kloosterman/projectdata/eyemem/MNI152_T1_3mm_brain.nii.gz');
% mri_mask = ft_read_mri('/Users/kloosterman/projectdata/eyemem/Standards/MNI152_T1_3mm_brain_mask.nii.gz')
mri_mask = ft_read_mri('/Users/kloosterman/projectdata/eyemem/Standards/old_masks(deprecated)/avg152_T1_gray_mask_90_3mm.nii.gz');
mri_mask.anatomy = round(mri_mask.anatomy); 
% ft_sourceplot([], mri, mri_mask); ft_sourceplot([],  mri_mask, mri)

region_of_interest = 'Hippocampus' % 'Hippocampus' 'Visual cortex'
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
% cfg=[];
% cfg.funparameter = 'tissue';
% cfg.anaparameter = 'anatomy';
% ft_sourceplot(cfg, mri_mask, mri);
% resample mask to 3 mm
cfg=[];
cfg.parameter = 'tissue';
mri_mask3mm = ft_sourceinterpolate(cfg, mri_mask, mri);

cfg=[];
cfg.funparameter = 'tissue';
cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, mri_mask3mm, mri);


%%
if ismac
  
  tmp = source;
  close all
  
  for itrial=2:5:25
    %     tmp.pow = squeeze(source.pow(itrial,:,:));
    %     tmp.powdimord = 'pos_time';
    %     tmp.pow = squeeze(source.pow(itrial,:,3))';
    % tmp.pow = squeeze(mean(source.pow(:,itrial,:),3));
    tmp.pow = squeeze(std(source.pow(:,itrial,:),0,3));
    tmp.powdimord = 'pos';
    cfg=[];
    cfg.funparameter = 'pow';
    cfg.method = 'ortho'; % slice ortho glassbrain vertex
    load colormap_jetlightgray.mat
    %     cfg.funcolormap = cmap(129:end,:);
    cfg.funcolormap = cmap;
    %         cfg.funcolorlim = 'zeromax';
          cfg.funcolorlim = 'maxabs';
    % cfg.funcolorlim = [-0.015 0.015];% 'zeromax'
    cfg.location = [25 11 30];
    cfg.locationcoordinates = 'voxel';
    %     cfg.colorbar = 'yes'
    %     ft_sourceplot(cfg, tmp, anat)
    ft_sourceplot(cfg, tmp, mri_mask);
  end
end

%% load all subjects, compute SD over time, append

sources = {};
cd(fullfile(fileparts(infile)))
subjects = dir(fullfile(fileparts(infile), 'source_sub*'));
labels = [];
load participantinfo.mat % TODO make this reliable

for isub = 1:length(subjects) % TODO leave out subject under investigation
  subj = subjects(isub).name(8:end-4);
  subjinfo = Participants(Participants.participant_id == subj, :);     % give different outfolder for OA and YA
  if subjinfo.group == 'old';    continue;   end

  disp(isub)
  % get eye data for labels
  eyedat = load(['eye' subjects(isub).name(7:end)]);
  label_subj = zscore(eyedat.trialinfo.HMAX_fix_lookregion_mean);
  labels = [labels label_subj];

  tmp = load(subjects(isub).name);
  tmp.pow = squeeze(std(tmp.pow, 0, 3));% compute SD over time

  tmp.pow = zscore(tmp.pow,0,2); % normalize each feaure separately across subject's trials
  
  tmp.powdimord = 'pos_rpt';
  % tmp.inside = logical(mri_mask.anatomy(:)); % use visual cortex mask
  tmp.inside = logical(mri_mask3mm.tissue(:)); % use visual cortex mask
  
  tmp.freq = 1;
  plotit=0;
  if plotit
    for itrial=2:5:20
      tmp2 = tmp; 
      tmp2.pow = tmp2.pow(:,itrial);
      cfg=[];
      cfg.funparameter = 'pow';
      cfg.method = 'ortho'; % slice ortho glassbrain vertex
      load colormap_jetlightgray.mat
      cfg.funcolormap = cmap;
      cfg.funcolorlim = 'maxabs';
      cfg.location = [25 11 30];
      cfg.locationcoordinates = 'voxel';
      ft_sourceplot(cfg, tmp2, mri); %, mri_mask
    end
  end
  sources{end+1} = tmp;
end

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