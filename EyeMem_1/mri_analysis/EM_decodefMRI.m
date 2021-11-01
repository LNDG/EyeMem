function accuracy = EM_decodefMRI(cfg)

%  TODO only nans: load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/LOTO/ERPremoved/YA/source_sub-50.mat
% if nargin == 0
%   load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/LOTO/ERPremoved/YA/source_sub-65.mat')
% end

infile = cfg.infile;
outfile = cfg.outfile;

disp(infile)
load(infile)

%%
if ismac
  
  source = src.source(1);
  tmp = source;
  close all
  
  for itrial=2:5:25
    %     tmp.pow = squeeze(source.pow(itrial,:,:));
    %     tmp.powdimord = 'pos_time';
    %     tmp.pow = squeeze(source.pow(itrial,:,3))';
    tmp.pow = squeeze(mean(source.pow(itrial,:,3:end),3))';
    tmp.powdimord = 'pos';
    cfg=[];
    cfg.funparameter = 'pow';
    cfg.method = 'ortho'; % slice ortho glassbrain vertex
    load colormap_jetlightgray.mat
    %     cfg.funcolormap = cmap(129:end,:);
    cfg.funcolormap = cmap;
    %         cfg.funcolorlim = 'zeromax';
    %       cfg.funcolorlim = 'maxabs';
    cfg.funcolorlim = [-0.015 0.015];% 'zeromax'
    cfg.location = [25 11 30];
    cfg.locationcoordinates = 'voxel';
    %     cfg.colorbar = 'yes'
    %     ft_sourceplot(cfg, tmp, anat)
    ft_sourceplot(cfg, tmp);
  end
end

%% decode image complexity from fMRI entropy
% To give a concrete code example, consider the ?faces vs. houses? experiment.
% For each trial, the BOLD response has been recorded for all voxels. This yields a
% [samples x voxels] data matrix for one subject, where the samples correspond to
% trials and the voxels serve as features. The matrix is denoted as X. Each trial
% corresponds to either a ?face? or a ?house? stimulus. This is encoded in a vector
% of class labels, denoted as clabel, that contains 1?s and 2?s (?face? = 1, ?house? = 2).
% Then the following piece of code performs 10-fold cross-validation with 2 repetitions.
% LDA is used as classifier and area under the ROC curve (AUC) is calculated as a classification metric.

% get data into trials X voxels
source = src.source(2);
accuracy = [];
for itim = 1:length(source.time)
  tmp = source;
  % tmp.pow = squeeze(source.pow(:,source.inside,3));
  % tmp.pow = squeeze(mean(source.pow(:,source.inside,3:end),3));
  tmp.pow = source.pow(:,source.inside,itim);
  tmp.powdimord = 'rpt_pos';
  X = tmp.pow;
  
  % get HMAX bins and run decoding
  disp 'make bins of trials based on hmax'
  [sortHMAX, sortinds] = sort(source.trialinfo(:,10));  %hmax in 10, ascending, trial inds
  nbins = 2;
  [~,~,clabel] = histcounts(sortinds, nbins);
  cfg = [];
  cfg.model = 'lda';
  cfg.metric = 'accuracy'; %auc accuracy
  cfg.cv = 'kfold'; % leaveout kfold
  cfg.k = 10;
  cfg.repeat = 20;
  accuracy(1,itim) = mv_classify(cfg, X, clabel)
end
disp(accuracy)
%%
disp(outfile)
save(outfile, 'accuracy')