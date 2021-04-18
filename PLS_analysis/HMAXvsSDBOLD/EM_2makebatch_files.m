function EM_2makebatch_files()
% Create matrix template for each individual containing condition info
% i.e. block onsets + lengths, we choose a fixed block length here

% Hack: To get a balanced matrix of condition onsets etc. -1 is entered for
% conditions that were not present for the run. Note that PLS give a
% warning for those, but will automatically drop them from processing.

% Note: For '2131' and '2237' only the first two runs are available.

%% paths
pn.root     = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem';
pn.trialinfo   = fullfile(pn.root, 'preproc/eye'); % eyedata.mat: location of info regarding run trialinfo matrices
pn.BOLDpath = fullfile(pn.root, 'variability/5TRspertrial');
pn.PLSfiles = fullfile(pn.root, 'variability/PLS_HMAXvsSDBOLD'); mkdir(pn.PLSfiles);

%% fill with necessary information
% N = 44 YA + 53 OA;
IDs={'sub-09', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-61', 'sub-62', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-71', 'sub-72', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-100', 'sub-101'}; % 'sub-70', 'sub-73', 

% runs = {'Run1'; 'Run2'; 'Run3'; 'Run4'};
% runs_boldname = {'run-1'; 'run-2'; 'run-3'; 'run-4'};
conditions = {};
for icond = 1:30
  conditions{end+1} = sprintf('hmax%d', icond);
end

for indID = 1:numel(IDs)
  % generate batch_file
  batch_file = [];
  batch_file.across_run = 1;
  batch_file.single_subj = 0;
  batch_file.is_analysis = 0;
  batch_file.prefix = ['task_', IDs{indID}];
  batch_file.cond_name = conditions;
  batch_file.brain_region = -10000; %0
  batch_file.ref_scan_onset = zeros(1,numel(conditions));
  batch_file.num_ref_scan = ones(1,numel(conditions));

  numberOfRuns = 1;
  for indRun = 1:numberOfRuns
    batch_file.data_files{1,indRun} = fullfile(pn.BOLDpath, IDs{indID}, ...
      'beta_series', sprintf('%s_sd_alltrials_beta_series.nii',  IDs{indID})); % ~750 volumes
    if exist(batch_file.data_files{1,indRun})
%       fprintf('found nii file %s\n', batch_file.data_files{1,indRun})
    else
      fprintf('nii file %s not found\n', batch_file.data_files{1,indRun})
      continue
    end
%     nvols = 750; % in principle 150 trials 5 betas per trial, but could be less if not enough time after 
    nvols = get_nii_frame(batch_file.data_files{1,indRun});    
    fprintf('%d vols\n', nvols)
%     onsets = transpose(1:5:nvols-5); % exclude last trial: not modelled well
    onsets = transpose(1:5:750); 
    onsets(end) = NaN; % exclude last trial: may be not modelled well in GLM

    % load trialinfo
    load(fullfile(pn.trialinfo, sprintf('eye_%s.mat', IDs{indID})))
    hmax = data.trialinfo(:,10);     % hmax per trial are in col 10
    hmax = hmax(1:length(onsets)); % only take the ones that have data
    % sort onsets by hmax low to high
    [hmaxsorted, hmaxsorted_inds] = sort(hmax);
    onsets_sorted = onsets(hmaxsorted_inds);
    
    % dimord cond_onsets, so 30 by 5
    onsets_sorted_rs = reshape(onsets_sorted, [5 30])';
    batch_file.block_onsets{1,indRun} = onsets_sorted_rs-1; % TODO check! keep in mind to indicate x-1 for PLS (unix style)
    batch_file.block_length{1,indRun} = zeros(size(onsets_sorted_rs)) + 5;    

    % Julian code
%     for icond = 1:numel(conditions)
%       % keep in mind to indicate x-1 for PLS (unix style)
%       numOnsets = numel(find(Regressors(:,4)==icond));
%       batch_file.block_onsets{1,indRun}(icond,1:numOnsets) = find(Regressors(:,4)==icond)-1;
%       % block length should be constant at approx. 125 volumes excl. feedback (find(Regressors(:,4)))
%       batch_file.block_length{1,indRun}(icond,1:numOnsets) = repmat(125,1,numOnsets);
%     end
  end
  % save individual configuration
  save( fullfile(pn.PLSfiles, sprintf('%s_task_PLS_info.mat', IDs{indID} )), 'batch_file')
end
