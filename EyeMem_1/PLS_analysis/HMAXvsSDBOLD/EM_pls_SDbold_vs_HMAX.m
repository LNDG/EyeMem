function EM_pls_SDbold_vs_HMAX(cfg)
% run task PLS on 5-trial binned SDbold vs Hmax. 30 conditions
% TODO turn into data format that PLS toolbox likes:
% examples:
% sesdat = load('/Users/kloosterman/gridmaster2012/LNDG/EyeMem_old/PLS/SDdatamats/EYEMEM009_BfMRIsessiondata.mat') % behav PLS?
% sesdat2 = load('/Volumes/FB-LIP/Projects/COBRA/data/mri+PETSharp/B_analyses/Paper1_2018/PLS_nback/SD_2mm/SD_C219_BfMRIsessiondata.mat') % task PLS?
% close all

sourcefile = cfg.sourcefile;
outfile_source = cfg.outfile_source;
outfile_sesdat = cfg.outfile_sesdat;
nbins = cfg.nbins;
behavfile = cfg.behavfile ;
PLStype = cfg.PLStype;
PLSbehav = cfg.PLSbehav;
subj = cfg.subj;
BOLDvar_measure = cfg.BOLDvar_measure;
binsubtract = cfg.binsubtract;
HMAXfolder = cfg.HMAXfolder;
eyefile = cfg.eyefile;
PREIN = cfg.PREIN;
PREOUT = cfg.PREOUT;
gazespecificHMAX= cfg.gazespecificHMAX;
fitcoeff = cfg.fitcoeff;
bintype = cfg.bintype;
nTRpertrial = cfg.nTRpertrial;
inducedortotalSD = cfg.inducedortotalSD;

disp(sourcefile)
source = load(sourcefile); % source comes out

disp 'only take GM voxels common in all subjects'
load(fullfile(PREIN, 'common_coords.mat'), 'common_coords');
source.inside = common_coords;
source.pow(~source.inside,:,:) = 0;
% source.pow = source.pow(source.inside,:,:);
source.time = [0 1 2 3 4]; % Breaks for classic LSS 1 beta per trial
source.cfg = [];
% source.trialinfo(:,end+1) = 1:150; % number trials to keep track

load(eyefile, 'trialinfo')
if size(trialinfo,1) ~= size(source.trialinfo,1)
  error('brain and behavior trials do not match')
end
source.trialinfo = trialinfo; % TODO check if ntrials match

disp 'Remove last trial since it is only zeros FIXME for non-gazespecific'
cfg = [];
cfg.trials = 1:149;
source = ft_selectdata(cfg, source); % TODO are trial eye and source still matched??

% Between or within-trial variability: subtract within trial mean beta weight per trial
switch inducedortotalSD
  case 'induced' % 
    source.pow = source.pow - mean(source.pow,3);    
  case 'within_trial' % within trial variability   
    source.pow = source.pow - source.pow(:,:,1);    
  case 'between_trial' % between trial variability     
    source.pow = repmat(source.pow(:,:,1), [1 1 5]);    
  case 'evoked'   
    source.pow = mean(source.pow,3); % average within each trial to isolate across trl var
    source.pow = repmat(source.pow, [1 1 5]); % tile it, strictly not necessary
  case 'total_pow' 
    disp 'All power considered'
end

% validtrials = find(source.trialinfo(:,end));
% validtrials = 1:ntrials;

%%
switch gazespecificHMAX
  case 'non-gazespecific' % bin based on overall HMAX, take SD over 5 trials    
    % sort onsets based on hmax TODO run for HMAX C2    
    [sortHMAX, sortinds] = sort(source.trialinfo.HMAX_fix);  %hmax in 10, ascending, trial inds - 10 is c1median
  case 'gaze-specific'
    %     [sortHMAX, sortinds] = sort(data.trialinfo(validtrials,17));  % gazelocked hmax in 17
    %     [sortHMAX, sortinds] = sort(source.trialinfo.HMAX_fix);  % gazelocked hmax in 17
    bin_variable = source.trialinfo.HMAX_fix_lookregion_mean;
end
cfg=[];
cfg.trials = bin_variable > 0;     % only take trials with valid HMAX_fix_lookregion
source = ft_selectdata(cfg, source);

[sortHMAX, sortinds] = sort(bin_variable(bin_variable > 0));

ntrials = size(source.trialinfo,1);
ntrlperbin = ntrials / nbins; % each subject has max 149 trials
[bininds, binedges] = discretize(1:ntrials, nbins);
bininds = sortrows([sortinds, bininds']);
bininds = bininds(:,2);
bininds(isnan(sortHMAX)) = NaN; % set outliers to nan
binedges = [sortHMAX(1); sortHMAX(binedges(2:end-1)); sortHMAX(end)];

disp 'make bins of trials based on hmax'
source_bin = source; 
source_bin.pow = nan(size(source_bin.pow,1), nbins);
source_bin.powdimord = 'pos_freq'; % freq is hmax condition
for ibin = 1:nbins
  seldat = source.pow(source.inside, bininds==ibin, :); %seldat = source.pow(:,cond_bins(:,ibin),:);   
  seldat =  seldat(:,:);
%   hmaxperbin(ibin,1) = mean(hmax_at_fix_trl(bininds==ibin));
  
  switch BOLDvar_measure
    case 'mean'
      source_bin.pow(source_bin.inside,ibin) = mean(seldat,2);
    case 'std'
      source_bin.pow(source_bin.inside,ibin) = std(seldat,0,2);
    case 'iqr'
      source_bin.pow(source_bin.inside,ibin) = iqr(seldat,2); % take IQR across 5 trials, 5 TR's each
    case 'mse'
      data = [];   % make data struct  % Required fields:  %   - time, trial, label
      for itrial = 1:size(source.pow, 2)
        data.time{itrial} = source.time;
        data.trial{itrial} = squeeze(source.pow(source.inside,itrial,:));
      end
      data.label = [];
      for ichan = 1:length(data.trial{1})
        data.label{ichan,1} = sprintf('%d', ichan);
      end
      
      cfg=[];
      cfg.trials = bininds==ibin;
      cfg.toi = 2;
      cfg.timwin = 5;
      cfg.timescales = 1;
      cfg.filtmethod = 'no';
      cfg.recompute_r = 'perscale_toi_sp';
      cfg.coarsegrainmethod = 'pointavg';
      mse = ft_entropyanalysis(cfg, data);      
      source_bin.pow(source_bin.inside,ibin) = mse.sampen;
  end
  source_bin.freq(ibin) = nanmean(binedges(ibin:ibin+1)); % use freq field for HMAX bin_No
end
source = source_bin;


%% make PLS sesssiondata structure, prepare important fields
% 1. standard stuff: TODO turn into function?
tmp=[];
% load(fullfile(PREIN, 'common_coords.mat'));
tmp.st_coords = find(common_coords); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
tmp.session_info.datamat_prefix = subj; %[subj '_' pattern];%stores common      %datamat prefix
switch PLStype
  case 'taskPLS'
    tmp.st_datamat = transpose(source.pow(tmp.st_coords,:)); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
    tmp.behavdata = [];
    tmp.behavname = {};%no behav names or data just yet, but set up the %field anyway...
  case 'behavPLSvsSDT'
    if isnumeric(binsubtract)
      tmp.st_datamat = transpose(source.pow(tmp.st_coords, binsubtract(1))) - transpose(source.pow(tmp.st_coords, binsubtract(2))); % highest - lowest BOLD variability
    else
      dat = transpose(source.pow(tmp.st_coords, :));
      for i = 1:length(dat)
        fit = polyfit( transpose(1:nbins), dat(:,i), 1 );
        tmp.st_datamat(1,i) = fit(1);
      end
    end
    tmp.behavname = {PLSbehav};%no behav names or data just yet, but set up the %field anyway...
    load(behavfile); % behavior comes out
    subjind = behavior.participants.participant_id == subj;
    tmp.behavdata = mean(behavior.(PLSbehav)(subjind,2,6)); % 2 is test, 6 is cond average
  case 'behavPLSvsDDM'
    if isnumeric(binsubtract)
      if numel(binsubtract) == 2
%       tmp.st_datamat = (transpose(source.pow(tmp.st_coords, binsubtract(1))) - transpose(source.pow(tmp.st_coords, binsubtract(2)))) ./ ...
%         transpose(source.pow(tmp.st_coords, binsubtract(2))); % highest - lowest BOLD variability psc
      tmp.st_datamat = (transpose(source.pow(tmp.st_coords, binsubtract(1))) - transpose(source.pow(tmp.st_coords, binsubtract(2)))); % highest - lowest BOLD variability
      elseif numel(binsubtract) == 1
        tmp.st_datamat = transpose(source.pow(tmp.st_coords, binsubtract(1)));
      end
    elseif strcmp(binsubtract, 'corrHmaxoverbins')
      disp 'get hmax_at_fix_trl per bin and correlate'
      dat = transpose(source.pow(tmp.st_coords, :));
      for i = 1:length(dat)
        tmp.st_datamat(1,i) = corr(hmaxperbin, dat(:,i), 'type', 'Pearson'); %figure; scatter(hmaxperbin, dat(:,i))
      end
    else
      dat = transpose(source.pow(tmp.st_coords, :));
      for i = 1:length(dat)
        fit = polyfit( transpose(1:nbins), dat(:,i), 1 );
        tmp.st_datamat(1,i) = fit(fitcoeff); % fit is in descending powers
      end
    end
    tmp.behavname = {PLSbehav};%no behav names or data just yet, but set up the %field anyway...
    load(behavfile); % behav comes out
    subjind = behavior.participants.participant_id == subj;
    tmp.behavdata = behavior.(PLSbehav)(subjind,2); % 2 is test phase
  case 'behavPLS_sdboldvsHmaxbins'
    dat = transpose(source.pow(tmp.st_coords, :));
    binvals = [binedges(1:end-1) binedges(2:end)];
    binvals = mean(binvals,2); % the mean val in each bin. TODO take median of trials in each bin?
    tmp.st_datamat(1,:) = corr(dat, binvals, 'type', 'Spearman'); % spearman accounts for nonlinearities in hmax vals
    tmp.behavname = {PLSbehav};%no behav names or data just yet, but set up the %field anyway...
    load(behavfile); % behav comes out
    subjind = behavior.participants.participant_id == subj;
    tmp.behavdata = behavior.(PLSbehav)(subjind,2); % 2 is test phase, Add memory performance
end
tmp.st_evt_list = 1:size(tmp.st_datamat,1);%how many conditions?
% tmp.st_sessionFile = [pls_dir subj '_' pattern '_BfMRIsession.mat'];
tmp.unequal_subj = 0;

% session_info struct, for bookkeeping I believe:
%       a = load([workdir '/charisaPLS/' subj '_BfMRIsession.mat']);%this      %loads a particular session file.
tmp.session_info.pls_data_path = fileparts(outfile_source); %this will add info to 'a' for   tracking purposes, and for commands below.
tmp.session_info.file_pattern = '*nii'; % NEEDED?
[pathstr, name]= fileparts(sourcefile);
tmp.session_info.run(1).data_path = pathstr; %[workdir subj];
tmp.session_info.run(1).data_files = {[name '.nii']}; % CHECK should be nifti? e.g. /Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/5TRspertrial/sub-77/beta_series
tmp.session_info.run(1).file_pattern = tmp.session_info.run(1).data_files; %['smooth_CommonAtlas_mc_tshift_cond' num2str(run) '.nii*'];

tmp.session_info.num_conditions = length(tmp.st_evt_list);
tmp.session_info.condition = {};
for i = 1:tmp.session_info.num_conditions
  tmp.session_info.condition{end+1} = sprintf('HMAX%d', i);
end
% tmp.session_info.condition_baseline = {[0 1]  [0 1]  [0 1]  [0 1]  [0 1]  [0 1]};
% % needed? silly baseline
% tmp.session_info.num_conditions0 =
% tmp.session_info.condition0 =
% tmp.session_info.condition_baseline0 =

tmp.session_info.num_runs = 1;
% bookkeeping of events:
[pathstr, name]= fileparts(source.hdr.fspec);
tmp.session_info.run.num_scans = source.hdr.nframes;
tmp.session_info.run.data_path = pathstr;
tmp.session_info.run.data_files = { [name '.nii'] };
tmp.session_info.run.file_pattern = [name '.nii'];
% tmp.session_info.run.blk_onsets = num2cell(cond_bins, 2)'; % =  {[9×1 double]  [9×1 double]  [9×1 double]  [2×1 double]} vectors with TR onsets
% tmp.session_info.run.blk_length = num2cell(ones(size(cond_bins)), 2)'; % =  {[9×1 double]  [9×1 double]  [9×1 double]  [2×1 double]} vectors with TR onsets% =  {[9×1 double]  [9×1 double]  [9×1 double]  [2×1 double]} vectors with lengths pe onset
for ibin = 1:nbins
  tmp.session_info.run.blk_onsets{ibin} = find(bininds==ibin);
  tmp.session_info.run.blk_length{ibin} = ones(size(tmp.session_info.run.blk_onsets{ibin}));
end

% all hardcoded stuff:
tmp.create_ver = '5.0807231';%not sure where value comes from; check.
tmp.st_win_size = 1;%PLS command...must check.

[i,j,k] = ind2sub(source.dim, find(source.pos(:,1)==0 & source.pos(:,2)==0 & source.pos(:,3)==0));
tmp.st_origin = [i j k]; % sets the origin...PLS default?

tmp.st_voxel_size = [source.hdr.xsize source.hdr.ysize source.hdr.zsize]; %vox size defined
tmp.st_dims = [ source.dim(1:2) 1 source.dim(3)]; % [40 48 1 34];%dimensions of matrix defined for PLS WHY 1 for dim3??
tmp.SingleSubject = 0;%not a single subj design
tmp.normalize_volume_mean = 0;
tmp.create_datamat_info.brain_mask_file = ''; % NEEDED??
% tmp.create_datamat_info.run_idx = run_idx; % NOT NEEDED it seems, bijv [1:7];%seven runs here
tmp.create_datamat_info.brain_coord_thresh = 0;%no ROIs here...check
tmp.create_datamat_info.consider_all_voxels_as_brain = 0;%we don't...
tmp.create_datamat_info.num_skipped_scans = 0;%no skipped scans
tmp.create_datamat_info.ignore_slices = [];%no slices to ignore, so %keep all
tmp.create_datamat_info.normalize_volume_mean = 0;%no normalization %here...
% tmp.create_datamat_info.normalize_with_baseline = 1;%yes, normalize %to baseline (but changes below, so must override?)
tmp.create_datamat_info.normalize_with_baseline = 0;%yes, normalize %to baseline (but changes below, so must override?)
% tmp.create_datamat_info.merge_across_runs = 1;%says to merge across %runs, but we don't do this below anyway...is overridden...??
tmp.create_datamat_info.merge_across_runs = 0;%says to merge across %runs, but we don't do this below anyway...is overridden...??
tmp.create_datamat_info.single_subject_analysis = 0;%not single subj

if strcmp(gazespecificHMAX, 'gaze-specific')
%   tmp.fixdur_keep = fixdur_keep;
%   tmp.hmax_at_fix_keep = hmax_at_fix_keep;
%   tmp.hmax_at_fix_trl = hmax_at_fix_trl; % to plot extracted saliency per trial
%   tmp.hmax_meanperbin = source_bin.freq;
  tmp.fixdur_keep = source.trialinfo.fixdur_mean;
  tmp.hmax_at_fix_keep = source.trialinfo.HMAX_fix;
%   tmp.hmax_at_fix_trl = hmax_at_fix_trl; % to plot extracted saliency per trial
  tmp.hmax_meanperbin = source_bin.freq;
end

disp(outfile_source)
save(outfile_source, 'source')

% save PLS sesdat
disp(outfile_sesdat)
save(outfile_sesdat, '-struct','tmp')

out = fullfile(fileparts(fileparts(outfile_sesdat)), [subj '_BfMRIsessiondata.mat']);
save(out , '-struct','tmp')
% % save([subj '_' pattern '_BfMRIdatamat.mat'], '-struct','tmp'); %,'-mat'

if ismac
  tmp = source_bin;
  tmp.powdimord = 'pos';
  %       vol=300;
  %   tmp.anatomy = tmp.anatomy(:,:,:,vol);
  %     tmp.pow = mean(tmp.pow(:,:),2);
  tmp.pow = tmp.pow(:,3);
  cfg=[];
  cfg.method = 'ortho'; % slice ortho glassbrain vertex
  cfg.funparameter = 'pow';
  cfg.funcolorlim = 'zeromax';% [-300 300];
  ft_sourceplot(cfg, tmp)
end

