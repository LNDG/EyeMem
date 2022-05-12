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
removeoutliers = cfg.removeoutliers;
Z_thresh = cfg.Z_thresh;
nbins = cfg.nbins;
do_kstest = cfg.do_kstest;
behavfile = cfg.behavfile ;
PLStype = cfg.PLStype;
subj = cfg.subj;
BOLDvar_measure = cfg.BOLDvar_measure;
binsubtract = cfg.binsubtract;
HMAXfolder = cfg.HMAXfolder;
eyefile = cfg.eyefile;
PREIN = cfg.PREIN;
gazespecificHMAX= cfg.gazespecificHMAX;
fitcoeff = cfg.fitcoeff;
bintype = cfg.bintype;

disp(sourcefile)
source = load(sourcefile); % source comes out

switch gazespecificHMAX
  case 'non-gazespecific' % bin based on overall HMAX, take SD over 5 trials    
    % sort onsets based on hmax TODO run for HMAX C2
    [sortHMAX, sortinds] = sort(source.trialinfo(:,10));  %hmax in 10, ascending, trial inds - 10 is c1median
  case 'gaze-specific'
    % load HMAX file
    hmaxlist=dir(fullfile(HMAXfolder, '*.mat' ));
    hmaxdat = {};
    for ih = 1:length(hmaxlist)
      load(fullfile(hmaxlist(ih).folder, hmaxlist(ih).name));
      hmaxdat{ih} = hmaxout;
    end
    
    load(eyefile)
    disp 'selecting only 0-5 s viewing time'
    cfg=[];
    cfg.latency = [0 5];
    data = ft_selectdata(cfg, data);
    ntrials = length(data.trial);
    hmax_at_fix_trl = nan(ntrials,1);
    
    for itrial = 1:ntrials
      % get HMAX data of pic shown
      catind = data.trialinfo(itrial, 2); %
      picno = data.trialinfo(itrial, 3); %
      picind = hmaxdat{catind}.picno == picno;
      curhmax = hmaxdat{catind}.c1(:,:,picind); % DONE add pic data to HMAX struct
      picdat = hmaxdat{catind}.picdat(:,:,picind);
      %     figure; imagesc(curhmax);
      
      % get fixation on and offsets: do based on non saccade episodes
      fixations = data.trial{itrial}(6,:)<0.1; % makes it 1 during fixation
      fixations(data.trial{itrial}(5,:) == 1) = 0; % blinks
      fixations(1) = 0; % so we get a fixation start at begin
      fixations(end) = 0; % so we get a fixation end at end
      fixtrig = [ find(diff(fixations) == 1)' find(diff(fixations) == -1)'];
      disp 'Drop first (trial starts with central fixation) and last fixation'
      fixtrig = fixtrig(2:end-1,:);
      
      % get XY coords of fixations: average XY within fixations
      nfix = size(fixtrig,1);
      fixloc = NaN(nfix,2);
      fixdur = NaN(nfix,1);
      cur_res = size(picdat); % resolution of the pics
      desiredres = size(curhmax);
      %     fixmap = zeros([size(curhmax) nfix]);
      fixloc_newres = NaN(nfix,2);
      disp 'resample gaze to curhmax resolution'
      xshift = (1024-640)/2;
      yshift = (768-480)/2;
      for ifix = 1:nfix
        fixinds = fixtrig(ifix,1):fixtrig(ifix,2);
        fixdur(ifix,1) = length(fixinds);
        fixloc(ifix,:) = round(mean(data.trial{itrial}(2:3,fixinds),2)); % Xgaze and Ygaze in chan2 and 3
        %       disp 'shift fixations, account for pic not fullscreen in scanner, but in the middle'
        fixloc(ifix,1) = fixloc(ifix,1)-xshift;
        fixloc(ifix,2) = fixloc(ifix,2)-yshift;
        % convert XY coords to resolution of HMAX
        fixloc_newres(ifix,:) = round(fixloc(ifix,:) ./ cur_res .* desiredres);
      end
      disp 'Drop fixations outside picture'
      %     validfix = NaN(size(fixloc,1),2);
      %     validfix(:,1) = fixloc(:,1) > 0 & fixloc(:,1) < 640;
      %     validfix(:,2) = fixloc(:,2) > 0 & fixloc(:,2) < 480;
      %     fixloc = fixloc(all(validfix,2),:);
      %     fixloc_newres = fixloc_newres(all(validfix,2),:); % also apply to resampled fix locations
      
      validfix = NaN(size(fixloc_newres,1),2);
      validfix(:,1) = fixloc_newres(:,1) > 0 & fixloc_newres(:,1) < desiredres(2);
      validfix(:,2) = fixloc_newres(:,2) > 0 & fixloc_newres(:,2) < desiredres(1);
      fixloc_newres = fixloc_newres(all(validfix,2),:); % also apply to resampled fix locations
      fixloc = fixloc(all(validfix,2),:);
      fixdur = fixdur(all(validfix,2));
      
      plotit=0;
      if ismac && plotit
        figure; hold on
        % The default EyeLink coordinates are those of a 1024 by 768 VGA display, with (0, 0) at the top left.
        imagesc(picdat) % for plotting transpose and flipud??
        ax=gca; ax.YDir = 'reverse';
        %       scatter(data.trial{itrial}(2,:)-xshift, data.trial{itrial}(3,:)-yshift, 'k'); hold on
        scatter(data.trial{itrial}(2,fixations)-xshift, data.trial{itrial}(3,fixations)-yshift, 'g'); hold on
        scatter(fixloc(:,1),fixloc(:,2), 'r', 'filled'); % already shifted
        %       xlim([0 1024]); ylim([0 768]); box on
        xlim([0 640]); ylim([0 480]); box on
        title(nfix)
        
        disp 'plot hmax and gaze'
        figure;
        imagesc(curhmax); hold on
        scatter(fixloc_newres(:,1),fixloc_newres(:,2), 'r', 'filled'); % already shifted
        ax=gca; ax.YDir = 'reverse';
        
        curhmax2 = curhmax;
        for ifix = 1:nfix
          curhmax2(fixloc_newres(ifix,2), fixloc_newres(ifix,1)) = 1;
        end
        figure; imagesc(curhmax2); hold on
        
        %       sc = scatter(fixind_newres(:,1), fixind_newres(:,2), 'filled');
        %       sc.SizeData = 50;
        %       sc.CData = [1 0 0];
      end
      
      disp 'get c1 HMAX vals at fixation locations'
      nfix = size(fixloc_newres,1);
      hmax_at_fix=NaN(nfix,1);
      for ifix = 1:nfix
        hmax_at_fix(ifix,1) = curhmax(fixloc_newres(ifix,2), fixloc_newres(ifix,1)); % Note the flip: Yaxis in dim1 (rows), Xaxis in dim2 (columns): scatter and plot need x,y, with indexing it's the other way around
      end
      if isempty(hmax_at_fix)
        continue
      end
      
      disp 'average over HMAX vals to get 1 val per trial'
      weightedmean = 0;
      if weightedmean == 1
        fixdur = fixdur / sum(fixdur);
        hmax_at_fix_trl(itrial,:) = sum((hmax_at_fix .* fixdur)) ;
      else
        hmax_at_fix_trl(itrial,:) = mean(hmax_at_fix);
      end
      
    end
    
    fprintf('%d trials without fixations found: ', sum(isnan(hmax_at_fix_trl)))
    if sum(isnan(hmax_at_fix_trl)) > 25
      warning('More than 50: skipping subject')
      return
    end

    switch bintype
      case 'fixednbins'
        % continue with sorting
        [sortHMAX, sortinds] = sort(hmax_at_fix_trl);  % gaze-specific HMAX values
        %         [sortHMAXold, sortindsold] = sort(source.trialinfo(:,10));
        ntrlperbin = ntrials / nbins; % each subject has 150 trials

        bininds = repmat(1:nbins, ntrlperbin, 1);
        bininds = bininds(sortinds(:));
        binedges = sortHMAX([1:ntrlperbin:ntrials ntrials])

      case 'uniformbinwidth'
        disp 'drop HMAX outlier trials'
%         figure; histogram(hmax_at_fix_trl, 100)
        [~,TF]=rmoutliers(hmax_at_fix_trl);
        fprintf('%d HMAX outliers found\n', sum(TF))
        hmax_at_fix_trl(TF) = NaN; % set outliers to nan
        [ntrlperbin,binedges,bininds] = histcounts(hmax_at_fix_trl,nbins);
        %         figure; histogram(hmax_at_fix_trl, 100)
    end
    
    disp 'TODO Does the rank change compared to using picture-averaged HMAX?'
%     figure; hold on; axis square; box on
%     %   scatter(sortHMAX, sortHMAXold)
%     scatter(sortinds, sortindsold)
%     [r,p] = corr(sortinds, sortindsold);
%     title(sprintf('%s, r = %1.2f, p = %1.3f', subj, r, p))
%     lsline
%     xlabel('Gaze-specific HMAX rank')
%     ylabel('Picture-averaged HMAX rank')
  otherwise
    error('Unknown gazespecificHMAX')
end

% fixed bin width, variable n trials per bin
disp 'make bins of trials based on hmax: fixed bin width, variable n trials per bin'
% ntrlperbin = 150 / nbins; % each subject has 150 trials
% cond_bins = reshape(sortinds, ntrlperbin, nbins); % dimord: TR trials cond
% hmax_bins = reshape(sortHMAX, ntrlperbin, nbins); % dimord: TR trials cond
source_bin = source; % binned
source_bin.pow = nan(size(source_bin.pow,1), nbins);
source_bin.powdimord = 'pos_freq'; % freq is hmax condition
if do_kstest; f = figure; f.Position = [  744          -9        1654        1059]; end
for ibin = 1:nbins
%   seldat = source.pow(source.inside,cond_bins(:,ibin),:); %seldat = source.pow(:,cond_bins(:,ibin),:); 
  seldat = source.pow(source.inside, bininds==ibin, :); %seldat = source.pow(:,cond_bins(:,ibin),:); 
  seldat =  seldat(:,:);
  if removeoutliers
    Z = zscore(seldat,1,2);
    %   Z = Z(Z~=0);  figure; histogram(Z(:))
    seldat_outliers = seldat(Z < -Z_thresh | Z > Z_thresh);
    source_bin.perc_outliers(ibin,1) = (numel(seldat_outliers) / numel(seldat(seldat~=0)))*100;
    seldat(Z < -Z_thresh | Z > Z_thresh) = NaN;
  end
  if do_kstest
    disp 'Kolmogorov-Smirnov test for normality'
    load(fullfile(PREIN, 'common_coords.mat'));
    
    % zscore each voxel
    %     ksdat = zscore(seldat(source.inside,:), 0, 2 );    % kstest tests for a standard normal distribution by default
    ksdat = seldat(common_coords,:);    % kstest tests for a standard normal distribution by default
    for ivox = 1:10
      ksdatsel = ksdat(ivox,:);
      %     ksdat = ksdat(:);    %ksdat = ksdat(1e5:2e5);
      %     [~,p] = kstest(ksdat);
      subplot(3,4,ivox);
      histogram(ksdatsel,15); %, 'Normalization', 'probability'
      %     cdfplot(ksdat);
      %     hold on
      %     x_values = linspace(min(ksdat),max(ksdat));
      %     plot(x_values,normcdf(x_values,0,1),'r-')
      %     legend('Empirical CDF','Standard Normal CDF','Location','best')
      %     xlim([-3 3])
      %     xlabel('GLM beta weight (Z-score)');
      %     ylabel('Cumulative frequency')
      %     title(sprintf('bin %d, kstest p = %g', ibin, p))
      %       xlim([-2.5e3 2.5e3])
      title(sprintf('bin %d', ibin))
      xlabel('GLM beta weight');
      ylabel('Frequency')
    end
  end
  switch BOLDvar_measure
    case 'nanstd'
      source_bin.pow(source_bin.inside,ibin) = nanstd(seldat,1,2); % take SD across 5 trials, 5 TR's each
    case 'iqr'
      source_bin.pow(source_bin.inside,ibin) = iqr(seldat,2); % take IQR across 5 trials, 5 TR's each
  end
%   source_bin.freq(ibin) = nanmean(hmax_bins(:,ibin)); % use freq field for HMAX bin_No
  source_bin.freq(ibin) = nanmean(binedges(ibin:ibin+1)); % use freq field for HMAX bin_No
end
source = source_bin;

%% old with fixed N trls per bin
% disp 'make bins of trials based on hmax'
% ntrlperbin = 150 / nbins; % each subject has 150 trials
% cond_bins = reshape(sortinds, ntrlperbin, nbins); % dimord: TR trials cond
% hmax_bins = reshape(sortHMAX, ntrlperbin, nbins); % dimord: TR trials cond
% source_bin = source; % binned
% source_bin.pow = nan(size(source_bin.pow,1), nbins);
% source_bin.powdimord = 'pos_freq'; % freq is hmax condition
% if do_kstest; f = figure; f.Position = [  744          -9        1654        1059]; end
% for ibin = 1:nbins
%   seldat = source.pow(source.inside,cond_bins(:,ibin),:); %seldat = source.pow(:,cond_bins(:,ibin),:); 
%   seldat =  seldat(:,:);
%   if removeoutliers
%     Z = zscore(seldat,1,2);
%     %   Z = Z(Z~=0);  figure; histogram(Z(:))
%     seldat_outliers = seldat(Z < -Z_thresh | Z > Z_thresh);
%     source_bin.perc_outliers(ibin,1) = (numel(seldat_outliers) / numel(seldat(seldat~=0)))*100;
%     seldat(Z < -Z_thresh | Z > Z_thresh) = NaN;
%   end
%   if do_kstest
%     disp 'Kolmogorov-Smirnov test for normality'
%     load(fullfile(PREIN, 'common_coords.mat'));
%     
%     % zscore each voxel
%     %     ksdat = zscore(seldat(source.inside,:), 0, 2 );    % kstest tests for a standard normal distribution by default
%     ksdat = seldat(common_coords,:);    % kstest tests for a standard normal distribution by default
%     for ivox = 1:10
%       ksdatsel = ksdat(ivox,:);
%       %     ksdat = ksdat(:);    %ksdat = ksdat(1e5:2e5);
%       %     [~,p] = kstest(ksdat);
%       subplot(3,4,ivox);
%       histogram(ksdatsel,15); %, 'Normalization', 'probability'
%       %     cdfplot(ksdat);
%       %     hold on
%       %     x_values = linspace(min(ksdat),max(ksdat));
%       %     plot(x_values,normcdf(x_values,0,1),'r-')
%       %     legend('Empirical CDF','Standard Normal CDF','Location','best')
%       %     xlim([-3 3])
%       %     xlabel('GLM beta weight (Z-score)');
%       %     ylabel('Cumulative frequency')
%       %     title(sprintf('bin %d, kstest p = %g', ibin, p))
%       %       xlim([-2.5e3 2.5e3])
%       title(sprintf('bin %d', ibin))
%       xlabel('GLM beta weight');
%       ylabel('Frequency')
%     end
%   end
%   switch BOLDvar_measure
%     case 'nanstd'
%       source_bin.pow(source_bin.inside,ibin) = nanstd(seldat,1,2); % take SD across 5 trials, 5 TR's each
%     case 'iqr'
%       source_bin.pow(source_bin.inside,ibin) = iqr(seldat,2); % take IQR across 5 trials, 5 TR's each
%   end
%   source_bin.freq(ibin) = nanmean(hmax_bins(:,ibin)); % use freq field for HMAX bin_No
% end
% source = source_bin;

%% make PLS sesssiondata structure, prepare important fields
% 1. standard stuff: TODO turn into function?
tmp=[];
% tmp.st_coords = find(source.inside); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
load(fullfile(PREIN, 'common_coords.mat'));
tmp.st_coords = find(common_coords); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
switch PLStype
  case 'taskPLS'
    tmp.st_datamat = transpose(source.pow(tmp.st_coords,:)); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
    tmp.session_info.datamat_prefix = 'SDbold_vs_HMAX'; %[subj '_' pattern];%stores common      %datamat prefix
    tmp.behavdata = [];
    tmp.behavname = {};%no behav names or data just yet, but set up the %field anyway...
  case 'behavPLSvsdprime'
    if isnumeric(binsubtract)
      tmp.st_datamat = transpose(source.pow(tmp.st_coords, binsubtract(1))) - transpose(source.pow(tmp.st_coords, binsubtract(2))); % highest - lowest BOLD variability
    else
      dat = transpose(source.pow(tmp.st_coords, :));
      for i = 1:length(dat)
        fit = polyfit( transpose(1:nbins), dat(:,i), 1 );
        tmp.st_datamat(1,i) = fit(1);
      end
    end
    tmp.session_info.datamat_prefix = 'SDboldHMAX_vs_dprime'; %[subj '_' pattern];%stores common      %datamat prefix
    tmp.behavname = {'dprime'};%no behav names or data just yet, but set up the %field anyway...
    load(behavfile); % behavior comes out
    subjind = behavior.participants.participant_id == subj;
    tmp.behavdata = mean(behav.test(subjind).dprime);
  case 'behavPLSvsDDM'
    if isnumeric(binsubtract)
      tmp.st_datamat = transpose(source.pow(tmp.st_coords, binsubtract(1))) - transpose(source.pow(tmp.st_coords, binsubtract(2))); % highest - lowest BOLD variability
    else
      dat = transpose(source.pow(tmp.st_coords, :));
      for i = 1:length(dat)
        fit = polyfit( transpose(1:nbins), dat(:,i), 1 );
        tmp.st_datamat(1,i) = fit(fitcoeff); % fit is in descending powers
      end
    end
    tmp.session_info.datamat_prefix = 'SDboldHMAX_vs_DDM'; %[subj '_' pattern];%stores common      %datamat prefix
    tmp.behavname = {'drift'};%no behav names or data just yet, but set up the %field anyway...
    load(behavfile); % behav comes out
    subjind = behavior.participants.participant_id == subj;
    tmp.behavdata = mean(behavior.ddmNiels.v(subjind));
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

disp(outfile_source)
save(outfile_source, 'source')

% save PLS sesdat
disp(outfile_sesdat)
save(outfile_sesdat, '-struct','tmp')
% % save([subj '_' pattern '_BfMRIdatamat.mat'], '-struct','tmp'); %,'-mat'


%   if ismac
%     tmp = source_bin;
%     tmp.powdimord = 'pos';
%     %       vol=300;
%     %   tmp.anatomy = tmp.anatomy(:,:,:,vol);
%     tmp.pow = nanmean(tmp.pow(:,:),2);
%     cfg=[];
%     cfg.method = 'ortho'; % slice ortho glassbrain vertex
%     cfg.funparameter = 'pow';
%     cfg.funcolorlim = 'zeromax';% [-300 300];
%     ft_sourceplot(cfg, tmp)
%   end

% %   source.time = 1:5;
% %   source.freq = 1;
% allsource(isub) = source_bin;
% % end
