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
source.time = 1:nTRpertrial;
source.cfg = [];
source.trialinfo(:,end+1) = 1:150; % number trials to keep track

source_ori = source;
% remove evoked response: subtract within trial mean beta weight per trial
switch inducedortotalSD
  case 'induced'
    
    source.pow = source.pow - mean(source.pow,3);
    
    plotit = 0;
    if ismac && plotit
      %   source.pow = mean(source.pow(:,5:5:end),2);
      source.pow = std(source.pow(:,:),0,2);
      source_ori.pow = std(source_ori.pow(:,:),0,2);
      source.powdimord = 'pos';
      source_ori.powdimord = 'pos';
      
      cfg=[];
      cfg.funparameter = 'pow';
      cfg.method = 'ortho'; % slice ortho glassbrain vertex
      load colormap_jetlightgray.mat
      cfg.funcolormap = cmap;
      ft_sourceplot(cfg, source)
      ft_sourceplot(cfg, source_ori)
    end
end

%% inspection of fMRI, not used in further processing
plotit = 0;
if ismac && plotit
    cfg = [];
    cfg.funparameter  = 'pow';
    cfg.maskparameter = cfg.funparameter;
  %   cfg.maskparameter = ;
  %   cfg.colorlim      = [-3 3]; % or 'maxabs'
  %   cfg.opacitymap    = 'vdown';
  %   cfg.opacitylim    = [-3 3]; % or 'maxabs'
    ft_sourceplot(cfg, source)
%   plotdat = squeeze(mean(source.pow(source.inside, :,1)));
  plotdat = source.pow(source.inside, :,1);
  figure; imagesc(plotdat); colorbar
  figure; plot(plotdat)
  figure
  for i=1:nTRpertrial
    subplot(2,3,i)
    plot(plotdat(:,i))
    % figure; plot(mean(plotdat,2))
  end
  source2plot = source;
  source2plot.pow = source2plot.pow(source2plot.inside,:,:);
  source2plot.pos = source2plot.pos(source2plot.inside,:);
  cfg=[];
  cfg.method = 'summary';
  sourcedata = ft_rejectvisual(cfg, source2plot)
  %% Detect outlier trials and reject: var taken per voxel, only 5 data points per trial.... still looks meaningful 
  cfg = [];
  cfg.cov_cut    = [0, 98]; % not used with zscorecut
  cfg.badtrs     = [];
  cfg.bad_trials = [];
  cfg.method = 'maxmin_perct'; % zscorecut (abs(min)+1 threshold) or maxmin_perct (original)
  [selecttrials, cfg] = EM_ft_varcut3(sourcedata, cfg, ismac); %https://dx.doi.org/10.1101/795799
end

remove_artf_eegstyle = 0;
if remove_artf_eegstyle
  disp 'Remove trials with max var taken over voxels, var computed per voxel over time1. 3 Zscores'
  powvar = zscore(max(var(source.pow,1,3)));
  zscorelim = 3;
  if ismac
    figure; scatter(1:150,powvar);
    line([0 150], [zscorelim zscorelim])
  end
  cfg = [];
  cfg.trials = find(powvar < zscorelim);
  fprintf('%d trials removed with zscore > %d\n',  size(source.pow,2) - length(cfg.trials), zscorelim)
  source = ft_selectdata(cfg, source);
  
  % disp 'remove trials with 0 (sometimes happens for last trial)'
  examplevoxel = squeeze(source.pow(find(source.inside,1,'first'),:,:));
  cfg=[];
  cfg.trials = find(~any(examplevoxel==0, 2));
  fprintf('%d trials removed with zeros in them\n', size(source.pow,2) - length(cfg.trials))
  source = ft_selectdata(cfg, source);
end

ntrials = size(source.trialinfo,1);
validtrials = find(source.trialinfo(:,end));

%%
switch gazespecificHMAX
  case 'non-gazespecific' % bin based on overall HMAX, take SD over 5 trials    
    % sort onsets based on hmax TODO run for HMAX C2    
    [sortHMAX, sortinds] = sort(source.trialinfo(validtrials,10));  %hmax in 10, ascending, trial inds - 10 is c1median
    ntrlperbin = ntrials / nbins; % each subject has 150 trials
    
    bininds = repmat(1:nbins, ntrlperbin, 1);
    %         bininds = bininds(:);
    %         bininds = bininds(sortinds); % does this reorder? NO just indexing
    bininds = sortrows([sortinds, bininds(:)]);
    bininds = bininds(:,2);
    bininds(isnan(sortHMAX)) = NaN; % set outliers to nan
    
    binedges = sortHMAX([1:ntrlperbin:ntrials ntrials]);
    
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
    
    disp 'remove trials rejected based on BOLD outliers'
    cfg=[];
    cfg.trials = validtrials;
    data = ft_selectdata(cfg, data);
    
%     ntrials = length(data.trial);
    hmax_at_fix_trl = nan(ntrials,1);
    hmax_at_fix_keep = [];
    fixdur_keep = [];

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
      %       disp 'Drop first (trial starts with central fixation) and last fixation'
      %       fixtrig = fixtrig(2:end-1,:);
      
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
        text(fixloc(:,1),fixloc(:,2), string(1:size(fixloc,1))); % already shifted
        
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
        hmax_at_fix_trl(itrial,:) = median(hmax_at_fix);   
      end
      
      % keep hmax and fix dur values to correlate: YA better track
      % complexity, i.e. get it better?
      hmax_at_fix_keep = [hmax_at_fix_keep; hmax_at_fix];
      fixdur_keep =      [fixdur_keep; fixdur];
      
    end
%     figure; scatter(hmax_at_fix_keep, fixdur_keep); lsline; title(corr(hmax_at_fix_keep, fixdur_keep))
% %     dat = [log(hmax_at_fix_keep) log(fixdur_keep)];
% %     dat = dat(all(~isinf(dat),2),:);
% %     figure; scatter(dat(:,1), dat(:,2)); lsline; title(corr(dat(:,1), dat(:,2)))

    
    fprintf('%d trials without fixations found: ', sum(isnan(hmax_at_fix_trl)))
    if sum(isnan(hmax_at_fix_trl)) > 25
      warning('More than 50: skipping subject')
      return
    end

    dropoutliers=0;
    if dropoutliers
      disp 'drop HMAX outlier trials'
      %         figure; histogram(hmax_at_fix_trl, 100)
      [~,TF]=rmoutliers(hmax_at_fix_trl);
      fprintf('%d HMAX outliers found\n', sum(TF))
      hmax_at_fix_trl(TF) = NaN; % set outliers to nan
    end
    
    switch bintype
      case 'fixednbins'
        % continue with sorting
        [sortHMAX, sortinds] = sort(hmax_at_fix_trl);  % gaze-specific HMAX values
        %         [sortHMAXold, sortindsold] = sort(source.trialinfo(:,10));
 
        bininds = repmat(transpose(1:nbins), ceil(ntrials/nbins), 1);
        bininds = sort(bininds(1:ntrials));
        
        ntrlperbin = floor(ntrials / nbins); % each subject has 150 trials
% 
%         bininds = repmat(1:nbins, ntrlperbin, 1);
% %         bininds = bininds(:);        
        bininds = sortrows([sortinds, bininds(:)]);
        bininds = bininds(:,2);
        bininds(isnan(hmax_at_fix_trl)) = NaN; % set outliers to nan
        
        binedges = sortHMAX([1:ntrlperbin:ntrials ntrials]);

      case 'uniformbinwidth'
        [ntrlperbin,binedges,bininds] = histcounts(hmax_at_fix_trl,nbins);
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
TFkeep=0;
for ibin = 1:nbins
%   seldat = source.pow(source.inside, cond_bins(:,ibin),:); %seldat = source.pow(:,cond_bins(:,ibin),:); 
  seldat = source.pow(source.inside, bininds==ibin, :); %seldat = source.pow(:,cond_bins(:,ibin),:);   
  seldat =  seldat(:,:);
  hmaxperbin(ibin,1) = mean(hmax_at_fix_trl(bininds==ibin));
  
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
    case 'std'
      disp 'removing outliers and taking SD'
      %       source_bin.pow(source_bin.inside,ibin) = std(seldat,1,2); % take SD across 5 trials, 5 TR's each
      
      inside_ind = find(source_bin.inside);
      %       seldat = source.pow(source.inside, :, :); %seldat = source.pow(:,cond_bins(:,ibin),:);
      %       seldat = seldat(:,:);
      for i = 1:size(seldat,1)  % 1:5000:size(seldat,1)
        %         close all
        %         f=figure;  f.Position = [        1000         997        1335         341];
        %         subplot(1,2,1); plot(seldat(i,:))
        
        rmboldoutliers = 0;
        if rmboldoutliers
          [seldat_clean, TF] = rmoutliers(seldat(i,:)); % , 'mean' is > 3 SD's from the mean
          %         subplot(1,2,2); plot(seldat_clean); xlim([1 150])
        else
          seldat_clean = seldat(i,:);
          TF=0;
        end
        source_bin.pow(inside_ind(i),ibin) = std(seldat_clean); % take SD across 5 trials, 5 TR's each
        TFkeep = TFkeep + length(find(TF));
      end
    case 'iqr'
      source_bin.pow(source_bin.inside,ibin) = iqr(seldat,2); % take IQR across 5 trials, 5 TR's each
    case 'mse'
      disp 'compute MSE'
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
      cfg.toi = 2.5;
      cfg.timwin = 5;
      cfg.timescales = 1;
      cfg.filtmethod = 'no';
      cfg.recompute_r = 'perscale_toi_sp';
      cfg.coarsegrainmethod = 'pointavg';
      mse = ft_entropyanalysis(cfg, data);
      
      inside_ind = find(source_bin.inside);
      source_bin.pow(inside_ind,ibin) = mse.sampen;
  end
  %   source_bin.freq(ibin) = nanmean(hmax_bins(:,ibin)); % use freq field for HMAX bin_No
  source_bin.freq(ibin) = nanmean(binedges(ibin:ibin+1)); % use freq field for HMAX bin_No
  source_bin.perc_BOLDremoved = (TFkeep / numel(seldat)) * 100;
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
%       tmp.st_datamat = (transpose(source.pow(tmp.st_coords, binsubtract(1))) - transpose(source.pow(tmp.st_coords, binsubtract(2)))) ./ ...
%         transpose(source.pow(tmp.st_coords, binsubtract(2))) .* 100; % highest - lowest BOLD variability
      tmp.st_datamat = (transpose(source.pow(tmp.st_coords, binsubtract(1))) - transpose(source.pow(tmp.st_coords, binsubtract(2)))); % highest - lowest BOLD variability
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
tmp.perc_BOLDremoved = source_bin.perc_BOLDremoved;

if strcmp(gazespecificHMAX, 'gaze-specific')
  tmp.fixdur_keep = fixdur_keep;
  tmp.hmax_at_fix_keep = hmax_at_fix_keep;
  tmp.hmax_at_fix_trl = hmax_at_fix_trl; % to plot extracted saliency per trial
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

% %   source.time = 1:5;
% %   source.freq = 1;
% allsource(isub) = source_bin;
% % end
