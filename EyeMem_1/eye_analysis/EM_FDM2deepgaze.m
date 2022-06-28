function [maps] = EM_FDM2deepgaze()
% load in FDM and Deepgaze, correlate maps to get a measure of how well
% subjects track the Deepgaze model
% supports both Deepgaze and Hmax

if ismac
  %   basepath = '/Users/kloosterman/gridmaster2012/kloosterman/';
  basepath = '/Users/kloosterman/Dropbox/tardis_code/';
else
  %   basepath = '/home/mpib/kloosterman/'; %/mnt/beegfs/home/
  basepath = '/mnt/beegfs/home/kloosterman/'; % to avoid ft path problems
end

addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'npy-matlab')));
PREINeye = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/eye';
PREINdeepgaze = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/DeepGazeII';
PREINhmax = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/D_paradigm/stimuli_640x480/hmax';

PREOUT = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/eye';

% saliencymodel = 'deepgaze';
saliencymodel = 'hmax';

nbins_x = 3;
nbins_y = 3;
omit_centerAOI = 1; 

% nbins_x = 832-192;
% nbins_y = 624-144;
% time = 0.5:1:5;
% timewin = 1;

% works well:
time = 0.25:0.25:4.75;
timewin = 0.5;

% % quick:
% time = 1.25:2.5:4.75;
% timewin = 2.5;

% no timebins:
% time = 2.5;
% timewin = 5;

agegroups = {'YA' 'OA'};

%% 1. load deepgaze maps and downsample
% deepgaze = 
%           dat: [5×30×3×3 double]
%        dimord: 'cond_pic_xpos_ypos'
%     picnolist: [5×30 double]

cond_names = {'fractals' 'landscapes' 'naturals1' 'streets1' 'streets2'}; % numbering
deepgaze=[];
deepgaze.dat = NaN( 5, 30, nbins_x, nbins_y);
deepgaze.dimord = 'cond_pic_xpos_ypos';

switch saliencymodel
  case 'deepgaze'
    cd(PREINdeepgaze)
    for icond=1:5
      piclist = dir(sprintf('outfile_%s*.npy', cond_names{icond}));
      for ipic = 1:length(piclist)
        
        temp = readNPY(piclist(ipic).name)';
        %       figure; imagesc(temp')
        deepgaze.dat(icond,ipic, :,:) = imresize(temp, [nbins_x, nbins_y]);
        strtok = tokenize(piclist(ipic).name, '_');
        strtok = tokenize(strtok{3}, '.');
        deepgaze.picnolist(icond,ipic) = str2double(strtok{1}); % pic nr in filename
      end
    end
  case 'hmax'
    cd(PREINhmax)
    hmaxlist=dir(fullfile(PREINhmax, '*.mat' ));
    hmaxdat = {};
    for icond = 1:length(hmaxlist)
      load(fullfile(hmaxlist(icond).folder, hmaxlist(icond).name));
      for ipic = 1:size(hmaxout.c1, 3)
%               figure; imagesc(hmaxout.c1(:,:,ipic))
%               figure; imagesc(hmaxout.picdat(:,:,ipic))
        deepgaze.dat(icond,ipic, :,:) = imresize(hmaxout.c1(:,:,ipic)', [nbins_x, nbins_y]);
        deepgaze.picnolist(icond,ipic) = hmaxout.picno(ipic);
      end
    end
end

%% 2. FDMs: Load timelocks, make % time spent in each AOI per age group

fixdens_age = cell(2,1);
for iage = 1:2
  cd(fullfile(PREINeye, agegroups{iage}))
  subjlist = dir('eye*.mat');
  nsub = length(subjlist);
  for isub = 1:nsub %
    disp(isub)
    load(subjlist(isub).name)
    
    cfg = [];
    cfg.keeptrials = 'yes';
    timelock = ft_timelockanalysis(cfg,data);

    fixdens = [];
    fixdens.map = NaN(size(timelock.trialinfo,1), 5, nbins_x, nbins_y);    
    %        (x_edges = np.linspace(192, 832, num=456),y_edges = np.linspace(144, 624, num=342))
    fixdens.x_edges = linspace(192, 832, nbins_x+1); % +1 bc edges
    fixdens.y_edges = linspace(144, 624, nbins_y+1);
    fixdens.dimord = 'rpt_time_xpos_ypos';
    fixdens.time = time;

    fixdens.trialinfo = timelock.trialinfo;
    runinfo = ft_findcfg(data.cfg, 'runinfo');
    fixdens.behavior = runinfo.behavior;
    
    disp 'binarize blinks and saccade chans'
    for ichan=5:6 % blinks and saccades
      blink_smp = timelock.trial(:,ichan,:);
      blink_smp(blink_smp > 0.5) = 1;
      blink_smp(blink_smp < 0.5) = 0;
      timelock.trial(:,ichan,:) = blink_smp;
    end
    
    for itoi = 1:length(time)
      cfg=[];
      %     cfg.latency = [0 5]; % only select viewing period, make timeresolved
      cfg.latency = [time(itoi)-timewin/2 time(itoi)+timewin/2]; % only select viewing period, make timeresolved
      datsel = ft_selectdata(cfg, timelock);
      for itrial = 1:size(datsel.trial,1)
        %       disp 'remove b and s samples'
        validsmp = squeeze(datsel.trial(itrial,5,:) == 0 & datsel.trial(itrial,6,:) == 0);
        gazedat = transpose(squeeze(datsel.trial(itrial,2:3,validsmp)));
        
        if isempty(gazedat)
          warning('No gaze data found')
        elseif numel(gazedat) < 30
          warning('<30 samples remain')
          continue
        end
        [FDM, x1, y1, binX, binY] = histcounts2(gazedat(:,1), gazedat(:,2), fixdens.x_edges, fixdens.y_edges);
        %       FDM = flipud(FDM); %The default EyeLink coordinates are those of a 1024 by 768 VGA display, with (0, 0) at the top left
        %       FDM = fliplr(FDM); %The default EyeLink coordinates are those of a 1024 by 768 VGA display, with (0, 0) at the top left
        %       disp 'TODO track Nfixations outside of picture'
        %             figure; subplot(2,2,1); scatter(gazedat(:,1),gazedat(:,2)); xlim([192, 832]); ylim([144, 624])
        %             subplot(2,2,2); imagesc([192, 832], [144, 624], FDM'); ax=gca; ax.YDir = 'Normal'; % transpose for imagesc bc it puts dim1 in columns
        %             disp 'TODO check if origin eyelink and pic match'
        
        %       fixdens.map(itrial,:,:) = (FDM / length(gazedat)) * 100;
        fixdens.map(itrial,itoi,:,:) = (FDM / length(gazedat)) * 100;
        
      end
    end % itoi
    %     fixdens_age{iage}{isub} = fixdens;
    if all(isnan(fixdens.map(:)))
      disp('subject dropped')
      continue
    end
    fixdens_age{iage}{end+1} = fixdens;
  end
end


%% 3. match FDM and Deepgaze, correlate fixed effects
incbins = true(nbins_x, nbins_y);
if omit_centerAOI & nbins_x == 3 & nbins_y == 3;
  incbins(2,2) = false;
end
incbins = incbins(:);

maps = [];
for iage = 1:2  
  maps(iage).deepgaze = [];
  maps(iage).fixation = [];
  maps(iage).dimord = 'subj_rpt_pos';
  maps(iage).omit_centerAOI = omit_centerAOI;
  maps(iage).nbins_x = nbins_x;
  maps(iage).nbins_y = nbins_y;
  maps(iage).x_edges = fixdens.x_edges;
  maps(iage).y_edges = fixdens.y_edges;
  maps(iage).agegroup = agegroups{iage};
  maps(iage).time = time;
  maps(iage).saliencymodel = saliencymodel;
  for isub = 1:length(fixdens_age{iage})
    fixdens = fixdens_age{iage}{isub};

    for itrial=1:size(fixdens.map,1)
      % find out pic watched
      trialinfo = fixdens.trialinfo(itrial,:);
      condind = trialinfo(2);
      picwatched = trialinfo(3);
      
      picind = deepgaze.picnolist(condind,:) == picwatched;
      maps(iage).deepgaze(isub,itrial,:) = deepgaze.dat(condind,picind,incbins);      % collapse over location here
      
      for itoi=1:length(fixdens.time)
          maps(iage).fixation(isub,itoi,itrial,:) = squeeze(fixdens.map(itrial,itoi,incbins));
      end
    end
    
    % add behavior
    maps(iage).behavior(isub,1) = nanmean(fixdens.behavior.study.dprime);
    maps(iage).behavior(isub,2) = nanmean(fixdens.behavior.study.RT);
    maps(iage).behavior(isub,3) = sum(fixdens.trialinfo(:,13))/length(fixdens.trialinfo); % prop remembered
    remtrls = logical(fixdens.trialinfo(:,13));
    maps(iage).behavior(isub,4) = nanmean(fixdens.trialinfo(remtrls,14)); % RT remembered old trials
    maps(iage).behavior(isub,5) = nanmean(fixdens.trialinfo(:,14)); % RT all old trials, remembered and not remembered
    
    % correlate
    for itoi=1:length(fixdens.time)
      corrdat = [squeeze(maps(iage).fixation(isub,itoi,:)) maps(iage).deepgaze(isub,:)'];
      corrdat = corrdat(~isnan(corrdat(:,1)),:);
      if isempty(corrdat)
        warning('No data remain!')
        maps(iage).corr(isub,itoi) = NaN; % to do inspect subj 45 OA
        maps(iage).corrdat{isub,itoi} = [NaN NaN];
        continue
      end
%       corrtype = 'Spearman'; % Spearman
      corrtype = 'Pearson'; % Spearman
      maps(iage).corr(isub,itoi) = corr( corrdat(:,1), corrdat(:,2), 'type', corrtype); % Spearman
      maps(iage).corrdat{isub,itoi} = corrdat;
    end
    maps(iage).corrtype = corrtype;
  end
end

%% 4. corr FDM to Deepgaze match to memory 

for iage = 1:2
  maps(iage).corr2behav = NaN(size(maps(iage).behavior,2), length(fixdens.time));
  for itoi=1:length(fixdens.time)
    % get behavoi, corr across subj
    for ibehav = 1:size(maps(iage).behavior,2)
      corrdat = [maps(iage).corr(:,itoi), maps(iage).behavior(:,ibehav)];
      corrdat = corrdat(~isnan(corrdat(:,1)),:);
      corrdat = corrdat(~isinf(corrdat(:,1)),:);
      corrdat = corrdat(~isinf(corrdat(:,2)),:); 
      disp 'TODO fix time bins loop'
      maps(iage).corr2behav(ibehav,itoi) = corr( corrdat(:,1), corrdat(:,2), 'type', corrtype );
      maps(iage).corrdatbehav{ibehav,itoi} = corrdat;
    end
  end
end

out = fullfile(PREOUT, ['maps_gaze_vs_' saliencymodel '.mat']);
disp(out)
save(out, 'maps')
disp 'done'
