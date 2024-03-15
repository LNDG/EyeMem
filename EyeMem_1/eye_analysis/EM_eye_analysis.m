function EM_eye_analysis(cfg)

if ismac
  edf2asc = '';
else
  edf2asc = '/mnt/beegfs/home/LNDG/EyeMem/tools/custom_tools/eyelink/linux/edf2asc'; %'/home/mpib/LNDG/EyeMem/tools/custom_tools/eyelink/linux/edf2asc';
end

PREIN =   cfg.PREIN;
PREOUT =  cfg.PREOUT;
subjno =  cfg.subjno;
edflist = cfg.edflist;
outfile = cfg.outfile;

disp(subjno)
cd(PREIN)

% fixation HMAX analysis
disp 'load HMAX file'
HMAXfolder = fullfile(fileparts(PREIN), 'D_paradigm', 'stimuli_640x480', 'hmax');
hmaxlist=dir(fullfile(HMAXfolder, '*.mat' ));
hmaxdat = {};
for ih = 1:length(hmaxlist)
  load(fullfile(hmaxlist(ih).folder, hmaxlist(ih).name));
  hmaxdat{ih} = hmaxout;
end

alldata = {};
for irun = 1:length(edflist)
  [~,eyename] = fileparts( edflist(irun).name );
  filename_asc = [eyename '.asc'];
  disp(filename_asc)
  filename_mat = [eyename '.mat'];
  if exist(filename_mat)
    disp 'Loading mat file'; disp(filename_mat)
    load(filename_mat)
  else
    disp 'convert edf to asc'
    if ~exist(filename_asc)
      system(sprintf('%s -y %s', edf2asc, edflist(irun).name )); %% convert edf to asc, overwrite
    end
    disp('preprocess eye data')
    cfg = [];
    cfg.dataset          = filename_asc;
    cfg.montage.tra      = eye(4);
    cfg.montage.labelorg = {'1', '2', '3', '4'};
    cfg.montage.labelnew = {'EYE_TIMESTAMP', 'EYE_HORIZONTAL', 'EYE_VERTICAL', 'EYE_DIAMETER'};
    data = ft_preprocessing(cfg);
    disp(filename_mat)
    save(filename_mat, 'data')
  end
  
  plotit = 0;
  if ismac && plotit
    cfg=[];    cfg.channel =[2 3]; ft_databrowser(cfg, data)
    figure; plot(data.trial{1}(2,:), data.trial{1}(3,:))
  end
  
  disp('interpolate blinks, add blinks, saccades and fixations as chans')
  hdr = ft_read_header(filename_asc); %, 'headerformat', 'eyelink_asc');
  data = interpolate_blinks(hdr, data); % TODO add microsaccades? getting channel 5 and 6
  
  if ismac && plotit
    cfg=[];    cfg.channel =[2 3]; ft_databrowser(cfg, data)
    figure; plot(data.trial{1}(2,:), data.trial{1}(3,:))
  end
  
  event = get_event(hdr); % subfunction below
  
  disp('get runinfo')
  [runinfo] = get_runinfo(event, subjno); % subfunction below
  
  if ismac && plotit
    cfg = [];
    cfg.event = event;
    cfg.channel = 'EYE_TIMESTAMP';
    cfg.preproc.demean = 'yes';
    ft_databrowser(cfg, data)
  end
  
  %   disp('low pass filter eye data') %only for pupils
  %   cfg = [];
  %   cfg.lpfilter = 'yes';
  %   cfg.lpfreq = 6; % cf de gee 2019 biorxiv
  %   cfg.lpfiltord = 3;
  %   cfg.channel = 'EYE_DIAMETER';
  %   data_pupil = ft_preprocessing(cfg, data);
  %   data.trial{1}(4,:) = data_pupil.trial{1}; % put filtered data back in with other chans
  
  disp('define trials')
  cfg=[];
  cfg.headerfile = filename_asc;
  cfg.runinfo = runinfo;
  cfg.trialdef.trg = 'stim';
  cfg.trialdef.begtim = -3;
  cfg.trialdef.endtim = 10;
  cfg.event = event;
  cfg.fsample = data.fsample;
  cfg.trialfun = 'EM_sortTrials_Marija';
  cfg = ft_definetrial(cfg); % make trl matrix  
  if all(cfg.trl(1:3) == [0 0 0]); disp 'TODO sort out resting state'; continue; end    
  data = ft_redefinetrial(cfg, data); %make trials
  
  data.trialinfo = array2table(data.trialinfo, 'VariableNames', ...
    {'exp_phase','pic_category', 'picno_study', 'picno_quiz', 'pic_quiz_oldnew', ...
    'response', 'resp_button', 'RT_quizpic', 'fMRI_trigger', 'HMAX_C1', 'HMAX_C2', 'runno', ...
    'mem_at_test', 'RT_at_test', 'trial_counter' });
  
  ntrials = size(data.trialinfo,1);
  varTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
  varNames = ["fixdur_mean", "fixdur_std", "fix_drift", "sacc_dur_mean", "sacc_dur_std", 'sacc_distance', ...
    "HMAX_fix", "HMAX_fix_weighted", "HMAX_fix_lookregion_mean", "HMAX_fix_lookregion_max", "HMAX_fix_lookregion_std", "Microsacc_count", "Microsacc_velocity"];
  sz = [ntrials length(varNames)];
  eyeinfo = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
  
  if ismac && plotit
    cfg2=[];
    cfg2.channel = 'EYE_SACCADES';  %{'EYE_TIMESTAMP'  'EYE_HORIZONTAL'  'EYE_VERTICAL'  'EYE_DIAMETER'  'EYE_BLINKS'  'EYE_SACCADES'}
    ft_databrowser(cfg2, data) %for plotting
    timelock = ft_timelockanalysis([], data);
    timelock.avg = timelock.var; % to look at variability
    ft_databrowser(cfg2, timelock)
  end
  
  disp 'shift gaze to remove gray boundary'
  xshift = -(1024-640)/2;  yshift = -(768-480)/2;
  chan_shifts = [0 xshift yshift 0 0 0 0]'; % Xgaze and Ygaze in chan2 and 3
  data.trial = cellfun(@(x) x + chan_shifts, data.trial, 'Uni',0) ; 

  %% Trial loop
  plotit = 1;
  if plotit
    close all
    f = figure; f.Position =[  1          58        1920         919 ];
  end
  
  for itrial = 1:ntrials
    cfg=[];
    cfg.trials = itrial; % select current trial
    cfg.latency = [0 5]; % select only viewing period
    data_trial = ft_selectdata(cfg, data);
    data_trial.sampleinfo = [1 length(data_trial.trial{1})];
    
    disp 'Detect fixations'
    fixation_detection = 'DM'; % EyeLink triggers or ft_detect_movement DM
    switch fixation_detection % make trl matrix for fixations as trials
      case 'EL_fix'
        fix_bool = logical(data_trial.trial{1}(7,:))';
      case 'EL_sacc'
        fix_bool = not(logical(data_trial.trial{1}(6,:)))';
      case 'DM' % detect fixations using Engbert and Kliegl method
        cfg=[];
        cfg.method = 'velocity2D';
        cfg.channel = { 'EYE_HORIZONTAL'  'EYE_VERTICAL' };
        cfg.velocity2D = [];
        cfg.velocity2D.mindur = 10; % minimum *saccade* duration, can be short
        cfg.velocity2D.velthres = 4; % lower = lower fixation time
        cfg.velocity2D.kernel   = [ones(1,16) zeros(1,8) -ones(1,16)].*(data_trial.fsample/6);% vector 1 x nsamples, kernel to compute velocity (default = [1 1 0 -1 -1].*(data.fsample/6);
        [~, movement] = ft_detect_movement(cfg, data_trial);
        fix_bool = true(5001,1);
        for im = 1:size(movement,1)
          fix_bool(movement(im,1):movement(im,2)) = false; % fixation is false at saccades
        end
    end
    fix_bool([1 end]) = false; % to have start and end fixation
    fixonsets  = find(diff(fix_bool) == 1);
    fixoffsets = find(diff(fix_bool) == -1);
    if numel(fixonsets) < 2;    disp('No fixations found');    continue;       end
    
    disp 'Make "trials" from fixations'
    cfg=[];
    cfg.trl = [fixonsets(:) fixoffsets(:) zeros(length(fixoffsets(:)),1)];
    data_fix = ft_redefinetrial(cfg, data_trial);    

    disp 'reject fixations with blinks'
    cfg=[];
    cfg.trials = cellfun(@(x) not(any(x(5,:))), data_fix.trial); % not any blinks
    data_fix = ft_selectdata(cfg, data_fix);
    if isempty(data_fix.trial);  disp('No fixations without blinks found');    continue;       end
        
    disp 'reject fixations outside picture'
    fixloc = cellfun(@(x) round(mean(x(2:3,:),2)), data_fix.trial, 'Uni',0) ; % Xgaze and Ygaze in chan2 and 3
    fixloc = [fixloc{:}]';
    cfg=[];
    cfg.trials = fixloc(:,1) > 0 & fixloc(:,1) < 640 & fixloc(:,2) > 0 & fixloc(:,2) < 480;
    data_fix = ft_selectdata(cfg, data_fix);
    if isempty(data_fix.trial);  disp('No fixations within picture without blinks found');    continue;       end
    
    fixdur = cellfun(@(x) size(x,2), data_fix.trial)';
    eyeinfo.fixdur_mean(itrial,1) = mean(fixdur); % fixation duration average
    eyeinfo.fixdur_std(itrial,1) = std(fixdur); % fixation duration SD
    % Fixation drift, see Clark et al. PNAS 2022 TODO diffusion constant
    eyeinfo.fix_drift(itrial,1) = mean(cellfun(@(x) std(x(2,:),0,2), data_fix.trial)); % x y 
    eyeinfo.fix_drift(itrial,2) = mean(cellfun(@(x) std(x(3,:),0,2), data_fix.trial)); % y
    
    disp 'Make "trials" from saccades'
    cfg=[];
    cfg.trl = [fixoffsets(1:end-1,1) fixonsets(2:end,1) zeros(length(fixoffsets(1:end-1)),1)]; % note the flip
    data_sacc = ft_redefinetrial(cfg, data_trial);    
    if isempty(data_sacc.trial);  disp('No saccades found');    continue;       end

    disp 'reject saccades with blinks'
    cfg=[];
    cfg.trials = cellfun(@(x) not(any(x(5,:))), data_sacc.trial); % not any blinks
    data_sacc = ft_selectdata(cfg, data_sacc);
    if isempty(data_sacc.trial);  disp('No saccades without blinks found');    continue;       end

    saccdur = cellfun(@(x) size(x,2), data_sacc.trial);
    eyeinfo.sacc_dur_mean(itrial,1) = mean(saccdur); % saccade duration average in ms
    eyeinfo.sacc_dur_std(itrial,1) = std(saccdur); % saccade duration SD
    % get Eucledian distance travelled: x1 to x2, y1 to y2: 
    euclides = @(coords) sqrt( (max(coords(2,:)) - min(coords(2,:)))^2 + ...
      (max(coords(3,:)) - min(coords(3,:)))^2 ); % TODO check with pdist2
    eyeinfo.sacc_distance(itrial,1) = mean(cellfun(euclides, data_sacc.trial)); % mean saccade distance per trial
    
    %% Fixation HMAX analysis
    disp 'get HMAX data of pic shown'
    catind = data.trialinfo.pic_category(itrial, 1); %
    picno = data.trialinfo.picno_study(itrial, 1); %
    picind = hmaxdat{catind}.picno == picno;
    curhmax = hmaxdat{catind}.c1(:,:,picind);
    picdat = hmaxdat{catind}.picdat(:,:,picind);
    %     figure; imagesc(curhmax);
    
    disp 'resample gaze to HMAX map resolution'
    fixloc = cellfun(@(x) round(mean(x(2:3,:),2)), data_fix.trial, 'Uni',0) ; % Xgaze and Ygaze in chan2 and 3
    fixloc = [fixloc{:}]';
    cur_res = size(picdat);     desiredres = size(curhmax); %resolutions of pic and hmax
    fixloc_newres = round(fixloc ./ cur_res .* desiredres);       % convert XY coords to resolution of HMAX,
    fixloc_newres(fixloc_newres==0) = 1; % if rounding makes location 0
    disp 'get c1 HMAX vals at fixation locations' % TODO get several pixels, loop?
    hmax_at_fix = diag(curhmax(fixloc_newres(:,2), fixloc_newres(:,1))); % Note the flip: Yaxis in dim1 (rows), Xaxis in dim2 (columns): scatter and plot need x,y, with indexing it's the other way around. diag bc all combinations of indices are returned
    
    %     lookregion = [1,1,1;1,1,1;1,1,1]; % 1 pix around fixloc
    lookregion = [...
      0,1,1,1,0;
      1,1,1,1,1;
      1,1,1,1,1;
      1,1,1,1,1;
      0,1,1,1,0;]; % 3 pix around fixloc
    hmax_at_fix_lookregion = [];
    for ifix = 1:size(fixloc_newres,1)
      M = zeros(size(curhmax));
      M(fixloc_newres(ifix,2), fixloc_newres(ifix,1)) = 1; % location
      tmp = curhmax(conv2(M,lookregion,'same') > 0);
      hmax_at_fix_lookregion.mean(ifix,1) = mean(tmp); % take mean within lookregion
      hmax_at_fix_lookregion.max(ifix,1) = max(tmp); % take max within lookregion
      hmax_at_fix_lookregion.std(ifix,1) = std(tmp); % take std within lookregion
    end
    %     disp 'average over HMAX vals to get 1 val per trial'
    eyeinfo.HMAX_fix(itrial,1) = mean(hmax_at_fix);
    eyeinfo.HMAX_fix_weighted(itrial,1) = sum((hmax_at_fix .*  (fixdur / sum(fixdur)))) ;
    eyeinfo.HMAX_fix_lookregion_mean(itrial,1) = mean(hmax_at_fix_lookregion.mean);
    eyeinfo.HMAX_fix_lookregion_max(itrial,1) = mean(hmax_at_fix_lookregion.max);
    eyeinfo.HMAX_fix_lookregion_std(itrial,1) = mean(hmax_at_fix_lookregion.std);

    plotit=0;
    if ismac && plotit
      %%
      figure; subplot(2,1,1); hold on
      imagesc(picdat) % for plotting transpose and flipud??
      ax=gca; ax.YDir = 'reverse';       % The default EyeLink coordinates are those of a 1024 by 768 VGA display, with (0, 0) at the top left.
      scatter(fixloc(:,1),fixloc(:,2), 'r', 'filled'); % already shifted
      
      text(fixloc(:,1),fixloc(:,2), string(1:size(fixloc,1))); % already shifted
      subplot(2,1,2); hold on %       disp 'plot hmax and gaze'
      imagesc(curhmax);
      scatter(fixloc_newres(:,1),fixloc_newres(:,2), 'r', 'filled'); % already shifted
      ax=gca; ax.YDir = 'reverse';
      %
      curhmax2 = curhmax;
      for ifix = 1:length(fixloc_newres)
        curhmax2(fixloc_newres(ifix,2), fixloc_newres(ifix,1)) = 1; % indexing so 2,1
      end
      figure; imagesc(curhmax2); hold on
    end
        
    %% Detect microsaccades within fixations
    % TODO plot MS in fixations, also velocity?
    disp 'Detect microsaccades within fixations'
    cfg=[];
    cfg.method = 'velocity2D';
    cfg.channel = { 'EYE_HORIZONTAL'  'EYE_VERTICAL' };
    cfg.velocity2D = [];
    cfg.velocity2D.mindur = 9; % Gao et al: 12 ms, Port et al (visual search) 9 ms
    cfg.velocity2D.velthres = 6; % SDs from median velocity? 6 default
    cfg.velocity2D.kernel = [ones(1,8) zeros(1,4) -ones(1,8)].*(data_trial.fsample/6);% vector 1 x nsamples, kernel to compute velocity (default = [1 1 0 -1 -1].*(data.fsample/6);
    [~, movement] = ft_detect_movement(cfg, data_fix);
    if isempty(movement);      movement = [0 0 0];    end
    eyeinfo.Microsacc_count(itrial,1) = size(movement,1); % Microsaccade count in 17
    eyeinfo.Microsacc_velocity(itrial,1) = mean(movement(:,3)); % Microsaccade peak velocity average in 18: % NK edit ft_detect_movement to get peak velocity     
    
    plotit=1;
    if plotit
%       subplot(5,6,itrial);
      xdat = data_trial.trial{1}(2,:); xdat(~fix_bool) = NaN;
      ydat = data_trial.trial{1}(3,:); ydat(~fix_bool) = NaN;
      plot(xdat, ydat); hold on
      xdat = data_trial.trial{1}(2,:); xdat(fix_bool) = NaN;
      ydat = data_trial.trial{1}(3,:); ydat(fix_bool) = NaN;
      plot(xdat, ydat); hold on
      set(gca,'Ydir','reverse'); title(itrial)
      % plot microsaccades
      if all(movement>0)
        plot(data_trial.trial{1}(2,movement(:,1)), data_trial.trial{1}(3,movement(:,1)), 'x', 'MarkerSize', 20)
        plot(data_trial.trial{1}(2,movement(:,2)), data_trial.trial{1}(3,movement(:,2)), 'x', 'MarkerSize', 20)
      end
      if itrial==1;      legend({'Fixations', 'Saccades', 'MS start', 'MS end'}); end
      
      %       % plot single fixations TODO add MS in red
      %       figure;
      %       for ifix = 1:length(data_fix.trial)
      %         subplot(6,6,ifix);
      %         plot(data_fix.trial{ifix}(2,:), data_fix.trial{ifix}(3,:)); hold on
      %       end
      
    end
    disp('plot it?')
  end % itrial
  if  plotit
    mkdir(fullfile(PREOUT, 'MSplots'))
    saveas(f, fullfile(PREOUT, 'MSplots', sprintf('%s', eyename)), 'png')
  end
  
  data.trialinfo = [data.trialinfo eyeinfo]; % combine trialinfos
  
  %   disp 'down sample eye data'
  %   cfg=[];
  %   cfg.resample = 'yes';
  %   cfg.resamplefs = 100;
  %   cfg.detrend = 'no';
  %   alldata{end+1} = ft_resampledata(cfg, data);
  
  alldata{end+1} = data;
  
end % irun

cfg=[];
cfg.keepsampleinfo = 'no';
data = ft_appenddata(cfg, alldata{:});
clear alldata

disp(outfile)
save(outfile, '-struct', 'data')

if ismac && plotit
  cfg2=[];
  cfg2.channel = 'EYE_DIAMETER';  %{'EYE_TIMESTAMP'  'EYE_HORIZONTAL'  'EYE_VERTICAL'  'EYE_DIAMETER'  'EYE_BLINKS'  'EYE_SACCADES'}
  ft_databrowser(cfg2, data)
  timelock = ft_timelockanalysis([], data);
  timelock.avg = timelock.var; % to look at variability
  ft_databrowser(cfg2, timelock)
end

end

%% helper functions

function [event] = get_event(hdr)
% create event structure
msg = hdr.orig.msg;
dat = hdr.orig.dat;

evcell = cell(length(msg),1);
event=struct('type', evcell, 'sample', evcell, 'value', evcell, 'offset', evcell, 'duration', evcell );
for i=1:length(msg)
  strtok = tokenize(msg{i});
  event(i).type = strtok{3};
  smpstamp = find(dat(1,:) == str2double(strtok{2})); % find sample index of trigger in ascii dat
  if ~isempty(smpstamp)
    event(i).sample = smpstamp(1);
  else
    event(i).sample = nan; % give nan if no data was recorded at msg time
  end
  event(i).value = msg{i}; % event(i).value = [strtok{3:end}]; %trigger value: e.g. ResponseMIBOff
end
end

function [runinfo] = get_runinfo(event, subjno)
% return runinfo with subject specific info, also loading in partc age, sex
% etc
% e.g. MSG	28982913 Subj9 Hand2 studyph naturals1 Run1 Start

% get the subj14 string to extract info
%temp = tokenize(char( {event( find(strfind(['Subj' cfg.subjno] ,{event.value})) ).value} ));
temp =tokenize(char({event( find(strcmp(sprintf('Subj%d', subjno),{event.type})) ).value}));
if any(cellfun(@(x) strcmp(x, 'RestingState'), temp ))
  disp 'TODO: make trl resting state'
  trl = [0 0 0];
  return
end
runinfo=[];
%runinfo.subjno = str2num(cfg.subjno); % gives err; because we set cfg.subjno to 14 manually?
runinfo.subjno = subjno;
runinfo.hand = str2num(temp{4}(5));
runinfo.phase = temp{5}(1:end-2);
runinfo.runno = str2num(temp{7}(4));
runinfo.category = temp{6};
% get subject info from file
if ismac
  load('participantinfo.mat') % /Users/kloosterman/Dropbox/tardis_code/MATLAB/eyemem_analysis/participantinfo/participantinfo.mat
else
  load('/mnt/beegfs/home/LNDG/EyeMem/analysis_scripts/participantinfo/participantinfo.mat') %
end
runinfo.subjID = Participants.participant_id(subjno);
runinfo.age = Participants.age(subjno);
runinfo.agegroup = Participants.group(subjno);
runinfo.sex = Participants.gender(subjno);
runinfo.weight = Participants.weight(subjno);
% % load study and test data and add
% if ismac
%   load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
% else
%   load /mnt/beegfs/home/LNDG/EyeMem/data/behavior/Eyemem_behavior.mat
% end
% runinfo.behavior.study = behav.study(subjno);
% runinfo.behavior.test = behav.test(subjno);

end



    %% ori code to get HMAX at fixations
%     % ORI: get fixation on and offsets: do based on non saccade episodes 
%     fixations = data.trial{itrial}(6,:)<0.1; % makes it 1 during fixation
%     %       fixations = logical(data.trial{itrial}(7,:)); % makes it 1 during fixation
%     fixations(data.trial{itrial}(5,:) == 1) = 0; % blinks
%     fixations(1) = 0; % so we get a fixation start at begin
%     fixations(end) = 0; % so we get a fixation end at end
%     fixtrig = [ find(diff(fixations) == 1)' find(diff(fixations) == -1)'];
%     %       disp 'Drop first (trial starts with central fixation) and last fixation'
%     %       fixtrig = fixtrig(2:end-1,:);
%     
%     % get XY coords of fixations: average XY within fixations
%     nfix = size(fixtrig,1);
%     fixloc = NaN(nfix,2);
%     fixdur = NaN(nfix,1);
%     cur_res = size(picdat); % resolution of the pics
%     desiredres = size(curhmax);
%     %     fixmap = zeros([size(curhmax) nfix]);
%     fixloc_newres = NaN(nfix,2);
%     disp 'resample gaze to curhmax resolution'
%     xshift = (1024-640)/2;
%     yshift = (768-480)/2;
%     for ifix = 1:nfix
%       fixinds = fixtrig(ifix,1):fixtrig(ifix,2);
%       fixdur(ifix,1) = length(fixinds);
%       fixloc(ifix,:) = round(mean(data.trial{itrial}(2:3,fixinds),2)); % Xgaze and Ygaze in chan2 and 3
%       %       disp 'shift fixations, account for pic not fullscreen in scanner, but in the middle'
%       fixloc(ifix,1) = fixloc(ifix,1)-xshift;
%       fixloc(ifix,2) = fixloc(ifix,2)-yshift;
%       % convert XY coords to resolution of HMAX
%       fixloc_newres(ifix,:) = round(fixloc(ifix,:) ./ cur_res .* desiredres);
%     end
%     disp 'Drop fixations outside picture'
%     
%     validfix = NaN(size(fixloc_newres,1),2);
%     validfix(:,1) = fixloc_newres(:,1) > 0 & fixloc_newres(:,1) < desiredres(2);
%     validfix(:,2) = fixloc_newres(:,2) > 0 & fixloc_newres(:,2) < desiredres(1);
%     fixloc_newres = fixloc_newres(all(validfix,2),:); % also apply to resampled fix locations
%     fixloc = fixloc(all(validfix,2),:);
%     fixdur = fixdur(all(validfix,2));
%     nfix =length(fixdur);
%     
%     plotit=1;
%     if ismac && plotit
%       figure; hold on
%       % The default EyeLink coordinates are those of a 1024 by 768 VGA display, with (0, 0) at the top left.
%       imagesc(picdat) % for plotting transpose and flipud??
%       ax=gca; ax.YDir = 'reverse';
%       %       scatter(data.trial{itrial}(2,:)-xshift, data.trial{itrial}(3,:)-yshift, 'k'); hold on
%       scatter(data.trial{itrial}(2,fixations)-xshift, data.trial{itrial}(3,fixations)-yshift, 'g'); hold on
%       scatter(fixloc(:,1),fixloc(:,2), 'r', 'filled'); % already shifted
%       text(fixloc(:,1),fixloc(:,2), string(1:size(fixloc,1))); % already shifted
%       
%       %       xlim([0 1024]); ylim([0 768]); box on
%       xlim([0 640]); ylim([0 480]); box on
%       title(nfix)
%       
%       disp 'plot hmax and gaze'
%       figure;
%       imagesc(curhmax); hold on
%       scatter(fixloc_newres(:,1),fixloc_newres(:,2), 'r', 'filled'); % already shifted
%       ax=gca; ax.YDir = 'reverse';
%       
%       curhmax2 = curhmax;
%       for ifix = 1:nfix
%         curhmax2(fixloc_newres(ifix,2), fixloc_newres(ifix,1)) = 1;
%       end
%       figure; imagesc(curhmax2); hold on
%       
%       %       sc = scatter(fixind_newres(:,1), fixind_newres(:,2), 'filled');
%       %       sc.SizeData = 50;
%       %       sc.CData = [1 0 0];
%     end
%     
%     disp 'get c1 HMAX vals at fixation locations'
%     nfix = size(fixloc_newres,1);
%     hmax_at_fix=NaN(nfix,1);
%     for ifix = 1:nfix
%       hmax_at_fix(ifix,1) = curhmax(fixloc_newres(ifix,2), fixloc_newres(ifix,1)); % Note the flip: Yaxis in dim1 (rows), Xaxis in dim2 (columns): scatter and plot need x,y, with indexing it's the other way around
%     end
%     if isempty(hmax_at_fix)
%       continue
%     end
%     
%     disp 'average over HMAX vals to get 1 val per trial'
%     weightedmean = 0;
%     if weightedmean == 1
%       fixdur = fixdur / sum(fixdur);
%       hmax_at_fix_trl(itrial,:) = sum((hmax_at_fix .* fixdur)) ;
%     else
%       hmax_at_fix_trl(itrial,:) = mean(hmax_at_fix);
%     end
%     
%     % keep hmax and fix dur values to correlate: YA better track
%     % complexity, i.e. get it better?
%     hmax_at_fix_keep = [hmax_at_fix_keep; hmax_at_fix];
%     fixdur_keep =      [fixdur_keep; fixdur];
%     
%     %     end
%     
%     
    % TODO put in trialinfo
    %% END TO: eye_analysis

