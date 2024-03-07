function [data, datainfo] = EM_eye_analysis(cfg)

plotit = 1;

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

datainfo = [];
datainfo.sampgaps = []; % 6 runs

alldata = {}; data_trial_cropped = {};
for irun = 1:length(edflist)
  disp 'convert edf to asc'
  [~,eyename] = fileparts( edflist(irun).name );
  filename_eye = sprintf('%s.asc', eyename);
  disp(filename_eye)
  if ~exist(filename_eye)
    system(sprintf('%s -y %s', edf2asc, edflist(irun).name )); %% convert edf to asc, overwrite
  end
  disp('preprocess eye data')
  cfg = [];
  cfg.dataset          = filename_eye;
  cfg.montage.tra      = eye(4);
  cfg.montage.labelorg = {'1', '2', '3', '4'};
  cfg.montage.labelnew = {'EYE_TIMESTAMP', 'EYE_HORIZONTAL', 'EYE_VERTICAL', 'EYE_DIAMETER'};
  data = ft_preprocessing(cfg); 
  
  if ismac && plotit
    cfg=[];    cfg.channel =[2 3]; ft_databrowser(cfg, data)
    figure; plot(data.trial{1}(2,:), data.trial{1}(3,:))
  end
  
  datainfo.sampgaps = [datainfo.sampgaps unique(diff(data.trial{1}(1,:)))];
  if data.fsample ~= 1000
    disp(datainfo.sampgaps)
    error('data is missing from recording! gaps are (in samples):')
% %     continue % TODO fix this
%     % put nans were missing data is    
%     nandata = nan(4, datainfo.sampgaps(2));
%     hiccupstart = find(diff(data.trial{1}(1,:)) > 1);
%     newdata = [data.trial{1}(:,1:hiccupstart-1) nandata data.trial{1}(:,hiccupstart:end)];
%     newdata(1,:) = (1:size(newdata,2)) + data.trial{1}(1,1); % new sample channel
%     data.time{1} = 0:1/data.fsample:((size(newdata,2)-1)/data.fsample);
%     data.trial{1} = newdata;
%     data.fsample = 1000; % 3 runs affected
  end
  
  disp('interpolate blinks, add blinks, saccades and fixations as chans')
  hdr = ft_read_header(filename_eye); %, 'headerformat', 'eyelink_asc');
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
  
  disp('low pass filter eye data') %only for pupils
  cfg = [];
  cfg.lpfilter = 'yes';
  cfg.lpfreq = 6; % cf de gee 2019 biorxiv
  cfg.lpfiltord = 3;
  cfg.channel = 'EYE_DIAMETER';
  data_pupil = ft_preprocessing(cfg, data);
  data.trial{1}(4,:) = data_pupil.trial{1}; % put filtered data back in with other chans
  
  disp('define trials')
  cfg=[];
  cfg.headerfile = filename_eye;
  cfg.runinfo = runinfo;
  cfg.trialdef.trg = 'stim';
  cfg.trialdef.begtim = -3;
  cfg.trialdef.endtim = 10;
  cfg.event = event;
  cfg.fsample = data.fsample;
  cfg.trialfun = 'EM_sortTrials_Marija';
  cfg = ft_definetrial(cfg); % make trl matrix
  
  if all(cfg.trl(1:3) == [0 0 0]); disp 'TODO sort out resting state'; continue; end
  
  %make trials
  data = ft_redefinetrial(cfg, data);
  
  if ismac && plotit
    cfg2=[];
    cfg2.channel = 'EYE_SACCADES';  %{'EYE_TIMESTAMP'  'EYE_HORIZONTAL'  'EYE_VERTICAL'  'EYE_DIAMETER'  'EYE_BLINKS'  'EYE_SACCADES'}
    ft_databrowser(cfg2, data) %for plotting
    timelock = ft_timelockanalysis([], data);
    timelock.avg = timelock.var; % to look at variability
    ft_databrowser(cfg2, timelock)
  end
  
  % Select only viewing period
  cfg=[];
  cfg.latency = [0 5];
  data_viewing = ft_selectdata(cfg, data);
  
  disp 'detect microsaccades within viewing period, per trial'
  %% Detect fixations in each trial, detect microsaccades within each fixation
  % TODO plot MS in fixations, also velocity?
  if  plotit
    close all
    f = figure; f.Position =[  1          58        1920         919 ];
  end
  microsaccades = []; 
  for itrial = 1:length(data_viewing.trial)
    cfg=[];
    cfg.trials = itrial;
    data_trial = ft_selectdata(cfg,data_viewing);
    data_trial.sampleinfo = [1 length(data_trial.trial{1})];
    
    disp 'Detect fixations'
    fixation_detection = 'DM'; % EyeLink triggers or ft_detect_movement DM
    switch fixation_detection % make trl matrix for fixations as trials
      case 'EL'
        data_trial.trial{1}(7,[1 end]) = 0; % start and end with no fix        
        fixonsets = find(diff( data_trial.trial{1}(7,:)) == 1); % -100
        fixoffsets = find(diff( data_trial.trial{1}(7,:)) == -1); % +100
        trl = [fixonsets(:) fixoffsets(:) zeros(length(fixoffsets(:)),1)];
        trl = trl(2:end,:); % remove 1st fixation because imposed above
      case 'DM' % detect fixations using Engbert and Kliegl method
        cfg=[];
        cfg.method = 'velocity2D';
        cfg.channel = { 'EYE_HORIZONTAL'  'EYE_VERTICAL' };
        cfg.velocity2D = [];
        cfg.velocity2D.mindur = 10; % minimum *saccade* duration, can be short
        cfg.velocity2D.velthres = 30; % lower = lower fixation time
        cfg.velocity2D.kernel   = [ones(1,16) zeros(1,8) -ones(1,16)].*(data_trial.fsample/6);% vector 1 x nsamples, kernel to compute velocity (default = [1 1 0 -1 -1].*(data.fsample/6);
        [~, movement] = ft_detect_movement(cfg, data_trial);
        fix_bool = true(5001,1);
        for im = 1:size(movement,1)
          fix_bool(movement(im,1):movement(im,2)) = false; % fixation is false at saccades
        end
        fix_bool([1 end]) = false; % to have start and end fixation
        fixonsets = find(diff(fix_bool) == 1); % -100
        fixoffsets = find(diff(fix_bool) == -1); % +100
        trl = [fixonsets(:) fixoffsets(:) zeros(length(fixoffsets(:)),1)];
        trl = trl(2:end,:); % remove 1st fixation because imposed above
        data_trial.trial{1}(7,:) = fix_bool;
        % TODO put fixation bool also in data?
    end

    if isempty(trl)
      microsaccades.movement{itrial} = [NaN NaN NaN];
      microsaccades.count(itrial,:) = 0;
      microsaccades.velocity(itrial,:) = NaN; % NK edit ft_detect_movement to get peak velocity
      continue
    end
    
    % keep 0-5 s trials with new fixation channel
    data_trial_cropped{end+1} = data_trial;
    
    disp 'Make "trials" from fixations'
    cfg=[]; 
    cfg.trl = trl;
    data_fix = ft_redefinetrial(cfg, data_trial);
    
    disp 'Detect microsaccades within fixations'
    cfg=[];
    cfg.method = 'velocity2D';
    cfg.channel = { 'EYE_HORIZONTAL'  'EYE_VERTICAL' };
    cfg.velocity2D = [];
    cfg.velocity2D.mindur = 9; % Gao et al: 12 ms, Port et al (visual search) 9 ms
    cfg.velocity2D.velthres = 6; % SDs from median velocity? 6 default
    cfg.velocity2D.kernel = [ones(1,8) zeros(1,4) -ones(1,8)].*(data_trial.fsample/6);% vector 1 x nsamples, kernel to compute velocity (default = [1 1 0 -1 -1].*(data.fsample/6);
    [~, movement] = ft_detect_movement(cfg, data_fix);
    
    if isempty(movement)
      movement = [0 0 0];
    end
    microsaccades.movement{itrial} = movement;
    microsaccades.count(itrial,:) = size(movement,1);
    microsaccades.velocity(itrial,:) = mean(movement(:,3)); % NK edit ft_detect_movement to get peak velocity
    if  plotit
      subplot(5,6,itrial);
      xdat = data_trial.trial{1}(2,:); xdat(~fix_bool) = NaN;
      ydat = data_trial.trial{1}(3,:); ydat(~fix_bool) = NaN;
      plot(xdat, ydat); hold on
      xdat = data_trial.trial{1}(2,:); xdat(fix_bool) = NaN;
      ydat = data_trial.trial{1}(3,:); ydat(fix_bool) = NaN;
      plot(xdat, ydat); hold on
      set(gca,'Ydir','reverse'); title(itrial)
      % plot microsaccades
      if all(movement>0)
        plot(data_viewing.trial{itrial}(2,movement(:,1)), data_viewing.trial{itrial}(3,movement(:,1)), 'x', 'MarkerSize', 20)
        plot(data_viewing.trial{itrial}(2,movement(:,2)), data_viewing.trial{itrial}(3,movement(:,2)), 'x', 'MarkerSize', 20)
      end
      if itrial==1;      legend({'Fixations', 'Saccades', 'MS start', 'MS end'}); end
            
%       % plot single fixations TODO add MS in red
%       figure;
%       for ifix = 1:length(data_fix.trial)
%         subplot(6,6,ifix);
%         plot(data_fix.trial{ifix}(2,:), data_fix.trial{ifix}(3,:)); hold on
%       end

    end
  end
  if  plotit
    mkdir(fullfile(PREOUT, 'MSplots'))
    saveas(f, fullfile(PREOUT, 'MSplots', sprintf('%s', eyename)), 'png')
  end
  
  %%
  data.trialinfo(:,15) = microsaccades.count;
  data.trialinfo(:,16) = microsaccades.velocity;
  
  disp 'down sample eye data'
  cfg=[];
  cfg.resample = 'yes';
  cfg.resamplefs = 100;
  cfg.detrend = 'no';
  %look at data here
  data = ft_resampledata(cfg, data);
  
  alldata{end+1} = data;
  
end

cfg=[];
cfg.keepsampleinfo = 'no';
data = ft_appenddata(cfg, alldata{:});
data_trial_cropped = ft_appenddata(cfg, data_trial_cropped{:});
clear alldata

disp(outfile)
save(outfile, 'data', 'data_trial_cropped')

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


