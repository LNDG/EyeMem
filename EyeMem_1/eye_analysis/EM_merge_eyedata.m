function [timelock] = EM_merge_eyedata()
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

eyedata=[];
load('/Users/kloosterman/Documents/GitHub/EyeMem/EyeMem_1/participantinfo/participantinfo.mat', 'Participants')
eyedata.participants = Participants;

PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsDDM/v/linearfit_fitcoeff1';
cd(fullfile(PREIN))
SUBJ = {};
subjlist = dir('sub*_BfMRIsessiondata.mat');
for isub=1:length(subjlist)
  tmp = tokenize(subjlist(isub).name, '_');
  SUBJ(isub) = tmp(1) ;
end
SUBJ = sort(SUBJ);

PREINeye = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/eye';

PREOUT = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/eye';
agegroups = {'YA' 'OA'};

%% load eye data

for iage = 1:2
  cd(fullfile(PREINeye, agegroups{iage}))
  subjlist = dir('eye*.mat');
  nsub = length(subjlist);
  temp = {};

  for isub = 1:nsub
    
    disp(isub)
    if ismember(subjlist(isub).name(5:end-4), SUBJ)
      load(subjlist(isub).name)
      
      cfg = [];
      cfg.keeptrials = 'no';
      temp{end+1} = ft_timelockanalysis(cfg,data);
    end
  end
  cfg=[];
  cfg.keepindividual = 'yes';
  cfg.parameter = 'avg';
  timelock{iage} = ft_timelockgrandaverage(cfg, temp{:})
end


out = fullfile(PREOUT, ['eye_timelock.mat']);
disp(out)
save(out, 'timelock')
disp 'done'
