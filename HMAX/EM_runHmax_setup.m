function EM_runHmax_setup()
% run on the cluster

if ismac  
%   basepath = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/study_information/D_paradigm/';
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/D_paradigm';
  
  backend = 'local';
  compile = 'no';
else
%   basepath = '/home/mpib/LNDG/EyeMem/study_information/D_paradigm/';
  basepath = '/home/mpib/kloosterman/projectdata/eyemem/D_paradigm';
  %     backend = 'slurm';
  backend = 'slurm';
  %     backend = 'local';
  compile = 'no';
end
timreq = 1000; %in minutes per run
memreq = 4000; % in MB

PREIN = fullfile(basepath, 'stimuli_640x480');
PREOUT = fullfile(basepath, 'stimuli_640x480');
mkdir(PREOUT)

overwrite = 1;

cd(basepath)
load(fullfile(basepath, 'scanner_stimuli.mat'))
scanner_stimuli

%make cells for each subject, to analyze in parallel
cfg = [];      cfglist = {};
for icond = 1:length(scanner_stimuli.conditions)
  
  cfg.pictures_path = fullfile(PREIN, scanner_stimuli.conditions{icond});
%   cfg.saveFolder = fullfile(PREIN, scanner_stimuli.conditions{icond});
  cfg.saveFolder = fullfile(PREIN, 'hmax');
  cfg.condition = scanner_stimuli.conditions{icond};
  cfg.piclist = {scanner_stimuli.picnames{:,icond}}';
  
  cfglist = [cfglist cfg];
end
mkdir(cfg.saveFolder)
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running EM_runHmax for %d cfgs\n', length(cfglist))

if strcmp(backend, 'slurm')
  options = '-D. -c4'; % --gres=gpu:1
else
  options =  '-l nodes=1:ppn=3'; % torque %-q testing or gpu
end

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');
if strcmp(compile, 'yes')
  fun2run = qsubcompile(@EM_runHmax, 'toolbox', {'signal', 'stats'}); %
  %   fun2run = qsubcompile({@MEG2afc_preproc @sortTrials_MEGhh_2afc @interpolate_blinks}, ...
  %       'executable', 'run_kloosterman_master_p10908_b18.sh'); % compiled function
else
  fun2run = @EM_runHmax;
end

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist)
  return
end

qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
  'StopOnError', false, 'backend', backend, 'options', options);
