% function EM_eye_analysis_setup()
% run from runMIBmeg_analysis

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem';
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
  load('participantinfo.mat')
else
  basepath = '/mnt/beegfs/home/kloosterman/projectdata/eyemem/'; %'/home/mpib/LNDG/EyeMem/data/'; %yesno or 2afc
  backend = 'slurm';
  %     backend = 'torque';
%   backend = 'local';
  compile = 'no';
  load('/mnt/beegfs/home/kloosterman/GitHub/EyeMem/EyeMem_1/participantinfo/participantinfo.mat')
end

timreq = 180; %in minutes per run
memreq = 5000; % in MB

PREIN = fullfile(basepath, 'data'); 
PREOUT = fullfile(basepath, 'preproc', 'eye');
mkdir(fullfile(PREOUT, 'YA'))
mkdir(fullfile(PREOUT, 'OA'))

overwrite = 1;

SUBJ= [9:101]; % TODO specify further?
SUBJ= [20]; % TODO specify further?

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;
cfglist = {};

for isub = 1:length(SUBJ)
  
  edflist = dir(fullfile(PREIN, sprintf('S%dp1*.edf', SUBJ(isub))));
  if isempty(edflist)
    fprintf('S%dp1*.edf not found\n', SUBJ(isub))
    continue
  end
  cfg.subjno = SUBJ(isub);
  cfg.edflist = edflist;
  if Participants.group(SUBJ(isub)) == 'young'
    agegroup = 'YA';
  else
    agegroup = 'OA';
  end    
  if SUBJ(isub) < 10
    cfg.outfile = fullfile(PREOUT, agegroup, sprintf('eye_sub-0%d.mat', SUBJ(isub)));
  else
    cfg.outfile = fullfile(PREOUT, agegroup, sprintf('eye_sub-%d.mat', SUBJ(isub)));
  end
  if ~exist(cfg.outfile, 'file') || overwrite
    cfglist{end+1} = cfg;
  end
end

% cfglist = cfglist(2)
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm') 
    options = ' --cpus-per-task=1 '; % --gres=gpu:1  -D. -c1 
else
    options =  '-l nodes=1:ppn=2'; % torque %-q testing or gpu
end

%setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');
if strcmp(compile, 'yes')
  % this appenrently not used on slurm
  ft_hastoolbox('ctf', 1); % for loading ctf data
  ft_hastoolbox('eeglab', 1); % for ica
  fun2run = qsubcompile({@EM_eye_analysis @EM_sorttrials}, 'toolbox', {'signal', 'stats'}); %
  %   fun2run = qsubcompile({@EM_eye_analysis @sortTrials_MEGhh_2afc @interpolate_blinks}, ...
  %       'executable', 'run_kloosterman_master_p10908_b18.sh'); % compiled function
else
  % this is used on slurm!
  fun2run = @EM_eye_analysis;
end

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist);
else
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', false, 'UniformOutput', true, 'backend', backend, 'options', options);
end

