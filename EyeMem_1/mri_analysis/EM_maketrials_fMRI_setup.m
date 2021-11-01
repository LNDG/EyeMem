function EM_maketrials_fMRI_setup()
% run from runMIBmeg_analysis
% select vols of interest from raw fMRI data, e.g. pic viewing TRs, and put
% in output file. Add time shift to account for HRF

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/'; %yesno or 2afc
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem/'; %yesno or 2afc
  backend = 'slurm';
%       backend = 'torque';
%     backend = 'local';
  compile = 'no';
end
timreq = 60; %in minutes per run
memreq = 15000; % in MB
% memreq = 12000; % in MB torque

PREIN = fullfile(basepath, 'preproc', 'mri');
PREINeye = fullfile(basepath, 'preproc', 'eye');

% PREOUT = fullfile(basepath, 'variability', 'VOIsel');
PREOUT = fullfile(basepath, 'preproc', 'mri_trials');

% eventshift = 5;
% latency = [1 5];
latency = [-2 11]; %until peak in many areas
% latency = 'all'; 

overwrite = 1;

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;
cfglist = {};

agedirs = {'YA' 'OA'};
mkdir(fullfile(PREOUT, agedirs{1}))
mkdir(fullfile(PREOUT, agedirs{2}))

for iage = 1:2
  cd(fullfile(PREIN, agedirs{iage}))
  
  % SUBJ= [9, 11:59, 61:69, 71,72, 74:101]; % TODO specify further?
  SUBJ = dir('sub-*');
%   if ismac
%     SUBJ = dir('sub-38*');
%   end
  
  for isub = 1:length(SUBJ)
    disp(SUBJ(isub).name)
    mrifile = dir( fullfile(PREIN, agedirs{iage}, SUBJ(isub).name, '*.nii') );
    if isempty(mrifile)
      fprintf('mri %s not found\n', SUBJ(isub).name);
      continue
    end
    cfg.mrifile = fullfile(mrifile.folder, mrifile.name);
    cfg.eyefile = fullfile(PREINeye, agedirs{iage}, sprintf('eye_%s.mat', SUBJ(isub).name));
    cfg.outfile = fullfile(PREOUT, agedirs{iage}, sprintf('source_%s.mat', SUBJ(isub).name));
    
    cfg.subjno = SUBJ(isub).name;
%     cfg.eventshift = eventshift; % shift by X s
    cfg.latency = latency;
    
    if ~exist(cfg.outfile, 'file') || overwrite
      cfglist{end+1} = cfg;
    end
  end
end

% cfglist = cfglist(2)
% cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm')
%   options = '-D. -c1'; % --gres=gpu:1 --partition gpu
  options = '-D. -c1 --gres=gpu:1 --partition gpu'; % 
else
  options =  '-l nodes=1:ppn=1'; % torque %-q testing or gpu
end

fun2run = @EM_maketrials_fMRI;

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist);
else
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
end

