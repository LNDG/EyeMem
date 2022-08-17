function EM_get_spm_residuals_setup()
% run from runMIBmeg_analysis

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
  %   backend = 'local';
  compile = 'no';
end
timreq = 2500; %in minutes per run
memreq = 10000; % in MB

nTRpertrial = 1; % 1 for classic LSS
PREIN = fullfile(basepath, 'variability2', sprintf('%dTRspertrial', nTRpertrial));
PREINeye = fullfile(basepath, 'preproc', 'eye');

PREOUT = fullfile(PREIN, 'ftsource', 'residuals');

mkdir(PREOUT)

overwrite = 1;

SUBJ= [9, 11:59, 61:69, 71,72, 74:101]; % TODO specify further?

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;
cfg.nTRpertrial = nTRpertrial;
cfglist = {};

for isub = 1:length(SUBJ)
  if SUBJ(isub) == 9
    subjstr = sprintf('sub-0%d', SUBJ(isub));
  else
    subjstr = sprintf('sub-%d', SUBJ(isub));
  end
  mripath = fullfile(PREIN, subjstr, 'spm');
  if SUBJ(isub) < 10
    cfg.eyefile = fullfile(PREINeye, sprintf('eye_sub-0%d.mat', SUBJ(isub)));
  else
    cfg.eyefile = fullfile(PREINeye, sprintf('eye_sub-%d.mat', SUBJ(isub)));
  end
  
  cfg.subjstr = subjstr;
  cfg.subjno = SUBJ(isub);
  cfg.mripath = mripath;
  cfg.outfile = fullfile(PREIN, subjstr, 'beta_series', sprintf('%s_sd_alltrials_residuals.nii', subjstr));
  if ~exist(cfg.outfile, 'file') || overwrite
    cfglist{end+1} = cfg;
  end
end

% cfglist = cfglist(2)
% cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm')
  options = '-D. -c1'; % --gres=gpu:1
else
  options =  '-l nodes=1:ppn=2'; % torque %-q testing or gpu
end

fun2run = @EM_get_spm_residuals;

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist);
else
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', 0, 'UniformOutput', true, 'backend', backend, 'options', options);
end

