function EM_runvartbx_matlabbatch_setup()
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

fun2run = @EM_runvartbx_matlabbatch;
% fun2run = @EM_getBOLDvar_LOTO;

% timreq = 60; %in minutes per run
timreq = 130*60; % for LOTO
memreq = 20000; % in MB

PREIN = fullfile(basepath, 'variability2/5TRspertrial/jobs');

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;

SUBJ = dir('sub-*sdmodel.mat');

cfglist = {};
for isub = 1:length(SUBJ)
  disp(SUBJ(isub).name)
  cfg.batchfile = fullfile(PREIN, SUBJ(isub).name);
  cfglist{end+1} = cfg;
end

% cfglist = cfglist(1)
% cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm')
%   options = '-D. -c1'; % --gres=gpu:1 --partition gpu
  options = '-D. -c2 --partition long'; % --gres=gpu:1
else
  options =  '-l nodes=1:ppn=1'; % torque %-q testing or gpu
  setenv('TORQUEHOME', 'yes')
end

mkdir('~/qsub'); cd('~/qsub');

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist);
else
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
end

