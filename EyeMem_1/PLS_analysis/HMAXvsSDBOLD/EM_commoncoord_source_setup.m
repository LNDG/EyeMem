function EM_commoncoord_source_setup()
% Can just run it locally, very fast

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/'; %yesno or 2afc
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem/'; %yesno or 2afc
%   backend = 'slurm';
%       backend = 'torque';
    backend = 'local';
  compile = 'no';
end
timreq = 10; %in minutes per run
memreq = 2000; % in MB

nTRpertrial = 1; % 1 for classic LSS
PREIN = fullfile(basepath, 'variability2', sprintf('%dTRspertrial', nTRpertrial), 'ftsource');
% PREIN  = fullfile(basepath, 'variability2', 'ftsource');
cd(PREIN)

SUBJ= [9, 11:59, 61:69, 71,72, 74:101]; % TODO specify further?

cfglist = {};
cfglist{1}.subjlist = dir('source*.mat');
cfglist{1}.PREIN = PREIN;

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm')
  options = '-D. -c1'; % --gres=gpu:1
else
  options =  '-l nodes=1:ppn=2'; % torque %-q testing or gpu
end

fun2run = @EM_commoncoord_source;

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist);
else
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
end

