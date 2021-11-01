function EM_pls_BOLDvsdprime_setup()
% PLS-preprocess data and collect all subj on the cluster, cast them to local machine

fun2run = @EM_pls_BOLDvsdprime;

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
else
  basepath = '/mnt/beegfs/home/kloosterman/projectdata/eyemem/';
  %   backend = 'slurm';
  backend = 'torque';
  %   backend = 'local';
  compile = 'no';
end
timreq = 10; %in minutes per run
memreq = 4000; % in MB

% analysis settings
load participantinfo.mat % TODO make this reliable

PREIN = fullfile(basepath, 'variability', 'VOIsel');
% PREINeye = fullfile(basepath, 'preproc', 'eye');


overwrite = 1;


%make cells for each subject, to analyze in parallel
cfg = [];
cfg.behavfile =  fullfile(basepath, 'preproc/behavior', 'Eyemem_behavior.mat');

cfglist = {};
agedirs = {'YA' 'OA'};
for iage = 1:2
  cfg.PREIN = fullfile(PREIN, agedirs{iage})

  cd(cfg.PREIN)
  PREOUT = fullfile(cfg.PREIN, 'PLS'); %
  mkdir(PREOUT)
  subjlist = dir('source*.mat');
  for isub = 1:length(subjlist)
    tok = tokenize(subjlist(isub).name, '_');
    [~, subj] = fileparts(tok{2});
    cfg.subj = subj;
    
    cfg.PREOUT = PREOUT;
    
    cfg.sourcefile = fullfile(cfg.PREIN, subjlist(isub).name);
    
    %     cfg.outfile_source = fullfile(PREOUT, 'source', subjlist(isub).name); % binned source data to source folder
    
    cfg.outfile_sesdat = fullfile(PREOUT, [subj '_BfMRIsessiondata.mat']); % PLS _BfMRIsessiondata
    
    if ~exist(cfg.outfile_sesdat, 'file') || overwrite
      cfglist{end+1} = cfg;
    end
   
  end
end
% cfglist = cfglist(2)
% cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm')
  options = '-D. -c1'; % --gres=gpu:1
else
  options =  '-l nodes=1:ppn=1'; % torque %-q testing or gpu
end

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist);
else
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
end

