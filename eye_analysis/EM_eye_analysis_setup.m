% function EM_eye_analysis_setup()
% run from runMIBmeg_analysis

if ismac
  basepath = '/Volumes/LNDG/Projects/EyeMem/Laura/data'; %yesno or 2afc %path for where the data is
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
  load('/Users/kloosterman/Dropbox/tardis_code/MATLAB/eyemem_analysis/participantinfo/participantinfo.mat')
else
  basepath = '/mnt/beegfs/home/LNDG/EyeMem/data/'; %'/home/mpib/LNDG/EyeMem/data/'; %yesno or 2afc
  backend = 'slurm';
  %     backend = 'torque';
  %backend = 'local';
  compile = 'no';
  load('/home/mpib/kloosterman/MATLAB/eyemem_analysis/participantinfo/participantinfo.mat')
end

timreq = 60; %in minutes per run
memreq = 5000; % in MB

PREIN = fullfile(basepath, 'data'); 
PREOUT = fullfile(basepath, 'preproc', 'eye');
mkdir(fullfile(PREOUT, 'YA'))
mkdir(fullfile(PREOUT, 'OA'))

overwrite = 1;

SUBJ= [9:101]; % TODO specify further?

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
% cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm') 
    options = '-D. -c1 '; % --gres=gpu:1
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
  [data, datainfo] = cellfun(fun2run, cfglist);
  
  if ismac
    timelock = cellfun(@(x) ft_timelockanalysis([], x),  num2cell(data));
    timelock = num2cell(timelock);
    timelockavg = ft_timelockgrandaverage([],  timelock{:});
    
    cfg=[];
    cfg.channel = 'EYE_SACCADES'; %{'EYE_TIMESTAMP'  'EYE_HORIZONTAL'  'EYE_VERTICAL'  'EYE_DIAMETER'  'EYE_BLINKS'  'EYE_SACCADES'}
    cfg.preproc.demean = 'yes';
    ft_databrowser(cfg, timelockavg)
  end
% if it's running on slurm
else
  [data, datainfo] = qsubcellfun(fun2run, cfglist, 'memreq', memreq, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
 
  %disp(data)
  %save(fullfile(PREOUT, 'eyedata'), 'data', 'datainfo')
  
  % %     jobidarray = cellfun(@qsubfeval, repmat({@EM_eye_analysis}, size(cfglist)), cfglist, 'memreq', memreq, 'timreq', timreq*60, 'stack', 1, ... )
  % %       'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options)
  %   for i = 1:length(cfglist)
  % %     jobidarray{i} = qsubfeval(@myfunction, inputarg{i});
  %     jobidarray{i} = qsubfeval(@EM_eye_analysis, cfglist{i}, 'memreq', memreq, 'timreq', timreq*60, 'stack', 1, ... )
  %       'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
  %   end
  %   save jobidarray.mat jobidarray
  % exit
  %
  % % start MATLAB again
  % load jobidarray.mat
  % outputarg = cellfun(@qsubget, jobidarray, 'UniformOutput', false);
  
  %   jobidarray = qsubfeval(@EM_eye_analysis,  cfglist{1}, 'memreq', memreq, 'timreq', timreq*60, 'stack', 1, ... )
  %     'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options)
  
end
