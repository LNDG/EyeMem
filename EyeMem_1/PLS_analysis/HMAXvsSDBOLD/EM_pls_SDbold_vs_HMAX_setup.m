function EM_pls_SDbold_vs_HMAX_setup()
% PLS-preprocess data and collect all subj on the cluster, cast them to local machine
close all

fun2run = @EM_pls_SDbold_vs_HMAX;

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/'; %yesno or 2afc
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem/'; 
  backend = 'slurm';
%   backend = 'torque';
  %   backend = 'local';
  compile = 'no';
end

timreq = 10; %in minutes per run
memreq = 2000; % in MB

% analysis settings
PLStype = 'taskPLS'; % behavPLSvsdprime taskPLS behavPLSvsDDM
nbins = 7; % no of bins used for Hmax binning
BOLDvar_measure = 'iqr'; % iqr, nanstd
bintype = 'uniformbinwidth';
removeoutliers = false;
Z_thresh = 3; % if removeoutliers
do_kstest = 0;
gazespecificHMAX = 'gaze-specific'; % gaze-specific  non-gazespecific

% binsubtract = [5 1]; % which bins to subtract: % [5 1] is bin5-bin1 ONLY behavPLSvsdprime
% binsubtract = [5 1; 4 1; 3 1; 5 3; 4 3]; % which bins to subtract: % [5 1] is bin5-bin1
binsubtract = 'linearfit';
fitcoeff = 1; % fit in descending powers: 1 = slope, 2 = intercept, for behavpls

load participantinfo.mat % TODO make this reliable

PREIN = fullfile(basepath, 'variability', 'ftsource');
% PREINeye = fullfile(basepath, 'preproc', 'eye');
HMAXfolder = fullfile(basepath, 'D_paradigm', 'stimuli_640x480', 'hmax');

overwrite = 1;
    
cd(PREIN);
subjlist = dir('source*.mat');

% SUBJ= [9, 11:59, 61:69, 71,72, 74:101]; % TODO specify further?

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.behavfile =  fullfile(basepath, 'preproc/behavior', 'Eyemem_behavior.mat');
cfg.PLStype = PLStype;
cfg.nbins = nbins;
cfg.removeoutliers = removeoutliers;
cfg.Z_thresh = Z_thresh;
cfg.do_kstest = do_kstest;
cfg.BOLDvar_measure = BOLDvar_measure;
cfg.gazespecificHMAX = gazespecificHMAX;
cfg.fitcoeff = fitcoeff;
cfg.bintype = bintype;

cfglist = {};

for isub = 1:length(subjlist)
  tok = tokenize(subjlist(isub).name, '_');
  [~, subj] = fileparts(tok{2});
  cfg.subj = subj;
  for ibinsubtr = 1:size(binsubtract,1)
    cfg.binsubtract = binsubtract(ibinsubtr,:);
    if isnumeric(binsubtract)
      binsubtractfolder = sprintf('bin%d-bin%d', binsubtract(ibinsubtr,:));
    else
      binsubtractfolder = sprintf('linearfit');
    end
    BOLDvar_binsfolder = sprintf('%s_%dbins', cfg.BOLDvar_measure, nbins);
    
    agefolder = Participants(Participants.participant_id == subj, :);     % give different outfolder for OA and YA
    
    if contains(PLStype, 'behav')
      PREOUT = fullfile(basepath, 'variability', 'ftsource', PLStype, BOLDvar_binsfolder, bintype, sprintf('%s_fitcoeff%d', binsubtractfolder, fitcoeff), gazespecificHMAX, char(agefolder.group)); % 'SDbold_vs_HMAX'
    else
      PREOUT = fullfile(basepath, 'variability', 'ftsource', PLStype, BOLDvar_binsfolder, bintype, gazespecificHMAX, char(agefolder.group)); % 'SDbold_vs_HMAX'
    end
    mkdir(PREOUT)
    mkdir(fullfile( PREOUT, 'source' ))
    cfg.PREOUT = PREOUT;

    eyefolder = char(agefolder.group);
    eyefolder = [upper(eyefolder(1)) 'A']; % I mean agefolder OA vs YA
    cfg.eyefile = fullfile(basepath, 'preproc', 'eye', eyefolder, ['eye_' subj '.mat']);
    
    cfg.sourcefile = fullfile(PREIN, subjlist(isub).name);

    cfg.outfile_source = fullfile(PREOUT, 'source', subjlist(isub).name); % binned source data to source folder
    
    cfg.outfile_sesdat = fullfile(PREOUT, [subj '_BfMRIsessiondata.mat']); % PLS _BfMRIsessiondata
    
    cfg.HMAXfolder = HMAXfolder;
    if ~exist(cfg.outfile_sesdat, 'file') || overwrite
      cfglist{end+1} = cfg;
    end
    if ~exist(fileparts(cfg.outfile_sesdat))
      mkdir(fileparts(cfg.outfile_sesdat))
    end
    
  end
end
% cfglist = cfglist(2)
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

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
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 5, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
end

