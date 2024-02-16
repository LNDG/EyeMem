function EM_pls_SDbold_vs_HMAX_setup(PLSbehav)
% PLS-preprocess data and collect all subj on the cluster, cast them to local machine

if nargin==0
  PLSbehav = '';
end

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
memreq = 10000; % in MB

% analysis settings
PLStype = 'taskPLS'; 
% PLStype = 'behavPLSvsDDM'; % behavPLS_sdboldvsHmaxbins
% PLStype = 'behavPLS_sdboldvsHmaxbins'; % 

% PLStype = 'behavPLSvsSDT';
% PLSbehav = 'dprime';
% PLSbehav = 'criterion';

nbins = 5; % no of bins used for Hmax binning, 750 samples
gazespecificHMAX = 'gaze-specific';
% gazespecificHMAX = 'non-gazespecific'; 
BOLDvar_measure = 'std'; % iqr, std mse
bintype = 'fixednbins';   %fixednbins   uniformbinwidth
% bintype = 'uniformbinwidth';   %fixednbins   uniformbinwidth
removeoutliers = false;
Z_thresh = 3; % if removeoutliers NOT USED
do_kstest = 0; % NOT USED
inducedortotalSD = 'induced';

% binsubtract = [5 1]; % Also possible in psc which bins to subtract: % [5 1] is bin5-bin1 ONLY behavPLSvsdprime
binsubtract = nbins; % just 1 number at binno = no subtraction
% binsubtract = [5 1; 4 1; 3 1; 5 3; 4 3]; % which bins to subtract: % [5 1] is bin5-bin1
% % binsubtract = 'linearfit';  
fitcoeff = 1; % fit in descending powers: 1 = slope, 2 = intercept, for behavpls
% binsubtract = 'corrHmaxoverbins';  % corr across bins Hmax vs SDbold

load participantinfo.mat % TODO make this reliable

nTRpertrial = 5; % 1 for classic LSS
PREIN = fullfile(basepath, 'variability2', sprintf('%dTRspertrial', nTRpertrial), 'ftsource');
% PREIN = fullfile(basepath, 'variability2', 'ftsource');
% PREINeye = fullfile(basepath, 'preproc', 'eye');
HMAXfolder = fullfile(basepath, 'D_paradigm', 'stimuli_640x480', 'hmax');

overwrite = 1;
    
cd(PREIN);
subjlist = dir('source_sub*.mat');

% SUBJ= [9, 11:59, 61:69, 71,72, 74:101]; % TODO specify further?

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.behavfile =  fullfile(basepath, 'preproc/behavior', 'Eyemem_behavior.mat');
cfg.PLStype = PLStype;
cfg.PLSbehav = PLSbehav;
cfg.nbins = nbins;
cfg.removeoutliers = removeoutliers;
cfg.Z_thresh = Z_thresh;
cfg.do_kstest = do_kstest;
cfg.BOLDvar_measure = BOLDvar_measure;
cfg.gazespecificHMAX = gazespecificHMAX;
cfg.fitcoeff = fitcoeff;
cfg.bintype = bintype;
cfg.nTRpertrial = nTRpertrial;
cfg.inducedortotalSD = inducedortotalSD;

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
      binsubtractfolder = sprintf(binsubtract);
    end
    BOLDvar_binsfolder = sprintf('%s_%dbins', cfg.BOLDvar_measure, nbins);
    
    agefolder = Participants(Participants.participant_id == subj, :);     % give different outfolder for OA and YA
    
    if contains(PLStype, 'behav')
      PREOUT = fullfile(PREIN, inducedortotalSD, BOLDvar_binsfolder, bintype, PLStype, PLSbehav, gazespecificHMAX, sprintf('%s_fitcoeff%d', binsubtractfolder, fitcoeff), char(agefolder.group)); % 'SDbold_vs_HMAX'
%       PREOUT = fullfile(PREIN, BOLDvar_binsfolder, bintype, PLStype, PLSbehav, gazespecificHMAX, sprintf('%s_fitcoeff%d', binsubtractfolder, fitcoeff), char(agefolder.group)); % 'SDbold_vs_HMAX'
    else
      PREOUT = fullfile(PREIN, inducedortotalSD, BOLDvar_binsfolder, bintype, PLStype, gazespecificHMAX, char(agefolder.group)); % 'SDbold_vs_HMAX'
    end
    disp(PREOUT)
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
% cfglist = cfglist(66)
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm')
  options = '-D. -c1'; % --gres=gpu:1
else
  options =  '-l nodes=1:ppn=1'; % torque %-q testing or gpu
end

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); if ~ismac; cd('~/qsub'); end

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist);
else
  qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 5, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
end

