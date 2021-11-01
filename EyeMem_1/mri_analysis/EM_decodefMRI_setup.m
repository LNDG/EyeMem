function [accuracy] = EM_decodefMRI_setup()
% run from runMIBmeg_analysis
% Decode HMAX bin from fMRI data

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

% fun2run = @EM_getBOLDvar;
fun2run = @EM_decodefMRI;

% timreq = 60; %in minutes per run
timreq = 5; 
memreq = 8000; % in MB
% memreq = 12000; % in MB torque

PREIN = fullfile(basepath, 'variability', 'LOTO');
PREOUT = PREIN;

% eventshift = 5;
% latency = [1 5];
% latency = [-2 11]; %until peak in many areas
% latency = 'all'; 

overwrite = 1;

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.removeERP = 'yes';
if strcmp(cfg.removeERP, 'yes')
  PREIN = fullfile(PREOUT, 'ERPremoved');
end
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;

cfglist = {};

agedirs = {'YA' 'OA'};
outname = '1_decodingrpara'; % decode HMAX from single trial MSE
mkdir(fullfile(PREOUT, agedirs{1}), outname)
mkdir(fullfile(PREOUT, agedirs{2}), outname)

for iage = 1:2
  cd(fullfile(PREIN, agedirs{iage}))
  
  % SUBJ= [9, 11:59, 61:69, 71,72, 74:101]; % TODO specify further?
  SUBJ = dir('source_sub-*');
%   if ismac
% %     SUBJ = dir('sub-38*');
%     SUBJ = dir('source_sub-09*');
%   end
  
  for isub = 1:length(SUBJ)
    disp(SUBJ(isub).name)
%     mrifile = dir( fullfile(PREIN, agedirs{iage}, SUBJ(isub).name, '*.nii') );
%     if isempty(mrifile)
%       fprintf('mri %s not found\n', SUBJ(isub).name);
%       continue
%     end
%     cfg.mrifile = fullfile(mrifile.folder, mrifile.name);
%     cfg.eyefile = fullfile(PREINeye, agedirs{iage}, sprintf('eye_%s.mat', SUBJ(isub).name));
    cfg.infile = fullfile(PREIN, agedirs{iage}, SUBJ(isub).name);
    cfg.outfile = fullfile(PREOUT, agedirs{iage}, outname, SUBJ(isub).name);
    
    cfg.subjno = SUBJ(isub).name;
%     cfg.eventshift = eventshift; % shift by X s
%     cfg.latency = latency;
    
%     if ~exist(cfg.outfile, 'file') || overwrite
      cfglist{end+1} = cfg;
%     end
  end
end

% cfglist = cfglist(1)
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running %s for %d cfgs\n', mfilename, length(cfglist))

if strcmp(backend, 'slurm')
%   options = '-D. -c1'; % --gres=gpu:1 --partition gpu
  options = '-D. -c1 --partition quick'; % 
else
  options =  '-l nodes=1:ppn=1'; % torque %-q testing or gpu
end

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');

if strcmp(backend, 'local')
  accuracy = cellfun(fun2run, cfglist);
else
  accuracy = qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
    'StopOnError', true, 'UniformOutput', true, 'backend', backend, 'options', options);
end
disp(accuracy)
try
save(fullfile(PREOUT, 'accuracy.mat'), 'accuracy')
catch
end

