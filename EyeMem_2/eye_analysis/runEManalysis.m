%restoredefaultpath
if ismac
  %   basepath = '/Users/kloosterman/gridmaster2012/kloosterman/';
  basepath = '/Users/terlau/';
else
  %   basepath = '/home/mpib/kloosterman/'; %/mnt/beegfs/home/
  %basepath = '/mnt/beegfs/home/kloosterman/'; % to avoid ft path problems
  basepath = '/home/mpib/LNDG/EyeMem/'; % to avoid ft path problems '/mnt/beegfs/home/LNDG/EyeMem/'
end


%addpath(genpath(fullfile(basepath, 'MATLAB', 'eyemem_analysis')));
%addpath(fullfile(basepath, 'MATLAB', 'tools', 'hmaxMatlab'));
addpath(fullfile(basepath, 'tools', 'fieldtrip')) % cloned on 13 09 19
%addpath(fullfile(basepath, 'MATLAB', 'tools', 'spm12'))
% addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'pls')))
%addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'PLS_rank'))) %
% only when running Spearman corr behavior PLS

% addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'pls_mat_in')))
% addpath(genpath('/Volumes/FB-LIP/Projects/StateSwitch/dynamic/data/mri/task/analyses/B4_PLS_preproc2/T_tools/pls'))
addpath(genpath(fullfile(basepath, 'tools/custom_tools/plotting')))

%addpath(fullfile(basepath, 'MATLAB', 'tools', 'NIFTI_toolbox'))
% addpath(fullfile(basepath, 'MATLAB', 'tools', 'qsub_tardis')) % qsub_tardis_slurmpreview
addpath(fullfile(basepath, 'tools', 'qsub_tardis_slurmpreview'))% qsub_tardis_slurmpreview
% addpath(genpath('/Volumes/LNDG/Programs_Tools_Scripts/data_processing_repo/PLS_repo/PLS_toolbox_modifications/LNDG2018_OnlyTXT/Pls'))


%addpath(fullfile(basepath, 'MATLAB', 'tools', 'hmaxMatlab'));
addpath(fullfile(basepath,'tools', 'fieldtrip')) 

addpath(fullfile(basepath,'tools', 'custom_tools'))% interpolateblinks
addpath(fullfile(basepath, 'LNDG', 'EyeMem', 'EyeMem_2', 'eye_analysis'))
addpath(fullfile(basepath, 'LNDG'))



ft_defaults

 
%% set up PLS analysis for HMAX vs SDBOLD
% 1 make model txt file for for pls toolbox
%EM_1makemodeltxtfile(behavoi) % TODO behavoi with HMAX vals. output: batch_files mats
% 2 gather info for pls toolbox
%EM_2makebatch_files() % make cond / data info for each subj
% 3 prepare input data for pls toolbox
%EM_3makeMEANbold_datamat() % runs batch_plsgui fun on batch_files
% 4 substitute SD for MEAN values
%EM_4makeSDbold_datamat() 
% 5. Start the analysis in MATLAB with: _batch\_plsgui('nameOfYourMATfile.mat')_
 
%batch_plsgui('modeltxtfile.txt') %% TODO which file???

%% run Hmax on stim pics
% EM_runHmax_setup()

%% eye pupil and gaze analysis
EM_eye_analysis_setup()
% %% behavioral analysis (txt files)
% % behav = EM_analysebehavioral()
%

%% set up vartbx for single-TR GLM
% % run 1 and 2 on tardis! TODO move all this from usuer kloosterman to LNDG
% % 1 make SPM model mat
% create_1stlevel_singletrial()
% % 2 make vartbx model mat
% create_varTBX_model_allstimuli_5sres_fortardis()
% % 3 run shell script at master node
% unix('runVarTbx_eyemem_allstimuli_5sres.sh')

%% Put fMRI data in fieldtrip source format
% EM_mri_to_ftsource_setup()
% 
%% set up PLS analysis for HMAX vs SDBOLD fieldtrip source way
% EM_pls_SDbold_vs_HMAX_setup()

%%  set up PLS analysis  SDbold YA vs OA 
%EM_pls_OAvsYA_setup()

% % make model specification txt file and run the PLS analysis
% EM_runPLSanalysis(PLSfolder)


%% Compare FDMs to deepgaze maps
% [maps] = EM_FDM2deepgaze();
% EM_FDM2deepgaze_plot(maps);

%% OLD remove soon
% %% set up PLS analysis for HMAX vs SDBOLD Julian way
% % 1 make model txt file for for pls toolbox
% EM_1makemodeltxtfile(behavoi) % TODO behavoi with HMAX vals. output: batch_files mats
% % 2 gather info for pls toolbox
% EM_2makebatch_files() % make cond / data info for each subj
% % 3 prepare input data for pls toolbox
% EM_3makeMEANbold_datamat() % runs batch_plsgui fun on batch_files
% % 4 substitute SD for MEAN values
% EM_5makeSDbold_datamat()
% % 5. Start the analysis in MATLAB with: _batch\_plsgui('nameOfYourMATfile.mat')_
%
% batch_plsgui('modeltxtfile.txt') %% TODO which file???
%
%
%
