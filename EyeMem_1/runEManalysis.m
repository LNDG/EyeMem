restoredefaultpath
if ismac
  %   basepath = '/Users/kloosterman/gridmaster2012/kloosterman/';
  basepath = '/Users/kloosterman/Dropbox/tardis_code/';
  addpath(genpath('/Users/kloosterman/Documents/GitHub/EyeMem/EyeMem_1'));
else
%     basepath = '/home/mpib/kloosterman/'; %/mnt/beegfs/home/
  basepath = '/mnt/beegfs/home/kloosterman/'; % to avoid ft path problems
  addpath(genpath('/mnt/beegfs/home/kloosterman/GitHub/EyeMem/EyeMem_1'));  
end

% addpath(genpath(fullfile(basepath, 'MATLAB', 'eyemem_analysis')));
addpath(fullfile(basepath, 'MATLAB', 'tools', 'hmaxMatlab'));
addpath(fullfile(basepath, 'MATLAB', 'tools', 'fieldtrip')) % cloned on 13 09 19
addpath(fullfile(basepath, 'MATLAB', 'tools', 'spm12'))
% addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'pls')))
addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'PLS_rank'))) %
% only when running Spearman corr behavior PLS

addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'mmse'))) %
addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'MVPA-Light/'))) %
startup_MVPA_Light

% addpath(genpath(fullfile(basepath, 'MATLAB', 'tools', 'pls_mat_in')))
% addpath(genpath('/Volumes/FB-LIP/Projects/StateSwitch/dynamic/data/mri/task/analyses/B4_PLS_preproc2/T_tools/pls'))
addpath(genpath(fullfile(basepath, 'MATLAB', 'tools/custom_tools/plotting')))

addpath(genpath(fullfile('/Users/kloosterman/Documents/GitHub/plotting-tools')))


addpath(fullfile(basepath, 'MATLAB', 'tools', 'custom_tools')) % interpolateblinks
addpath(fullfile(basepath, 'MATLAB', 'tools', 'NIFTI_toolbox'))
ft_defaults
% addpath(fullfile(basepath, 'MATLAB', 'tools', 'qsub_tardis')) % qsub_tardis_slurmpreview
addpath(fullfile(basepath, 'MATLAB', 'tools', 'qsub_tardis_slurmpreview'))% qsub_tardis_slurmpreview
% addpath(genpath('/Volumes/LNDG/Programs_Tools_Scripts/data_processing_repo/PLS_repo/PLS_toolbox_modifications/LNDG2018_OnlyTXT/Pls'))

%% behavioral analysis (txt files)
% % behavior = EM_analysebehavioral()
% if ismac
%   load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat')
% end
% EM_plotbehavioral(behavior)
%% run Hmax on stim pics, also save original pic and full HMAX maps
% EM_runHmax_setup()

%% eye pupil and gaze analysis
% EM_eye_analysis_setup()   % EM_sortTrials_Marija

% timelock=EM_merge_eyedata()

%% fMRI ANALYSIS:
%% set up vartbx for single-TR GLM (LEast squares, single)
% % run 1 and 2 on tardis!
% % 1 make SPM model mat
% create_1stlevel_singletrial()
% % 2 make vartbx model mat
% create_varTBX_model_allstimuli_5sres_fortardis()
% % 3 run shell script at master node
% unix('runVarTbx_eyemem_allstimuli_5sres.sh')

% EM_runvartbx_matlabbatch_setup % run with qsubcellfun!

%% Put fMRI beta_series in fieldtrip source format
% EM_mri_to_ftsource_setup()
% EM_commoncoord_source_setup() % determine voxels present in all subjects

%% set up PLS analysis for HMAX vs SDBOLD fieldtrip source way
% EM_pls_SDbold_vs_HMAX_setup() % check settings in function!!

EM_pls_SDbold_vs_HMAX_setup('v')
EM_pls_SDbold_vs_HMAX_setup('a')
EM_pls_SDbold_vs_HMAX_setup('t')
EM_pls_SDbold_vs_HMAX_setup('z')
EM_pls_SDbold_vs_HMAX_setup('dc')

%%  set up PLS analysis  SDbold YA vs OA 
% EM_pls_OAvsYA_setup()

% % make model specification txt file and run the PLS analysis
% EM_runPLSanalysis(PLSfolder)

% plot results of gaze-specific HMAX vs IQRbold
% EM_plot_gazeHMAXvsSDbold

% split pls result into age groups, test interaction
% EM_plsgroupeffectANOVA()

%% NO GLM APPROACH: Cut out trials in continuous fMRI data
% EM_maketrials_fMRI_setup()

%% compute mean, std, and mse across trials, also LOTO
% EM_getBOLDvar_setup()

%% decode HMAX from single trial MSE 
% EM_decodefMRI_setup()

%% merge all subject voxel time courses in a struct array
% source = EM_merge_source();


%% Compare FDMs to deepgaze maps
% [maps] = EM_FDM2deepgaze();
% or 
% load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/eye/maps_gaze_vs_hmax.mat', 'maps')
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
