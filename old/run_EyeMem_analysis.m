restoredefaultpath
if ismac
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB'; % '/mnt/homes/home022/nkloost1';
else
    basepath = '/home/mpib/kloosterman/MATLAB/'; % '/mnt/homes/home022/nkloost1';
end
addpath(genpath(fullfile(basepath, 'EyeMem_analysis')))
% addpath(fullfile(basepath, 'tools', 'fieldtrip-20150803'))
addpath(fullfile(basepath, 'tools', 'fieldtrip-20170611'))
addpath(fullfile(basepath, 'tools', 'qsub_tardis'))
addpath(fullfile(basepath, 'tools', 'gencode'))
addpath(genpath(fullfile(basepath, 'tools', 'plotting')))
addpath(fullfile(basepath, 'tools', 'custom_tools', 'stats')) 
addpath(fullfile(basepath, 'tools', 'custom_tools', 'plotting')) 

ft_defaults

addpath(fullfile(basepath, 'EyeMem_analysis'))

% EM_convertMRIdata
% EM_convertEYEdata
% EM_copyBEHAVdata
% EM_copyPHYSIOdata


%% subject_batch = 
subject_batch = EM_generate_subject_batch();
% % % subject_batch =
% % % 1x101 struct array with fields:
% % %     SUBJ
% % %     dateofbirth
% % %     agegroup
% % %     responsemapping
% % %     scandate
% % %     gender
% % %     handedness


behav = EM_analysebehavioral(subject_batch);

%% plot behavior
EM_plotbehavioral(subject_batch, behav)

%% create trialinfo and evfiles for pls analysis
Eyemem_sorttrials



