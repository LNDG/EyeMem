function EM_3makeMEANbold_datamat()

% Create the mean datamats based on the info mats.

% pn.plstoolbox   = [pn.root, 'analyses/B4_PLS_preproc2/T_tools/'];  addpath(genpath(pn.plstoolbox));

pn.root     = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem';
pn.trialinfo   = fullfile(pn.root, 'preproc/eye'); % eyedata.mat: location of info regarding run trialinfo matrices
pn.BOLDpath = fullfile(pn.root, 'variability/5TRspertrial');
pn.PLSfiles = fullfile(pn.root, 'variability/PLS_HMAXvsSDBOLD'); 

cd(pn.PLSfiles); % files will be output in the CD directory

% N = 43 YA + 53 OA;
IDs={'sub-09', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-61', 'sub-62', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-71', 'sub-72', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-100', 'sub-101'}; % 'sub-70', 'sub-73', 

for indID = 1:numel(IDs)
    disp(['Processing subject ', IDs{indID}, '.']);
%     BATCHPATH = fullfile(pn.PLSfiles, IDs{indID}, '_task_PLS_info.mat';
    
    BATCHPATH = fullfile(pn.PLSfiles, sprintf('%s_task_PLS_info.mat', IDs{indID} ));
    
    fprintf('Creating PLS_data_mat for subject %s \n', IDs{indID});
    batch_plsgui(BATCHPATH);
    fprintf('Finished creating PLS_data_mat for subject %s \n', IDs{indID});
end
