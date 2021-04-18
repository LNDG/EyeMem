function replace_mean_with_SD

% This function loads mean_datamat files created with PLSgui (from sd_niftis 
% that were generated with VarTbx). The current content of the
% datamat files is not meaningful and they are only used for providing the
% correct structure input for later analysis with PLSgui. The actual data
% is taken from sd_nifti files.

% created by MT, adapted from JQK

pn.root = '/home/mpib/LNDG/EyeMem/PLS/';
pn.meanFiles = [ pn.root, 'mean_datamats'];
pn.sdFiles = [pn.root, 'sd_datamats/'];
pn.sdniftis = '/home/mpib/LNDG/EyeMem/data/mri/task/variability/B_data/';
pn.tools = ['/home/mpib/LNDG/toolboxes/preprocessing_tools/']; addpath(genpath(pn.tools));

IDs={'sub-66'};%{'sub-09', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-61', 'sub-62', 'sub-64', 'sub-65', 'sub-67', 'sub-68', 'sub-69', 'sub-70', 'sub-71', 'sub-72', 'sub-73', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-99', 'sub-100', 'sub-101'}; % List all subjectsSessionID={}; % If no sessions, leave completely empty as a 0x0 cell

% including subjects sub-10 sub-60 sub-63 sub-66
%IDs={'sub-09', 'sub-10', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-60', 'sub-61', 'sub-62', 'sub-63', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-70', 'sub-71', 'sub-72', 'sub-73', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-99', 'sub-100', 'sub-101'}; % List all subjectsSessionID={}; % If no sessions, leave completely empty as a 0x0 cell

% % load info on run numbers
% load('/home/mpib/LNDG/EyeMem/scripts/PLS-scripts/conditionnum.mat');

numConds_raw = 5; % number of conditions may have changed in PLS input structure during preproc

for indID = 1:numel(IDs)
    disp(['Processing subject ', IDs{indID}, '.']);
    
%     % determine number of runs
%     numConds_raw = numel(find(~isnan(conditionnum.data(:,60))));
    
    % load subject's sessiondata file
    a = load([pn.meanFiles, '/', 'mean_', IDs{indID}, '_BfMRIsessiondata.mat']);

    % load common coordinates
    load('/home/mpib/LNDG/EyeMem/scripts/PLS-scripts/coords_EVAL.mat');
    
    a.session_info.condition{3} = 'naturals';
    conditions=a.session_info.condition(1:numConds_raw);
    a = rmfield(a,'st_datamat');
    a = rmfield(a,'st_coords');

    %replace fields with correct info.
    a.session_info.datamat_prefix   = (['SD_', IDs{indID}]); % SD PLS file name; _Bfmirsessiondata will be automatically appended!
    a.st_coords = final_coords;     % constrain analysis to shared non-zero GM voxels
    a.pls_data_path = pn.sdFiles;

    % load subject values for corresponding coordinates
    for indCond = 1:numConds_raw
        fname = [pn.sdniftis, IDs{indID} '/' conditions{indCond} '/SDNIFTI/', IDs{indID},'_sd_',conditions{indCond},'_sd.nii'];
        img = double(S_load_nii_2d(fname)); clear fname;
        img = img(final_coords,:); % restrict to final_coords
        a.st_datamat(indCond,:) = img;
        clear img;
    end
    a.st_evt_list = a.st_evt_list(1:size(a.st_datamat,1));
    a.num_subj_cond = a.num_subj_cond(1:size(a.st_datamat,1));
    
    a.session_info.num_conditions = size(a.st_datamat,1);
    a.session_info.condition = a.session_info.condition(1:size(a.st_datamat,1));
    a.session_info.condition_baseline = a.session_info.condition_baseline(1:size(a.st_datamat,1));
    a.session_info.num_conditions0 = size(a.st_datamat,1);
    a.session_info.condition0 = a.session_info.condition0(1:size(a.st_datamat,1));
    a.session_info.condition_baseline0 = a.session_info.condition_baseline0(1:size(a.st_datamat,1));
    a.session_info.run = 1:4;
    
    save([pn.sdFiles, a.session_info.datamat_prefix ,'_BfMRIsessiondata.mat'],'-struct','a','-mat');
    disp ([IDs{indID} ' done!'])
end