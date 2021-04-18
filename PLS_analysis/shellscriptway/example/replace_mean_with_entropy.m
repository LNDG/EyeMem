function replace_mean_with_entropy

% This function loads mean_datamat files created with PLSgui (from sd_niftis 
% that were generated with VarTbx). The current content of the
% datamat files is not meaningful and they are only used for providing the
% correct structure input for later analysis with PLSgui. 
% The data used for subsitutinh mean_datamats is taken from BOLD_entropy files.

% created by MT, adapted from JQK

pn.root = '/home/mpib/LNDG/EyeMem/PLS/conditions/';
pn.meanFiles = [ pn.root, 'mean_datamats'];
pn.entrFiles = [pn.root, 'entropy_datamats/'];
pn.Bentropy = '/home/mpib/LNDG/EyeMem/data/mri/task/entropy/B_data/';
pn.tools = ['/home/mpib/LNDG/toolboxes/preprocessing_tools/']; addpath(genpath(pn.tools));

IDs={'sub-09','sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-61', 'sub-62', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-70', 'sub-71', 'sub-72', 'sub-73', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-99', 'sub-100', 'sub-101'}; % List all subjectsSessionID={}; % If no sessions, leave completely empty as a 0x0 cell
subnum = [9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101];

% load info on run numbers
load('/home/mpib/LNDG/EyeMem/scripts/PLS-scripts/conditions/conditionnum.mat');

% load common coordinates
load('/home/mpib/LNDG/EyeMem/scripts/PLS-scripts/conditions/coords_EVAL.mat');

numConds_raw = 5; % number of conditions may have changed in PLS input structure during preproc

for indID = 1:numel(IDs)
    disp(['Processing subject ', IDs{indID}, '.']);
    
    % get run order for condition (BOLD entropy)
    run_nums = [conditionnum.data(:,subnum(indID))]';
      
    % load subject's sessiondata file
    a = load([pn.meanFiles, '/', 'mean_', IDs{indID}, '_BfMRIsessiondata.mat']);
   
    a.session_info.condition{3} = 'naturals';
    conditions=a.session_info.condition(1:numConds_raw);
    a = rmfield(a,'st_datamat');
    a = rmfield(a,'st_coords');

    %replace fields with correct info.
    a.session_info.datamat_prefix   = (['entr_', IDs{indID}]); % SD PLS file name; _Bfmirsessiondata will be automatically appended!
    a.st_coords = final_coords;     % constrain analysis to shared non-zero GM voxels
    a.pls_data_path = pn.entrFiles;

    % load subject values for corresponding coordinates
    for indCond = run_nums
        
        % condition order (for SD BOLD)     
        cond_num = find(run_nums == indCond);
        
        % load BOLD_entropy run according to condition (starting with fractals, aso.)
        fname = [pn.Bentropy, IDs{indID} '/' IDs{indID} '_run-0' num2str(indCond) '_BOLD_entropy_Sload.mat']; 
        load(fname); clear fname;
        img = BOLD_entropy.sampen(:,4);
        img = img(final_coords,:); % restrict to final_coords
        a.st_datamat(cond_num,:) = img;
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
    a.session_info.run = 1:5;
    
    save([pn.entrFiles, a.session_info.datamat_prefix ,'_BfMRIsessiondata.mat'],'-struct','a','-mat');
    disp ([IDs{indID} ' done!'])
end