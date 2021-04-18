function make_meanbold_datamat
%% 2nd pre-step for PLS.
% Input: 1)SubjectID
% Create meanbold PLS files from text batch files.

% For compilation, use:
% /opt/matlab/R2014b/bin/mcc -m make_meanbold_datamat -a /home/mpib/LNDG/toolboxes/PLS/pls

%% Subject List 
SubjectID={'sub-09', 'subj-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-61', 'sub-62', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-70', 'sub-71', 'sub-72', 'sub-73', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-100', 'sub-101'}; 
% List all subjectsSessionID={}; % If no sessions, leave completely empty as a 0x0 cell
 % If running serially, list all subjects

%% Paths if not running compiled
% addpath(genpath('/home/mpib/LNDG/toolboxes/PLS/pls'));
% BASEPATH=('/home/mpib/LNDG/EyeMem/');
% MEANPATH=([BASEPATH 'PLS/stim_vs_ITI/mean_datamats/']);
% 
addpath(genpath('/Users/kloosterman/gridmaster2012/LNDG/toolboxes/PLS/pls')); 
BASEPATH=('/Users/kloosterman/gridmaster2012/LNDG/EyeMem/');
MEANPATH=([BASEPATH 'data/PLS/HMAXvsSDBOLD/mean_datamats/']); 

%% Create mean-BOLD PLS files
for i=1:length(SubjectID)
    TextBatch=([SubjectID{i} '_meanbold_batch.txt']); % Text batch filename
    cd (MEANPATH) % The toolbox outputs our mean-BOLD PLS files in the current directory
    fprintf('Creating mean-BOLD PLS datamat for subject %s \n', SubjectID{i});
    batch_plsgui([ MEANPATH '/' TextBatch ]);
    fprintf('Finished creating PLS meanbold_datamat \n');
    
end
