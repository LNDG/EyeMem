function create_varTBX_model_allstimuli_5sres_fortardis()
% This function creates an spm job specification mat file that is used as input to VarTbx.
% for testing,: 
% spm_jobman('run', 'sub-53_5TRspertrial_sdmodel.mat')
% % Excluded subjcets because missing runs: 10, 60, 63
% SubjectID_BIDS = {'sub-09','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19','sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28','sub-29','sub-30','sub-31','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-38','sub-39','sub-40','sub-41','sub-42','sub-43','sub-44','sub-45','sub-46','sub-47','sub-48','sub-49','sub-50','sub-51','sub-52','sub-53','sub-54','sub-55','sub-56','sub-57','sub-58','sub-59','sub-61','sub-62','sub-64','sub-65','sub-66','sub-67','sub-68','sub-69','sub-70','sub-71','sub-72','sub-73','sub-74','sub-75','sub-76','sub-77','sub-78','sub-79','sub-80','sub-81','sub-82','sub-83','sub-84','sub-85','sub-86','sub-87','sub-88','sub-89','sub-90','sub-91','sub-92','sub-93','sub-94','sub-95','sub-96','sub-97','sub-98','sub-99','sub-100','sub-101'};

% run on tardis!

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem';
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem';
end

PREIN_mri = fullfile(basepath, 'preproc', 'mri');

% analysis_name = 'GLM_TRwise';
% PREOUT = fullfile(basepath, analysis_name); % keep track of output folder here

analysis_name = '1TRspertrial';
PREOUT = fullfile(basepath, 'variability2', analysis_name); % keep track of output folder here
jobdir = fullfile(PREOUT, 'jobs');
cd(jobdir)
 
subjlist = [];
subjlist.YA = dir(fullfile(PREIN_mri, 'YA', 'sub-*'));
subjlist.OA = dir(fullfile(PREIN_mri, 'OA', 'sub-*'));
ageleg = {'YA' 'OA'};

for iage=1:2
  cursubjlist = subjlist.(ageleg{iage});
  nsub = length(cursubjlist);
  
  for isub = 1:nsub
    disp(cursubjlist(isub).name)
    if strcmp(cursubjlist(isub).name, 'sub-10');
      disp('TODO sort out subj10')
      continue
    end
    
    % make directory for nifti file output
    outdir = fullfile(PREOUT, cursubjlist(isub).name);
    %   mkdir (outdir)
    
    % specify matlabbatch
    matlabbatch{1}.spm.tools.variability.modeltype = 'lss';
    matlabbatch{1}.spm.tools.variability.modelmat = {fullfile(jobdir, sprintf('%s_%s_model.mat', cursubjlist(isub).name, analysis_name)) };  %{fullfile(jobdir, 'SPM.mat')}; %% directory of 1st-level model;
    matlabbatch{1}.spm.tools.variability.metric = 'sd';
    matlabbatch{1}.spm.tools.variability.resultprefix = [cursubjlist(isub).name '_sd'];
    matlabbatch{1}.spm.tools.variability.resultdir = {outdir};
    
%     % following is removed by SPM, hardcoded savebeta
%     matlabbatch{1}.spm.tools.variability.savebeta = true; 
    
    save(fullfile(jobdir, sprintf('%s_%s_sdmodel.mat', cursubjlist(isub).name, analysis_name)), 'matlabbatch');
    
  end %isub
  
end