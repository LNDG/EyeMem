function create_1stlevel_singletrial()
% This function creates a 1st level model specification mat file per
% subject OVER ALL CONDITIONS, which is called in the mat file specifying the job to be used as input to VarTbx.

% excluded subjects: EYEMEM001 - EYEMEM008
% TODO additionally excluded becase not all runs present: 10, 60, 63
% SubjectID = { 'EYEMEM009', 'EYEMEM011','EYEMEM012', 'EYEMEM013', 'EYEMEM014', 'EYEMEM015', 'EYEMEM016', 'EYEMEM017', 'EYEMEM018', 'EYEMEM019', 'EYEMEM020', 'EYEMEM021', 'EYEMEM022', 'EYEMEM023', 'EYEMEM024', 'EYEMEM025', 'EYEMEM026', 'EYEMEM027', 'EYEMEM028', 'EYEMEM029', 'EYEMEM030', 'EYEMEM031', 'EYEMEM032', 'EYEMEM033', 'EYEMEM034', 'EYEMEM035', 'EYEMEM036', 'EYEMEM037', 'EYEMEM038', 'EYEMEM039', 'EYEMEM040', 'EYEMEM041', 'EYEMEM042', 'EYEMEM043', 'EYEMEM044', 'EYEMEM045', 'EYEMEM046', 'EYEMEM047', 'EYEMEM048', 'EYEMEM049', 'EYEMEM050', 'EYEMEM051', 'EYEMEM052', 'EYEMEM053', 'EYEMEM054', 'EYEMEM055', 'EYEMEM056', 'EYEMEM057', 'EYEMEM058', 'EYEMEM059', 'EYEMEM061', 'EYEMEM062', 'EYEMEM064', 'EYEMEM065', 'EYEMEM066', 'EYEMEM067', 'EYEMEM068', 'EYEMEM069', 'EYEMEM070', 'EYEMEM071', 'EYEMEM072', 'EYEMEM073', 'EYEMEM074', 'EYEMEM075', 'EYEMEM076', 'EYEMEM077', 'EYEMEM078', 'EYEMEM079', 'EYEMEM080', 'EYEMEM081', 'EYEMEM082', 'EYEMEM083', 'EYEMEM084', 'EYEMEM085', 'EYEMEM086', 'EYEMEM087', 'EYEMEM088', 'EYEMEM089', 'EYEMEM090', 'EYEMEM091', 'EYEMEM092', 'EYEMEM093', 'EYEMEM094', 'EYEMEM095', 'EYEMEM096', 'EYEMEM097', 'EYEMEM098', 'EYEMEM099', 'EYEMEM100', 'EYEMEM101'};
% SubjectID_BIDS = {'sub-09','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19','sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28','sub-29','sub-30','sub-31','sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-38','sub-39','sub-40','sub-41','sub-42','sub-43','sub-44','sub-45','sub-46','sub-47','sub-48','sub-49','sub-50','sub-51','sub-52','sub-53','sub-54','sub-55','sub-56','sub-57','sub-58','sub-59','sub-61','sub-62','sub-64','sub-65','sub-66','sub-67','sub-68','sub-69','sub-70','sub-71','sub-72','sub-73','sub-74','sub-75','sub-76','sub-77','sub-78','sub-79','sub-80','sub-81','sub-82','sub-83','sub-84','sub-85','sub-86','sub-87','sub-88','sub-89','sub-90','sub-91','sub-92','sub-93','sub-94','sub-95','sub-96','sub-97','sub-98','sub-99','sub-100','sub-101'};

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem';
else
  basepath = '/home/mpib/kloosterman/projectdata/eyemem';
end
PREIN_mri = fullfile(basepath, 'preproc', 'mri');
PREIN_eye = fullfile(basepath, 'preproc', 'eye');

analysis_name = '1TRspertrial';  %1TRspertrial  5TRspertrial
PREOUT = fullfile(basepath, 'variability2', analysis_name); % keep track of output folder here
% analysis_name = 'GLM_TRwise';
% PREOUT = fullfile(basepath, analysis_name); % keep track of output folder here
jobdir = fullfile(PREOUT, 'jobs');
mkdir(jobdir)

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
    %-----------------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(PREOUT, ageleg{iage}, cursubjlist(isub).name)};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1; % TR
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16; % microtime resolution
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % microtime onset
    nvols = 2310;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {};
    for ivol = 1:nvols
      matlabbatch{1}.spm.stats.fmri_spec.sess.scans{end+1,1} = fullfile(PREIN_mri, ageleg{iage}, cursubjlist(isub).name, ...
        sprintf('%s_task-eyemem_allruns_bold_feat_detrended_highpassed_denoised_MNI_demeaned.nii,%d', cursubjlist(isub).name, ivol));
    end
    
    load(fullfile(PREIN_eye, ageleg{iage}, ['eye_' cursubjlist(isub).name '.mat'])); % load eye struct data in ft format
    
    % update triggers so they match the concatenated data
    nvolsperrun = 462;  % 2310/462 = 5;
    onsetshift = (data.trialinfo(:,12)-1) .* nvolsperrun;  % runno-1 * 462 = shifts due to concatenation
    piconsets = data.trialinfo(:,9) + onsetshift;
    
%     % fsl txt files test
%     ev = table([piconsets zeros(size(piconsets))+5 ones(size(piconsets)) ])
%     writetable(ev, '~/Desktop/evstim.txt', 'Delimiter', ' ')
    
    switch analysis_name
      case '5TRspertrial'
        % make onsets for each TR in trial
        onsets = repmat(piconsets, 1,5);
        onsets = onsets + [0 1 2 3 4]; % dimord: trials TR
        duration = ones(length(onsets),1);
      case '1TRspertrial'
        onsets = piconsets; % onsets as is
        duration = zeros(length(onsets),1)+5; % 5 s duration
    end
        
    %% Specify conditions single trial
    removeonsets_near_run_end = true;  % remove onsets within 20 s from run end, due to modelling issues

    onsets = sort(onsets(:)); % order in time and put together
    if removeonsets_near_run_end
      onsets = onsets(onsets < nvols-20);
      duration = duration(1:length(onsets));
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'alltrials';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = onsets;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = duration;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    
    switch analysis_name
      case '1TRspertrial'
        matlabbatch{1}.spm.stats.fmri_est.spmmat = '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % set to 1 to keep residuals
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    end
%     % add patch task as "condition"
%     patchonsets = piconsets + 6.5; % patch comes 6.5 s after pic
%     matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'nuisance_patch';
%     matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = patchonsets;
%     matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = ones(length(patchonsets),1) + 2; % leave on for X s
%     matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
%     matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
%     matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;

%     % this treats each trial as a single condition
%     for itrial = 1:size(onsets, 1)
%       matlabbatch{1}.spm.stats.fmri_spec.sess.cond(itrial).name = sprintf('trial%d', itrial) ;
%       
%       curonsets = onsets(itrial,:); % flatten matrix and sort to be sure
%       if removeonsets_near_run_end
%         curonsets = curonsets(curonsets < nvols-20);
%       end
%       matlabbatch{1}.spm.stats.fmri_spec.sess.cond(itrial).onset = curonsets(:);
%       matlabbatch{1}.spm.stats.fmri_spec.sess.cond(itrial).duration = ones(length(curonsets),1); % [5; 5; 5];
%       matlabbatch{1}.spm.stats.fmri_spec.sess.cond(itrial).tmod = 0;
%       matlabbatch{1}.spm.stats.fmri_spec.sess.cond(itrial).pmod = struct('name', {}, 'param', {}, 'poly', {});
%       matlabbatch{1}.spm.stats.fmri_spec.sess.cond(itrial).orth = 1;
%     end
    
    % is this all needed?
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % spm_jobman('run', matlabbatch)
    save( fullfile(jobdir, sprintf('%s_%s_model.mat', cursubjlist(isub).name, analysis_name)) , 'matlabbatch');
  end % isub
end % iage


