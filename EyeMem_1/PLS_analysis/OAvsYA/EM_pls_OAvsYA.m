function EM_pls_OAvsYA(cfg)
% run task PLS on 5-trial binned SDbold vs Hmax. 30 conditions
% TODO turn into data format that PLS toolbox likes:
% examples:
% sesdat = load('/Users/kloosterman/gridmaster2012/LNDG/EyeMem_old/PLS/SDdatamats/EYEMEM009_BfMRIsessiondata.mat') % behav PLS?
% sesdat2 = load('/Volumes/FB-LIP/Projects/COBRA/data/mri+PETSharp/B_analyses/Paper1_2018/PLS_nback/SD_2mm/SD_C219_BfMRIsessiondata.mat') % task PLS?

sourcefile = cfg.sourcefile;
outfile_source = cfg.outfile_source;
outfile_sesdat = cfg.outfile_sesdat;
removeoutliers = cfg.removeoutliers;
Z_thresh = cfg.Z_thresh;
nbins = cfg.nbins;
do_kstest = cfg.do_kstest;
behavfile = cfg.behavfile ;
PLStype = cfg.PLStype;
subj = cfg.subj;
BOLDvar_measure = cfg.BOLDvar_measure;
binsubtract = cfg.binsubtract;

disp(sourcefile)
source = load(sourcefile); % source comes out
if ismac;  try
    cfg=[];
    cfg.method = 'ortho'; % slice
    cfg.funparameter = 'pow';
    ft_sourceplot(cfg, source)
  end
end


switch BOLDvar_measure
  case 'nanstd'
%     source.pow = nanstd(seldat,1,2); % take SD across 5 trials, 5 TR's each
    source.pow = nanstd(source.pow(:,:),0,2); % take SD across all trials and times
    source.powdimord = 'pos';
  case 'iqr'
    source.pow = iqr(source.pow(:,:),2); % take IQR across all trials and times
%     source.pow(:,ibin) = iqr(seldat,2); % take IQR across 5 trials, 5 TR's each
end

%% make PLS sesssiondata structure, prepare important fields
% 1. standard stuff: TODO turn into function?
tmp=[];
% tmp.st_coords = find(source.inside); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
load(fullfile(cfg.PREIN, 'common_coords.mat'));
tmp.st_coords = find(common_coords); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
switch PLStype
  case 'taskPLS'
    tmp.st_datamat = transpose(source.pow(tmp.st_coords,:)); %these are defined above as a mask of voxels %to use. All others are excluded (e.g., in a GM mask). % cond x voxels
    tmp.session_info.datamat_prefix = 'SDbold_vs_HMAX'; %[subj '_' pattern];%stores common      %datamat prefix
    tmp.behavdata = [];
    tmp.behavname = {};%no behav names or data just yet, but set up the %field anyway...
end
tmp.st_evt_list = 1:size(tmp.st_datamat,1);%how many conditions?
% tmp.st_sessionFile = [pls_dir subj '_' pattern '_BfMRIsession.mat'];
tmp.unequal_subj = 0;

% session_info struct, for bookkeeping I believe:
%       a = load([workdir '/charisaPLS/' subj '_BfMRIsession.mat']);%this      %loads a particular session file.
tmp.session_info.pls_data_path = fileparts(outfile_source); %this will add info to 'a' for   tracking purposes, and for commands below.
tmp.session_info.file_pattern = '*nii'; % NEEDED?
[pathstr, name]= fileparts(sourcefile);
tmp.session_info.run(1).data_path = pathstr; %[workdir subj];
tmp.session_info.run(1).data_files = {[name '.nii']}; % CHECK should be nifti? e.g. /Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/5TRspertrial/sub-77/beta_series
tmp.session_info.run(1).file_pattern = tmp.session_info.run(1).data_files; %['smooth_CommonAtlas_mc_tshift_cond' num2str(run) '.nii*'];

tmp.session_info.num_conditions = length(tmp.st_evt_list);
tmp.session_info.condition = {};
for i = 1:tmp.session_info.num_conditions
  tmp.session_info.condition{end+1} = sprintf('HMAX%d', i); % OA vs YA??
end

tmp.session_info.num_runs = 1;
% bookkeeping of events:
[pathstr, name]= fileparts(source.hdr.fspec);
tmp.session_info.run.num_scans = source.hdr.nframes;
tmp.session_info.run.data_path = pathstr;
tmp.session_info.run.data_files = { [name '.nii'] };
tmp.session_info.run.file_pattern = [name '.nii'];
tmp.session_info.run.blk_onsets = []; %num2cell(cond_bins, 2)'; % =  {[9�1 double]  [9�1 double]  [9�1 double]  [2�1 double]} vectors with TR onsets
tmp.session_info.run.blk_length = []; %num2cell(ones(size(cond_bins)), 2)'; % =  {[9�1 double]  [9�1 double]  [9�1 double]  [2�1 double]} vectors with TR onsets% =  {[9�1 double]  [9�1 double]  [9�1 double]  [2�1 double]} vectors with lengths pe onset

% all hardcoded stuff:
tmp.create_ver = '5.0807231';%not sure where value comes from; check.
tmp.st_win_size = 1;%PLS command...must check.

[i,j,k] = ind2sub(source.dim, find(source.pos(:,1)==0 & source.pos(:,2)==0 & source.pos(:,3)==0));
tmp.st_origin = [i j k]; % sets the origin...PLS default?

tmp.st_voxel_size = [source.hdr.xsize source.hdr.ysize source.hdr.zsize]; %vox size defined
tmp.st_dims = [ source.dim(1:2) 1 source.dim(3)]; % [40 48 1 34];%dimensions of matrix defined for PLS WHY 1 for dim3??
tmp.SingleSubject = 0;%not a single subj design
tmp.normalize_volume_mean = 0;
tmp.create_datamat_info.brain_mask_file = ''; % NEEDED??
% tmp.create_datamat_info.run_idx = run_idx; % NOT NEEDED it seems, bijv [1:7];%seven runs here
tmp.create_datamat_info.brain_coord_thresh = 0;%no ROIs here...check
tmp.create_datamat_info.consider_all_voxels_as_brain = 0;%we don't...
tmp.create_datamat_info.num_skipped_scans = 0;%no skipped scans
tmp.create_datamat_info.ignore_slices = [];%no slices to ignore, so %keep all
tmp.create_datamat_info.normalize_volume_mean = 0;%no normalization %here...
% tmp.create_datamat_info.normalize_with_baseline = 1;%yes, normalize %to baseline (but changes below, so must override?)
tmp.create_datamat_info.normalize_with_baseline = 0;%yes, normalize %to baseline (but changes below, so must override?)
% tmp.create_datamat_info.merge_across_runs = 1;%says to merge across %runs, but we don't do this below anyway...is overridden...??
tmp.create_datamat_info.merge_across_runs = 0;%says to merge across %runs, but we don't do this below anyway...is overridden...??
tmp.create_datamat_info.single_subject_analysis = 0;%not single subj

if ismac
  cfg=[];
  cfg.method = 'ortho'; % slice
  cfg.funparameter = 'pow';
  ft_sourceplot(cfg, source)
end
disp(outfile_source)
save(outfile_source, 'source')

% save PLS sesdat
disp(outfile_sesdat)
save(outfile_sesdat, '-struct','tmp')
% % save([subj '_' pattern '_BfMRIdatamat.mat'], '-struct','tmp'); %,'-mat'


%   if ismac
%     tmp = source_bin;
%     tmp.powdimord = 'pos';
%     %       vol=300;
%     %   tmp.anatomy = tmp.anatomy(:,:,:,vol);
%     tmp.pow = nanmean(tmp.pow(:,:),2);
%     cfg=[];
%     cfg.method = 'ortho'; % slice ortho glassbrain vertex
%     cfg.funparameter = 'pow';
%     cfg.funcolorlim = 'zeromax';% [-300 300];
%     ft_sourceplot(cfg, tmp)
%   end

% %   source.time = 1:5;
% %   source.freq = 1;
% allsource(isub) = source_bin;
% % end
