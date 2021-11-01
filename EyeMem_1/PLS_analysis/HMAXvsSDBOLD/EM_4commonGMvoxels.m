function S1D_commonGMvoxels_summary_STSWD()
% Identify common GM coordinates to constrain analyses to.
%   nii_coords                    | non-NaN voxels across subjects
%   GM_coords                     | MNI GM mask voxels
%   final_coords                  | non-NaN MNI GM voxels
%   nii_nonZero                   | non-NaN, non-zero-power voxels across subjects
%   final_coords_withoutZero      | non-NaN, non-zero-power GM voxels across subjects

% PLS requires exlusively non-NaN voxels. GM voxels and non-Zero-power is
% optional, however may be recommended. Note that by focussing on
% non-Zero-power voxels, each subject contributes an identical amount of
% samples (N) to the analysis. Non-Zero-power voxels are usually located
% outside the brain and may inter-individually differ in number depending
% on the coregistration with MNI.

% Folder 'VoxelOverlap' has to be created manually and GM MNI mask has to
% be shifted into this directory.

% 171220 | adapted from SW, JR by JQK
% 180223 | adapted for STSWD
% 180322 | created this version, which uses the individual summary images
%           rather than loading the raw values.
% 180801 | adapted for StateSwitch task

% Note: For '2131' and '2237' only the first two runs are available. Skip.

%% paths & setup

%     pn.root    = '/Volumes/LNDG/Projects/StateSwitch/dynamic/data/mri/task/';
%     pn.overlapFolder    = [pn.root, 'analyses/B4_PLS_preproc2/B_data/VoxelOverlap/'];
%     pn.standards = [pn.root, 'analyses/B4_PLS_preproc2/B_data/A_standards/'];
%     pn.SDfiles          = [pn.root, 'analyses/B4_PLS_preproc2/B_data/X1_IndividualSD/'];
%     addpath(genpath([pn.root,'analyses/B4_PLS_preproc2/T_tools/NIFTI_toolbox/']));
%     addpath(genpath([pn.root,'analyses/B4_PLS_preproc2/T_tools/preprocessing_tools/']));

pn.root     = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem';
pn.overlapFolder    = fullfile(pn.root, 'variability/VoxelOverlap');
pn.standards = fullfile(pn.root, 'variability/A_standards');
pn.SDfiles          = fullfile(pn.root, 'variability/PLS_HMAXvsSDBOLD/SDfiles');
mkdir(pn.SDfiles)

% N = 43 YA + 53 OA;
IDs={'sub-09', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-61', 'sub-62', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-71', 'sub-72', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-100', 'sub-101'}; % 'sub-70', 'sub-73',

RunID={'alltrials'};

%% get grey matter coordinates based on MNI GM mask

GM_coords = load_nii(fullfile(pn.standards, 'tissuepriors/avg152T1_gray_MNI_3mm.nii')); % JQK: based on FSl tissue probability (thresholded @.35)
GM_coords = reshape(GM_coords.img, [],1);
GM_coords = find(GM_coords);

%% identify non-NaN voxels & non-zero power voxels

nii_coords=1:64*76*64 ; % 3mm
nii_nonZero=nii_coords;

for indID = 1:numel(IDs)
  for indRun = 1:numel(RunID)
    NIINAME=(['SD_', IDs{indID}, '_', RunID{indRun}]);
    coords=load_nii([pn.SDfiles NIINAME '.nii']);
    coords=reshape(coords.img,[],coords.hdr.dime.dim(5));
    coords(1:100,:) = NaN;
    coords_noNaN = find(~isnan(coords));
    nii_coords=intersect(nii_coords,coords_noNaN);
    % find non-zero power voxels across subjects
    coords_nonZero = find(~isnan(coords) & coords~=0);
    nii_nonZero=intersect(nii_nonZero,coords_nonZero);
    clear coords
    disp(['Done with ' RunID{indRun}])
  end
  disp (['Done with ' IDs{indID}])
end

%% indices of non-NaN GM voxels

final_coords = intersect(nii_coords, GM_coords);                            % non-NaN GM voxels across subjects
final_coords_withoutZero = intersect(nii_nonZero, GM_coords);               % non-NaN & non-zero GM voxels across subjects

%% save coordinate mat

save([pn.overlapFolder, 'coords_N',num2str(numel(IDs)),'.mat'] ,'nii_coords', 'final_coords', 'GM_coords', 'final_coords_withoutZero', 'nii_nonZero');

%% save niftis

tempNii = load_nii([pn.standards, 'tissuepriors/avg152T1_gray_MNI_3mm.nii.gz']);

nii_coords=zeros(64,76,64) ; % 3mm
nii_coords(final_coords_withoutZero) = 1;
tempNii.img = nii_coords;
save_nii(tempNii,[pn.overlapFolder, 'coords_nozero_N',num2str(numel(IDs)),'.nii'])

nii_coords=zeros(64,76,64) ; % 3mm
nii_coords(final_coords) = 1;
tempNii.img = nii_coords;
save_nii(tempNii,[pn.overlapFolder, 'coords_noNaN_N',num2str(numel(IDs)),'.nii'])

end
