function make_commoncoords(MNI)
%% 3rd pre-step for PLS analysis
% Input: 1)SubjectID
% Create GM masked common coordiantes for all subject/conditions/runs

% Best to run as an interactive job rather than as a compilation. % For compilation, use:
% /opt/matlab/R2014b/bin/mcc -m make_commoncoord -a /home/mpib/LNDG/toolboxes/preprocessing_tools

%% Settings

% MNI Size
try
    MNI = str2num(MNI); % input should be double, not string
catch
end

% ID Lists
SubjectID={'sub-09', 'sub-10', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-60', 'sub-61', 'sub-62', 'sub-63', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-70', 'sub-71', 'sub-72', 'sub-73', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-99', 'sub-100', 'sub-101'}; % List all subjectsSessionID={}; % If no sessions, leave completely empty as a 0x0 cell
RunID={'fractals', 'landscapes', 'naturals', 'streets1', 'streets2'}; % List all run/condition images

% Preprocessing input infortmation
PreprocPipe={'variability'}; % Location of preprocessing pipeline used for analysis
PreprocSuffix={'sd'}; % Suffix of preprocessing steps performed on data, ex: '_FEAT_detrend_highpass_denoised_2mm_MNI'

%% Initiate paths
ROOTPATH=('/home/mpib/LNDG/');
BASEPATH=([ROOTPATH 'EyeMem/']);
DATAPATH=([BASEPATH 'data/']);
SAVEPATH=([BASEPATH, 'scripts/PLS-scripts/coords_EVAL.mat']);
GRAYPATH=([ROOTPATH '/Standards/binary_masks/avg152_T1_gray_mask_90_' num2str(MNI) 'mm_binary.nii']); % Binary Gray matter mask

% Comment out if running compiled
addpath(genpath([ROOTPATH '/toolboxes/NIfTI_20140122']));
addpath(genpath([ROOTPATH '/toolboxes/preprocessing_tools']));

%% Initiate coordinate size for NIfTI files
if MNI >= 1
    nii_coords= 1:floor(182/MNI)*floor(218/MNI)*floor(182/MNI) ; % 1MM = 182*218*182
elseif MNI == 0.5
    nii_coords= 1:364*436*364;
end

%% Find common coordiantes for all subjects/runs
for i = 1:numel(SubjectID)
%     if size(SessionID)==0
%         SessionID=1;
%     end  
%     for s=1:length(SessionID)
%         if iscellstr(SessionID)>=1
%             SessionName=[SessionID{s}, '_'];
%             SessionFolder=[SessionID{s}, '/'];
%         else
%             SessionName=[];
%             SessionFolder=[];
%         end 
        for j = 1:numel(RunID)
            %% Set image name/path
            NIIPATH=([DATAPATH 'mri/task/' PreprocPipe{1} '/B_data/' SubjectID{i} '/' RunID{j} '/SDNIFTI/']);
            NIINAME=([SubjectID{i} '_sd_' RunID{j} '_' PreprocSuffix{1}]);
            
            %% Find common coords
            fname=([ NIIPATH NIINAME '.nii']);
            if exist(fname)
            [ coords ] = double(S_load_nii_2d( fname ));
            %coords=load_nii([ NIIPATH NIINAME '.nii']);
            %coords=reshape(coords.img,[],coords.hdr.dime.dim(5));
            coords(isnan(coords))=0; % Necessary to avoid NaNs, which cause problems for PLS analysis
            coords = find(coords(:,1));
            
            % Add to common coordinates
            nii_coords=intersect(nii_coords,coords);
            end
            
            clear coords
        end
%     end
    disp (['Finished with subject ' SubjectID{i}])
end

%% Gray Matter Mask
% Find coordiantes of gray matter mask
GM_coords = load_nii(GRAYPATH);
GM_coords = reshape(GM_coords.img, [],1);
GM_coords = find(GM_coords);

% Perform gray matter masking of common coordinates
final_coords = intersect(nii_coords, GM_coords);

%% Export coordinates
save ( SAVEPATH ,'nii_coords', 'final_coords', 'GM_coords');

end
