function EM_4makeSDbold_datamat()
%% 2nd pre-step for PLS.
% fill matrix with standard deviations for various conditions, runs and blocks
% IDs{indID} - subject ID (string)

% NOTE: This script is expected to be run off-tardis for subjects for which
%% paths
pn.root     = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem';
pn.trialinfo   = fullfile(pn.root, 'preproc/eye'); % eyedata.mat: location of info regarding run trialinfo matrices
pn.BOLDpath = fullfile(pn.root, 'variability/5TRspertrial');
pn.PLSfiles = fullfile(pn.root, 'variability/PLS_HMAXvsSDBOLD'); mkdir(pn.PLSfiles);
% pn.root     = '/Volumes/LNDG/Projects/StateSwitch/dynamic/data/mri/task/';
% pn.tools	= [pn.root, 'analyses/B4_PLS_preproc2/T_tools/'];  addpath(genpath(pn.tools));
% pn.data     = [pn.root, 'analyses/B4_PLS_preproc2/B_data/'];
% pn.PLSpath  = [pn.data,'SD_STSWD_task_v2/'];
% pn.maskPath = [pn.data, 'VoxelOverlap/'];
% pn.pcaPath  = [pn.data,'PCAdim/'];

%% fill with necessary information
% N = 44 YA + 53 OA;
IDs={'sub-09', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29', 'sub-30', 'sub-31', 'sub-32', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40', 'sub-41', 'sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-46', 'sub-47', 'sub-48', 'sub-49', 'sub-50', 'sub-51', 'sub-52', 'sub-53', 'sub-54', 'sub-55', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-61', 'sub-62', 'sub-64', 'sub-65', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-71', 'sub-72', 'sub-74', 'sub-75', 'sub-76', 'sub-77', 'sub-78', 'sub-79', 'sub-80', 'sub-81', 'sub-82', 'sub-83', 'sub-84', 'sub-85', 'sub-86', 'sub-87', 'sub-88', 'sub-89', 'sub-90', 'sub-91', 'sub-92', 'sub-93', 'sub-94', 'sub-95', 'sub-96', 'sub-97', 'sub-98', 'sub-100', 'sub-101'}; % 'sub-70', 'sub-73', 

numConds_raw = 1; % number of conditions may have changed in PLS input structure during preproc

for indID = 1:numel(IDs)
    disp(['Processing subject ', IDs{indID}, '.']);

%% create the SDBOLD datamats

    % load subject's sessiondata file
    a = load(fullfile(pn.PLSfiles, sprintf('task_%s_BfMRIsessiondata.mat', IDs{indID})));

    % TODO load common coordinates
    load([pn.maskPath, 'coords_N95.mat'], 'final_coords_withoutZero');
    final_coords = final_coords_withoutZero;
    
    conditions=a.session_info.condition(1:numConds_raw);
    a = rmfield(a,'st_datamat');
    a = rmfield(a,'st_coords');

    %replace fields with correct info.
    a.session_info.datamat_prefix   = (['SD_', IDs{indID}, '_v6']); % SD PLS file name; _Bfmirsessiondata will be automatically appended!
    a.st_coords = final_coords;     % constrain analysis to shared non-zero GM voxels
    a.pls_data_path = pn.PLSpath;

    % initialize this subject's datamat
    a.st_datamat = NaN(numel(conditions),numel(final_coords)); %(indCond voxel)

    % intialize indCond specific scan count for populating cond_data
    clear count cond_data block_scan;
    for indCond = 1:numel(conditions)
        count{indCond} = 0;
    end

    %% Header
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % within each indBlock express each scan as deviation from indBlock's
    % temporal mean.Concatenate all these deviation values into one
    % long condition specific set of scans that were normalized for
    % indBlock-to-indBlock differences in signal strength. In the end calculate
    % stdev across all normalized indCond scans
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % for each condition identify its scans within  runs
    % and prepare where to put indCond specific normalized data
    for indCond = 1:numel(conditions)
        tot_num_scans = 0;
        for indRun = 1:a.session_info.num_runs
            if a.session_info.run(indRun).num_scans==0
                disp(['No data for ' IDs{indID} ' in ', num2str(indRun)]);
                continue;
            end
            onsets = a.session_info.run(indRun).blk_onsets{indCond}+1; % +1 is because we need matlab indexing convention (damain)
            lengths = a.session_info.run(indRun).blk_length{indCond};
            for indBlock = 1:numel(onsets)
                block_scans{indCond}{indRun}{indBlock} = onsets(indBlock)-1+[1:lengths(indBlock)];
                this_length = lengths(indBlock);
                if max(block_scans{indCond}{indRun}{indBlock}>a.session_info.run(indRun).num_scans)
                    disp(['Problem: ' IDs{indID} ' something wrong in ', num2str(indRun)]);
                    block_scans{indCond}{indRun}{indBlock} = intersect(block_scans{indCond}{indRun}{indBlock},[1:a.session_info.run(indRun).num_scans]);
                    this_length = numel(block_scans{indCond}{indRun}{indBlock});
                end
                tot_num_scans = tot_num_scans + this_length;
            end
        end
        % create empty matrix with dimensions coords (rows) by total # of scans (columns).
        cond_data{indCond} = NaN(numel(final_coords),tot_num_scans);
    end

    %% Load NIfTI file

    for indRun = 1:a.session_info.num_runs

        %% Load NIfTI
        % load nifti file for this run; %(x by y by z by time)
        % check for nii or nii.gz, unzip .gz and reshape to 2D

        fname = [a.session_info.run(indRun).data_path '/' a.session_info.run(indRun).data_files{:}];
        
        if ~exist(fname) || isempty(a.session_info.run(indRun).data_files)
            warning(['File not available: ', IDs{indID}, ' Run ', num2str(indRun)]);
            continue;
        end
        
        % If using preprocessing tools
        [img] = double(S_load_nii_2d( fname ));
        img = img(final_coords,:); %restrict to final_coords
        
        %% replace DVARS volumes with NaN across voxels
        % these volumes will naturally be dropped --> hardcore censoring
        
        load([pn.maskPath, 'S2_DVARS.mat'], 'DVARSout');
        DVARS_idx = find(strcmp(DVARSout.IDs, IDs{indID}));
        Volumes2Censor = DVARSout.SignificantDVARS{DVARS_idx,indRun,1};
        
        img(:,Volumes2Censor) = NaN;

        %% Now, proceed with creating SD datamat...

        for indCond = 1:numel(conditions)
            for indBlock = 1:numel(block_scans{indCond}{indRun})
                block_data = img(:,block_scans{indCond}{indRun}{indBlock});% (vox time)
                % normalize block_data to global indBlock mean = 100.
                block_data = 100*block_data/nanmean(nanmean(block_data));
                % temporal mean of this indBlock
                block_mean = nanmean(block_data,2); % (vox) - this should be 100
                % express scans in this indBlock as deviations from block_mean
                % and append to cond_data
                good_vox = find(block_mean);
                for t = 1:size(block_data,2)
                    count{indCond} = count{indCond} + 1;
                    cond_data{indCond}(good_vox,count{indCond}) = block_data(good_vox,t) - block_mean(good_vox);%must decide about perc change option here!!??
                end
            end
        end
        clear img;
    end

    %% Calculate standard deviation for each condition across runs
    for indCond = 1:numel(conditions)
        a.st_datamat(indCond,:) = squeeze(nanstd(cond_data{indCond},[],2))';
    end
    % Save SD datamat
    save([pn.PLSpath, a.session_info.datamat_prefix '.mat'],'-struct','a','-mat');
            
%     %% Calculate PCA dimensionality for each condition across runs
%     
%     for indCond = 1:numel(conditions)
% 
%         tmpDataIn = cond_data{indCond};
%         tmpDataIn(:,isnan(squeeze(mean(tmpDataIn)))) = [];
%         
%         %% Dimensions
%         tic;[coeff.RAW{indCond}, scores{indCond}, ~, ~, EXPLAINED{indCond}] = pca(tmpDataIn', 'VariableWeights','variance', 'Centered', true); toc;   % spatial PCA using correlation matrix
%         % Rows of X correspond to observations and columns correspond to variables
% 
%         %% Extracting components explaining min. 90% variance
%         criterion = 90;
% 
%         cumulExplained = cumsum(EXPLAINED{indCond});
%         Dimensions(1, indCond) = find(round(cumulExplained)>=criterion,1,'first');
% 
% %         %% standardize coeff scores to gain a comparable correlation matrix
% %         % to standardize, the scores (principal component scores) will be 
% %         % correlated with coeff.weights columnwise. I.e. each voxels TS is 
% %         % correlated with coeff.weights
% %         % for spatial PCA matrix has to be transposed 884*171922 -> 171922*884
% % 
% %         try
% %             coeff.STAND{indCond} = corr(scores{indCond}, tmpDataIn');
% %         catch ME
% %             disp (ME.message)
% %         end
% % 
% %         %% create coefficient/eigenvector matrix of the rotated solution
% %         % with the exact number of dimenions for the subject previously calculated
% %         try
% %             coeff.STAND_ROT{indCond}=rotatefactors(coeff.STAND{indCond}(:, 1:Dimensions(1, indCond)));%create rotated versions of standardized coeffs above. Default is varimax... 
% %         catch ME
% %             disp (ME.message)
% %         end
% 
%     end % condition
%     % save PCA dimensionality
%     save([pn.pcaPath, IDs{indID}, 'PCAcorr_dim_spatial.mat'], 'Dimensions', 'EXPLAINED'); %note that EXPLAINED here is from typical, unrotated solution. 

    %% clean up

    disp ([IDs{indID} ' done!'])
    
end