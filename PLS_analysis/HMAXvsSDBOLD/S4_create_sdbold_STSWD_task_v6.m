function S4_create_sdbold_STSWD_task_v6()
%% 2nd pre-step for PLS.
% fill matrix with standard deviations for various conditions, runs and blocks
% IDs{indID} - subject ID (string)

% NOTE: This script is expected to be run off-tardis for subjects for which
% the analysis had to be repeated.

% 171220 | adapted by JQK; use only non-zero-power GM voxels
% 180202 | adapted for repeats
% 180302 | non-zero voxels lost occipital cortex, use non-NaN voxels
% 181212 | fixed run encoding

% v4:   revert to use the non-zero version
% v6:   build on v4
%       exclude volumes with significant DVARS contamination
%       include calculation of spatial PCA scores

% If data is missing (e.g. a run), command window issues a warning, but
% continues with the next available data.
    
pn.root     = '/Volumes/LNDG/Projects/StateSwitch/dynamic/data/mri/task/';
pn.tools	= [pn.root, 'analyses/B4_PLS_preproc2/T_tools/'];  addpath(genpath(pn.tools));
pn.data     = [pn.root, 'analyses/B4_PLS_preproc2/B_data/'];
pn.PLSpath  = [pn.data,'SD_STSWD_task_v2/'];
pn.maskPath = [pn.data, 'VoxelOverlap/'];
pn.pcaPath  = [pn.data,'PCAdim/'];

% N = 44 YA + 51 OA; 2131, 2237 dropped due to missing runs
IDs = {'1117';'1118';'1120';'1124';'1125';'1126';'1131';'1132';'1135';'1136';...
    '1151';'1160';'1164';'1167';'1169';'1172';'1173';'1178';'1182';'1214';'1215';...
    '1216';'1219';'1223';'1227';'1228';'1233';'1234';'1237';'1239';'1240';'1243';...
    '1245';'1247';'1250';'1252';'1257';'1261';'1265';'1266';'1268';'1270';'1276';'1281';...
    '2104';'2107';'2108';'2112';'2118';'2120';'2121';'2123';'2125';'2129';'2130';...
    '2132';'2133';'2134';'2135';'2139';'2140';'2145';'2147';'2149';'2157';...
    '2160';'2201';'2202';'2203';'2205';'2206';'2209';'2210';'2211';'2213';'2214';...
    '2215';'2216';'2217';'2219';'2222';'2224';'2226';'2227';'2236';'2238';...
    '2241';'2244';'2246';'2248';'2250';'2251';'2252';'2258';'2261'};

numConds_raw = 4; % number of conditions may have changed in PLS input structure during preproc

for indID = 1:numel(IDs)
    disp(['Processing subject ', IDs{indID}, '.']);

%% create the SDBOLD datamats

    % load subject's sessiondata file
    a = load([pn.PLSpath, 'task_', IDs{indID}, '_BfMRIsessiondata.mat']);

    % load common coordinates
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