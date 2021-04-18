% Create matrix template for each individual containing condition info
% i.e. block onsets + lengths, we choose a fixed block length here

% Hack: To get a balanced matrix of condition onsets etc. -1 is entered for
% conditions that were not present for the run. Note that PLS give a
% warning for those, but will automatically drop them from processing.

% Note: For '2131' and '2237' only the first two runs are available.

%% paths

pn.root     = '/Volumes/LNDG/Projects/StateSwitch/dynamic/data/mri/task/';
pn.template = [pn.root, 'analyses/B4_PLS_preproc2/T_tools/batch_subject_datamat_from_matrix/'];
pn.timing   = [pn.root, 'analyses/A_extractDesignTiming/B_data/A_regressors/']; % location of info regarding run timing matrices
pn.BOLDpath = [pn.root, 'analyses/B4_PLS_preproc2/B_data/BOLDin/'];
pn.PLSfiles = [pn.root, 'analyses/B4_PLS_preproc2/B_data/SD_STSWD_task_v2/']; mkdir(pn.PLSfiles);

%% fill with necessary information

% N = 44 YA + 53 OA;
IDs = {'1117';'1118';'1120';'1124';'1125';'1126';'1131';'1132';'1135';'1136';...
    '1151';'1160';'1164';'1167';'1169';'1172';'1173';'1178';'1182';'1214';'1215';...
    '1216';'1219';'1223';'1227';'1228';'1233';'1234';'1237';'1239';'1240';'1243';...
    '1245';'1247';'1250';'1252';'1257';'1261';'1265';'1266';'1268';'1270';'1276';'1281';...
    '2104';'2107';'2108';'2112';'2118';'2120';'2121';'2123';'2125';'2129';'2130';...
    '2131';'2132';'2133';'2134';'2135';'2139';'2140';'2145';'2147';'2149';'2157';...
    '2160';'2201';'2202';'2203';'2205';'2206';'2209';'2210';'2211';'2213';'2214';...
    '2215';'2216';'2217';'2219';'2222';'2224';'2226';'2227';'2236';'2237';'2238';...
    '2241';'2244';'2246';'2248';'2250';'2251';'2252';'2258';'2261'};

IDs = {'1126'};

runs = {'Run1'; 'Run2'; 'Run3'; 'Run4'};
runs_boldname = {'run-1'; 'run-2'; 'run-3'; 'run-4'};
conditions = {'dim1'; 'dim2'; 'dim3'; 'dim4'};

for indID = 1:numel(IDs)
    % load example template
    load([pn.template, 'task_BfMRI_batchfile.mat'], 'batch_file');
    batch_file.is_analysis = 0;
    batch_file.prefix = ['task_', IDs{indID}];
    batch_file.cond_name = conditions;
    batch_file.brain_region = 0;
    batch_file.ref_scan_onset = repmat(0,1,numel(conditions));
    batch_file.num_ref_scan = repmat(1,1,numel(conditions));
    if strcmp(IDs{indID}, '2131') || strcmp(IDs{indID}, '2237')
        numberOfRuns = 2;
    else
        numberOfRuns = numel(runs);
    end
    for indRun = 1:numberOfRuns
        if strcmp(IDs{indID},'1126')
            batch_file.data_files{1,indRun} = [pn.BOLDpath, IDs{indID}, '_',runs_boldname{indRun},'_feat_detrended_highpassed_denoised_nlreg_2009c_3mm.nii'];
        else
            batch_file.data_files{1,indRun} = [pn.BOLDpath, IDs{indID}, '_',runs_boldname{indRun},'_feat_detrended_highpassed_denoised_t2nlreg_2009c_3mm.nii'];
        end
        % load timing information
        load([pn.timing, IDs{indID}, '_',runs{indRun},'_regressors.mat'], 'Regressors', 'RegressorInfo');
        % onsets of block dimensionality are encoded the fourth column
        % Note: blocks were pseudo-randomized such that two blocks of each dim-condition appeared in each run. 
        batch_file.block_onsets{1,indRun} = repmat(-1,numel(conditions),2);
        batch_file.block_length{1,indRun} = repmat(-1,numel(conditions),2);
        for indCond = 1:numel(conditions)
            % keep in mind to indicate x-1 for PLS (unix style)
            numOnsets = numel(find(Regressors(:,4)==indCond));
            batch_file.block_onsets{1,indRun}(indCond,1:numOnsets) = find(Regressors(:,4)==indCond)-1;
            % block length should be constant at approx. 125 volumes excl. feedback (find(Regressors(:,4)))
            batch_file.block_length{1,indRun}(indCond,1:numOnsets) = repmat(125,1,numOnsets);
        end
    end
    % save individual configuration
    save([pn.PLSfiles, IDs{indID}, '_task_PLS_info.mat'], 'batch_file')
end
