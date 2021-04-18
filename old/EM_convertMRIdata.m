function EM_convertdata( )
% Convert raw data and make file structure

overwrite = 0;

inpath = '/Volumes/LNDG/Eye_mem/data';
% outpath = '/Users/kloosterman/gridmaster2012/projectdata/EyeMem/data';
outpath = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data';
dcm2nii = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/toolbox/dcm2nii';

mkdir(outpath)

%% MRI data
PREIN_mri = fullfile(inpath, 'mri');

cd(PREIN_mri)
subj_list = dir(fullfile(PREIN_mri, 'EYEMEM*'));

for isub = 1:length(subj_list)
    %prepare output destinations
    subj_id = subj_list(isub).name(1:9);
    outpath_subj = fullfile(outpath, subj_id);
    mkdir(outpath_subj)    
    % go into dest folder, make dir structure
    cd(outpath_subj)

    %% convert eye edf to asc
    mkdir('eye')
    cd(fullfile(inpath, 'eye'))
    
    subj_no = str2num(subj_list(isub).name(7:9));
    runlist = dir(sprintf('S%dp*', subj_no));
    for irun = 1:length(runlist)
        
        [asc, data, event] = EM_eyelink_to_fieldtrip(runlist(irun).name);
        
        % TODO save data 
    end
    
    
    mkdir('behavior')
    mkdir('physio')

       
    %% Convert MRI data
    mkdir('mri')
    % go to raw data 
    cd(fullfile(PREIN_mri, subj_list(isub).name))
    datafolder = dir('LIP*');
    cd(datafolder.name)    
    
    %convert the EPI's
%     runlist = dir('EP2D*');
    runlist = dir('T1_MPRAGE_SAG_MD_1MM_00*');
%     runlist = dir('GRE_FIELD_MAPPING_00*');
    for irun = 1:length(runlist)
    
        outdir_split = strsplit(runlist(irun).name, '_');
%         outdir = lower(outdir_split{3}); % for epi
        outdir = lower(outdir_split{1}); % for t1
        
        PREOUT = fullfile(outpath_subj, 'mri', outdir);
        mkdir(PREOUT)
        
        % count n volumes
        volume_list = dir(runlist(irun).name);
        if length(volume_list) > 602 % exported twice to same folder
            disp('Run exported twice!')
%             continue
        end
        
        if isempty(dir(fullfile(PREOUT, '*nii.gz'))) || overwrite %convert
%         nii_list = dir(fullfile(PREOUT, '*nii.gz'));
%         if length(nii_list) < 2 || overwrite %convert
            fprintf('\n\n\nConverting %s to %s\n\n\n', runlist(irun).name, PREOUT )
            unix(sprintf('%s -o %s %s', dcm2nii, PREOUT, runlist(irun).name));
        else
            disp('Nifti exists, skipping . . .')
        end
        
    end


end


