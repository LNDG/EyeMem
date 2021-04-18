function EM_copyBEHAVdata( )
% Convert raw data and make file structure

% overwrite = 0;

if ismac
    root_folder = '/Users/kloosterman/gridmaster2012';
else
    root_folder = '/home/mpib';
end
PREIN = fullfile(root_folder, 'LNDG/EyeMem/behavior_raw' );
PREOUT = fullfile(root_folder, 'LNDG/EyeMem/data' );

cd(PREIN)
subj_list = dir(fullfile(PREIN, 'EYEMEM*'));

for isub = 1:length(subj_list)
    
    cd(fullfile(PREOUT, subj_list(isub).name, 'behavior')) % go to destination

    sourcefile = fullfile(PREIN, subj_list(isub).name, '*.txt');
    fprintf('Moving %s to %s\n', sourcefile, pwd)
    copyfile(sourcefile)

end
cd(PREIN)

