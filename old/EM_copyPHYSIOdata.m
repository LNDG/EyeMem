function EM_copyPHYSIOdata( )
% Convert raw data and make file structure

% overwrite = 0;

if ismac
    root_folder = '/Users/kloosterman/gridmaster2012';
else
    root_folder = '/home/mpib';
end
PREIN = fullfile(root_folder, 'LNDG/EyeMem/physio' );
% PREIN = '/Volumes/FB-LIP/Eye_mem/data/physio' ;
PREOUT = fullfile(root_folder, 'LNDG/EyeMem/data' );

cd(PREIN)
subj_list = dir(fullfile(PREIN, 'EyeMem*'));

for isub = 1:length(subj_list)
    
    cd(fullfile(PREOUT, upper(subj_list(isub).name), 'physio')) % go to destination

    sourcefile = fullfile(PREIN, subj_list(isub).name, '*');
%     unix(sprintf('sudo chown kloosterman %s', sourcefile))
    fprintf('Moving %s to %s\n', sourcefile, pwd)
%     try
    copyfile(sourcefile)
%     catch
%         fprintf('%s not copied', sourcefile)
%     end
end
cd(PREIN)

