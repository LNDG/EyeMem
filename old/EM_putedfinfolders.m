
% move edf files into folders with name EYEMEM#

PREIN = '/Volumes/FB-LIP/Eye_mem/data/eye';
cd(PREIN)

for isub = 1:101
    
    runlist = dir(sprintf('S%dp*', isub));
    
    if isub < 10
        dirname = sprintf('EYEMEM00%d', isub);
    elseif isub < 100
        dirname = sprintf('EYEMEM0%d', isub);
    else
        dirname = sprintf('EYEMEM%d', isub);
    end
    mkdir(dirname)
    cd(dirname)
    
    for irun = 1:length(runlist)
        sourcefile = fullfile(PREIN, runlist(irun).name );
        fprintf('Moving %s to %s\n', sourcefile, pwd)
        movefile(sourcefile)
    end
    
    cd ..
    
end