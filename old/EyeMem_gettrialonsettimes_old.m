function EyeMem_gettrialonsettimes()
% read in eyelink event structure and parse messages. Output:
% trial struct array with fields
% stimpic, startstim, endstim, quizpic, startquiz, endquiz, response, correct, RT


% datapath = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data/EYEMEM004/preproc2/fractals';
datapath = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data';  % /EYEMEM004/preproc2/fractals
cd(datapath)
cond_names = {'fractals' 'landscapes' 'naturals1' 'streets1' 'streets2'};

for isub = 1:101
    if isub < 10
        SUBJ = sprintf('EYEMEM00%d', isub);
    elseif isub < 100
        SUBJ = sprintf('EYEMEM0%d', isub);
    else
        SUBJ = sprintf('EYEMEM%d', isub);
    end
    fprintf('SUBJ: %s   ', SUBJ)
    cd(fullfile(datapath, SUBJ, 'eye'))
    runs = dir('run*');
    % get events and match them to the conditions
    allevent = cell(5,1);
    for irun = 1:length(runs)
        cd(runs(irun).name)
        load(sprintf('%s_%s_event.mat', SUBJ, runs(irun).name ));
        
        % get trial fields
        trialind = find(strcmp('Trial', {event.type}));
        val = {event.value};
        trialvals = val(trialind)';
        
        
        % find out which condition
        ind = find(ismember({event.type}, sprintf('Subj%d', isub)));
        condstr = tokenize(event(ind).value);
        condstr = condstr{6};
        [~,cond_ind] = ismember(condstr, cond_names);
        allevent{cond_ind} = event;
        cd ..
    end
    fprintf('%d runs found\n', irun)
    
    for icond = 1:5
        %     load('/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data/EYEMEM020/eye/run1/EYEMEM020_run1_event.mat')
        event = allevent{icond};
        if isempty(event); continue; end
        
        %unique({event.type})
        trialind = find(strcmp('Trial', {event.type}));
        val = {event.value};
        trialvals = val(trialind)';
        
        %% get timestamps in col 2
        trialtimes=[]; ctr=0;
        for i = 1:length(trialvals)
            temp = tokenize(trialvals{i}, '	');
            if contains(temp{2}, 'Start') & ~contains(temp{2}, 'quiz', 'IgnoreCase', true )
                ctr=ctr+1;
                trialtimes(ctr,:) = str2num(strtok(temp{2})); % , '	'
            end
        end
        
        %% find out scan start
        % 'MSG	5102637 Run1 fMRI scanning Start time 5077.28'
        beginind = find(contains({event.value}, 'scanning Start time'));
        
        scanstartval = val(beginind);
        temp = tokenize(scanstartval{:}, '	');
        scanstart = str2num(strtok(temp{2})); % , '	'
        
        %% get trial times wrt to scan onset
        trialonsets = (trialtimes - scanstart) /1000;
        trialonsets = trialonsets - 12; % 12 vols discarded
        outpath = fullfile(datapath, SUBJ, 'preproc2', cond_names{icond}, 'evfiles');
        mkdir(outpath);  
        outpathfile = fullfile(outpath, 'trialonsets.txt');
        fid = fopen(outpathfile, 'w');
        fprintf(fid, '%6.2f\n', trialonsets);
        fclose(fid);
        
    end
end
