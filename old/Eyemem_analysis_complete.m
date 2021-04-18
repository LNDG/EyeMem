function Eyemem_analysis_complete()

% run edf data analysis with fieldtrip
cd('~'); %clear all; close all;

matlabpath = '/Users/kloosterman/Dropbox/MATLAB/';
addpath(fullfile(matlabpath, 'toolbox/fieldtrip-20150803')) %inc JJ edit ft_artifact_zvalue
ft_defaults

basepath = '/Users/kloosterman/Dropbox/PROJECTS/Memory_Eyetracking/Pilotstudie/' ;
scriptsfolder = fullfile(basepath, 'analysis');
addpath(genpath(scriptsfolder))
edf2asc = fullfile(basepath, 'analysis', 'edf2asc');

batchlists = {
% %     'batch_eyemem_S9' % S11 onwards with fixed stims TODO use other categories for S1-10

% Full subjects (10 runs):
%     'batch_eyemem_S13'
%     'batch_eyemem_S14'
%     'batch_eyemem_S15'
%     'batch_eyemem_S16'
%     'batch_eyemem_S20'
%     'batch_eyemem_S21'
    'batch_eyemem_S22'
%     'batch_eyemem_S23'
%     'batch_eyemem_S24'
%     'batch_eyemem_S25'
};

% trigger = 'resp';  %one at time!
trigger = 'stim';
% trigger = 'fixation'; % baseline fixation is in separate file

begtim = -2;
endtim = 2;

if strcmp(trigger, 'resp')
    artfcrittoilim = [-1 1]; % trials w blinks in this window are rejected
elseif strcmp(trigger, 'stim')
    artfcrittoilim = [-1 1]; % trials w blinks in this window are rejected
end

overwrite = 0;

disp(batchlists);
disp(trigger);
% disp(['overwrite ' overwrite])
tic
%% ET analysis
fprintf('Running eye mem pilot analysis . . .\n\n')

cfg2=[]; % visualization purposes to inspect data
cfg2.channel = 'p';
cfg2.method = 'channel';

for ibatch = 1:length(batchlists) % batch is list of runs
    
	batch=[];
    eval(batchlists{ibatch}); %load in batchlist file, batch, PREOUT and PREIN come out
    disp(batchlists{ibatch})
    PREIN = fullfile(basepath, 'data', 'edf', PRE);
    cd(PREIN)
    PREOUT = fullfile(basepath, 'data', 'preproc', PRE);
    mkdir(PREOUT)
    rundata = {};
    for irun=1:length(batch)

        if isempty(batch(irun).dataset) % || (exist([studyoutfile '.mat'], 'file') || (exist([testoutfile '.mat'], 'file') && ~overwrite) % % in case batch entries are commented out, or if the output file exists and no overwrite enabled
            warning('Batch is empty, output mat file exists or overwrite switched off!')
            continue 
        end       
        
        [~,infile] = fileparts(batch(irun).dataset);
        if ~exist([infile  '.asc'], 'file') %convert edf to ascii
            unix(sprintf('%s -sg %s', edf2asc, infile)); %-y
        end
        if ~exist([infile  '_ascdat.mat'], 'file') % && overwrite
            disp('Reading asc . . .') % read asc file into matlab
            ascdat = read_eyelink_ascNK( [infile  '.asc'] );
            save([infile '_ascdat.mat'], 'ascdat')            

            % create event structure
            evcell = cell(length(ascdat.msg),1);
            event=struct('type', evcell, 'sample', evcell, 'value', evcell, 'offset', evcell, 'duration', evcell );
            for i=1:length(ascdat.msg)
                strtok = tokenize(ascdat.msg{i});
                event(i).type = strtok{3};
                smpstamp = find(ascdat.dat(1,:) == str2double(strtok{2})); % find sample index of trigger in ascii dat
                if ~isempty(smpstamp)
                    event(i).sample = smpstamp(1);
                else
                    event(i).sample = nan; % give nan if no data was recorded at msg time
                end
                event(i).value = ascdat.msg{i}; % event(i).value = [strtok{3:end}]; %trigger value: e.g. ResponseMIBOff
            end
            save([infile '_event.mat'], 'event')
        else
            load([infile '_ascdat.mat'])
            load([infile '_event.mat'])
        end
        
        
        %make data struct
        data = [];
        data.label ={'h'; 'v'; 'p'};
        data.trial = {ascdat.dat(2:end,:)};
        data.fsample = ascdat.fsample;
        data.time = {0:1/data.fsample:length(ascdat.dat(1,:))/data.fsample-1/data.fsample};
        data.sampleinfo = [1 length(ascdat.dat(1,:))];
        
        %                                 ft_rejectvisual(cfg2,data)
        
        % Interpolate blinks (linear) interpolate across multiple blinks
        if length(ascdat.eblink) > 0
            arg1 = repmat({'%*s%*s%d%d'}, length(ascdat.eblink), 1);
            blinktimes = cellfun(@sscanf, ascdat.eblink, arg1, 'UniformOutput', false); % parse blinktimes from ascdat
            blinktimes = cell2mat(cellfun(@transpose, blinktimes, 'UniformOutput', false)); %transpose and turn into matrix
            timestamps = ascdat.dat(1,:);
            if ascdat.fsample == 250
                blinksmp = arrayfun(@(x) find(timestamps > x-3 & timestamps < x+3 , 1,'first'), blinktimes, 'UniformOutput', true ); %find sample indices of blinktimes in timestamps
            else
                blinksmp = arrayfun(@(x) find(timestamps == x, 1,'first'), blinktimes, 'UniformOutput', true ); %find sample indices of blinktimes in timestamps
            end
            blinksmp(:,1) = blinksmp(:,1) - 0.1*ascdat.fsample; % some extra padding: 0.1 cf de gee etal
            blinksmp(:,2) = blinksmp(:,2) + 0.1*ascdat.fsample;
            blinksmp = blinksmp(find(blinksmp(:,1) > 0),:);
            % remove blinks during which no data was recorded
            [blinksin,~] = find(blinksmp(:,2) < size(data.time{1},2));
            blinksmp = blinksmp(blinksin,:);
            
            curblink=0;
            nblinks = size(blinksmp,1);
            disp('Interpolating blinks')
            for iblink=1:size(blinksmp,1)
                if iblink ~= curblink+1; % continue until next series of blinks
                    continue
                end
                curblink = iblink;
                if curblink+1 < nblinks
                    while blinksmp(curblink+1,1) - blinksmp(curblink,2) < 0.5*data.fsample % if succeeding blinks are < 0.5 s away from each other, interpolate at once
                        disp('Blink series detected')
                        %                         disp(iblink)
                        %                         disp(blinksmp(curblink+1,1) - blinksmp(curblink,2))
                        curblink = curblink+1; % if yes, take blink end of succeeding blink
                        if curblink+1 > nblinks
                            break % no more successive nblinks
                        end
                    end
                end
                data.trial{1}(1,blinksmp(iblink,1):blinksmp(curblink,2)) = linspace(data.trial{1}(1,blinksmp(iblink,1)), data.trial{1}(1,blinksmp(curblink,2)), length(blinksmp(iblink,1):blinksmp(curblink,2))); % interpolate x
                data.trial{1}(2,blinksmp(iblink,1):blinksmp(curblink,2)) = linspace(data.trial{1}(2,blinksmp(iblink,1)), data.trial{1}(2,blinksmp(curblink,2)), length(blinksmp(iblink,1):blinksmp(curblink,2))); % interpolate y
                data.trial{1}(3,blinksmp(iblink,1):blinksmp(curblink,2)) = linspace(data.trial{1}(3,blinksmp(iblink,1)), data.trial{1}(3,blinksmp(curblink,2)), length(blinksmp(iblink,1):blinksmp(curblink,2))); % interpolate p
            end
            %                                             ft_rejectvisual(cfg2,data);
            % figure; plot(data.time{1}, data.trial{1})
            % figure; plot(data.trial{1})
        end
        
        % define trials
        cfg=[];
        cfg.headerfile = fullfile(PREIN, batch(irun).dataset);
        cfg.subjno = batch(irun).subj(2:end);
        cfg.trialdef.trg = trigger; %stim or resp
        cfg.trialdef.begtim = begtim; 
        cfg.trialdef.endtim = endtim;
        cfg.event = event;
        cfg.fsample = ascdat.fsample;
        cfg.trialfun = 'Eyemem_ET_analysis_sortTrials';
        if strcmp(batch(irun).type, 'test') 
            infile = fullfile(PREOUT,  sprintf('%s_%s_%s_data.mat', batch(irun).subj, 'study', trigger));
            study = load(infile); % data from study, needed for oldnew correctness in test            
            cfg.studytrialinfo = study.data.trialinfo;
        end
        cfg = ft_definetrial(cfg); % make trl matrix
                
        %make trials
        rundata{batch(irun).exp} = ft_redefinetrial(cfg, data); 
                
        if irun == 5 || irun == 10 % append all runs per phase, save
            
            data = ft_appenddata([], rundata{:});
            
            outfile = fullfile(PREOUT,  sprintf('%s_%s_%s_data.mat', batch(irun).subj, batch(irun).type, trigger));
            [pathstr, name] = fileparts(outfile);
            fprintf('Saving %s to...\n %s\n', name, pathstr)
            save(outfile, 'data');

        end
    end % irun

end % ibatch
toc








%%
%         % do timelock per condition, save
%         cfg=[];
%         cfg.keeptrials = 'yes';
%         timelock = ft_timelockanalysis(cfg,data);

%         fixdata=data; % keep unfiltered data for fixation error rejection

%         % LOWpass filter data
%         cfg=[];
%         cfg.lpfilter = 'yes';
%         cfg.lpfreq = 4; 
%         cfg.lpfiltord = 3;
%         data = ft_preprocessing(cfg, data);
        %                             ft_rejectvisual(cfg2,contdata)
