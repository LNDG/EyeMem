function trlall = Eyemem_sorttrials()
% read in eyelink event structure and parse messages. Output: trialinfo as
% in fieldtrip
% columns:
% 1 stim start
% 2 stim end
% 3 quiz start
% 4 quiz end
% 5 condition (1:5)
% 6 trialno
% 7 stimpicno
% 8 quizpicno
% 9 response (yes: 1, no 0, omission nan)
% 10 correct (1 or 0)
% 11 RT in sec
% 12 hand (response mapping)
% 13 age group 1 young, 2 old
% 14 gender m 1, f 2

% OR
% trial struct array with fields
% stimpic, startstim, endstim, quizpic, startquiz, endquiz, response, correct, RT

ncols = 14;

% datapath = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data/EYEMEM004/preproc2/fractals';
datapath = '/Users/kloosterman/gridmaster2012/LNDG/EyeMem/data';  % /EYEMEM004/preproc2/fractals
cd(datapath)
cond_names = {'fractals' 'landscapes' 'naturals1' 'streets1' 'streets2'};

batch = EM_generate_subject_batch(); % info per subject

trlall = nan(30,ncols,101,5); % run trial subj col

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
    for irun = 1:length(runs)
        cd(runs(irun).name)
        load(sprintf('%s_%s_event.mat', SUBJ, runs(irun).name ));
        
        % get relevant trial MSG
        trialind = find(strcmp('Trial', {event.type}));
        val = {event.value};
        trialvals = val(trialind)';
        [~,trialvals] = strtok(trialvals);
        trl = nan(30, ncols); % trials, cols
        for i = 1:length(trialvals)
            
            msg = tokenize(trialvals{i}); %  ''    '4554310'    'Trial'    '1'    'Stim'    'landscapes'    'study'    'image262'    'Start'
            
            if any(contains(msg, 'during')) % stim msg
                %'Button'    'press'    'during'    'study'    'pic'
                continue 
            end            
            trialno = str2double(msg{4});
            trl(trialno, 6) = trialno;
            if any(contains(msg, 'Stim')) % stim msg
                if any(contains(msg, 'Start'))
                    if any(contains(msg, 'quiz'))
                        trl(trialno, 3) = str2double(msg{2}); % subtract start time later
                        if contains(msg{9}, 'image') % stimpicno
                            trl(trialno, 8) = str2double(msg{9}(6:end)); 
                        else
                            trl(trialno, 8) = str2double(msg{9}); % stimpicno
                        end
                    else % study pic
                        trl(trialno, 1) = str2double(msg{2});
                        [Lia,cond_ind] = ismember(msg, cond_names);
                        trl(trialno, 5) = cond_ind(Lia); % cond no
                        if contains(msg{8}, 'image')
                            trl(trialno, 7) = str2double(msg{8}(6:end)); % stimpicno
                        else
                            trl(trialno, 7) = str2double(msg{8}); % stimpicno
                        end
                    end
                elseif any(contains(msg, 'End'))
                    if any(contains(msg, 'Quiz'))
                        trl(trialno, 4) = str2double(msg{2}); % subtract start time later
                    else % study pic
                        trl(trialno, 2) = str2double(msg{2});
                    end
                end
            else % response MSG
                % response mapping
                if any(contains(val, 'Hand1'))
                    trl(trialno, 12) = 1; 
                else
                    trl(trialno, 12) = 2; 
                end
                %                 leftbutton = 'z'; % Yellow left % button broke couple of times
                %                 rightbutton = 'g'; % Green right
                % response: same or different
                % fprintf('%s', msg{6})
                if ~isnan(trl(trialno, 9))
                    continue % first bp counts: response already given this trial
                end
                if trl(trialno, 12) == 1 % Hand1: z = same pic, g = different pic
                    if strcmp(msg{6}, 'z') || strcmp(msg{6}, 'b')
                        trl(trialno, 9) = 1; % same
                    else
                        trl(trialno, 9) = 2; % different
                    end
                else
                    if strcmp(msg{6}, 'z') ||  strcmp(msg{6}, 'b')
                        trl(trialno, 9) = 2; % different
                    else
                        trl(trialno, 9) = 1; % same
                    end
                end    
                % correct
                if trl(trialno,7) == trl(trialno,8) && trl(trialno,9) == 1 || ...
                        trl(trialno,7) ~= trl(trialno,8) && trl(trialno,9) == 2
                    trl(trialno, 10) = 1; % correct
                else
                    trl(trialno, 10) = 0; % incorrect
                end                
                trl(trialno, 11) = str2double(msg{12}); % RT
                
            end
        end
        
        %% find out scan start
        % 'MSG	5102637 Run1 fMRI scanning Start time 5077.28'
        beginind = find(contains({event.value}, 'scanning Start time'));
        scanstartval = val(beginind);
        msg = tokenize(scanstartval{:}, '	');
        scanstart = str2num(strtok(msg{2})); % , '	'
        trl(:,1:4) = (trl(:,1:4) - scanstart - 12000) /1000; % 12 TR's discarded 
        trl(:,13) = strcmp(batch(isub).agegroup, 'old') + 1;
        trl(:,14) = strcmp(batch(isub).gender, 'f') + 1;
        
        % put in big mat
        trlall(:,:,isub, trl(1,5)) = trl;
        cd ..
    end
    fprintf('%d runs found\n', irun)
    
    for icond = 1:5
        trl = trlall(:,:,isub, icond);
        if isnan(trl(1)); continue; end
                
        outpath = fullfile(datapath, SUBJ, 'preproc2', cond_names{icond}, 'evfiles');
        if ~exist(outpath, 'dir');        mkdir(outpath); end

        % start times all trials
        trialonsets = trl(:,1);
        
        outpathfile = fullfile(outpath, 'trialonsets.txt');
        fid = fopen(outpathfile, 'w');
        fprintf(fid, '%6.2f\n', trialonsets);
        fclose(fid);

        % start times corrects
        ind = trl(:,10) == 1;
        trialonsets = trl(ind,1);
        
        outpathfile = fullfile(outpath, 'trialonsets_correct.txt');
        fid = fopen(outpathfile, 'w');
        fprintf(fid, '%6.2f\n', trialonsets);
        fclose(fid);
        
        % start times incorrects and omissions
        ind = trl(:,10) ~= 1;
        trialonsets = trl(ind,1);
        
        outpathfile = fullfile(outpath, 'trialonsets_incorrect.txt');
        fid = fopen(outpathfile, 'w');
        fprintf(fid, '%6.2f\n', trialonsets);
        fclose(fid);
    end
end
