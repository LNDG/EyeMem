function [behav] = EM_analysebehavioral
%Read in behavior text files, compute d', RT, etc.

% load('participantinfo.mat') % 
load /Users/kloosterman/Dropbox/tardis_code/MATLAB/eyemem_analysis/participantinfo/participantinfo.mat

category_labels = {'fractals'	'landscapes'	'naturals1'	'streets1'	'streets2'}; %1-5

exp_phases = {'study' 'test'};
ntrials_per_run = [30 60];

SUBJ= [9:101]; % TODO specify further?

behav=[];
% for iphase = 1:2
%   behav.(exp_phases{iphase}).propcorrect = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).dprime = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).criterion = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).RT = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).RTsd = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).RTsd2 = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).omissions = nan(length(SUBJ), 6); % SUBJ, categories
%   
% %   behav.(exp_phases{iphase}).driftrate = nan(length(SUBJ), 6); % SUBJ, categories
% %   behav.(exp_phases{iphase}).boundsep = nan(length(SUBJ), 6); % SUBJ, categories
% %   behav.(exp_phases{iphase}).nondectime = nan(length(SUBJ), 6); % SUBJ, categories
%   
%   behav.(exp_phases{iphase}).RT_hits = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).RT_misses = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).RT_fas = nan(length(SUBJ), 6); % SUBJ, categories
%   behav.(exp_phases{iphase}).RT_crs = nan(length(SUBJ), 6); % SUBJ, categories
% end

if ismac
  PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/preproc/behavior/raw';
  PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/preproc/behavior';
else
%   basepath = '/home/mpib';
end
% PREOUT = fullfile(basepath, 'LNDG/EyeMem/plots' );
cd(PREIN)

ddm_dat{1} = []; ddm_dat{2} = [];
for isub = 1:length(SUBJ)
  
  studied_pics = []; % to keep track of pics seen during study (the same for all subjects)
  for icat = 1:5
    studied_pics.(category_labels{icat}) = {};
  end
  n_omissions = 0;
  
  singletrial{1} = []; singletrial{2} = [];
  for iphase = 1:2 % study, test
    
    %     if isempty(subject_batch(isub).SUBJ)
    %       fprintf('SUBJ %d dropped\n', isub)
    %       continue
    %     end
    
    txtfile = fullfile(PREIN, sprintf('S%d_%s_log.txt', SUBJ(isub), exp_phases{iphase} ));
    disp(txtfile)
    %     txtfile = dir(sprintf('*%s*.txt', exp_phases{iphase})); % PREIN, subj_list(isub).name,
    
    fid = fopen(txtfile, 'rt');
    if fid == -1;
      fprintf('not found:   %s\n', txtfile)
      continue
    end
    tline = fgetl(fid);
    strtok = strsplit(tline, ' ');
    subNo = str2double(strtok{1});%                             fprintf(datafilepointer,'%i %i %s %i %i %s %s %s %s %i %i\n', ...

    ac_accum = 0; resp=[]; rt = []; rts_correct{isub} = [];
    rt_hits = [];  rt_misses = []; rt_crs = []; rt_fas = [];
    hits=0; misses=0; fas=0; crs=0;
    n_absent = 0; n_present = 0;
    behav.(exp_phases{iphase})(subNo).p_repeatbalanced = NaN(5,2);
    while ~feof(fid)
      %             i = i + 1;
      %             ft_progress(i/n_lines);
      
      tline = fgetl(fid);
      strtok = strsplit(tline, ' ');
      
%       subNo = str2double(strtok{1});%                             fprintf(datafilepointer,'%i %i %s %i %i %s %s %s %s %i %i\n', ...
      hand  = str2double(strtok{2});
      phaselabel  = strtok{3};
      irun = str2double(strtok{4});
      itrial = str2double(strtok{5});
      resptmp = strtok{6};
      category = strtok{7};
      view_file = strtok{8};
      test_file = strtok{9}; % nan for test phase
      ac = str2double(strtok{10});
      rt(itrial,:) = str2double(strtok{11}) / 1000; % convert to s
      
      if rt(itrial,:) < 0 % subj 60 has 1 negative rt
        rt(itrial,:) = NaN;
      end
      %               leftbutton = 'z'; % Yellow left
      %         rightbutton = 'g'; % Green right
      if strcmp(resptmp, 'z') | strcmp(resptmp, 'LeftGUI')
        resp(itrial,1) = 1;
      elseif strcmp(resptmp, 'g') | strcmp(resptmp, 'RightGUI')
        resp(itrial,1) = 2;
      else
        resp(itrial,1) = NaN;
      end
      %             if iphase == 2 % counted from stim offset, only affects
      %             nondec time
      %                 rt(itrial,:) = rt(itrial,:) - 0.75;
      %             end
      if isnan(resp)
        n_omissions = n_omissions + 1;
      end
      
      ac_accum = ac_accum + ac;
      
      if iphase == 1 % remember which pics were shown during study
        studied_pics.(category) = [studied_pics.(category); {view_file}];
        target_present = strcmp(view_file, test_file); % target present?
      else
        target_present = ismember(view_file, studied_pics.(category)); % find out if view_file was shown during study
      end
      
      % get stim pic nr from view_file string
      if contains(view_file, 'image') % for landscapes and streets2 format is imageXXX
        view_file = view_file(6:end);
      end
      temp = tokenize(view_file, '.');
      picno = str2double(temp{1});
      
      if target_present
        n_present = n_present + 1;
        if ac == 1 % correct, meaning it was a Yes
          hits = hits + 1;
          rts_correct{isub} = [rts_correct{isub}; rt(itrial,:)];
          rt_hits = [rt_hits; rt(itrial,:)];
        else %if ~ac % incorrect
          misses = misses + 1;
          rt_misses = [rt_misses; rt(itrial,:)];
        end
      else % target absent
        n_absent = n_absent + 1;
        if ac % correct, meaning it was a No
          crs = crs + 1;
          rts_correct{isub} = [rts_correct{isub}; rt(itrial,:)];
          rt_crs = [rt_crs; rt(itrial,:)];
        elseif ~ac % incorrect
          fas = fas + 1;
          rt_fas = [rt_fas; rt(itrial,:)];
        end
      end
      
      icat = find(strcmp(category_labels, category));
      
      % For HDDM: Put subjid, category, stim, ac, rt agegroup
      if Participants.group(SUBJ(isub)) == 'young'
        ageind = 0;
      else
        ageind = 1;
      end
      ddm_dat{iphase} = [ddm_dat{iphase}; isub-1  icat target_present resp(itrial,:) ac rt(itrial,:) ageind];
      singletrial{iphase} = [singletrial{iphase}; icat target_present resp(itrial,:) ac rt(itrial,:) picno];
      
      if itrial == ntrials_per_run(iphase) % calculate things per run
        
        if n_omissions < 4 % only allow runs with < 10 % missed responses
%           behav.(exp_phases{iphase}).propcorrect(subNo, icat) = ac_accum/itrial; % TODO omit first trial(s)
          behav.(exp_phases{iphase})(subNo).propcorrect(icat) = ac_accum/itrial; % TODO omit first trial(s)
          
          Hitrate = hits/n_present;
          if Hitrate == 1; Hitrate = 0.95; end
          if Hitrate == 0; continue; end % wrong response pressed or st
          
          FArate = fas/n_absent;
          if FArate == 0; FArate = 0.05; end
          
          behav.(exp_phases{iphase})(subNo).dprime(icat) = norminv(Hitrate) - norminv(FArate); % TODO omit first trial(s)
          behav.(exp_phases{iphase})(subNo).criterion(icat) = -0.5 * (norminv(Hitrate) + norminv(FArate)); % TODO omit first trial(s)
          
          behav.(exp_phases{iphase})(subNo).RT(icat) = nanmean(rt);
          behav.(exp_phases{iphase})(subNo).RTsd(icat) = nanstd(rt);
          behav.(exp_phases{iphase})(subNo).RTsd2(icat) = nanstd(rt) / nanmean(rt);
          
          behav.(exp_phases{iphase})(subNo).RT_hits(icat) = nanmean(rt_hits);
          behav.(exp_phases{iphase})(subNo).RT_misses(icat) = nanmean(rt_misses);
          behav.(exp_phases{iphase})(subNo).RT_fas(icat) = nanmean(rt_fas);
          behav.(exp_phases{iphase})(subNo).RT_crs(icat) = nanmean(rt_crs);
          
%           p_repeatbalanced(irun,1) = sum(diff(button) == 0 & button(2:end,:) == 1) / (sum(button(2:end)==1));
%           p_repeatbalanced(irun,2) = sum(diff(button) == 0 & button(2:end,:) == 2) / (sum(button(2:end)==2));
          behav.(exp_phases{iphase})(subNo).p_repeatbalanced(icat,1) = sum(diff(resp) == 0 & resp(2:end,:) == 1) / (sum(resp(2:end)==1));
          behav.(exp_phases{iphase})(subNo).p_repeatbalanced(icat,2) = sum(diff(resp) == 0 & resp(2:end,:) == 2) / (sum(resp(2:end)==2));
          
          
        end
        
        behav.(exp_phases{iphase})(subNo).omissions(icat) = n_omissions;
        
        ac_accum = 0; rt = [];
        hits=0; misses=0; fas=0; crs=0;
        n_absent = 0; n_present = 0;     n_omissions = 0;
        
      end
    end
    behav.(exp_phases{iphase})(subNo).singletrial = singletrial{iphase};
    behav.(exp_phases{iphase})(subNo).singletrialleg = 'condition target_present response accuracy RT picno';
    fclose(fid);
    disp 'compute history bias'
    

  end
end

for iphase = 1:2 % study, test
%   behav.(exp_phases{iphase}).propcorrect(:,6) = nanmean(behav.(exp_phases{iphase}).propcorrect(:,1:5) ,2);
%   behav.(exp_phases{iphase}).dprime(:,6) = nanmean(behav.(exp_phases{iphase}).dprime(:,1:5) ,2);
%   behav.(exp_phases{iphase}).criterion(:,6) = nanmean(behav.(exp_phases{iphase}).criterion(:,1:5) ,2);
%   
%   behav.(exp_phases{iphase}).RT(:,6) = nanmean(behav.(exp_phases{iphase}).RT(:,1:5) ,2);
%   behav.(exp_phases{iphase}).RTsd(:,6) = nanmean(behav.(exp_phases{iphase}).RTsd(:,1:5) ,2);
%   behav.(exp_phases{iphase}).RTsd2(:,6) = nanmean(behav.(exp_phases{iphase}).RTsd2(:,1:5) ,2);
%   
%   behav.(exp_phases{iphase}).omissions(:,6) = nansum(behav.(exp_phases{iphase}).omissions(:,1:5),2);
%   
%   behav.(exp_phases{iphase}).RT_hits(:,6) = nanmean(behav.(exp_phases{iphase}).RT_hits(:,1:5) ,2);
%   behav.(exp_phases{iphase}).RT_misses(:,6) = nanmean(behav.(exp_phases{iphase}).RT_misses(:,1:5) ,2);
%   behav.(exp_phases{iphase}).RT_fas(:,6) = nanmean(behav.(exp_phases{iphase}).RT_fas(:,1:5) ,2);
%   behav.(exp_phases{iphase}).RT_crs(:,6) = nanmean(behav.(exp_phases{iphase}).RT_crs(:,1:5) ,2);
  
  ddm_dat{iphase} = ddm_dat{iphase}(~isnan(ddm_dat{iphase}(:,5)),:);
  % save ddm_dat to csv:             % For HDDM: Put subjid, category, stim, ac, rt
  csv_file = sprintf('/Users/kloosterman/Dropbox/PROJECTS/EyeMem/HDDM/EyeMem_hddm_%s.csv', exp_phases{iphase});
  fid = fopen(csv_file, 'w');
  fprintf(fid, 'subj_idx,category,stim,response,accuracy,rt,age\n');
  dlmwrite(csv_file , ddm_dat{iphase},'delimiter',',','-append');
end

% % TODO Fit EZ DDM and get drift rate, boundary sep and non dec time
% for isub = 1:length(SUBJ)
%   for iphase = 1:2 % study, test
% %     if isempty(SUBJ(isub).SUBJ)
% %       fprintf('SUBJ %d dropped\n', isub)
% %       continue
% %     end
%     fprintf('%s EZ DDM fitting . . .\n', SUBJ(isub))
%     
%     propCorrect = behav.(exp_phases{iphase}).propcorrect(isub,6);
%     if propCorrect == 0 || propCorrect == 0.5 || propCorrect == 1
%       warning('Propcorrect is off')
%     end
%     unix_string = sprintf('R CMD BATCH --no-timing --no-save --no-restore ''--args Pc=%g VRT=%g MRT=%g'' /Users/kloosterman/Dropbox/PROJECTS/EyeMem/DDM_EZ/runvaTer.R /Users/kloosterman/Dropbox/PROJECTS/EyeMem/DDM_EZ/runvaTer.out', ...
%       propCorrect, var(rts_correct{isub}), mean(rts_correct{isub}) );
%     test = unix(unix_string);
%     ddmpars = load('/Users/kloosterman/Dropbox/PROJECTS/EyeMem/DDM_EZ/DDMfit.txt');
%     behav.(exp_phases{iphase}).driftrate(isub,6) = ddmpars(1);
%     behav.(exp_phases{iphase}).boundsep(isub,6) = ddmpars(2);
%     behav.(exp_phases{iphase}).nondectime(isub,6) = ddmpars(3);
%     
%   end
% end

behav.participants = Participants;
% behav.agegroup = transpose(strcmp({SUBJ.agegroup}, 'old' ) + 1); % young is 1, old is 2
% dropped_subj = cellfun(@isempty, {SUBJ.agegroup});
% behav.agegroup(find(dropped_subj)) = nan;
% todo save behav

disp('Saving')
disp(fullfile(PREOUT, 'Eyemem_behavior.mat'))
save(fullfile(PREOUT, 'Eyemem_behavior.mat'), 'behav')
