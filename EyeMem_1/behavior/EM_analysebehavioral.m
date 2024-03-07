function [behavior] = EM_analysebehavioral()
%Read in behavior text files, compute d', RT, etc.
% TODO organize MEG2afc way, track age with Participants table
% dimord subj_phase_cond

behavior=[];
load('/Users/kloosterman/Documents/GitHub/EyeMem/EyeMem_1/participantinfo/participantinfo.mat', 'Participants')
behavior.participants = Participants;

category_labels = {'fractals'	'landscapes'	'naturals1'	'streets1'	'streets2'}; %1-5
exp_phases = {'study' 'test'};
ntrials_per_run = [30 60];

SUBJ= [9:101]; % TODO specify further?
% SUBJ = [69:70];
% find out which SUBJ are still in the mix in fMRI
% PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsdprime/std_3bins/fixednbins/linearfit_fitcoeff1/gaze-specific';
% cd(fullfile(PREIN))
% SUBJ = [];
% subjlist = dir('sub*_BfMRIsessiondata.mat');
% for isub=1:length(subjlist)
%   tmp = tokenize(subjlist(isub).name, '_');
%   SUBJ(isub,1) = find(behavior.participants.participant_id == tmp{1});
% end
% SUBJ = sort(SUBJ);

behavior.SUBJ = SUBJ;
behavior.dimord = 'subj_phase_cond';

PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/raw';
behavior.PREOUT = fullfile('/Users/kloosterman/Dropbox/PROJECTS/EyeMem/plots');
cd(PREIN)

ddm_dat{1} = []; ddm_dat{2} = [];
for isub = 1:length(SUBJ)
  
  studied_pics = []; % to keep track of pics seen during study (the same for all subjects)
  for icond = 1:5
    studied_pics.(category_labels{icond}) = {};
  end
  n_omissions = 0;
  
  singletrial{1} = []; singletrial{2} = [];
  for iphase = 1:2 % study, test
        
    txtfile = fullfile(PREIN, sprintf('S%d_%s_log.txt', SUBJ(isub), exp_phases{iphase} ));
    disp(txtfile)    
    fid = fopen(txtfile, 'rt');
    if fid == -1
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
%     behavior.(exp_phases{iphase})(subNo).p_repeatbalanced = NaN(5,2);
    while ~feof(fid)      
      tline = fgetl(fid);
      strtok = strsplit(tline, ' ');
      
      %       subNo = str2double(strtok{1});%                             fprintf(datafilepointer,'%i %i %s %i %i %s %s %s %s %i %i\n', ...
      hand  =       str2double(strtok{2});
      phaselabel  = strtok{3};
      irun =        str2double(strtok{4});
      itrial =      str2double(strtok{5});
      resptmp =     strtok{6};
      category =    strtok{7};
      view_file =   strtok{8};
      test_file =   strtok{9}; % nan for test phase
      ac = str2double(strtok{10});
      rt(itrial,:) = str2double(strtok{11}) / 1000; % convert to s
      
      if rt(itrial,:) <= 0 % subj 60 has 1 negative rt
        rt(itrial,:) = NaN;
      end
      %               leftbutton = 'z'; % Yellow left
      %         rightbutton = 'g'; % Green right
      if hand == 1
        if strcmp(resptmp, 'z') | strcmp(resptmp, 'b') | strcmp(resptmp, 'LeftGUI')
          resp(itrial,1) = 2; % 1 = old, 2 = new
        elseif strcmp(resptmp, 'g') | strcmp(resptmp, 'RightGUI')
          resp(itrial,1) = 1;
        else
          resp(itrial,1) = NaN;
        end
      elseif hand == 2
        if strcmp(resptmp, 'z') | strcmp(resptmp, 'b') | strcmp(resptmp, 'LeftGUI')
          resp(itrial,1) = 1;
        elseif strcmp(resptmp, 'g') | strcmp(resptmp, 'RightGUI')
          resp(itrial,1) = 2;
        else
          resp(itrial,1) = NaN;
        end
      end
      %             if iphase == 2 % counted from stim offset, only affects
      %             nondec time
      %                 rt(itrial,:) = rt(itrial,:) - 0.75;
      %             end
      if isnan(resp(itrial,1))
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
      
      icond = find(strcmp(category_labels, category));
      
      % For HDDM: Put subjid, category, stim, ac, rt agegroup
      if Participants.group(SUBJ(isub)) == 'young'
        ageind = 0;
      else
        ageind = 1;
      end
      %       ddm_dat{iphase} = [ddm_dat{iphase}; isub-1  icond target_present resp(itrial,:) ac rt(itrial,:) ageind];
      ddm_dat{iphase} = [ddm_dat{iphase}; subNo icond target_present resp(itrial,:) ac rt(itrial,:) ageind];
      singletrial{iphase} = [singletrial{iphase}; icond target_present resp(itrial,:) ac rt(itrial,:) picno];
      
      if any(rt)==0
        warning('zero RT remains')
      end
      
      if itrial == ntrials_per_run(iphase) % calculate things per run
        
        if n_omissions < 100 % only allow runs with < 10 % missed responses
          behavior.propcorrect(subNo,iphase,icond) = ac_accum/itrial; % TODO omit first trial(s)
          
          Hitrate = hits/n_present;
          if Hitrate == 1; Hitrate = 0.95; end
          if Hitrate == 0;
            warning('0 Hitrate')
            continue; 
          end % wrong response pressed or st
          
          FArate = fas/n_absent;
          if FArate == 0; FArate = 0.05; end
          
          behavior.dprime(subNo,iphase,icond) = norminv(Hitrate) - norminv(FArate); % TODO omit first trial(s)
%           if  behavior.dprime(subNo,iphase,icond) == 0
          if  subNo == 27
            warning('stop')
          end
          
          behavior.criterion(subNo,iphase,icond) = -0.5 * (norminv(Hitrate) + norminv(FArate)); % TODO omit first trial(s)
          behavior.propNo(subNo,iphase,icond) = 1-(hits+fas)/(n_present+n_absent);
          
          behavior.RT(subNo,iphase,icond) = nanmean(rt);
          behavior.RTsingletrial{iphase}(subNo,:,icond) = rt;
          behavior.RTsd(subNo,iphase,icond) = nanstd(rt);
          behavior.RTsd2(subNo,iphase,icond) = nanstd(rt) / nanmean(rt);
          
          behavior.RT_hits(subNo,iphase,icond) = nanmean(rt_hits);
          behavior.RT_misses(subNo,iphase,icond) = nanmean(rt_misses);
          behavior.RT_fas(subNo,iphase,icond) = nanmean(rt_fas);
          behavior.RT_crs(subNo,iphase,icond) = nanmean(rt_crs);
          
%           %           p_repeatbalanced(irun,1) = sum(diff(button) == 0 & button(2:end,:) == 1) / (sum(button(2:end)==1));
%           %           p_repeatbalanced(irun,2) = sum(diff(button) == 0 & button(2:end,:) == 2) / (sum(button(2:end)==2));
%           behavior.p_repeatbalanced(icond,1) = sum(diff(resp) == 0 & resp(2:end,:) == 1) / (sum(resp(2:end)==1));
%           behavior.p_repeatbalanced(icond,2) = sum(diff(resp) == 0 & resp(2:end,:) == 2) / (sum(resp(2:end)==2));
          
        end
        
        behavior.omissions(subNo,iphase,icond) = n_omissions;
        
        ac_accum = 0; rt = [];
        hits=0; misses=0; fas=0; crs=0;
        n_absent = 0; n_present = 0;     n_omissions = 0;
        
      end
    end
    fclose(fid);
    behavior.singletrial{subNo,1} = singletrial;
  end
end
behavior.singletrialleg = 'condition target_present response accuracy RT picno';

disp 'average over conditions'
bmeas = {'dprime' 'criterion' 'propcorrect' 'RT' 'RTsd' 'RTsd2' 'RT_misses' 'RT_hits' 'RT_fas' 'RT_crs' 'omissions' 'propNo'}; % p_repeatbalanced dim5 is LR
for im = 1:length(bmeas)
  if strcmp(bmeas{im}, 'omissions')
    behavior.(bmeas{im})(:,:,6) = nansum(behavior.(bmeas{im}), 3); % avg over diff
  else
    behavior.(bmeas{im})(:,:,6) = nanmean(behavior.(bmeas{im}), 3); % avg over diff
  end
end

for iphase = 1:2 % study, test
  ddm_dat{iphase} = ddm_dat{iphase}(~isnan(ddm_dat{iphase}(:,5)),:);
  % save ddm_dat to csv:             % For HDDM: Put subjid, category, stim, ac, rt
  csv_file = sprintf('/Users/kloosterman/Dropbox/PROJECTS/EyeMem/HDDM/EyeMem_hddm_%s.csv', exp_phases{iphase});
  fid = fopen(csv_file, 'w');
  fprintf(fid, 'subj_idx,category,stim,response,accuracy,rt,age\n');
  dlmwrite(csv_file , ddm_dat{iphase},'delimiter',',','-append');
end

%% Add DDM data Laura (7 subjects missing)
% 
% disp 'load ddm fits'
% model_names = {'eyemem1_params_biasmodel_YA' 'eyemem1_params_biasmodel_OA' }; % chi_accuracy_basic_runs
% pars = {'a' 't' 'v' 'z' 'dc'};
% ddmpath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/HDDM/Laura';
% for imodel = 1:length(model_names)
%   model_name = model_names{imodel};
%   ddmfit = readtable(fullfile(ddmpath, [model_name '.csv'] ));
%   
%   behavior.(model_name) = [];
%   behavior.(model_name).dimord = 'subj';
%   
%   for ipar = 1:length(pars)
%     ddmstr = sprintf('%s_subj.',  pars{ipar});
%     line_inds = contains(ddmfit.Var1, ddmstr);
%     % get subjID from Var1 strings, to use as index
%     subj = ddmfit.Var1(line_inds);
%     subj = cellfun(@(x) tokenize(x, '.'), subj, 'UniformOutput', false);
%     subj = vertcat(subj{:});
%     subj = cellfun(@str2double, subj(:,2));
%     
%     temp = NaN(101,1);
%     temp(subj,1) = ddmfit.mean(line_inds);
%     behavior.(model_name).(pars{ipar}) = temp;
%     behavior.ddmLaura.(pars{ipar})(subj,1) = ddmfit.mean(line_inds);
% %     behavior.eyemem1_params_biasmodel.(pars{ipar})(behavior.eyemem1_params_biasmodel.(pars{ipar}) == 0) = NaN;
%   end
%   
%   if ~any(line_inds); warning('No data found'); continue;  end
% end
% 
% disp 'TODO remove zeros from eyemem1_params_biasmodel!'
% behavior.participants = Participants;

%% Add DDM data Niels (separate YA and OA)
% 
% disp 'load ddm fits'
% model_names = {'params_HDDMbias_YA' 'params_HDDMbias_OA' }; % chi_accuracy_basic_runs
% pars = {'a' 't' 'v' 'z' 'dc'};
% ddmpath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/HDDM/Niels';
% for imodel = 1:length(model_names)
%   model_name = model_names{imodel};
%   ddmfit = readtable(fullfile(ddmpath, [model_name '.csv'] ));
%   
%   behavior.(model_name) = [];
%   behavior.(model_name).dimord = 'subj';
%   
%   for ipar = 1:length(pars)
%     ddmstr = sprintf('%s_subj',  pars{ipar});
%     line_inds = contains(ddmfit.Var1, ddmstr);
%     % get subjID from Var1 strings, to use as index
%     subj = ddmfit.Var1(line_inds);
%     subj = cellfun(@(x) tokenize(x, '.'), subj, 'UniformOutput', false);
%     subj = vertcat(subj{:});
%     subj = cellfun(@str2double, subj(:,2));
%     
%     temp = NaN(101,1);
%     temp(subj,1) = ddmfit.mean(line_inds);
%     behavior.(model_name).(pars{ipar}) = temp;
%     behavior.params_HDDMbias_YAOA.(pars{ipar})(subj,1) = ddmfit.mean(line_inds);
% %     behavior.eyemem1_params_biasmodel.(pars{ipar})(behavior.eyemem1_params_biasmodel.(pars{ipar}) == 0) = NaN;
%   end
%   
%   if ~any(line_inds); warning('No data found'); continue;  end
% end
% 
% disp 'TODO remove zeros from eyemem1_params_biasmodel!'
% behavior.participants = Participants;



%% Add DDM data YAOA fit together, does not differ from fitting separately 

disp 'load ddm fits'
% model_names = {'params_HDDMbias_YAOA' }; % hddm, young and old together
model_names = {'params_HDDMbias_study' 'params_HDDMbias_test' }; % hddm, young and old together
pars = {'a' 't' 'v' 'z' 'dc'};
agegroup = {'YA' 'OA'};
ddmpath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/HDDM/Niels';
for imodel = 1:length(model_names)
  model_name = model_names{imodel};
  ddmfit = readtable(fullfile(ddmpath, [model_name '.csv'] ));
  
  for iage = 1:2
    behavior.(model_name).(agegroup{iage}) = [];
    behavior.(model_name).(agegroup{iage}).dimord = 'subj';
    
    for ipar = 1:length(pars)
      ddmstr = sprintf('%s_subj(%s).',  pars{ipar}, agegroup{iage});
      line_inds = contains(ddmfit.Var1, ddmstr);
      % get subjID from Var1 strings, to use as index
      subj = ddmfit.Var1(line_inds);
      subj = cellfun(@(x) tokenize(x, '.'), subj, 'UniformOutput', false);
      subj = vertcat(subj{:});
      subj = cellfun(@str2double, subj(:,2));
      
      temp = NaN(101,1);
      temp(subj,1) = ddmfit.mean(line_inds);
      behavior.(model_name).(agegroup{iage}).(pars{ipar}) = temp;
%       behavior.params_HDDMbias_YAOA.(pars{ipar})(subj,1) = ddmfit.mean(line_inds);
      behavior.ddmNiels.(pars{ipar})(subj,1) = ddmfit.mean(line_inds);
      behavior.(pars{ipar})(subj,imodel) = ddmfit.mean(line_inds); % set as primary ddm
    end
  end
  if ~any(line_inds); warning('No data found'); continue;  end
end
behavior.age_groups = {'young' 'old'};
behavior.phases = {'study' 'test'};
behavior.cond = category_labels;

disp 'TODO remove zeros from eyemem1_params_biasmodel!'

%%
disp('Saving')
PREOUT = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior';
disp(fullfile(PREOUT, 'Eyemem_behavior.mat'))
save(fullfile(PREOUT, 'Eyemem_behavior.mat'), 'behavior')

