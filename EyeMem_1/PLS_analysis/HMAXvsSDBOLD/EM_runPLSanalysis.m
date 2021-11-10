function EM_runPLSanalysis(analysisname)
% make model specification txt file and run the PLS analysis
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat

if nargin==0
%   analysisname = 'corrSDbold_vsdprime';
  %    analysisname = 'SDbold_OAvsYA_task'
  analysisname = 'SDbold_vs_HMAX';
end
switch analysisname
  case 'SDbold_OAvsYA_task'
    if contains(which('plsgui'), 'rank')
      warning 'PLS_rank does not work with task PLS, switching to regular PLS toolbox'
      rmpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/PLS_rank'))
      addpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/pls'))
    end

    PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/ftsource/taskPLS/OAvsYA_SD';
        
    disp 'Generate model txt file'    
    txtfilename = sprintf('%s_BfMRIanalysis.txt', analysisname);
    resultfilename = sprintf('%s_BfMRIresult.mat', analysisname);    
    pls_option = '1';    
    mean_type = '1';
    num_perm = '1000';
    num_split = '0';
    num_boot = '1000';
    boot_type = 'strat';
    clim = '95';
    save_data = '0';
    cormode = '0';
    selected_cond = []; %num2str(ones(1,5)); disp 'TODO get ncond somewhere'
    behavior_data = {};
    behavior_name = {};
    
    agegroups = {'young' 'old'};
       
    cd(fullfile(PREIN))
    subjlist = dir('*_BfMRIsessiondata.mat');
    id_list = cell(1,2);
    for isub=1:length(subjlist)
      tmp = tokenize(subjlist(isub).name, '_');
      subjind = behav.participants.participant_id == tmp{1};
      if behav.participants{subjind, 'group'} == 'old'
        age = 1;
      else
        age = 2;
      end
      id_list{age}{end+1} = tmp{1};
    end
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
    batch_plsgui(txtfilename)
    
  case 'corrSDbold_vsdprime'
    corrtype = 'Spearman'; %Spearman Pearson
    agegroup = 'OA'; % ALLsubj OA YA
    behavnames = {...
      %       {'study' 'dprime'};
      %       {'study' 'criterion'};
      %       {'study' 'RT'};
      %       {'study' 'RTsd'};
      {'test' 'dprime'};
%       {'test' 'RT'}
%       {'test' 'criterion'};
%       {'test' 'RTsd'}
%       {'test' 'p_repeatbalanced'};
      };

%     basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
%     PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/ftsource/taskPLS/OAvsYA_SD';
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
%     PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/VOIsel/YA/PLS';
    PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/VOIsel/OA/PLS';
        
    disp 'Generate model txt file'
    %     txtfilename = 'corrSDbold_vsRT_OA_BfMRIanalysis.txt';
    %     resultfilename = 'corrSDbold_vsRT_OA_BfMRIresult.mat';
    outname = analysisname;
    txtfilename = sprintf('%s%s_%s_BfMRIanalysis.txt', outname, agegroup, corrtype);
    resultfilename = sprintf('%s%s_%s_BfMRIresult.mat', outname, agegroup, corrtype);    
    pls_option = '3';    
    mean_type = '2';
    if strcmp(corrtype, 'Pearson' )
      cormode = '0'; % Pearson
    elseif strcmp(corrtype, 'Spearman' )
      cormode = '8'; % Spearman
    end

    num_perm = '1000';
    num_split = '0';
    num_boot = '1000';
    boot_type = 'strat';
    clim = '95';
    save_data = '0';
    selected_cond = []; %num2str(ones(1,5)); disp 'TODO get ncond somewhere'
    behavior_data = {};
    behavior_name = {};
    
%     agegroups = {'young' 'old'};
    % for iage = 1:2
    %   cd(fullfile(PREIN, agegroups{iage}))
    cd(fullfile(PREIN))
    subjlist = dir('*_BfMRIsessiondata.mat');
    id_list = cell(1,1);         behavior_data_keep = {}; behavior_name_keep = {};

    for ibehav = 1:size(behavnames,1)
      behavior_name = {[behavnames{ibehav}{:}]};
      behavior_name_keep = [behavior_name_keep [behavnames{ibehav}{:}]];
      behavior_data = {};
      
      for isub=1:length(subjlist)
        tmp = tokenize(subjlist(isub).name, '_');
        subjind = behav.participants.participant_id == tmp{1};
        if strcmp(agegroup, 'YA'); 
          agegroup = 'young'; 
        elseif strcmp(agegroup, 'OA'); 
          agegroup = 'old'; 
        end
        if behav.participants{subjind, 'group'} == agegroup || strcmp(agegroup, 'ALLsubj')
          behavoi = behav.(behavnames{ibehav}{1});
          behav_val = nanmean(behavoi(subjind).(behavnames{ibehav}{2}));
          %         behav_val
          %         if behav_val < 0.25 % dprime lower than that is fishy
          %           fprintf('behav_val < 0.25, skipping subject\n')
          %           continue
          %         end
          id_list{1}{end+1} = tmp{1};
          behavior_data{end+1,1} = num2str(behav_val); %
          behavior_data_keep{end+1,ibehav} = num2str(behav_val);
          
        end
      end
    end
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
    batch_plsgui(txtfilename)
    % end
    
  case 'SDbold_vs_HMAX'
    warning 'PLS_rank does not work with task PLS'
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
%     PREIN = fullfile(basepath, 'variability', 'ftsource', 'SDbold_vs_HMAX', 'old', '5bins');
%     PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/linearfit/young'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/linearfit/old'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/bin5-bin1/gazespecific'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_10bins/bin5-bin1/gazespecific'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_3bins/bin5-bin1/gazespecific'
        
    disp 'Generate model txt file'
    txtfilename = 'SDbold_vs_HMAX_gazespec_OAvsYA_BfMRIanalysis.txt';
    resultfilename = 'SDbold_vs_HMAX_gazespec_OAvsYA_BfMRIresult.mat';
    pls_option = '1';
    mean_type = '0';
    cormode = '0';
    num_perm = '1000';
    num_split = '0';
    num_boot = '1000';
    boot_type = 'strat';
    clim = '95';
    save_data = '0';
    selected_cond = num2str(ones(1,5)); disp 'TODO get ncond somewhere'
    behavior_data = {};
    behavior_name = {};
    
    agegroups = {'young' 'old'};
%     agegroups = {'young'};
    id_list = cell(length(agegroups),1);
    for iage = 1:length(agegroups)
      cd(fullfile(PREIN, agegroups{iage}))
%       cd(fullfile(PREIN))
      subjlist = dir('*_BfMRIsessiondata.mat');
      for isub=1:length(subjlist)
        tmp = tokenize(subjlist(isub).name, '_');
        id_list{iage}{end+1} = tmp{1};
      end
    end
    cd(PREIN)
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
    batch_plsgui(txtfilename)
    
    %%
%   case 'behavPLSvsdprime' % load behav for dprime
%     corrtype = 'Pearson';
%     basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
%     %     PREIN = fullfile(basepath, 'variability', 'ftsource', 'behavPLSvsdprime', 'iqr_5bins', 'bin5-bin1', 'young');
%     PREIN = fullfile(basepath, 'variability', 'ftsource', 'behavPLSvsdprime', 'iqr_5bins');
%     agegroups = {'young' 'old'};
%     
%     disp 'Generate model txt file'
%     pls_option = '3'; % behav PLS    
%     mean_type = '0';
%     if strcmp(corrtype, 'Pearson' )
%       cormode = '0'; % Pearson
%     elseif strcmp(corrtype, 'Spearman' )
%       cormode = '8'; % Spearman
%     end
%     num_perm = '1000';
%     num_split = '0';
%     num_boot = '1000';
%     boot_type = 'strat';
%     clim = '95';
%     save_data = '0';
%     selected_cond = num2str(ones(1,1)); disp 'TODO get ncond somewhere'
%     
%     load /Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
%     behavnames = {...
%       %       {'study' 'dprime'};
%       %       {'study' 'criterion'};
%       %       {'study' 'RT'};
%       %       {'study' 'RTsd'};
%       {'test' 'dprime'};
% %       {'test' 'RT'}
% %       {'test' 'criterion'};
% %       {'test' 'RTsd'}
% %       {'test' 'p_repeatbalanced'};
%       };
%     
%     cd(PREIN)
% %     runlist = dir('bin*');
%     runlist = dir('lin*');
%     for irun = 1:length(runlist) % analysis run
%       for iage = 1%:2
%         cd(fullfile(PREIN, runlist(irun).name, agegroups{iage}))
%         disp(pwd)
%         subjlist = dir('sub*_BfMRIsessiondata.mat');
%         behavior_data_keep = {}; behavior_name_keep = {};
%         for ibehav = 1:size(behavnames,1)
%           behavior_name = [behavnames{ibehav}{:}];
%           behavior_name_keep = [behavior_name_keep [behavnames{ibehav}{:}]];
%           behavior_data = {};
%           txtfilename = sprintf('Model_behavPLSvs_%s_%s_BfMRIanalysis.mat', behavior_name, corrtype) ; %  'Model_behavPLSvsdprime_BfMRIanalysis.txt';
%           resultfilename = sprintf('Model_behavPLSvs_%s_%s_BfMRIresult.mat', behavior_name, corrtype) ;
%           id_list = {};
%           for isub=1:length(subjlist)
%             tmp = tokenize(subjlist(isub).name, '_');
%             id_list{end+1} = tmp{1};
%             subjind = behav.participants.participant_id == tmp{1};
%             %       behavior_data{end+1,1} = num2str(mean(behav.test(subjind).dprime));
%             %         behavior_data{end+1,1} = num2str(mean(behav.test(subjind).dprime));
%             %         getfield(behav.(behavnames{ibehav}{1}), behavnames{ibehav}{2});
%             
%             behavoi = behav.(behavnames{ibehav}{1});
%             behavior_data{isub,1} = num2str(mean(behavoi(subjind).(behavnames{ibehav}{2}))); %
%             behavior_data_keep{isub,ibehav} = num2str(mean(behavoi(subjind).(behavnames{ibehav}{2}))); %
%           end
%           %           PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,{behavior_name})
%           %           batch_plsgui(txtfilename)
%         end
%         disp 'model with all behav at once'
%         txtfilename = sprintf('Model_behavPLSvs_%s_%s_BfMRIanalysis.mat', [behavior_name_keep{:}], corrtype) ; %  'Model_behavPLSvsdprime_BfMRIanalysis.txt';
%         resultfilename = sprintf('Model_behavPLSvs_%s_%s_BfMRIresult.mat',[behavior_name_keep{:}], corrtype) ;
%         PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,...
%           clim,save_data,selected_cond,behavior_data_keep,behavior_name_keep)
%         batch_plsgui(txtfilename)
%       end
%     end
end
