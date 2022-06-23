function EM_runPLSanalysis(analysisname)
% make model specification txt file and run the PLS analysis
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat

if nargin==0
  analysisname = 'corrSDbold'; % behav PLS vs DDM drift
%   analysisname = 'SDbold_vs_HMAX';  % task PLS
  %    analysisname = 'SDbold_OAvsYA_task' % overall BSV YA vs OA
end
%%
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
    num_perm = '100';
    num_split = '0';
    num_boot = '100';
    boot_type = 'strat';
    clim = '95';
    save_data = '0';
    cormode = '0';
    selected_cond = []; %num2str(ones(1,5)); disp 'TODO get ncond somewhere'
    behavior_data = {};
    behavior_name = {};
    
    agegroups = {'young' 'old'};
       
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
    batch_plsgui(txtfilename)
    
  case 'corrSDbold'
%%
corrtype = 'Pearson'; %Spearman Pearson
% corrtype = 'Spearman'; %Spearman Pearson
%     agegroup = 'OA'; % ALLsubj OA YA
    behavnames = {...
      %       {'study' 'dprime'};
      %       {'study' 'criterion'};
      %       {'study' 'RT'};
      %       {'study' 'RTsd'};
%       {'test' 'dprime'};
      %       {'test' 'RT'}
      %       {'test' 'criterion'};
      %       {'test' 'RTsd'}
      %       {'test' 'p_repeatbalanced'};
%       {'params_HDDMbias_YAOA' 'v'};
%       {'params_HDDMbias_YAOA' 'a'};
%       {'params_HDDMbias_YAOA' 't'};
%       {'params_HDDMbias_YAOA' 'dc'};
%       {'params_HDDMbias_YAOA' 'z'};
%       {'ddmLaura' 'v'};
      {'ddmNiels' 'v'};
%       {'eyemem1_params_biasmodel_Niels' 'v'};
      };
    
    %     basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
    %     PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/ftsource/taskPLS/OAvsYA_SD';
    basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/'; %yesno or 2afc
    %     PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/VOIsel/YA/PLS';
    %     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/VOIsel/OA/PLS';
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/linearfit/gaze-specific'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_3bins/linearfit/gaze-specific'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/uniformbinwidth/linearfit_fitcoeff1/gaze-specific/old';

    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/fixednbins/linearfit_fitcoeff1/gaze-specific/young';

    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/gaze-specific';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/gaze-specific/young';
    
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_3bins/fixednbins/linearfit_fitcoeff1/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_3bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/old'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_3bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_3bins/fixednbins/bin2-bin1_fitcoeff1/gaze-specific'
    
    %     agegroups = {'young' 'old'};
    
%     agegroups = {'young'};
    agegroups = {'old'};
%     agegroups = {''};

    disp 'Generate model txt file'
    %     txtfilename = 'corrSDbold_vsRT_OA_BfMRIanalysis.txt';
    %     resultfilename = 'corrSDbold_vsRT_OA_BfMRIresult.mat';
    outname = analysisname;
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
    
    %%% both groups code
    
    behavior_data =  cell(length(agegroups),1);
    
    id_list = cell(length(agegroups),1);
    ages=table('Size', [1 1], 'VariableTypes', ["string"]);
    for iage = 1:length(agegroups)
      cd(fullfile(PREIN, agegroups{iage}))
      subjlist = dir('sub*_BfMRIsessiondata.mat');
      for isub=1:length(subjlist)
        tmp = tokenize(subjlist(isub).name, '_');
        
        subjind = behavior.participants.participant_id == tmp{1};
        behav_valkeep = [];
        behavior_name = [];
        for ibehav = 1:length(behavnames)
          behavior_name = [behavior_name '  ' behavnames{ibehav}{:}];
          behavoi = behavior.(behavnames{ibehav}{1});
          behav_val = behavoi.(behavnames{ibehav}{2})(subjind,1);
          ages.Var1(isub,1) = behavior.participants.group(find(subjind),:);
          if behav_val == 0
            disp(agegroups{iage})
            warning('zero found!')
            disp(behavior.participants.participant_id(subjind))
            disp('dropping subject'); continue
          end
          behav_valkeep = [behav_valkeep behav_val];
        end

        % could also take from data itself, same result
%         load(subjlist(isub).name, 'behavdata')
%         behav_valkeep = behavdata;
%         behavior_name = 'drift';

%         if behav_val > 0
          behavior_data{iage}{end+1} = num2str(behav_valkeep); %
          id_list{iage}{end+1} = tmp{1};
%         end
      end
    end
    if numel(agegroups) == 2
      cd(PREIN)
    end
    
    outfilename = sprintf('%s_%s_%s_%d_%s', outname, [behavnames{:}{:}], [agegroups{:}], cellfun(@length, id_list), corrtype); %
    disp(outfilename)
    txtfilename = [ outfilename '_BfMRIanalysis.txt'];
    resultfilename = [ outfilename '_BfMRIresult.mat'];

    save ages ages
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
    batch_plsgui(txtfilename)
%%    
  case 'SDbold_vs_HMAX'
    if contains(which('plsgui'), 'rank')
      warning 'PLS_rank does not work with task PLS, switching to regular PLS toolbox'
      rmpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/PLS_rank'))
      addpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/pls'))
    end
%     warning 'PLS_rank does not work with task PLS'
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc
%     PREIN = fullfile(basepath, 'variability', 'ftsource', 'SDbold_vs_HMAX', 'old', '5bins');
%     PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/linearfit/young'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/linearfit/old'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/bin5-bin1/gazespecific'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_10bins/bin5-bin1/gazespecific'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_3bins/bin5-bin1/gazespecific'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/nanstd_3bins/bin5-bin1/gazespecific'
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/nanstd_5bins/bin5-bin1/gazespecific'
    
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/non-gazespecific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_3bins/non-gazespecific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_3bins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_10bins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_10bins/non-gazespecific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_25bins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_150bins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_75bins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/nanstd_3bins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/nanstd_5bins/uniformbinwidth/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/uniformbinwidth/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/uniformbinwidth/gaze-specific';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_9bins/uniformbinwidth/gaze-specific/old'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/gaze-specific/old'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/fixednbins/gaze-specific/young'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/fixednbins/gaze-specific/old'

    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/fixednbins/non-gazespecific/young'
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/iqr_5bins/fixednbins/non-gazespecific/old'
       
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/std_5bins/fixednbins/non-gazespecific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/std_5bins/fixednbins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/std_3bins/fixednbins/gaze-specific/young';
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/std_3bins/fixednbins/gaze-specific';
    
    disp 'Generate model txt file'
%     txtfilename = 'SDbold_vs_HMAX_gazespec_OAvsYA_BfMRIanalysis.txt';
%     resultfilename = 'SDbold_vs_HMAX_gazespec_OAvsYA_BfMRIresult.mat';
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
    % %     id_list = cell(1,1);
    for iage = 1:length(agegroups)
      cd(fullfile(PREIN, agegroups{iage}))
      subjlist = dir('sub*_BfMRIsessiondata.mat');
      outlier_perc_BOLDremoved = zeros(length(subjlist),1)
      for isub=1:length(subjlist)
        tmp = tokenize(subjlist(isub).name, '_');
        id_list{iage}{end+1} = tmp{1};
%         load(subjlist(isub).name, 'perc_BOLDremoved');
%         outlier_perc_BOLDremoved(isub) = perc_BOLDremoved;
      end
      %     disp 'avg BOLD removed in outlier removal:'
      %     mean(perc_BOLDremoved)
    end
    outname = analysisname;
    outfilename = sprintf('%s_%d', outname, cellfun(@length, id_list)); %
    disp(outfilename)
    txtfilename = [ outfilename '_BfMRIanalysis.txt'];
    resultfilename = [ outfilename '_BfMRIresult.mat'];
    
    cd(PREIN)
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
    batch_plsgui(txtfilename)
    
end
