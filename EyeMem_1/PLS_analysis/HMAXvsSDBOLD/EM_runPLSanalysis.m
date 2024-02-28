function EM_runPLSanalysis(analysisname)
% make model specification txt file and run the PLS analysis
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat

if nargin==0
%   analysisname = 'corrSDbold'; % behav PLS vs DDM drift
  analysisname = 'SDbold_vs_HMAX';  % task PLS
end
%%
switch analysisname
    case 'SDbold_vs_HMAX'
    if contains(which('plsgui'), 'rank')
      warning 'PLS_rank does not work with task PLS, switching to regular PLS toolbox'
      rmpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/PLS_rank'))
      addpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/pls'))
    end
    %     warning 'PLS_rank does not work with task PLS'
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/'; %yesno or 2afc    
    PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/mean_5bins/fixednbins/taskPLS/non-gazespecific';
    
    %     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_10bins/fixednbins/taskPLS/gaze-specific';
    
    %mse!
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/mse_5bins/fixednbins/taskPLS/gaze-specific';
%     % look at induced power, by subtracting mean per 5TR trial
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_3bins/fixednbins/taskPLS/gaze-specific'
    cd(PREIN)
    agegroups = {'young' 'old'};
%     agegroups = {'young'};
%     agegroups = {'old'};
%     agegroups = {''};

    disp 'Generate model txt file'
%     txtfilename = 'SDbold_vs_HMAX_gazespec_OAvsYA_BfMRIanalysis.txt';
%     resultfilename = 'SDbold_vs_HMAX_gazespec_OAvsYA_BfMRIresult.mat';
    pls_option = '1';
    mean_type = '0';
    cormode = '0';
    num_perm = '100';
    num_split = '0';
    num_boot = '100';
    boot_type = 'strat';
    clim = '95';
    save_data = '0';
    selected_cond = num2str(ones(1,5)); disp 'TODO get ncond somewhere'
    behavior_data = {};
    behavior_name = {};
    
    id_list = cell(length(agegroups),1);
    for iage = 1:length(agegroups)
      cd(fullfile(PREIN, agegroups{iage}))
      subjlist = dir('sub*_BfMRIsessiondata.mat');
      for isub=1:length(subjlist)
        tmp = tokenize(subjlist(isub).name, '_');
        id_list{iage}{end+1} = tmp{1};
      end
    end
    outname = analysisname;

    outfilename = sprintf('%s_%s_%d_%d', outname, [agegroups{:}], cellfun(@length, id_list)'); %
    disp(outfilename)
    txtfilename = [ outfilename '_BfMRIanalysis.txt'];
    resultfilename = [ outfilename '_BfMRIresult.mat'];
    
    cd(PREIN)
    PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
    batch_plsgui(txtfilename)

  case 'corrSDbold'
    %%
    corrtype = 'Pearson'; %Spearman Pearson
%     corrtype = 'Spearman'; %Spearman Pearson
    
    if strcmp(corrtype, 'Spearman')
      addpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/PLS_rank'))
      rmpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/pls'))
    end
    
    basepath = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/'; %yesno or 2afc
    % new dir structure naming
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/dprime/linearfit_fitcoeff1';
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/criterion/linearfit_fitcoeff1';
%     PREIN = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/RT/linearfit_fitcoeff1';
        
%     params = {'v' 'a' 't' 'dc' 'z'};
    params = {'v' };
%     gazetype = 'non-gazespecific';
    gazetype = 'gaze-specific';
%     PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/ftsource/std_3bins/fixednbins/behavPLSvsDDM';
    
%     PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/1TRspertrial/ftsource/std_3bins/fixednbins/behavPLSvsDDM'
%     PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/std_3bins/fixednbins/behavPLSvsDDM';
%     PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/totalSD/std_3bins/fixednbins/behavPLSvsDDM/';
%     % induced power: correlation gone
%     PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_3bins/fixednbins/behavPLSvsDDM'

% PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/std_3bins/fixednbins/behavPLSvsDDM'
% PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/std_5bins/fixednbins/behavPLS_sdboldvsHmaxbins/'
PRE = '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/mean_5bins/fixednbins/behavPLSvsDDM/'

%     agegroups = {'young' 'old'};
%     agegroups = {'young'};
%     agegroups = {'old'};
        agegroups = {''};
    
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
    
    num_perm = '100';
    num_split = '0';
    num_boot = '100';
    boot_type = 'strat';
    clim = '95';
    save_data = '0';
    selected_cond = []; %num2str(ones(1,5)); disp 'TODO get ncond somewhere'
    
    for iparam = 1:length(params)
%       PREIN = fullfile(PRE, params{iparam}, gazetype, 'linearfit_fitcoeff1');
%       PREIN = fullfile(PRE, params{iparam}, gazetype, 'corrHmaxoverbins_fitcoeff1');
      PREIN = fullfile(PRE, params{iparam}, gazetype, 'bin5-bin1_fitcoeff1');
      cd(PREIN)
      behavior_data =  cell(length(agegroups),1);
      id_list = cell(length(agegroups),1);
      ages=table('Size', [1 1], 'VariableTypes', ["string"]);
      for iage = 1:length(agegroups)
        cd(fullfile(PREIN, agegroups{iage}))
        subjlist = dir('sub*_BfMRIsessiondata.mat');
        for isub=1:length(subjlist)
          tmp = tokenize(subjlist(isub).name, '_');
          subjind = behavior.participants.participant_id == tmp{1};
          % take behav from data
          load(subjlist(isub).name, 'behavdata', 'behavname')
          if isnan(behavdata)
            disp(subjlist(isub).name)
            continue
          end
          behav_valkeep = behavdata;
          behavior_name = behavname{:};
          
          ages.Var1(isub,1) = behavior.participants.group(find(subjind),:);
          
          %         if behav_val > 0
          behavior_data{iage}{end+1} = num2str(behav_valkeep); %
          id_list{iage}{end+1} = tmp{1};
          %         end
        end
      end
      if numel(agegroups) == 2
        cd(PREIN)
      end
      
      outfilename = sprintf('%s_%s_%s_%d_%d_%s', outname, behavior_name, [agegroups{:}], cellfun(@length, id_list), corrtype); %
      disp(outfilename)
      txtfilename = [ outfilename '_BfMRIanalysis.txt'];
      resultfilename = [ outfilename '_BfMRIresult.mat'];
      
      save ages ages
      save id_list id_list
      PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)
      batch_plsgui(txtfilename)
    end

    
end
