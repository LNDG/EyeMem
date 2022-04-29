function PLSmodeltxtfilegenerator(txtfilename,resultfilename,id_list,pls_option,mean_type,cormode,num_perm,num_split,num_boot,boot_type,clim,save_data,selected_cond,behavior_data,behavior_name)

%Native generation of PLS model txt files from variables and parameters of interest. Avoids manual manipulation of txt files!
%
%txtfilename: string-based name of the PLS model txt file you want to generate (e.g., 'ModelXX_BfMRIanalysis.txt')
%
%resultfilename: string-based name of model output (e.g., 'ModelXX_BfMRIresult.mat')
%
%id_list: cell array of string-based IDs in your preferred order (e.g., {'id1','id2'}).
%
%pls_option: PLS model type. 1. Mean-Centering PLS; 2. Non-Rotated Task PLS; 3. Regular Behav PLS, 4. Multiblock PLS, 5. Non-Rotated Behav PLS, 6.Non-Rotated Multiblock PLS (e.g., '3')
%
%mean_type: Mean-Centering Type. 0. Remove group condition means from conditon means within each group, 1. Remove grand condition means from each group condition mean, 2. Remove grand mean over all subjects and conditions, 3. Remove all main effects by subtracting condition and group means. (e.g., '1')
%
%cormode: Correlation Mode. 0. Pearson correlation, 2. covariance, 4. cosine angle, 6. dot product (e.g., '0')
%
%num_perm: # permutations (e.g., '1000')
%
%num_split: # of "Natasha" split halfs (e.g., '0')
%
%num_boot: # bootstrap samples (e.g., '1000')
%
%boot_type: Bootstrap type. Either 'strat' or 'nonstrat' (e.g., 'strat')
%
%clim: Bootstrap confidence Level (e.g., '95')
%
%save_data: Set to '1' to save stacked datamat.
%
%selected_cond: binary string of conditions of interest. (e.g., if three conds and want only the first cond, then '1 0 0')
%
%behavior_data: cell array of string data for behaviours of interest in columnar format (e.g., single behavioral variable for 2 subjects = {'13';'65'}
%
%behavior_name: name of behaviours of interest (e.g., {'fluid intelligence','crystallized intelligence'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open txt file for writing
fid = fopen(txtfilename,'w');

%now write text lines for PLS model files
fprintf(fid,'\n');
fprintf(fid,'%%------------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  Result File Name Start  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%  Note: Result file must be listed first, and must follow the file\n');
fprintf(fid,'%%	 name format of xxxx_yyyyresult.mat, where xxxx stands for\n');
fprintf(fid,'%%	 "any result file name prefix" and yyyy stands for the name\n');
fprintf(fid,'%%	 of PLS module (either PET ERP fMRI BfMRI STRUCT or SmallFC).\n');
fprintf(fid,'%%	 File name is case sensitive on  Unix or Linux computers.\n');
fprintf(fid,'\n');
fprintf(fid,['result_file   ',resultfilename,'\n']);
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  Result File Name End  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%------------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  Group Section Start  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
for igroup = 1:length(id_list)
  %Append file suffix to ID list so don't have to include in id list
  %itself...
  id_list_append = {};
  for i = 1:length(id_list{igroup})
    id_list_append{i} = [id_list{igroup}{i},'_BfMRIsessiondata.mat'];
  end
  fprintf(fid,'%s ','group_files   ',id_list_append{:});
  fprintf(fid,'\n\n');
end
fprintf(fid,'\n');
fprintf(fid,'%% ... following above pattern for more groups\n');
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  Group Section End  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%------------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  PLS Section Start  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%  Notes:\n');
fprintf(fid,'%%    1. Mean-Centering PLS\n');
fprintf(fid,'%%    2. Non-Rotated Task PLS (please also fill out contrast data below)\n');
fprintf(fid,'%%    3. Regular Behav PLS (please also fill out behavior data & name below)\n');
fprintf(fid,'%%    4. Multiblock PLS (please also fill out behavior data & name below)\n');
fprintf(fid,'%%    5. Non-Rotated Behav PLS (please also fill out contrast data and\n');
fprintf(fid,'%%	behavior data & name below)\n');
fprintf(fid,'%%    6. Non-Rotated Multiblock PLS (please also fill out contrast data and\n');
fprintf(fid,'%%	behavior data & name below)\n');
fprintf(fid,'\n');
fprintf(fid,['pls   ',pls_option,'	%% PLS Option (between 1 to 6, see above notes)\n']);
fprintf(fid,'\n');
fprintf(fid,'%%  Mean-Centering Type:\n');
fprintf(fid,'%%    0. Remove group condition means from conditon means within each group\n');
fprintf(fid,'%%    1. Remove grand condition means from each group condition mean\n');
fprintf(fid,'%%    2. Remove grand mean over all subjects and conditions\n');
fprintf(fid,'%%    3. Remove all main effects by subtracting condition and group means\n');
fprintf(fid,'\n');
fprintf(fid,['mean_type     ',mean_type,'		%% Mean-Centering Type (between 0 to 3, see above)\n']);
fprintf(fid,'\n');
fprintf(fid,'%%  Correlation Mode:\n');
fprintf(fid,'%%    0. Pearson correlation\n');
fprintf(fid,'%%    2. covaraince\n');
fprintf(fid,'%%    4. cosine angle\n');
fprintf(fid,'%%    6. dot product\n');
fprintf(fid,'\n');
fprintf(fid,['cormode	',cormode,'		%% Correlation Mode (can be 0,2,4,6, see above)\n']);
fprintf(fid,'\n');
fprintf(fid,['num_perm	',num_perm,'		%% Number of Permutation\n']);
fprintf(fid,['num_split	',num_split,'		%% Natasha Perm Split Half\n']);
fprintf(fid,['num_boot	',num_boot,'		%% Number of Bootstrap\n']);
fprintf(fid,['boot_type	',boot_type,'		%% Either strat or nonstrat bootstrap type\n']);
fprintf(fid,['clim		',clim,'            %% Confidence Level for Behavior PLS\n']);
fprintf(fid,['save_data	',save_data,'		%% Set to 1 to save stacked datamat\n']);
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  PLS Section End  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%------------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  Condition Selection Start  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%  Notes: If you don''t need to deselect conditions, just leave\n');
fprintf(fid,'%%  "selected_cond" and "selected_bcond" to be commented.\n');
fprintf(fid,'\n');
fprintf(fid,'%%  First put k number of 1 after "selected_cond" keyword, where k is the\n');
fprintf(fid,'%%  number of conditions in sessiondata file. Then, replace with 0 for\n');
fprintf(fid,'%%  those conditions that you would like to deselect for any case except\n');
fprintf(fid,'%%  behavior block of multiblock PLS. e.g. If you have 3 conditions in\n');
fprintf(fid,'%%  sessiondata file, and you would like to deselect the 2nd condition,\n');
fprintf(fid,'%%  then you should enter 1 0 1 after selected_cond.\n');
fprintf(fid,'%%\n');
fprintf(fid,['selected_cond	',selected_cond,' \n']);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'	%%  Condition Selection End  %%\n');
fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%------------------------------------------------------------------------\n');
fprintf(fid,'\n');
if ~isempty(behavior_data)
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'	%%  Behavior Data Start  %%\n');
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'\n');
  fprintf(fid,'%%  Notes: only list selected conditions (selected_cond)\n');
  fprintf(fid,'\n');
  fprintf(fid,'%%behavedata\n');
  fprintf(fid,'\n');
  %Now repmat 'behavior_data' convention for PLS in front of each behavioural data line entry.
  label = {'behavior_data'};
  label_rep = repmat(label,size(behavior_data,1),1);
  combo = horzcat(label_rep,behavior_data)';%bring data together
  % s_rep = [repmat('%s ',1,size(combo,1)),'\n'];
  % %now write whole behavioral data cell array in single fprintf command and
  % %the on to rest of lines to write...
  % fprintf(fid,s_rep,combo{:});%%%%%
  
  for igroup = 1:length(id_list)
    %Append file suffix to ID list so don't have to include in id list
    %itself...
    id_list_append = {};
    for i = 1:length(id_list{igroup})
      %     id_list_append{i} = [id_list{igroup}{i},'_BfMRIsessiondata.mat'];
      fprintf(fid,'behavior_data %s\n', combo{2,igroup}{i});
    end
  end
  fprintf(fid,'\n\n');
  
  fprintf(fid,'\n');
  fprintf(fid,'%% ... following above pattern for more groups\n');
  fprintf(fid,'\n');
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'	%%  Behavior Data End  %%\n');
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'\n');
  fprintf(fid,'%%------------------------------------------------------------------------\n');
  fprintf(fid,'\n');
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'	%%  Behavior Name Start  %%\n');
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'\n');
  fprintf(fid,'%%  Numbers of Behavior Name should match the Behavior Data above\n');
  fprintf(fid,'\n');
  fprintf(fid,'%s ','behavior_name	',behavior_name);
  fprintf(fid,'\n');
  fprintf(fid,'\n');
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'	%%  Behavior Name End  %%\n');
  fprintf(fid,'	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'\n');
  fprintf(fid,'%%------------------------------------------------------------------------\n');
  fprintf(fid,'\n');
end

%now close file
fclose(fid);
