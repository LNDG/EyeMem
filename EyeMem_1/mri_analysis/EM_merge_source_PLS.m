function result = EM_merge_source_PLS(source)
% Behav PLS vs dprime across space and time.

% All datamats must be in the form of "subject in condition".
% rows subjects, voxels columns

result = cell(3,3);
for iage = 1%
  for imeas = 3
    inside = source.source(1).inside;
%     agegroup = source.ageleg{iage};
    if iage == 3      
      datamat_lst = vertcat(source.source(1:2, imeas, 1).pow);
      datamat_lst = datamat_lst(:,inside,:);
    else      
      datamat_lst = source.source(iage, imeas).pow(:,inside,:);
    end
    datamat_lst = datamat_lst(:,:);
    
    option = [];
    option.method = 3; % [1] | 2 | 3 | 4 | 5 | 6
    option.num_perm = 100; %( single non-negative integer )
    % option.is_struct = [0] | 1
    % option.num_split = 100 %( single non-negative integer )
    option.num_boot = 100; %500 % ( single non-negative integer )
    % option.clim = ( [95] single number between 0 and 100 )
    % option.bscan = ( subset of  1:num_cond )
    % option.stacked_designdata = ( 2-D numerical matrix )
    
    num_subj_lst = length(source.SUBJ{iage});
    behavdat = NaN(num_subj_lst, 1);
    for isub = 1:num_subj_lst
      subj_id = source.behavior.participants.participant_id == source.SUBJ{iage}{isub};
      behavdat(isub,1) =  nanmean(source.behavior.test(subj_id).dprime);
    end
    
    option.stacked_behavdata = behavdat; %( 2-D numerical matrix )
    
    % option.meancentering_type = [0] | 1 | 2 | 3
%     option.cormode = 8; % [0] | 2 | 4 | 6
    option.cormode = 8; % [0] | 2 | 4 | 6
    % option.boot_type = ['strat'] | 'nonstrat'
    
    num_cond = 1;
    result{iage, imeas} = pls_analysis({datamat_lst}, num_subj_lst, num_cond, option)
    
  end
end

% plot scatter of correlation

% integrate r-values using bootstrap ratios as mask

% 
