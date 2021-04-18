%% plot results of first EyeMem PLS analysis

result = load('/Users/kloosterman/gridmaster2012/LNDG/EyeMem/PLS/SDdatamats/EYEMEM_TaskPLS_1000P_500B_BfMRIresult.mat');

nsub = result.result.num_subj_lst;
ncond=result.result.num_conditions;

f = figure;
f.Position = [2117 642 905 276];
for iLV = 1:2
    
    LVdatY = reshape(result.result.usc(1:nsub(1)*ncond, iLV), nsub(1), ncond ); % dimord subj cond 
    LVdatO = reshape(result.result.usc(nsub(1)*ncond+1:end, iLV), nsub(2), ncond ); 
    
    meandat = [];
    meandat(1,:) = mean(LVdatY, 1); % young
    meandat(2,:) = mean(LVdatO, 1); % old
    
    semdat = [];
    semdat(1,:) = std(LVdatY, 0, 1) / sqrt(nsub(1)); % young
    semdat(2,:) = std(LVdatO, 0, 1) / sqrt(nsub(2)); % old

    subplot(1,2,iLV)
    barweb(meandat, semdat, 0.5, {'YOUNG', 'OLDER'});
    title(sprintf('LV %d', iLV))
    if iLV == 2; legend(result.cond_name, 'Location', 'NorthWest'); end
    ax = gca;
    ax.FontSize = 12;
end