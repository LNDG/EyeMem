function EM_regressage(result)
% account for age group while correlating e.g. SDbold vs dprime
% 1 load in data from all subjects
% 2 load in dprime data
% 3 partialcorr controlling for age group
% save residual SDbold and residual dprime for use in subsequent PLS

path = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/eyemem/variability/ftsource/taskPLS/OAvsYA_SD';
if nargin == 0
  load(fullfile(path, 'corrSDbold_vsdprime_ALLsubj_BfMRIresult.mat'))
end
PREOUT = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/eyemem_analysis/plots';
SAV = 0;
%%
r=[];
r(1) = corr(result.usc(1:43,1), result.stacked_behavdata(1:43), 'type', 'Spearman');
r(2) = corr(result.usc(44:end,1), result.stacked_behavdata(44:end), 'type', 'Spearman');
r(3) = corr(result.usc(:,1), result.stacked_behavdata, 'type', 'Spearman');

%%
close all
figure; hold on; axis square; box on; clear l s
s(1) = scatter(result.usc(1:43,1), result.stacked_behavdata(1:43), 'b');
l = lsline; 
s(2) = scatter(result.usc(44:end,1), result.stacked_behavdata(44:end), 'r');
l = lsline; 
s(3) = scatter(result.usc(:,1), result.stacked_behavdata(:), 'w.');
l = lsline; 
l(1).Color = 'b'; l(2).Color = 'r'; l(3).Color = 'k';

% legend(l, {'OA' 'YA' 'ALL'}); %legend boxoff
legend(l, {sprintf('OA r=%1.2f', r(1)), sprintf('YA r=%1.2f', r(2)), sprintf('ALL r=%1.2f', r(3))}, 'Location', 'Northwest' )
legend boxoff

% xlabel('SD-BOLD (ranked?)')
% ylabel('Dprime (ranked)')
xlabel('SD-BOLD')
ylabel('Dprime')
if SAV
  outfile = sprintf('corr_SDboldvsDprime.png')
  saveas(gcf, fullfile(PREOUT, outfile ))
  cd(PREOUT)
end

