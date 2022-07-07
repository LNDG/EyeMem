function EM_corrSDboldvsbehav_plot(result)
% TODO, plot nice scatters YA and OA with the group PLS model
if nargin==0
  load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/dprime/linearfit_fitcoeff1/corrSDbold_dprime_youngold_45_41_Pearson_BfMRIresult.mat')
end
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
% % load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific/ages.mat'
% load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/linearfit/gaze-specific/ages.mat'
% % load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/ages.mat'

load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/dprime/linearfit_fitcoeff1/ages.mat')
% load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_3bins/fixednbins/linearfit_fitcoeff1/gaze-specific/ages.mat

% age = behavior.participants.group;
brainscores = result.usc;
behav = result.stacked_behavdata;

t = table(brainscores, behav);

%% plot
agegroups= {'young' 'old'};
% agecolors = {'b' 'r'};
agecolors = [0 0 1; 1 0 0];
corrtype = 'Pearson'; % 'Spearman'; Pearson
f=figure; hold on;
clear l
for iage=1:2%:2
%   sc(iage)=scatter(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'MarkerFaceColor', agecolors(iage,:));
%   [r(iage),p(iage)]=corr(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'type', corrtype);
  sc(iage) = scatter(t.behav(t.Var3 == agegroups{iage}), t.brainscores(t.Var3 == agegroups{iage}), 'MarkerFaceColor', agecolors(iage,:));
  [r(iage),p(iage)] = corr(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'type', corrtype);
end
l = lsline;

l(1).Color = agecolors(2,:);
l(2).Color = agecolors(1,:);
legend(sc,agegroups, 'Location', 'Southeast')
title(sprintf('young: r = %1.2f, p = %1.3f\nold: r = %1.2f, p =  %1.3f', r(1), p(1), r(2), p(2)))
box on; axis square

ranked=0;
if ranked
  ylabel('brainscore (ranked)')
  xlabel('drift rate (ranked)')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold_ranked.pdf')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold_ranked.png')
else
  ylabel('brainscore')
  xlabel('drift rate')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold.pdf')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold.png')
end


%% ANOVA

anovan(brainscores, {ages.Var1 behav},'model',2, 'continuous', 2, 'varnames',{'age','drift rate'})

