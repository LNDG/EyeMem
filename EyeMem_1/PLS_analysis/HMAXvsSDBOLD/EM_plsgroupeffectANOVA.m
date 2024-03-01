function EM_plsgroupeffectANOVA(result)
% split pls result into age groups, test interaction
PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/EyeMem/plots';

if nargin==0
  %   result = load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific/corrSDbold_ddmNielsv__87_Spearman_BfMRIresult.mat')
  %   result = load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific/corrSDbold_ddmNielsv__87_Pearson_BfMRIresult.mat')
%     result = load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/linearfit/gaze-specific/corrSDbold_ddmNielsv__86_Pearson_BfMRIresult.mat')
%   result = load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_3bins/linearfit/gaze-specific/corrSDbold_ddmNielsv__86_Spearman_BfMRIresult.mat')
  %unibinwidth:
%   load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/uniformbinwidth/linearfit_fitcoeff1/gaze-specific/corrSDbold_ddmNielsv__86_Pearson_BfMRIresult.mat')
%   load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/corrSDbold_ddmNielsv__88_Pearson_BfMRIresult.mat')
%   load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/gaze-specific/corrSDbold_ddmNielsv__86_Pearson_BfMRIresult.mat')
%   load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/dprime/linearfit_fitcoeff1/corrSDbold_dprime__86_80_earson_BfMRIresult.mat')
% load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsDDM/v/linearfit_fitcoeff1/corrSDbold_v__86_83_pearman_BfMRIresult.mat'
% load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__85_80_earson_BfMRIresult.mat')
% load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1_psc/corrSDbold_v__86_80_earson_BfMRIresult.mat')

% load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__87_80_earson_BfMRIresult.mat')
% load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__85_80_earson_BfMRIresult.mat')

load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/corrSDbold_v__87_83_pearman_BfMRIresult.mat')

end
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
% % load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific/ages.mat'
% load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/linearfit/gaze-specific/ages.mat'
% % load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/ages.mat'

% load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/dprime/linearfit_fitcoeff1/ages.mat')
% load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_3bins/fixednbins/linearfit_fitcoeff1/gaze-specific/ages.mat
% load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsDDM/v/linearfit_fitcoeff1/ages.mat
% load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1_psc/ages.mat
load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/total_pow/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1/ages.mat')

% age = behavior.participants.group;
brainscores = result.usc;
behav = result.stacked_behavdata;

t = table(brainscores, behav, ages.Var1);

%% plot
agegroups= {'young' 'old'};
% agecolors = {'b' 'r'};
agecolors = [1 0 0; 0 0 1; 0 0 0];
corrtype = 'Pearson'; % 'Spearman'; Pearson
% corrtype = 'Spearman'; % ''; Pearson
f=figure; f.Position = [1000        1130         200         220];
hold on;
clear l
for iage=1:2%:2
  sc(iage) = scatter(t.behav(t.Var3 == agegroups{iage}), t.brainscores(t.Var3 == agegroups{iage}), 'MarkerFaceColor', agecolors(iage,:), 'MarkerEdgeColor', 'w', 'sizedata', 50);
  [r(iage),p(iage)] = corr(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'type', corrtype);
end
% sc(3) = scatter(t.behav, t.brainscores, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'sizedata', 1);
[r(3),p(3)] = corr(t.brainscores, t.behav, 'type', corrtype);
% [rho, pval]=partialcorr(result.usc(:,1), behavdata_lst{1}, agegroup)

l = lsline;

% l(1).Color = agecolors(3,:);
% l(2).Color = agecolors(2,:);
% l(3).Color = agecolors(1,:);
% l(1).Color = agecolors(3,:);
l(1).Color = agecolors(2,:);
l(2).Color = agecolors(1,:);

% agegroups= {'Young' 'Older' 'All'};
agegroups= {'Young' 'Older'};
legend(sc, agegroups, 'Location', 'Northwest'); legend boxon
title(sprintf('Young: r = %1.2f, p = %1.3f\nOlder: r = %1.2f, p =  %1.3f\nAll: r = %1.2f, p =  %1.3f', r(1), p(1), r(2), p(2), r(3), p(3)))
box on; axis square

ranked=0;
if ranked
  ylabel('brainscore (ranked)')
  xlabel('drift rate (ranked)')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold_ranked.pdf')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold_ranked.png')
else
  ylabel('bin 5-1 SDbold')
  xlabel('DDM drift')
%   xlim([0 1.3])
  saveas(f, fullfile(PREOUT, sprintf('group_behavpls_youngold_%s.pdf', corrtype )))
  saveas(f, fullfile(PREOUT, sprintf('group_behavpls_youngold_%s.eps', corrtype )), 'epsc')
  saveas(f, fullfile(PREOUT, sprintf('group_behavpls_youngold_%s.png', corrtype )))
%     saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold.png')
end


%% ANOVA

anovan(brainscores, {ages.Var1 behav},'model',2, 'continuous', 2, 'varnames',{'age','drift rate'})

% %% Doug meeting 29022024: Order not correct, file not stacked
% load('corrSDbold_v__86_80_earson_BfMRIresult.mat')
% agegroup = [zeros(44,1); ones(42,1)];
% [rho, pval]=partialcorr(result.usc(:,1), behavdata_lst{1}, agegroup)
% [rho, pval]=corr(result.usc(1:44,1), behavdata_lst{1}(1:44,:), 'type', 'Spearman')
% [rho, pval]=corr(result.usc(1:44,1), behavdata_lst{1}(1:44,:), 'type', 'Pearson')
% [rho, pval]=corr(result.usc(45:end,1), behavdata_lst{1}(45:end,:), 'type', 'Pearson')

