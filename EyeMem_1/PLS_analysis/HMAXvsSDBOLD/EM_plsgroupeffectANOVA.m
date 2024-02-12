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
load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1_psc/corrSDbold_v__86_80_earson_BfMRIresult.mat')

end
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
% % load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific/ages.mat'
% load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/iqr_5bins/linearfit/gaze-specific/ages.mat'
% % load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_5bins/fixednbins/linearfit_fitcoeff1/non-gazespecific/ages.mat'

% load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsSDT/dprime/linearfit_fitcoeff1/ages.mat')
% load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/std_3bins/fixednbins/linearfit_fitcoeff1/gaze-specific/ages.mat
% load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/std_3bins/fixednbins/behavPLSvsDDM/v/linearfit_fitcoeff1/ages.mat
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/variability2/5TRspertrial/ftsource/induced/std_5bins/fixednbins/behavPLSvsDDM/v/gaze-specific/bin5-bin1_fitcoeff1_psc/ages.mat

% age = behavior.participants.group;
brainscores = result.usc;
behav = result.stacked_behavdata;

t = table(brainscores, behav, ages.Var1);

%% plot
agegroups= {'young' 'old'};
% agecolors = {'b' 'r'};
agecolors = [1 0 0; 0 0 1;];
corrtype = 'Pearson'; % 'Spearman'; Pearson
f=figure; f.Position = [1000        1130         200         220];
hold on;
clear l
for iage=1:2%:2
%   sc(iage)=scatter(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'MarkerFaceColor', agecolors(iage,:));
%   [r(iage),p(iage)]=corr(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'type', corrtype);
  sc(iage) = scatter(t.behav(t.Var3 == agegroups{iage}), t.brainscores(t.Var3 == agegroups{iage}), 'MarkerFaceColor', agecolors(iage,:), 'MarkerEdgeColor', 'w', 'sizedata', 50);
  [r(iage),p(iage)] = corr(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'type', corrtype);
end
l = lsline;

l(1).Color = agecolors(2,:);
l(2).Color = agecolors(1,:);
agegroups= {'Young' 'Older'};
legend(sc, agegroups, 'Location', 'Northwest'); legend boxon
title(sprintf('Young: r = %1.2f, p = %1.3f\nOlder: r = %1.2f, p =  %1.3f\n', r(1), p(1), r(2), p(2)))
box on; axis square

ranked=0;
if ranked
  ylabel('brainscore (ranked)')
  xlabel('drift rate (ranked)')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold_ranked.pdf')
  %   saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold_ranked.png')
else
  ylabel('HMAX slope')
  xlabel('dprime')
  saveas(f, fullfile(PREOUT, 'group_behavpls_youngold.pdf'))
  saveas(f, fullfile(PREOUT, 'group_behavpls_youngold.eps'), 'epsc')
  saveas(f, fullfile(PREOUT, 'group_behavpls_youngold.png'))
%     saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold.png')
end


%% ANOVA

anovan(brainscores, {ages.Var1 behav},'model',2, 'continuous', 2, 'varnames',{'age','drift rate'})

