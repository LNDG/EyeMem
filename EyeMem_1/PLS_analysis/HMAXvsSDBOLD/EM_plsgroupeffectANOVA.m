function EM_plsgroupeffectANOVA(result)
% split pls result into age groups, test interaction
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific/ages.mat'

% age = behavior.participants.group;
brainscores = result.result.usc;
behav=result.result.stacked_behavdata;

t = table(brainscores, behav, ages.Var1);

%% plot
agegroups= {'young' 'old'};
% agecolors = {'b' 'r'};
agecolors = [0 0 1; 1 0 0];
corrtype = 'Pearson'; % 'Spearman'; Pearson
f=figure; hold on;
clear l
for iage=1:2%:2
  sc(iage)=scatter(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'MarkerFaceColor', agecolors(iage,:));
  [r(iage),p(iage)]=corr(t.brainscores(t.Var3 == agegroups{iage}), t.behav(t.Var3 == agegroups{iage}), 'type', corrtype);
end
l = lsline;

l(1).Color = agecolors(2,:);
l(2).Color = agecolors(1,:);
legend(sc,agegroups, 'Location', 'Southeast')
title(sprintf('young: r = %1.2f, p = %1.3f\nold: r = %1.2f, p =  %1.3f', r(1), p(1), r(2), p(2)))
box on; axis square
xlabel('brainscore (ranked)')
ylabel('dprime (ranked)')
  
saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold.pdf')
saveas(f, '/Users/kloosterman/gridmaster2012/projectdata/eyemem/plots/group_behavpls_youngold.png')



