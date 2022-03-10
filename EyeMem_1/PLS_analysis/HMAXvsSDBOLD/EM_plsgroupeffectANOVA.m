function EM_plsgroupeffectANOVA(result)
% split pls result into age groups, test interaction
load /Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat
load '/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/behavPLSvsDDM/nanstd_5bins/linearfit/gazespecific/ages.mat'

% age = behavior.participants.group;
brainscores = result.result.usc;
behav=result.result.stacked_behavdata;

t = table(brainscores, behav, ages )





