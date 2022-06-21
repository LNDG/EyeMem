
load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/std_3bins/fixednbins/gaze-specific/SDbold_vs_HMAX_44)__BfMRIresult.mat')

nsub = result.num_subj_lst;
ncond = result.num_conditions;
inds = [1 nsub(1); 
  nsub(1)+1  nsub(1)*2;
  nsub(1)+1  nsub(1)*2]

ind_start = linspace(1, nsub(1)*ncond, ncond+1)
ind_end   = linspace(nsub(1), nsub(1)*ncond, ncond)

lv=1;
inds=0;
iage=1;
plotdat = {};

ageinds = 0;
for iage =1:2
  
  uscdat = result.usc(ageinds,lv)
  for icond = 1:ncond
    curinds = inds+1:(inds+1+nsub(iage));
    plotdat{iage, icond} = uscdat(curinds,:)
    
    inds = curinds(end);
    
  end
end