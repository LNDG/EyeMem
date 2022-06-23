
load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/variability/ftsource/taskPLS/std_3bins/fixednbins/gaze-specific/SDbold_vs_HMAX_youngold_45_41_BfMRIresult.mat')

nsub = result.num_subj_lst;
ncond = result.num_conditions;
inds = [1 nsub(1);
  nsub(1)+1  nsub(1)*2;
  nsub(1)+1  nsub(1)*2]

ind_start = linspace(1, nsub(1)*ncond, ncond+1)
ind_end   = linspace(nsub(1), nsub(1)*ncond, ncond)

lv=1;

% age groups in different cells
uscdat{1} = result.usc(1:nsub(1)*ncond,lv);
uscdat{2} = result.usc(nsub(1)*ncond+1:end,lv);

plotdat = {};
ageinds = 0;
for iage =1:2
  inds=0;
  for icond = 1:ncond
    curinds = inds+1:(inds+nsub(iage));
    plotdat{icond, iage} = uscdat{iage}(curinds,:)
    
    inds = curinds(end);
    
  end
end

%%
flipit = 1;
if flipit
  plotdat = cellfun(@(x) x*-1, plotdat, 'UniformOutput', false);
end
ylims = [14.2e4 15.4e4; 11.6e4 12.8e4];
ylims = [1.1e5 2.1e5; 0.9e5 1.9e5];
titleg = {'YA' 'OA'};
figure;
for iage=1:2
  subplot(1,2,iage); hold on
  %   bar(cellfun(@mean, plotdat(:,iage)))
  plotSpread(plotdat(:,iage))
  xdat=nan(2,3);
  ydat=nan(2,3);
  for icond=1:3
    xdat(:,icond) = [icond-0.25 icond+0.25];
    ydat(:,icond) = [mean(plotdat{icond,iage}); mean(plotdat{icond,iage})];
  end
%   if flipit
%     ydat = ydat*-1;
%   end
  line(xdat, ydat, 'linewidth', 4)
  
  % plot single subj lines
  ydat = [plotdat{:,iage}];
%   if flipit
%     ydat = ydat*-1;
%   end
  xdat = repmat([1 2 3], size(ydat,1),1);
  line(xdat', ydat', 'linewidth', 1, 'Color', [0.5 0.5 0.5])
  
  
%   ylim(ylims(iage,:))
  title(titleg{iage})
  ylabel('Brain score')
end
saveas(gcf, 'brainscores.pdf')
saveas(gcf, 'brainscores.png')